#pragma once
#include <igl/min_quad_with_fixed.h>
#include <queue>
#include <set>

#include "geom_helpers.hpp"
#include "gcs_defs.hpp"
#include "print.hpp"

enum class SIMP_OP;
class Mapping_operation;

struct customLess
{
    bool operator()(const std::shared_ptr<std::pair<int, double>> l, const std::shared_ptr<std::pair<int, double>> r) const {
      if (l->second == r->second)
      {
        return l->first < r->first;
      }
      return l->second < r->second;
    }
};

struct iSimpData {
  // input (original vertices & faces)
  Eigen::Matrix<double, -1, 3> V;
  Eigen::Matrix<int, -1, 3> F;

  // edge lengths in the intrinsic domain
  Eigen::Matrix<double, -1, 1> L;

  // intrinsic simplification quantities
  // Vertex masses (M_minus, M_plus)
  Eigen::Matrix<double, -1, 2> M;
  // Error vectors t_minus, t_plus for positive/negative curvature in polar form (first coordinate: normalized angle in [0, 1], second coordinate: vector length)
  Eigen::Matrix<double, -1, 2> T_minus;
  Eigen::Matrix<double, -1, 2> T_plus;
  // NOTE: I needed to use std::set for the priority queue, because std::priority_queue only allows access to the highest priority element (And we also need to update others that are queued)
  std::set<std::shared_ptr<std::pair<int, double>>, customLess> Q;
  std::vector<std::pair<int, double>> Q_elems; // vertex index, intrinsic curvature error pairs used in the priority queue Q

  // geometries
  std::unique_ptr<gcs::ManifoldSurfaceMesh> inputMesh;
  std::unique_ptr<gcs::VertexPositionGeometry> inputGeometry;
  std::unique_ptr<gcs::ManifoldSurfaceMesh> intrinsicMesh;

  // edges that define the angle 0 of the tangent spaces, centered on a vertex
  Eigen::Matrix<int, -1, 1> tangent_space_reference_edges;

  // bijective mapping consisting of a sequence of atomic mapping operations
  std::vector<std::unique_ptr<Mapping_operation>> mapping;
  // points tracked through the bijective mapping from the fine to the coarse triangulation
  std::vector<BarycentricPoint> tracked_points;
  std::vector<BarycentricPoint> mapped_points;
  // mapping from point indices to their respective texture coordinates
  std::vector<TexCoord> tracked_texcoords;
  // TODO: implement (This is the inverse reference from tracked by triangle, which maps from triangles to points. This vector maps points to triangles.)
  // the index mapping of the barycentric points to their triangle in the intrinsic simplification
  // std::vector<int> tracked_to_triangle;
  // std::vector<int> mapped_to_triangle;
  // the index of the outer vector denotes to which triangle the inner list of barycentric point indices belongs to
  std::vector<std::vector<int>> tracked_by_triangle;
  std::vector<std::vector<int>> mapped_by_triangle;
  bool hasConverged;
};

// types of simplification operations
enum class SIMP_OP {
    V_FLATTEN = 0,
    E_FLIP = 1,
    V_REMOVAL = 2
};

class Mapping_operation {
    public:
        SIMP_OP op_type;

        virtual void map_to_coarse(std::vector<BarycentricPoint>& tracked_points, std::vector<std::vector<int>>& tracked_by_triangle)
        {
            std::cout << "ALERT! Called map_to_coarse of superclass!" << std::endl;
        }
};

// intrinsic flip:
// => save normalized normal and bias
// => save triangle corner vertices in the unfolding
// => in which triangle a barycentric coordinate falls depends on the linear equation (n^Tx + b >= 0), where b is the point in the unfolding corresponding to the barycentric coordinate of the point to be transformed.

// mapping for intrinsic edge flips
class Edge_Flip : public Mapping_operation {
  public:
    int Fk_idx;
    int Fl_idx;

    // mapping parameters
    int Fk_unique_v_idx;
    int Fl_unique_v_idx;
    Point2D v0, v1, v2, v3;
    Normal2D n;
    Bias b;

    // in case I want to write order permutations into tables to remove the switch case
    // std::array<std::array<Point2D, 3>, 3> Fk_v_order;
    // std::array<std::array<Point2D, 3>, 3> Fl_v_order;

    Edge_Flip(const int Fk_idx_, const int Fl_idx_, const Quad2D quad_, const int Fk_unique_v_idx_, const int Fl_unique_v_idx_)
    {
      op_type = SIMP_OP::E_FLIP;
      Fk_idx = Fk_idx_;
      Fl_idx = Fl_idx_;

      Fk_unique_v_idx = Fk_unique_v_idx_;
      Fl_unique_v_idx = Fl_unique_v_idx_;
      v0 = quad_.row(0); v1 = quad_.row(1); v2 = quad_.row(2); v3 = quad_.row(3);
      // printEigenMatrixXd("v0", v0);
      // printEigenMatrixXd("v1", v1);
      // printEigenMatrixXd("v2", v2);
      // printEigenMatrixXd("v3", v3);

      Vector2D P3P2 = v2 - v3;
      n = Normal2D(P3P2(1), -P3P2(0)).normalized();
      b = - n.dot(v2);
      // printEigenMatrixXd("v2", v2);
      // printEigenMatrixXd("v3", v3);
      // printEigenMatrixXd("p3p2", P3P2);
      // printEigenMatrixXd("n", n);
      // std::cout << "b: " << b << std::endl;
      // std::cout << "Fk_unique_idx = " << Fk_unique_v_idx;
      // std::cout << ", Fl_unique_idx = " << Fl_unique_v_idx << std::endl;
    }

    bool lies_in_Fk(const Point2D& p) {
      // std::cout << "Checking if point p: " + to_str(p) + " is inside Fk. n * p = " << n.dot(p) << " + b: " << b << std::endl;
      return n.dot(p) + b >= 0;  // TODO: check if we need >= 0 or <= 0
    }

    void map_to_coarse(std::vector<BarycentricPoint>& tracked_points, std::vector<std::vector<int>>& tracked_by_triangle)
    {
      // copy tracked_by_triangle index lists for Fk, Fl and clear them
      std::vector<int> Fk_points = tracked_by_triangle[Fk_idx];
      std::vector<int> Fl_points = tracked_by_triangle[Fl_idx];
      tracked_by_triangle[Fk_idx].clear();
      tracked_by_triangle[Fl_idx].clear();

      // TODO: remove prints
      // std::cout << "\n\nMapping to coarse!" << std::endl;
      // std::cout << "Fk_idx = " << Fk_idx;
      // std::cout << ", Fl_idx = " << Fl_idx << std::endl;
      // std::cout << "Planar unfolding of edge lengths:\n";
      // printEigenVector2d("v0", v0);
      // printEigenVector2d("v1", v1);
      // printEigenVector2d("v2", v2);
      // printEigenVector2d("v3", v3);
      // reinsert the barycentrically reexpressed vertices into the their correct new face Fk or Fl
      for (int p_idx : Fk_points)
      {
        Point2D p;
        // TODO: remove prints
        // std::cout << "Mapping Fk point " << to_str(tracked_points[p_idx]) << " through intrinsic edge flip." << std::endl;
        // std::cout << "Fk_unique_idx: " << Fk_unique_v_idx << std::endl;
        // TODO: create table and replace this switch case with a single to_explicit call with table lookups
        switch (Fk_unique_v_idx)
        {
          // This is what I thought it should be like: 
          case 0:
            // std::cout << "making explicit (201), barycentric point " << to_str(tracked_points[p_idx]) << " as convex combination of A: " << to_str(v2) << ", B: " << to_str(v0) << " and C: " << to_str(v1) << std::endl;
            p = to_explicit(tracked_points[p_idx], v2, v0, v1);
            break;
          case 1:
            // std::cout << "making explicit (120), barycentric point " << to_str(tracked_points[p_idx]) << " as convex combination of A: " << to_str(v1) << ", B: " << to_str(v2) << " and C: " << to_str(v0) << std::endl;
            p = to_explicit(tracked_points[p_idx], v1, v2, v0);
            break;
          default:
          case 2:
            // std::cout << "making explicit (012), barycentric point " << to_str(tracked_points[p_idx]) << " as convex combination of A: " << to_str(v0) << ", B: " << to_str(v1) << " and C: " << to_str(v2) << std::endl;
            p = to_explicit(tracked_points[p_idx], v0, v1, v2);
            break;
        }
        if(tracked_points[p_idx][1] < tracked_points[p_idx][0]) {
          // std::cout << "explicit point p: " << to_str(p) << ", from barycentric: " << to_str(tracked_points[p_idx]) << std::endl;
        }
        if (lies_in_Fk(p))
        {
          // std::cout << "p lies in Fk. Reexpressing p barycentrically: " << std::endl;
          tracked_points[p_idx] = to_barycentric(p, v2, v3, v1);
          tracked_by_triangle[Fk_idx].push_back(p_idx);
        } else {
          // std::cout << "p lies in Fl. Reexpressing p barycentrically: " << std::endl;
          tracked_points[p_idx] = to_barycentric(p, v3, v2, v0);
          tracked_by_triangle[Fl_idx].push_back(p_idx);
        }
      }
      for (int p_idx : Fl_points)
      {
        Point2D p;
        // TODO: create table and replace this switch case with a single to_explicit call with table lookups
        switch (Fl_unique_v_idx)
        {
          case 0:
            p = to_explicit(tracked_points[p_idx], v3, v1, v0);
            break;
          case 1:
            p = to_explicit(tracked_points[p_idx], v0, v3, v1);
            break;
          default:
          case 2:
            p = to_explicit(tracked_points[p_idx], v1, v0, v3);
            break;
        }
        if (lies_in_Fk(p))
        {
          tracked_points[p_idx] = to_barycentric(p, v2, v3, v1);
          tracked_by_triangle[Fk_idx].push_back(p_idx);
        } else {
          tracked_points[p_idx] = to_barycentric(p, v3, v2, v0);
          tracked_by_triangle[Fl_idx].push_back(p_idx);
        }
      }
    }
};

// mapping for vertex flattening
class Vertex_Flattening : public Mapping_operation {
  public:
    // mapping parameters
    // mapped faces
    std::vector<int> F_idxs;
    // indices of flattened vertex in mapped faces
    std::vector<int> v_idxs;
    // int v_idx;
    // scaling factor resulting from removal
    double v_u;

    Vertex_Flattening(const std::vector<int> F_idxs_, const std::vector<int> v_idxs_, const double v_u_)
    {
      op_type = SIMP_OP::V_FLATTEN;
      F_idxs = F_idxs_;
      v_idxs = v_idxs_;
      v_u = v_u_;
    }

    void map_to_coarse(std::vector<BarycentricPoint>& tracked_points, std::vector<std::vector<int>>& tracked_by_triangle)
    {
      for (int i=0; i<F_idxs.size(); i++)
      {
        int F_idx = F_idxs[i];
        int F_v_idx = v_idxs[i];
        for (int j=0; j<tracked_by_triangle[F_idx].size(); j++)
        {
          // scale correct barycentric coordinate & renormalize
          int p_idx = tracked_by_triangle[F_idx][j];
          BarycentricPoint& bp = tracked_points[p_idx];
          bp[F_v_idx] *= std::exp(v_u);
          double norm = bp.sum();
          bp /= norm;
        }
      }
    }
};

// mapping for vertex removal
class Vertex_Removal : public Mapping_operation {
  public:
    // mapping parameters
    std::vector<int> F_idxs;
    int F_res_idx;
    int F_res_permutation;

    int F_jk_shared_v_idx;
    int F_jl_shared_v_idx;
    int F_kl_shared_v_idx;
    Point2D vi, vj, vk, vl;


    // Vertex_Removal(const std::vector<int> F_idxs_, const int F_res_idx_, const Quad2D quad_, const int F_jk_shared_v_idx_, const int F_jl_shared_v_idx_, const int F_kl_shared_v_idx_)
    Vertex_Removal(const std::vector<int> F_idxs_, const int F_res_idx_, const Quad2D quad_, const std::vector<int> shared_idxs, const int F_res_permutation_)
    {
      op_type = SIMP_OP::V_REMOVAL;
      F_idxs = F_idxs_;
      F_res_idx = F_res_idx_;
      F_res_permutation = F_res_permutation_;

      // F_jk_shared_v_idx = F_jk_shared_v_idx_;
      // F_jl_shared_v_idx = F_jl_shared_v_idx_;
      // F_kl_shared_v_idx = F_kl_shared_v_idx_;
      F_jk_shared_v_idx = shared_idxs[0];
      F_jl_shared_v_idx = shared_idxs[1];
      F_kl_shared_v_idx = shared_idxs[2];
      vi = quad_.row(0); vj = quad_.row(1); vk = quad_.row(2); vl = quad_.row(3);
    }

    void map_to_coarse(std::vector<BarycentricPoint>& tracked_points, std::vector<std::vector<int>>& tracked_by_triangle)
    {
      // copy tracked_by_triangle index lists for Fk, Fl and clear them
      std::vector<int> Fjk_points = tracked_by_triangle[F_idxs[0]];
      std::vector<int> Fjl_points = tracked_by_triangle[F_idxs[1]];
      std::vector<int> Fkl_points = tracked_by_triangle[F_idxs[2]];
      tracked_by_triangle[F_idxs[0]].clear();
      tracked_by_triangle[F_idxs[1]].clear();
      tracked_by_triangle[F_idxs[2]].clear();

      for (int p_idx : Fjk_points)
      {
        Point2D p;
        // TODO: create table and replace this switch case with a single to_explicit call with table lookups
        switch (F_jk_shared_v_idx)
        {
          default:
          case 0:
            p = to_explicit(tracked_points[p_idx], vi, vj, vk);
            break;
          case 1:
            p = to_explicit(tracked_points[p_idx], vk, vi, vj);
            break;
          case 2:
            p = to_explicit(tracked_points[p_idx], vj, vk, vi);
            break;
        }

        switch(F_res_permutation)
        {
          default:
          case 0:
            tracked_points[p_idx] = to_barycentric(p, vj, vk, vl);
            break;
          case 1:
            tracked_points[p_idx] = to_barycentric(p, vl, vj, vk);
            break;
          case 2:
            tracked_points[p_idx] = to_barycentric(p, vk, vl, vj);
            break;
        }
        tracked_by_triangle[F_res_idx].push_back(p_idx);
      }

      for (int p_idx : Fjl_points)
      {
        Point2D p;
        // TODO: create table and replace this switch case with a single to_explicit call with table lookups
        switch (F_jl_shared_v_idx)
        {
          default:
          case 0:
            p = to_explicit(tracked_points[p_idx], vi, vl, vj);
            break;
          case 1:
            p = to_explicit(tracked_points[p_idx], vj, vi, vl);
            break;
          case 2:
            p = to_explicit(tracked_points[p_idx], vl, vj, vi);
            break;
        }

        switch(F_res_permutation)
        {
          default:
          case 0:
            tracked_points[p_idx] = to_barycentric(p, vj, vk, vl);
            break;
          case 1:
            tracked_points[p_idx] = to_barycentric(p, vl, vj, vk);
            break;
          case 2:
            tracked_points[p_idx] = to_barycentric(p, vk, vl, vj);
            break;
        }
        tracked_by_triangle[F_res_idx].push_back(p_idx);
      }

      for (int p_idx : Fkl_points)
      {
        Point2D p;
        // TODO: create table and replace this switch case with a single to_explicit call with table lookups
        switch (F_kl_shared_v_idx)
        {
          default:
          case 0:
            p = to_explicit(tracked_points[p_idx], vi, vk, vl);
            break;
          case 1:
            p = to_explicit(tracked_points[p_idx], vl, vi, vk);
            break;
          case 2:
            p = to_explicit(tracked_points[p_idx], vk, vl, vi);
            break;
        }

        switch(F_res_permutation)
        {
          default:
          case 0:
            tracked_points[p_idx] = to_barycentric(p, vj, vk, vl);
            break;
          case 1:
            tracked_points[p_idx] = to_barycentric(p, vl, vj, vk);
            break;
          case 2:
            tracked_points[p_idx] = to_barycentric(p, vk, vl, vj);
            break;
        }
        tracked_by_triangle[F_res_idx].push_back(p_idx);
      }
    }
};