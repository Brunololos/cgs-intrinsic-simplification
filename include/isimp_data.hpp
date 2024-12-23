#pragma once
#include <igl/min_quad_with_fixed.h>

#include "geom_helpers.hpp"
#include "gcs_defs.hpp"
#include "print.hpp"

enum class SIMP_OP;
class Mapping_operation;

struct iSimpData {
  // input (original vertices & faces)
  Eigen::Matrix<double, -1, 3> V;
  Eigen::Matrix<int, -1, 3> F;

  // edge lengths in the intrinsic domain
  Eigen::Matrix<double, -1, 1> L;

  // intrinsic simplification quantities
  Eigen::Matrix<double, -1, 2> Masses;
  Eigen::Matrix<double, -1, 3> Err_Vecs; // TODO: check if error vectors are supposed to be 2d intrinic or 3d extrinsic

  // geometries
  // std::unique_ptr<gcs::ManifoldSurfaceMesh> inputMesh;
  std::unique_ptr<gcs::ManifoldSurfaceMesh> inputMesh;
  std::unique_ptr<gcs::VertexPositionGeometry> inputGeometry;
  std::unique_ptr<gcs::ManifoldSurfaceMesh> intrinsicMesh;

  // bijective mapping consisting of a sequence of atomic mapping operations
  std::vector<std::unique_ptr<Mapping_operation>> mapping;
  // points tracked through the bijective mapping from the fine to the coarse triangulation
  std::vector<BarycentricPoint> tracked_points;
  std::vector<BarycentricPoint> mapped_points;
  // TODO: implement (This is the inverse reference from tracked by triangle, which maps from triangles to points. This vector maps points to triangles.)
  // the index mapping of the barycentric points to their triangle in the intrinsic simplification
  std::vector<int> tracked_to_triangle;
  std::vector<int> mapped_to_triangle;
  // the index of the outer vector denotes to which triangle the inner list of barycentric point indices belongs to
  std::vector<std::vector<int>> tracked_by_triangle;
  std::vector<std::vector<int>> mapped_by_triangle;
};

// types of simplification operations
enum class SIMP_OP {
    V_FLATTEN = 0,
    E_FLIP = 1,
    V_REMOVAL = 2
};

// TODO:
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
    int Fk_unique_v;
    int Fl_unique_v;
    Point2D v0, v1, v2, v3;
    Normal2D n;
    Bias b;

    // in case I want to write order permutations into tables to remove the switch case
    // std::array<std::array<Point2D, 3>, 3> Fk_v_order;
    // std::array<std::array<Point2D, 3>, 3> Fl_v_order;

    Edge_Flip(const int Fk_idx_, const int Fl_idx_, const Quad2D quad_, const int Fk_unique_v_, const int Fl_unique_v_)
    {
      op_type = SIMP_OP::E_FLIP;
      Fk_idx = Fk_idx_;
      Fl_idx = Fl_idx_;

      Fk_unique_v = Fk_unique_v_;
      Fl_unique_v = Fl_unique_v_;
      v0 = quad_.row(0); v1 = quad_.row(1); v2 = quad_.row(2); v3 = quad_.row(3);

      Vector2D P3P2 = v2 - v3;
      n = Normal2D(P3P2(1), -P3P2(0)).normalized();
      b = - n.dot(v2);
    }

    bool lies_in_Fk(const Point2D& p) {
       return n.dot(p) + b >= 0;  // TODO: check if we need >= 0 or <= 0
    }

    void map_to_coarse(std::vector<BarycentricPoint>& tracked_points, std::vector<std::vector<int>>& tracked_by_triangle)
    {
      // copy tracked_by_triangle index lists for Fk, Fl and clear them
      std::vector<int> Fk_points = tracked_by_triangle[Fk_idx];
      std::vector<int> Fl_points = tracked_by_triangle[Fl_idx];
      tracked_by_triangle[Fk_idx].clear();
      tracked_by_triangle[Fl_idx].clear();

      std::cout << "\n\nMapping to coarse!" << std::endl;
      std::cout << "Fk_idx = " << Fk_idx;
      std::cout << ", Fl_idx = " << Fl_idx << std::endl;
      std::cout << "Planar unfolding of edge lengths:\n";
      printEigenVector2d("v0", v0);
      printEigenVector2d("v1", v1);
      printEigenVector2d("v2", v2);
      printEigenVector2d("v3", v3);
      // reinsert the barycentrically reexpressed vertices into the their correct new face Fk or Fl
      for (int p_idx : Fk_points)
      {
        Point2D p;
        std::cout << "Mapping Fk point " << to_str(tracked_points[p_idx]) << " through intrinsic edge flip." << std::endl;
        std::cout << "Fk_unique_idx: " << Fk_unique_v << std::endl;
        // TODO: create table and replace this switch case with a single to_explicit call with table lookups
        switch (Fk_unique_v)
        {
          case 0:
            p = to_explicit(tracked_points[p_idx], v2, v0, v1);
            break;
          case 1:
            p = to_explicit(tracked_points[p_idx], v1, v2, v0);
            break;
          default:
          case 2:
            p = to_explicit(tracked_points[p_idx], v0, v1, v2);
            break;
        }
        std::cout << "explicit point p from barycentric: " << to_str(p) << std::endl;
        if (lies_in_Fk(p))
        {
          std::cout << "p lies in Fk. Reexpressing p barycentrically: " << std::endl;
          tracked_points[p_idx] = to_barycentric(p, v2, v3, v1);
          tracked_by_triangle[Fk_idx].push_back(p_idx);
        } else {
          std::cout << "p lies in Fl. Reexpressing p barycentrically: " << std::endl;
          tracked_points[p_idx] = to_barycentric(p, v3, v2, v0);
          tracked_by_triangle[Fl_idx].push_back(p_idx);
        }
      }
      for (int p_idx : Fl_points)
      {
        Point2D p;
        // TODO: create table and replace this switch case with a single to_explicit call with table lookups
        switch (Fl_unique_v)
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
    // TODO:
    // mapped faces
    std::vector<int> F_idxs;
    // resulting face
    int Fr_idx;
    // removed vertex
    int v_idx;
    // scaling factor resulting from removal
    double scaling;

    Vertex_Flattening(const std::vector<int> F_idxs_, const int Fr_idx_, const int v_idx_, const double scaling_)
    {
      F_idxs = F_idxs_;
      Fr_idx = Fr_idx_;
      v_idx = v_idx_;
      scaling = scaling_;
    }

    void map_to_coarse(std::vector<BarycentricPoint>& tracked_points, std::vector<std::vector<int>>& tracked_by_triangle)
    {
      for (int i=0; i<F_idxs.size(); i++)
      {
        int F_idx = F_idxs[i];
        for (int j=0; j<tracked_by_triangle[F_idx].size(); j++)
        {
          // TODO: scale correct barycentric coordinate and insert all datapoints into resulting face
        }
      }
    }
};

// mapping for vertex removal
class Vertex_Removal : public Mapping_operation {
  public:
    // mapping parameters
    // TODO:
    int F1_idx;
    void map_to_coarse(std::vector<BarycentricPoint>& tracked_points, std::vector<std::vector<int>>& tracked_by_triangle)
    {

    }
};