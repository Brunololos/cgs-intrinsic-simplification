#include "isimp.hpp"

#include <iostream>
#include <fstream>

#include <igl/parallel_for.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/polar_svd3x3.h>

#include <Eigen/Sparse>

typedef Eigen::SparseVector<double>::InnerIterator SVIter;

void initMesh(const Eigen::Matrix<double, -1, 3>& V, const Eigen::Matrix<int, -1, 3>& F, iSimpData& data){
  data.V = V;
  data.F = F;
  data.M = Eigen::MatrixXd::Zero(V.rows(), 2);
  data.T_minus = Eigen::MatrixXd::Zero(V.rows(), 2);
  data.T_plus = Eigen::MatrixXd::Zero(V.rows(), 2);
  data.Q = std::set<std::shared_ptr<std::pair<int, double>>, customLess>();
  data.Q_elems = std::vector<std::pair<int, double>>(V.rows());

  data.inputMesh.reset(new gcs::ManifoldSurfaceMesh(F));
  data.inputGeometry.reset(new gcs::VertexPositionGeometry(*data.inputMesh, V));

  // initialize intrinsic connectivity
  data.intrinsicMesh.reset(new gcs::ManifoldSurfaceMesh(F));
  data.tangent_space_reference_edges = Eigen::MatrixXi::Zero(V.rows(), 1);

  // initialize intrinsic edge lengths
  int nEdges = data.intrinsicMesh->nEdges();
  data.L = Eigen::VectorXd::Zero(nEdges);
  int i = 0;
  // std::cout << "Calculating L:\n" << std::endl; // TODO: remove
  for (gcs::Edge e : data.inputMesh->edges())
  {
    data.L(e.getIndex()) = data.inputGeometry->edgeLength(e);
    // /* TODO: remove */ std::cout << "edge(" << e.firstVertex().getIndex() << ", " << e.secondVertex().getIndex() << ") L(" << i << ") = " << data.L(i) << std::endl;
    // if(data.L(e.getIndex()) == 0.0) {
    //   std::cout << RED << "edge(" << e.firstVertex().getIndex() << ", " << e.secondVertex().getIndex() << ") L(" << i << ") = " << data.L(i) << RESET << std::endl;
    // }
    i++;
  }

  // init masses & error vectors
  for (gcs::Vertex v : data.intrinsicMesh->vertices())
  {
    // calculate curvature
    double curvature = gaussian_curvature(data, v);
    double K_minus = - std::min(curvature, 0.0);
    double K_plus = std::max(curvature, 0.0);
    data.M.row(v.getIndex()) = Eigen::Vector2d(K_minus, K_plus);

    for (gcs::Edge e : v.adjacentEdges())
    {
      data.tangent_space_reference_edges(v.getIndex(), 0) = e.getIndex();
      break;
    }
  }

  // insert vertices with intrinsic curvature error into priority queue
  for (gcs::Vertex v : data.intrinsicMesh->vertices())
  {
    double ice = intrinsic_curvature_error(data, v);
    std::pair<int, double> v_elem = std::pair<int, double>(v.getIndex(), ice);
    data.Q_elems[v.getIndex()] = v_elem;
    data.Q.insert(std::make_shared<std::pair<int, double>>(v_elem));
  }

  data.mapping = std::vector<std::unique_ptr<Mapping_operation>>();
  data.tracked_points = std::vector<BarycentricPoint>();
  data.tracked_texcoords = std::vector<TexCoord>();
  data.tracked_by_triangle = std::vector<std::vector<int>>(data.F.rows());
  data.hasConverged = false;
}