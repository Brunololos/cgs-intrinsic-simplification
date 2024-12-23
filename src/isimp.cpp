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

  data.inputMesh.reset(new gcs::ManifoldSurfaceMesh(F));
  data.inputGeometry.reset(new gcs::VertexPositionGeometry(*data.inputMesh, V));

  // initialize intrinsic connectivity
  data.intrinsicMesh.reset(new gcs::ManifoldSurfaceMesh(F));

  // initialize intrinsic edge lengths
  int nEdges = data.inputMesh->nEdges();
  data.L = Eigen::VectorXd::Zero(nEdges);
  int i = 0;
  std::cout << "Calculating L:\n"; // TODO: remove
  for (gcs::Edge e : data.inputMesh->edges())
  {
    data.L(e.getIndex()) = data.inputGeometry->edgeLength(e);
    /* TODO: remove */ std::cout << "edge(" << e.firstVertex().getIndex() << ", " << e.secondVertex().getIndex() << ") L(" << i << ") = " << data.L(i) << std::endl;
    i++;
  }

  data.mapping = std::vector<std::unique_ptr<Mapping_operation>>();
  data.tracked_points = std::vector<BarycentricPoint>();
  data.tracked_by_triangle = std::vector<std::vector<int>>(data.F.rows());
}