#pragma once

#include <Eigen/Core>
#include <Eigen/Dense>
#include <vector>
#include "geometrycentral/surface/signpost_intrinsic_triangulation.h"
#include <igl/min_quad_with_fixed.h>

#include "isimp_utils.hpp"

namespace gc = geometrycentral;
namespace gcs = gc::surface;

struct iSimpData {
  // input (original vertices & faces)
  Eigen::Matrix<double, -1, 3> V;
  Eigen::Matrix<int, -1, 3> F;

  // geometries
  std::unique_ptr<gcs::ManifoldSurfaceMesh> inputMesh;
  std::unique_ptr<gcs::VertexPositionGeometry> inputGeometry;
};
