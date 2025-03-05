#pragma once

#include "isimp_data.hpp"
#include "isimp_utils.hpp"

void initMesh(const Eigen::Matrix<double, -1, 3>& V, const Eigen::Matrix<int, -1, 3>& F, iSimpData& data);
bool iSimp_step(iSimpData& data);