#pragma once

#include "isimp_data.hpp"
#include "heat_diff_data.hpp"
#include "isimp_utils.hpp"

void initMesh(const Eigen::Matrix<double, -1, 3>& V, const Eigen::Matrix<int, -1, 3>& F, iSimpData& data);
bool iSimp_step(iSimpData& data);

// recover the intrinsic triangulation (connectivity & edge lengths) for a specific simplification state
// e.g. to calculate something on them (They are written into data.recoveredMesh & data.recovered_L)
void recover_intrinsics_at(const int mapping_idx, iSimpData& data);
void recover_heat_and_intrinsics_at(const int mapping_idx, heatDiffData& hdData, iSimpData& data);