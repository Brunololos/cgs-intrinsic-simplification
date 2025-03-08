#ifndef HEAT_DIFF_H
#define HEAT_DIFF_H

#include "geom_helpers.hpp"
#include "isimp_utils.hpp"

struct heatDiffData {
    Eigen::Matrix<double, -1, 1> initial_heat;
    Eigen::Matrix<double, -1, 1> final_heat;

    // Laplacian // TODO: store laplacian for multiple different simplification states
    Eigen::MatrixXd M;
    Eigen::MatrixXd L;
};

void initHeat(heatDiffData& hdData, const iSimpData& iSData);
void diffuse_heat(heatDiffData& hdData, const int steps);

// build laplacian with current intrinsic mesh
void build_cotan_mass(heatDiffData& hdData, const iSimpData& iSData);
void build_cotan_laplacian(heatDiffData& hdData, const iSimpData& iSData);

double barycentric_lumped_mass(const iSimpData& iSData, const int vertex_idx);
double cotan_weight(const iSimpData& iSData, const int edge_idx);

#endif