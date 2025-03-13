#ifndef HEAT_DIFF_H
#define HEAT_DIFF_H

#include "geom_helpers.hpp"
#include "isimp_utils.hpp"

enum class CircumcentricShape { ZERO, TRIANGLE, CHOPPED_TRIANGLE, RECTANGLE, CHOPPED_QUADRILATERAL };

struct heatDiffData {
    Eigen::Matrix<double, -1, 1> initial_heat;
    Eigen::Matrix<double, -1, 1> final_heat;

    // ordered list of vertex indices, selected as a subset from the vertices of the original mesh
    // NOTE: this is to ensure a consistent ordering & allows for calculating heat diff on a simplified mesh,
    // while being able to map the results back to the original mesh via the original vertex indices.
    std::vector<int> heat_sample_vertices;
    std::unordered_map<int, int> vertex_idcs_to_heat_idcs;

    // Laplacian // TODO: store laplacian for multiple different simplification states
    Eigen::MatrixXd M;
    Eigen::MatrixXd L;
};

void initHeatDiff(heatDiffData& hdData, const iSimpData& iSData);
void initHeat(heatDiffData& hdData, const iSimpData& iSData);
void build_heat_sample_indices(heatDiffData& hdData, const iSimpData& iSData);
void diffuse_heat(heatDiffData& hdData, const double stepsize, const int steps);

double get_heat(heatDiffData& hdData, const int vertex_idx);

// build laplacian with current intrinsic mesh
void build_cotan_mass(heatDiffData& hdData, const iSimpData& iSData);
void build_cotan_laplacian(heatDiffData& hdData, const iSimpData& iSData);

// laplacian helpers
double barycentric_lumped_mass(const iSimpData& iSData, const int vertex_idx);
double circumcentric_dual_mass(const iSimpData& iSData, const int vertex_idx);
double cotan_weight(const iSimpData& iSData, const int edge_idx);

// Verification
double total_mass(const heatDiffData& hdData);

#endif