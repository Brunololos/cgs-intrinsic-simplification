#ifndef HEAT_DIFF_H
#define HEAT_DIFF_H

#include <random>
#include "heat_diff_data.hpp"
#include "geom_helpers.hpp"
#include "isimp_utils.hpp"

enum class CircumcentricShape { ZERO, TRIANGLE, CHOPPED_TRIANGLE, RECTANGLE, CHOPPED_QUADRILATERAL };

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