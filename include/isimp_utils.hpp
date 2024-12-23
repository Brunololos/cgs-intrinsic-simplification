#ifndef ISIMP_UTILS_H
#define ISIMP_UTILS_H

#include <utility>

#include "tinyad_defs.hpp"
#include "print.hpp"

#include "isimp_data.hpp"
// register a point to be tracked through the simplification
void register_point(iSimpData& iSData, const BarycentricPoint point, const int face_idx);
void map_registered(iSimpData& iSData);
// TODO: implement simplification operations
// returns whether the edge flip could be performed/was successful
bool flip_intrinsic(iSimpData& iSData, const gcs::Edge edge);
void flatten_vertex(iSimpData& iSData, const int vertex_idx);
void remove_vertex(iSimpData& iSData, const int vertex_idx);

std::array<int, 5> order_quad_edge_indices(const gcs::Edge& edge, const int Fk_unique_v, const int Fl_unique_v);
std::array<int, 4> order_quad_vertex_indices(const gcs::Edge& edge);
std::array<int, 4> order_quad_vertex_indices(const gcs::Edge& edge, const int Fk_unique_v, const int Fl_unique_v);

// returns faceindex and value
std::pair<int, int> find_unique_vertex_index(const gcs::Face& F, const gcs::Edge& edge);

#endif