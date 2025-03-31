#ifndef ISIMP_UTILS_H
#define ISIMP_UTILS_H

#include <utility>
#include <unordered_map>

#include "tinyad_defs.hpp"
#include "print.hpp"

#include "isimp_data.hpp"
#include "heat_diff_data.hpp"

// register a point to be tracked through the simplification
void register_point(iSimpData& iSData, const BarycentricPoint point, const TexCoord texcoord, const int face_idx);
void map_registered(iSimpData& iSData, const int n_vertices);
void map_registered(iSimpData& iSData);
// DANGER! These functions expect a specific matching state of iSData & mapping_idx to work properly.
// void map_current_from(iSimpData& iSData, const int from_n_vertices);
// void map_current_from_to(iSimpData& iSData, const int from_n_vertices, const int n_vertices);
void map_current_from(iSimpData& iSData, const int from_mapping_idx);
void map_current_from_to(iSimpData& iSData, const int from_mapping_idx, const int to_mapping_idx);

// during the simplification heat is redistributed.
// for the result of the heat diffusion on simplified triangulations
// to be comparable to the result of heat diffusion on the fine triangulations
// we need to compute the reverse mapping.
void unredistribute_heat(iSimpData& iSData, Eigen::VectorXd& heat);
double circumcentric_dual_mass(const iSimpData& iSData, const int vertex_idx);

// simplification operations (returns whether the operation could be performed/was successful)
bool is_edge_flippable(const iSimpData& iSData, const Quad2D& quad, const gcs::Edge edge);
double compute_flattening_factor(const Eigen::VectorXd& L, const gcs::Vertex v, const bool verbose);
bool flip_intrinsic(iSimpData& iSData, const gcs::Edge edge, const bool verbose=false);
bool flatten_vertex(iSimpData& iSData, const int vertex_idx, const bool verbose=false);
bool remove_vertex(iSimpData& iSData, const int vertex_idx, const bool verbose=false);
void flip_to_delaunay(iSimpData& iSData, const bool verbose=false);
void flip_to_delaunay(iSimpData& iSData, const std::vector<gcs::Vertex>& vertices, const bool verbose=false);
bool flip_vertex_to_deg3(iSimpData& iSData, const int vertex_idx, const bool verbose=false);
void rotate_reference_edge(iSimpData& iSData, const gcs::Vertex& vertex);

PolarVector2D next_error_vector(const iSimpData& iSData, const double alpha, const gcs::Vertex vertex, const gcs::Vertex neighbor, bool for_T_minus);
double intrinsic_curvature_error(const iSimpData& iSData, const gcs::Vertex vertex, const bool verbose=false);

// replay operations for intrinsics recovery after initial simplification
void replay_intrinsic_flip(iSimpData& iSData, const int edge_idx);
void replay_vertex_flattening(iSimpData& iSData, const int vertex_idx);
void replay_vertex_flattening_with_heat(iSimpData& iSData, heatDiffData& hdData, const int vertex_idx);
void replay_vertex_removal(iSimpData& iSData, const int vertex_idx);
void replay_vertex_removal_with_heat(iSimpData& iSData, heatDiffData& hdData, const int vertex_idx);

std::array<int, 3> order_triangle_vertex_indices(const gcs::Face& face, const int F_first_v);
std::array<int, 3> order_triangle_edge_indices(const gcs::Face& face, const int F_first_v);

std::array<int, 5> order_quad_edge_indices(const gcs::Edge& edge, const int Fk_unique_v, const int Fl_unique_v);
std::array<int, 5> order_quad_edge_indices(const gcs::Edge& edge, const int Fk_unique_v, const int Fl_unique_v, const int vi_idx);
std::array<int, 4> order_quad_vertex_indices(const gcs::Edge& edge);
std::array<int, 4> order_quad_vertex_indices(const gcs::Edge& edge, const int Fk_unique_v, const int Fl_unique_v);
std::array<int, 4> order_quad_vertex_indices(const gcs::Edge& edge, const int Fk_unique_v, const int Fl_unique_v, const int vi_idx);
std::array<gcs::Face, 2> get_edge_faces(const gcs::Edge& edge);

// returns faceindex and value
std::pair<int, int> find_unique_vertex_index(const gcs::Face& F, const gcs::Edge& edge);
int find_vertex_face_index(const gcs::Face& F, const int vertex_idx);

double total_angle(const iSimpData& iSData, const gcs::Vertex& vertex);
double gaussian_curvature(const iSimpData& iSData, const gcs::Vertex& vertex);
double find_tangent_space_angle(const iSimpData& iSData, const gcs::Edge& edge, const gcs::Vertex& vertex);
PolarVector2D parallel_transport(const iSimpData& iSData, PolarVector2D tangent_vector, const gcs::Vertex origin, const gcs::Vertex target);

bool validate_intrinsic_edge_lengths(const iSimpData& iSData);

#endif