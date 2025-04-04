#include "isimp_utils.hpp"

void register_point(iSimpData &iSData, const BarycentricPoint point, const TexCoord texcoord, const int face_idx)
{
  int point_idx = iSData.tracked_points.size();
  iSData.tracked_by_triangle[face_idx].push_back(point_idx);
  iSData.tracked_points.push_back(point);
  iSData.tracked_texcoords.push_back(texcoord);
}

void map_registered(iSimpData &iSData)
{
  iSData.mapped_points = iSData.tracked_points;
  iSData.mapped_by_triangle = iSData.tracked_by_triangle;

  for (int i = 0; i < iSData.mapping.size(); i++)
  {
    iSData.mapping[i]->map_to_coarse(iSData.mapped_points, iSData.mapped_by_triangle);
  }
}

void map_registered(iSimpData &iSData, const int n_vertices)
{
  iSData.mapped_points = iSData.tracked_points;
  iSData.mapped_by_triangle = iSData.tracked_by_triangle;
  int remove_n_vertices = iSData.inputMesh->nVertices() - n_vertices;

  int removal_count = 0;
  for (int i = 0; i < iSData.mapping.size() && removal_count < remove_n_vertices; i++)
  {
    iSData.mapping[i]->map_to_coarse(iSData.mapped_points, iSData.mapped_by_triangle);
    if (iSData.mapping[i]->op_type == SIMP_OP::V_REMOVAL)
    {
      removal_count++;
    }
  }
}

// WATCH OUT! mapped_points & mapped_by_triangle are assumed to be initialized correctly.
// And the mapping index has to be chosen correctly as well.
// Use responsibly!
void map_current_from(iSimpData &iSData, const int from_mapping_idx)
{
  for (int i = from_mapping_idx; i < iSData.mapping.size(); i++)
  {
    iSData.mapping[i]->map_to_coarse(iSData.mapped_points, iSData.mapped_by_triangle);
  }
}

// WATCH OUT! mapped_points & mapped_by_triangle are assumed to be initialized correctly.
// And the mapping index has to be chosen correctly as well.
// Use responsibly!
void map_current_from_to(iSimpData &iSData, const int from_mapping_idx, const int to_mapping_idx)
{
  for (int i = from_mapping_idx; i < to_mapping_idx; i++)
  {
    iSData.mapping[i]->map_to_coarse(iSData.mapped_points, iSData.mapped_by_triangle);
  }
}

void unredistribute_heat(iSimpData &iSData, Eigen::VectorXd &heat)
{
  for (int i = iSData.heat_mapping.size() - 1; i >= 0; i--)
  {
    heatRedist hredist = iSData.heat_mapping[i];
    heat(hredist.vertex_idx) = 0.0;
    for (int i = 0; i < hredist.neighbor_idcs.rows(); i++)
    {
      int neighbor_idx = hredist.neighbor_idcs(i);
      int mass_delta = hredist.neighbor_mass_delta(i);
      double weight = hredist.neighbor_mass_delta(i) / (hredist.vertex_mass);
      heat(hredist.vertex_idx) += weight * heat(neighbor_idx);
    }
    for (int i = 0; i < hredist.neighbor_idcs.rows(); i++)
    {
      int neighbor_idx = hredist.neighbor_idcs(i);
      int mass_delta = hredist.neighbor_mass_delta(i);
      double weight = hredist.neighbor_mass_delta(i) / (hredist.neighbor_mass_delta(i) + hredist.neighbor_mass(i));
      heat(neighbor_idx) = (heat(neighbor_idx) - weight * heat(hredist.vertex_idx)) / (1 - weight);
    }
  }
}

double circumcentric_dual_mass(const iSimpData& iSData, const int vertex_idx)
{
    double total_area = 0.0;
    gcs::Vertex V = iSData.recoveredMesh->vertex(vertex_idx);
    double epsilon = 0.000001;

    for (gcs::Face F : V.adjacentFaces())
    {
        CircumcentricShape calculation_case = CircumcentricShape::ZERO;
        double area_contribution;
        std::array<int, 3> face_edges = order_triangle_edge_indices(F, vertex_idx);

        int e_ij = face_edges[0];
        int e_ik = face_edges[1];
        int e_jk = face_edges[2];
        double l_ij = iSData.recovered_L(e_ij);
        double l_ik = iSData.recovered_L(e_ik);
        double l_jk = iSData.recovered_L(e_jk);
        double midpoint_l_ij = l_ij/2.0;
        double midpoint_l_ik = l_ik/2.0;

        // Swap so the first edge is always the shorter one
        if (midpoint_l_ij > midpoint_l_ik)
        {
            double tmp_len = l_ij;
            l_ij = l_ik;
            l_ik = tmp_len;

            tmp_len = midpoint_l_ij;
            midpoint_l_ij = midpoint_l_ik;
            midpoint_l_ik = tmp_len;

            int tmp_e = e_ij;
            e_ij = e_ik;
            e_ik = tmp_e;
        }

        // The angle beta is only is reasonable, when considering a right triangle, i.e. in the TRIANGLE/CHOPPED TRIANGLE CASES
        double alpha = angle_i_from_lengths(l_ij, l_ik, l_jk, false, "circumcentric_dual_mass");
        double beta = M_PI/2.0 - alpha;

        // calc triangle side lengths
        double b = midpoint_l_ij;
        double a = tan(alpha)*b;
        double c = sqrt(a*a + b*b);

        // NOTE: If the angle is 0 or 180 degrees, the circumcentric area of p0 is 0.
        // NOTE: If the angle is greater than 0 degrees and less than 90 degrees, the circumcentric area of p0 can be calculated based on the area of the triangle between the midpoint1 and midpoint2 edges and the line containing midpoint1 with the direction of the normal of the midpoint1 edge.
        // NOTE: If the angle is 90 degrees, the circumcentric area of p0 can be calculated as the square of the length of the edge between p0 and midpoint1.
        // NOTE: If the angle is greater than 90 degrees and less than 180 degrees, the circumcentric area of p0 can be calculated based on the quadrilateral between the midpoint1 and midpoint2 edges and the two lines each containing either midpoint1 or midpoint2 with either the direction of the normal of the midpoint1 or midpoint2 edge.
        if (alpha < M_PI/2 - epsilon && c <= midpoint_l_ik)                     { calculation_case = CircumcentricShape::TRIANGLE; }
        if (alpha < M_PI/2 - epsilon && c > midpoint_l_ik)                      { calculation_case = CircumcentricShape::CHOPPED_TRIANGLE; }
        if (alpha >= M_PI/2 - epsilon && alpha <= M_PI/2 + epsilon)             { calculation_case = CircumcentricShape::RECTANGLE; }
        if (alpha > M_PI/2 + epsilon && alpha < M_PI - epsilon)                 { calculation_case = CircumcentricShape::CHOPPED_QUADRILATERAL; }

        switch(calculation_case)
        {
            default:
            case CircumcentricShape::ZERO:
                area_contribution = 0;
                break;
            case CircumcentricShape::TRIANGLE:
                area_contribution = a*b/2.0;
                break;
            case CircumcentricShape::CHOPPED_TRIANGLE:
            {
                // calc chopped off triangle side lengths
                double sa = c - midpoint_l_ik;
                double sb = tan(beta)*sa;
                area_contribution = (a*b - sa*sb)/2.0;
                break;
            }
            case CircumcentricShape::RECTANGLE:
                area_contribution = midpoint_l_ij*midpoint_l_ik;
                break;
            case CircumcentricShape::CHOPPED_QUADRILATERAL:
            {
                double triangle_area = 2.0*midpoint_l_ij*midpoint_l_ik*sin(alpha);

                double beta1 = angle_i_from_lengths(l_ij, l_jk, l_ik, false, "circumcentric_dual_mass");
                double beta2 = angle_i_from_lengths(l_ik, l_jk, l_ij, false, "circumcentric_dual_mass");
                double alpha1 = M_PI/2.0 - beta1;
                double alpha2 = M_PI/2.0 - beta2;
                double b1 = midpoint_l_ij / tan(alpha1);
                double b2 = midpoint_l_ik / tan(alpha2);
                double sta1 = abs(midpoint_l_ij*b1/2.0);
                double sta2 = abs(midpoint_l_ik*b2/2.0);

                area_contribution = triangle_area - sta1 - sta2;
                break;
            }
        }
        total_area += area_contribution;
    }

    return total_area;
}

// For an intrinsic edge to be (intrinsically) flippable two conditions have to hold:
//   1. the edge adjacent triangles have to form a convex quad
//   2. the edge endpoints have degree two or more
bool is_edge_flippable(const iSimpData &iSData, const Quad2D &quad, const gcs::Edge edge)
{
  gcs::Vertex v0 = edge.firstVertex();
  gcs::Vertex v1 = edge.secondVertex();

  // check if unfolding worked
  // (The unfolding can fail when (sometimes due to numerics) the edge lengths barely fail to satisfy the triangle inequality.
  // In such cases the angle calculation results in nan's which propagate to the unfolding.)
  if (quad.hasNaN())
  {
    return false;
  }

  // check vertex degrees
  if (v0.degree() < 2 || v1.degree() < 2)
  {
    return false;
  }
  // check if quad is convex
  if (!is_convex(quad))
  {
    return false;
  }
  return true;
}

double compute_flattening_factor(const Eigen::VectorXd &L, const gcs::Vertex v, const bool verbose)
{
  double THETA_HAT = (!v.isBoundary()) ? 2.0 * M_PI : M_PI;
  int vertex_idx = v.getIndex();
  int maximum_iterations = 5;
  int backtracking_iterations = 10;
  double threshold = 0.000001;
  double u_i = 0.0;

  // step 1000.0 for shock value. It is overridden in each iteration anyways. But it needs to be set to enter the loop.
  double step = 1000.0;
  double step_size = 1.0;
  int iteration = 0;
  while (std::abs(step_size * step) > threshold && iteration < maximum_iterations)
  {
    // Calculate newton step
    double enumerator = THETA_HAT;
    double denominator = 0.0;
    step_size = 1.0;
    for (gcs::Face F : v.adjacentFaces())
    {
      // get edge lengths
      std::array<int, 3> edge_indices = order_triangle_edge_indices(F, vertex_idx);
      double l_ij = std::exp(u_i / 2.0) * L[edge_indices[0]];
      double l_ik = std::exp(u_i / 2.0) * L[edge_indices[1]];
      double l_jk = L[edge_indices[2]];

      // calculate angles
      enumerator -= angle_i_from_lengths(l_ij, l_ik, l_jk, !verbose, "compute_flattening_factor");
      double theta_ij = angle_i_from_lengths(l_ik, l_jk, l_ij, !verbose, "compute_flattening_factor");
      double theta_ik = angle_i_from_lengths(l_ij, l_jk, l_ik, !verbose, "compute_flattening_factor");

      // calculate hessian
      if (theta_ij <= 0.0 || theta_ij >= M_PI)
      {
        denominator += 0.0;
      }
      else
      {
        denominator += 1.0 / std::tan(theta_ij);
      }
      if (theta_ik <= 0.0 || theta_ik >= M_PI)
      {
        denominator += 0.0;
      }
      else
      {
        denominator += 1.0 / std::tan(theta_ik);
      }
    }

    denominator /= 2.0;
    step = enumerator / denominator;

    // Backtracking
    bool found_valid_step = false;
    for (int i = 0; i < backtracking_iterations; i++)
    {
      // validate the result
      bool is_valid = true;
      for (gcs::Face F : v.adjacentFaces())
      {
        // get edge lengths
        std::array<int, 3> edge_indices = order_triangle_edge_indices(F, vertex_idx);
        double l_ij = std::exp((u_i - step_size * step) / 2.0) * L[edge_indices[0]];
        double l_ik = std::exp((u_i - step_size * step) / 2.0) * L[edge_indices[1]];
        double l_jk = L[edge_indices[2]];
        if (!satisfies_triangle_ineq(l_ij, l_ik, l_jk))
        {
          is_valid = false;
          break;
        }
      }
      if (!is_valid)
      {
        step_size = 0.5 * step_size;
      }
      else
      {
        found_valid_step = true;
        u_i = u_i - step_size * step;
        break;
      }
    }
    if (!found_valid_step)
    {
      return -1.0;
    }

    // safety check for NaN
    if (u_i != u_i)
    {
      if (verbose)
      {
        std::cout << dye("Encountered NaN!", RED) << std::endl;
      }
      return -1.0;
    }
    iteration++;
  }
  return u_i;
}

bool flip_intrinsic(iSimpData &iSData, const gcs::Edge edge, const bool verbose)
{
  std::array<gcs::Face, 2> faces = get_edge_faces(edge);
  gcs::Face Fk = faces[0];
  gcs::Face Fl = faces[1];

  std::pair<int, int> Fk_unique = find_unique_vertex_index(Fk, edge);
  std::pair<int, int> Fl_unique = find_unique_vertex_index(Fl, edge);

  // prevent flipping into self-edges
  if (Fk_unique.second == Fl_unique.second)
  {
    return false;
  }
  // get ordered intrinsic edge lengths from neighborhood
  std::array<int, 5> edge_indices = order_quad_edge_indices(edge, Fk_unique.second, Fl_unique.second);
  std::array<int, 4> vertex_indices = order_quad_vertex_indices(edge, Fk_unique.second, Fl_unique.second);
  double l_ij = iSData.L(edge_indices[0]);
  double l_ik = iSData.L(edge_indices[1]);
  double l_jk = iSData.L(edge_indices[2]);
  double l_il = iSData.L(edge_indices[3]);
  double l_jl = iSData.L(edge_indices[4]);
  Quad2D unfolded = unfold(l_ij, l_ik, l_jk, l_il, l_jl);

  // check if edge can be flipped
  if (!is_edge_flippable(iSData, unfolded, edge))
  {
    return false;
  }
  double l_flipped = flipped_edgelength(unfolded);
  // prevent flipping into degenerate triangles
  if (!satisfies_triangle_ineq(l_flipped, l_ik, l_il))
  {
    return false;
  }
  if (!satisfies_triangle_ineq(l_flipped, l_jk, l_jl))
  {
    return false;
  }
  // if (!is_non_degenerate(l_flipped, l_ik, l_il)) { return false; }
  // if (!is_non_degenerate(l_flipped, l_jk, l_jl)) { return false; }
  // update tangent space reference edges that would be flipped
  if (iSData.tangent_space_reference_edges[edge.firstVertex().getIndex()] == edge.getIndex())
  {
    rotate_reference_edge(iSData, edge.firstVertex());
  }
  if (iSData.tangent_space_reference_edges[edge.secondVertex().getIndex()] == edge.getIndex())
  {
    rotate_reference_edge(iSData, edge.secondVertex());
  }
  bool couldFlip = iSData.intrinsicMesh->flip(edge, false);
  if (!couldFlip)
  {
    return false;
  }

  // update edge length
  iSData.L(edge.getIndex()) = l_flipped;
  if (verbose)
  {
    bool result = validate_intrinsic_edge_lengths(iSData);
    if (!result)
    {
      std::cout << "occurred in flip_intrinsic for edge: " << edge.getIndex() << std::endl;
    }
  }
  // add intrinsic flip into mapping queue
  iSData.mapping.push_back(std::make_unique<Edge_Flip>(edge.getIndex(), Fk.getIndex(), Fl.getIndex(), unfolded, Fk_unique.first, Fl_unique.first));

  return true;
}

bool flatten_vertex(iSimpData &iSData, const int vertex_idx, const bool verbose)
{
  gcs::Vertex v = iSData.intrinsicMesh->vertex(vertex_idx);

  // TODO: if vertex is boundary vertex, check if it is incident on a boundary self-edge, or is a vertex of a self-face. If it is, it is skipped and flattening fails.
  std::unordered_map<int, double> L_opt = std::unordered_map<int, double>();
  for (gcs::Edge E : v.adjacentEdges())
  {
    int idx = E.getIndex();
    L_opt.insert({idx, iSData.L[idx]});
  }

  if (L_opt.size() == 2)
  {
    if (verbose)
    {
      std::cout << dye("Vertex to be flattened has only 2 incident faces!", RED) << std::endl;
    }
    return false;
  }
  for (gcs::Face F : v.adjacentFaces())
  {
    for (gcs::Face neighbor : F.adjacentFaces())
    {
      if (F.getIndex() == neighbor.getIndex())
      {
        if (verbose)
        {
          std::cout << dye("Found self-face during flattening!", RED) << std::endl;
        }
        return false;
      }
    }
  }

  double u_i = compute_flattening_factor(iSData.L, v, verbose);
  if (u_i == -1.0)
  {
    return false;
  }

  double eps = 0.000001;
  // if u_i is close to zero, the flattening is an identity operation (so we return a success, but don't insert a mapping operation into the mapping)
  // if(u_i < eps) { return true; } // NOTE: This supposed optimization lead to big errors in accuracy and NaN.

  // check if length's satisfy the triangle inequality
  for (gcs::Face F : v.adjacentFaces())
  {
    // get edge lengths
    std::array<int, 3> edge_indices = order_triangle_edge_indices(F, vertex_idx);
    double l_ij = std::exp(u_i / 2.0) * iSData.L[edge_indices[0]];
    double l_ik = std::exp(u_i / 2.0) * iSData.L[edge_indices[1]];
    double l_jk = iSData.L[edge_indices[2]];

    if (!satisfies_triangle_ineq(l_ij, l_ik, l_jk))
    {
      // TODO: perform edge flips to try to mitigate it
      if (verbose)
      {
        std::cout << dye("Violated triangle inequality", RED) << std::endl;
      }
      return false;
    }
  }

  // save previous curvatures for weighted mass redistribution
  std::unordered_map<int, double> redistribution_weights = std::unordered_map<int, double>();
  for (gcs::Vertex neighbor : v.adjacentVertices())
  {
    int neighbor_idx = neighbor.getIndex();
    redistribution_weights.insert({neighbor_idx, gaussian_curvature(iSData, neighbor)});
  }

  for (gcs::Edge E : v.adjacentEdges())
  {
    iSData.L[E.getIndex()] = std::exp(u_i / 2.0) * iSData.L[E.getIndex()];
  }
  if (verbose)
  {
    bool result = validate_intrinsic_edge_lengths(iSData);
    if (!result)
    {
      std::cout << "occurred in flatten_vertex" << std::endl;
    }
  }
  std::vector<int> faces = std::vector<int>();
  std::vector<int> face_vertex_indices = std::vector<int>();
  for (gcs::Face F : v.adjacentFaces())
  {
    faces.push_back(F.getIndex());
    face_vertex_indices.push_back(find_vertex_face_index(F, vertex_idx));
  }
  iSData.mapping.push_back(std::make_unique<Vertex_Flattening>(vertex_idx, faces, face_vertex_indices, u_i));

  // calculate weighting
  double total_curvature_change = 0.0;
  for (auto &pair : redistribution_weights)
  {
    gcs::Vertex neighbor = iSData.intrinsicMesh->vertex(pair.first);
    pair.second = std::abs(gaussian_curvature(iSData, neighbor) - pair.second);
    total_curvature_change += pair.second;
  }
  for (auto &pair : redistribution_weights)
  {
    // if already flat, distribute evenly
    if (total_curvature_change == 0.0)
    {
      pair.second = 1.0 / ((double)v.degree());
      continue;
    }
    pair.second = pair.second / total_curvature_change;
  }

  // redistribute mass & update error vectors
  for (gcs::Vertex neighbor : v.adjacentVertices())
  {
    int neighbor_idx = neighbor.getIndex();
    double alpha = redistribution_weights[neighbor_idx];
    iSData.T_minus.row(neighbor_idx) = next_error_vector(iSData, alpha, v, neighbor, true);
    iSData.T_minus.row(neighbor_idx) = next_error_vector(iSData, alpha, v, neighbor, false);

    iSData.M(neighbor_idx, 0) += alpha * iSData.M(vertex_idx, 0);
    iSData.M(neighbor_idx, 1) += alpha * iSData.M(vertex_idx, 1);
  }
  // clear flattened vertices values
  iSData.T_minus.row(vertex_idx) = PolarVector2D::Zero();
  iSData.T_plus.row(vertex_idx) = PolarVector2D::Zero();
  iSData.M(vertex_idx, 0) = 0.0;
  iSData.M(vertex_idx, 1) = 0.0;
  return true;
}

bool remove_vertex(iSimpData &iSData, const int vertex_idx, const bool verbose)
{
  // check if vertex is intrinsically flat (we can only remove it, when it is)
  gcs::Vertex vertex = iSData.intrinsicMesh->vertex(vertex_idx);
  double THETA_HAT = (!vertex.isBoundary()) ? 2.0 * M_PI : M_PI;
  if (vertex.isDead())
  {
    if (verbose)
    {
      std::cout << "Vertex has already been removed from mesh." << std::endl;
    }
    return false;
  }

  // NOTE: gcs::manifold_surface_mesh implementation does not allow boundary vertex removal
  if (vertex.isBoundary())
  {
    if (verbose)
    {
      std::cout << "Vertex is boundary vertex." << std::endl;
    }
    return false;
  }

  // check if vertex is degree 3
  if (vertex.degree() != 3)
  {
    return false;
  }

  int F_jk = -1;
  int F_jl = -1;
  int F_kl = -1;
  int F_res;
  int F_res_permutation;
  for (gcs::Face F : vertex.adjacentFaces())
  {

    if (F_jk == -1)
    {
      F_jk = F.getIndex();
      continue;
    }
    if (F_jl == -1)
    {
      F_jl = F.getIndex();
      continue;
    }
    if (F_kl == -1)
    {
      F_kl = F.getIndex();
      break;
    }
  }
  F_res = F_jk;

  gcs::Edge e_ij;
  for (gcs::Edge E1 : iSData.intrinsicMesh->face(F_jk).adjacentEdges())
  {
    for (gcs::Edge E2 : iSData.intrinsicMesh->face(F_jl).adjacentEdges())
    {
      if (E1.getIndex() == E2.getIndex())
      {
        e_ij = E1;
      }
    }
  }

  std::pair<int, int> Fk_unique = find_unique_vertex_index(iSData.intrinsicMesh->face(F_jk), e_ij);
  std::pair<int, int> Fl_unique = find_unique_vertex_index(iSData.intrinsicMesh->face(F_jl), e_ij);

  // get ordered intrinsic edge lengths from neighborhood
  std::array<int, 5> edge_indices = order_quad_edge_indices(e_ij, Fk_unique.second, Fl_unique.second, vertex_idx);
  double l_ij = iSData.L(edge_indices[0]);
  double l_ik = iSData.L(edge_indices[1]);
  double l_jk = iSData.L(edge_indices[2]);
  double l_il = iSData.L(edge_indices[3]);
  double l_jl = iSData.L(edge_indices[4]);
  Quad2D unfolded = unfold(l_ij, l_ik, l_jk, l_il, l_jl);

  // determining at which position in each face the shared vertex resides
  std::vector<int> shared_face_vertex_indices = std::vector<int>();
  int F_jk_shared = find_vertex_face_index(iSData.intrinsicMesh->face(F_jk), vertex_idx);
  shared_face_vertex_indices.push_back(F_jk_shared);

  int F_jl_shared = find_vertex_face_index(iSData.intrinsicMesh->face(F_jl), vertex_idx);
  shared_face_vertex_indices.push_back(F_jl_shared);
  int F_kl_shared = find_vertex_face_index(iSData.intrinsicMesh->face(F_kl), vertex_idx);
  shared_face_vertex_indices.push_back(F_kl_shared);

  // collect adjacent vertices
  std::vector<int> F_idxs = std::vector<int>();
  F_idxs.push_back(F_jk);
  F_idxs.push_back(F_jl);
  F_idxs.push_back(F_kl);
  std::array<int, 4> ordered_verts = order_quad_vertex_indices(e_ij, Fk_unique.second, Fl_unique.second, vertex_idx);

  // update reference edges for tangent spaces, if they are removed due to vertex removal
  for (gcs::Vertex V : iSData.intrinsicMesh->vertex(vertex_idx).adjacentVertices())
  {
    if (iSData.tangent_space_reference_edges[V.getIndex()] == iSData.intrinsicMesh->connectingEdge(vertex, V).getIndex())
    {
      rotate_reference_edge(iSData, V);
    }
  }

  // remove vertex from intrinsicMesh
  iSData.intrinsicMesh->removeVertex(vertex);
  F_res_permutation = (find_vertex_face_index(iSData.intrinsicMesh->face(F_res), vertex_idx) + 1 % 3);

  // insert vertex removal into mapping queue
  iSData.mapping.push_back(std::make_unique<Vertex_Removal>(vertex_idx, F_idxs, F_res, unfolded, shared_face_vertex_indices, F_res_permutation));
  return true;
}

// greedily flip edges, until delaunay
void flip_to_delaunay(iSimpData &iSData, const bool verbose)
{
  for (gcs::Edge E : iSData.intrinsicMesh->edges())
  {
    if (E.isBoundary())
    {
      continue;
    }
    std::array<gcs::Face, 2> faces = get_edge_faces(E);
    gcs::Face Fk = faces[0];
    gcs::Face Fl = faces[1];
    std::pair<int, int> Fk_unique = find_unique_vertex_index(Fk, E);
    std::pair<int, int> Fl_unique = find_unique_vertex_index(Fl, E);
    std::array<int, 5> edge_indices = order_quad_edge_indices(E, Fk_unique.second, Fl_unique.second);
    double l_ij = iSData.L(edge_indices[0]);
    double l_ik = iSData.L(edge_indices[1]);
    double l_jk = iSData.L(edge_indices[2]);
    double l_il = iSData.L(edge_indices[3]);
    double l_jl = iSData.L(edge_indices[4]);

    if (verbose)
    {
      if (!satisfies_triangle_ineq(l_ij, l_ik, l_jk))
      {
        std::cout << "flip_to_delaunay found triangle inequality violation!" << std::endl;
      }
      if (!satisfies_triangle_ineq(l_ij, l_il, l_jl))
      {
        std::cout << "flip_to_delaunay found triangle inequality violation!" << std::endl;
      }
    }
    // bool result = validate_intrinsic_edge_lengths(iSData);
    // if (!result) { std::cout << "occurred in flip_to_delaunay: " << std::endl; }

    double theta_k = angle_i_from_lengths(l_ik, l_jk, l_ij, !verbose, "flip_to_delaunay");
    double theta_l = angle_i_from_lengths(l_il, l_jl, l_ij, !verbose, "flip_to_delaunay");

    if (theta_k + theta_l > M_PI + 0.000001)
    {
      flip_intrinsic(iSData, E);
    }
  }
  if (verbose)
  {
    bool result = validate_intrinsic_edge_lengths(iSData);
    if (!result)
    {
      std::cout << "occurred in flip_to_delaunay: " << std::endl;
    }
  }
}

// greedily flip edges, until delaunay
// NOTE: Tested this, but yielded worse results than just going over all edges of the mesh
void flip_to_delaunay(iSimpData &iSData, const std::vector<gcs::Vertex> &vertices, const bool verbose)
{
  // gather edge indices a priori, so we don't get weird behavior due to concurrent modifications of the mesh during iteration.
  std::vector<int> edge_idcs = std::vector<int>();
  for (int i = 0; i < vertices.size(); i++)
  {
    gcs::Vertex V = vertices[i];
    for (gcs::Edge E : V.adjacentEdges())
    {
      edge_idcs.push_back(E.getIndex());
    }
  }

  for (int i = 0; i < edge_idcs.size(); i++)
  {
    gcs::Edge E = iSData.intrinsicMesh->edge(i);
    if (E.isBoundary())
    {
      continue;
    }
    std::array<gcs::Face, 2> faces = get_edge_faces(E);
    gcs::Face Fk = faces[0];
    gcs::Face Fl = faces[1];
    std::pair<int, int> Fk_unique = find_unique_vertex_index(Fk, E);
    std::pair<int, int> Fl_unique = find_unique_vertex_index(Fl, E);
    std::array<int, 5> edge_indices = order_quad_edge_indices(E, Fk_unique.second, Fl_unique.second);
    double l_ij = iSData.L(edge_indices[0]);
    double l_ik = iSData.L(edge_indices[1]);
    double l_jk = iSData.L(edge_indices[2]);
    double l_il = iSData.L(edge_indices[3]);
    double l_jl = iSData.L(edge_indices[4]);

    if (verbose)
    {
      if (!satisfies_triangle_ineq(l_ij, l_ik, l_jk))
      {
        std::cout << "flip_to_delaunay found triangle inequality violation!" << std::endl;
      }
      if (!satisfies_triangle_ineq(l_ij, l_il, l_jl))
      {
        std::cout << "flip_to_delaunay found triangle inequality violation!" << std::endl;
      }
    }

    double theta_k = angle_i_from_lengths(l_ik, l_jk, l_ij, !verbose, "flip_to_delaunay");
    double theta_l = angle_i_from_lengths(l_il, l_jl, l_ij, !verbose, "flip_to_delaunay");

    if (theta_k + theta_l > M_PI + 0.000001)
    {
      flip_intrinsic(iSData, E);
    }
  }
}

bool flip_vertex_to_deg3(iSimpData &iSData, const int vertex_idx, const bool verbose)
{
  // stack of performed edge flips in case we need to revert them
  std::vector<int> flips = std::vector<int>();

  gcs::Vertex V = iSData.intrinsicMesh->vertex(vertex_idx);
  if (V.degree() == 3)
  {
    return true;
  }

  bool found_flip = true;
  while (found_flip)
  {
    found_flip = false;
    for (gcs::Edge E : V.adjacentEdges())
    {
      bool could_flip = flip_intrinsic(iSData, E);
      if (could_flip)
      {
        // std::cout << "flipped edge: " << E << "(" << E.firstVertex().getIndex() << ", " << E.secondVertex().getIndex() << ")" << std::endl;
        flips.push_back(E.getIndex());
        found_flip = true;
        break;
      }
      // else {
      //   std::cout << "could not flip edge: " << E << "(" << E.firstVertex().getIndex() << ", " << E.secondVertex().getIndex() << ")" << std::endl;
      // }
    }
    // found valid edge flip sequence
    if (V.degree() == 3)
    {
      return true;
    }
  }

  // TODO: the edge flip reversion currently leads to some issues
  // It seems there is a case that (possibly due to numerics) can arise, where a first edge flip can be performed, but not the three edge flips to undo it.
  // Because in my reasoning a quad that can be flipped once, should be able to be flipped any number of times.
  // revert edge-flips
  // while(flips.size() > 0)
  // {
  //   gcs::Edge E = iSData.intrinsicMesh->edge(flips.back());
  //   flips.pop_back();

  //   // NOTE: we need to flip three times to undo one edge flip
  //   flip_intrinsic(iSData, E);
  //   flip_intrinsic(iSData, E);
  //   flip_intrinsic(iSData, E);

  //   // NOTE: we remove the total four edge flip operations on edge E from the mapping
  //   iSData.mapping.pop_back();
  //   iSData.mapping.pop_back();
  //   iSData.mapping.pop_back();
  //   iSData.mapping.pop_back();
  // }
  return false;
}

void rotate_reference_edge(iSimpData &iSData, const gcs::Vertex &vertex)
{
  if (vertex.degree() < 2)
  {
    return;
  }
  int vertex_idx = vertex.getIndex();
  int ref_edge_idx = iSData.tangent_space_reference_edges[vertex_idx];

  gcs::Edge next_ref_edge;
  for (gcs::Edge edge : vertex.adjacentEdges())
  {
    if (edge.getIndex() != ref_edge_idx)
    {
      next_ref_edge = edge;
      break;
    }
  }

  double next_ref_edge_angle = find_tangent_space_angle(iSData, next_ref_edge, vertex);

  // update tangent space values
  iSData.T_minus(vertex_idx, 0) -= next_ref_edge_angle;
  iSData.T_plus(vertex_idx, 0) -= next_ref_edge_angle;
  if (iSData.T_minus(vertex_idx, 0) < 0.0)
  {
    iSData.T_minus(vertex_idx, 0) += 1.0;
  }
  if (iSData.T_plus(vertex_idx, 0) < 0.0)
  {
    iSData.T_plus(vertex_idx, 0) += 1.0;
  }
  iSData.tangent_space_reference_edges[vertex_idx] = next_ref_edge.getIndex();
}

PolarVector2D next_error_vector(const iSimpData &iSData, const double alpha, const gcs::Vertex vertex, const gcs::Vertex neighbor, bool for_T_minus)
{
  int vertex_idx = vertex.getIndex();
  int neighbor_idx = neighbor.getIndex();
  // std::cout << "next_error_vector from vertex " << vertex_idx << " to vertex " << neighbor_idx << std::endl;
  PolarVector2D T;
  Vector2D nT;
  double mi;
  double mj;
  if (for_T_minus)
  {
    T = iSData.T_minus.row(vertex_idx);
    nT = to_cartesian(iSData.T_minus.row(neighbor_idx));
    mi = alpha * iSData.M(vertex_idx, 0);
    mj = iSData.M(neighbor_idx, 0);
  }
  else
  {
    T = iSData.T_plus.row(vertex_idx);
    nT = to_cartesian(iSData.T_plus.row(neighbor_idx));
    mi = alpha * iSData.M(vertex_idx, 1);
    mj = iSData.M(neighbor_idx, 1);
  }

  gcs::Edge connection = iSData.intrinsicMesh->connectingEdge(vertex, neighbor);
  if (connection.isBoundary())
  {
    std::cout << RED << "Chose dead connection!" << RESET << std::endl;
  }
  if (connection.firstVertex().isDead())
  {
    std::cout << RED << "Chose dead connection!" << RESET << std::endl;
    return PolarVector2D(0.0, 0.0);
  }
  if (connection.secondVertex().isDead())
  {
    std::cout << RED << "Chose dead connection!" << RESET << std::endl;
    return PolarVector2D(0.0, 0.0);
  }
  double edge_angle = find_tangent_space_angle(iSData, connection, neighbor);
  double edge_length = iSData.L[connection.getIndex()];
  Vector2D transport = to_cartesian({edge_angle, edge_length});
  Vector2D transported = to_cartesian(parallel_transport(iSData, T, vertex, neighbor));

  if (mi == 0.0 && mj == 0.0)
  {
    return PolarVector2D(0.0, 0.0);
  }
  PolarVector2D next_T = to_polar((alpha * mi * (transported + transport) + mj * nT) / (alpha * mi + mj));
  return next_T;
}

double intrinsic_curvature_error(const iSimpData &iSData, const gcs::Vertex vertex, const bool verbose)
{
  if (vertex.isBoundary())
  {
    return INFINITY;
  } // guard, because we chose to forgo boundary vertex simplification
  double ice = 0.0;
  int vertex_idx = vertex.getIndex();

  // TODO: if vertex is boundary vertex, check if it is incident on a boundary self-edge, or is a vertex of a self-face. If it is, it is skipped and flattening fails.
  std::unordered_map<int, double> L_opt = std::unordered_map<int, double>();
  for (gcs::Edge E : vertex.adjacentEdges())
  {
    int idx = E.getIndex();
    L_opt.insert({idx, iSData.L[idx]});
  }

  if (L_opt.size() == 2)
  {
    if (verbose)
    {
      std::cout << dye("Found 2-face vertex during flattening!", RED) << std::endl;
    }
    return INFINITY;
  }
  for (gcs::Face F : vertex.adjacentFaces())
  {
    for (gcs::Face neighbor : F.adjacentFaces())
    {
      if (F.getIndex() == neighbor.getIndex())
      {
        if (verbose)
        {
          std::cout << dye("Found self-face during flattening!", RED) << std::endl;
        }
        return INFINITY;
      }
    }
  }

  double u_i = compute_flattening_factor(iSData.L, vertex, false);
  if (u_i == -1)
  {
    return INFINITY;
  }
  double eps = 0.0000001;

  // check if length's satisfy the triangle inequality
  for (gcs::Face F : vertex.adjacentFaces())
  {
    // get edge lengths
    std::array<int, 3> edge_indices = order_triangle_edge_indices(F, vertex_idx);
    double l_ij = std::exp(u_i / 2.0) * iSData.L[edge_indices[0]];
    double l_ik = std::exp(u_i / 2.0) * iSData.L[edge_indices[1]];
    double l_jk = iSData.L[edge_indices[2]];

    if (!satisfies_triangle_ineq(l_ij, l_ik, l_jk))
    {
      // TODO: perform edge flips to try to mitigate it
      if (verbose)
      {
        std::cout << dye("Violated triangle inequality", RED) << std::endl;
      }
      return INFINITY;
    }
  }

  // save previous curvatures for weighted mass redistribution
  std::unordered_map<int, double> redistribution_weights = std::unordered_map<int, double>();
  for (gcs::Vertex neighbor : vertex.adjacentVertices())
  {
    int neighbor_idx = neighbor.getIndex();
    redistribution_weights.insert({neighbor_idx, gaussian_curvature(iSData, neighbor)});
  }

  // calculate weighting
  double total_curvature_change = 0.0;
  for (auto &pair : redistribution_weights)
  {
    gcs::Vertex neighbor = iSData.intrinsicMesh->vertex(pair.first);
    pair.second = std::abs(gaussian_curvature(iSData, neighbor) - pair.second);
    total_curvature_change += pair.second;
  }
  for (auto &pair : redistribution_weights)
  {
    // if already flat, distribute evenly
    if (total_curvature_change == 0.0)
    {
      pair.second = 1.0 / ((double)vertex.degree());
      continue;
    }
    pair.second = pair.second / total_curvature_change;
  }

  // redistribute mass & update error vectors
  PolarVector2D T_minus;
  PolarVector2D T_plus;
  for (gcs::Vertex neighbor : vertex.adjacentVertices())
  {
    int neighbor_idx = neighbor.getIndex();
    double alpha = redistribution_weights[neighbor_idx];
    T_minus = next_error_vector(iSData, alpha, vertex, neighbor, true);
    T_plus = next_error_vector(iSData, alpha, vertex, neighbor, false);

    double m_minus = iSData.M(neighbor_idx, 0) + alpha * iSData.M(vertex_idx, 0);
    double m_plus = iSData.M(neighbor_idx, 1) + alpha * iSData.M(vertex_idx, 1);
    ice += m_minus * T_minus[1];
    ice += m_plus * T_plus[1];
  }
  return ice;
}

// orders triangle vertex indices of triangle ijk by a vertex i
// index 0: vertex i
// index 1: vertex j of ijk-triangle
// index 2: vertex k of ijk-triangle
std::array<int, 3> order_triangle_vertex_indices(const gcs::Face &face, const int F_first_v)
{
  std::array<int, 3> indices = std::array<int, 3>();

  int first_idx = find_vertex_face_index(face, F_first_v);

  int i = 0;
  for (gcs::Vertex V : face.adjacentVertices())
  {
    indices[(i + (3 - first_idx)) % 3] = V.getIndex();
    i++;
  }

  return indices;
}

void replay_intrinsic_flip(iSimpData &iSData, const int edge_idx)
{
  gcs::Edge edge = iSData.recoveredMesh->edge(edge_idx);
  std::array<gcs::Face, 2> faces = get_edge_faces(edge);
  gcs::Face Fk = faces[0];
  gcs::Face Fl = faces[1];

  std::pair<int, int> Fk_unique = find_unique_vertex_index(Fk, edge);
  std::pair<int, int> Fl_unique = find_unique_vertex_index(Fl, edge);
  // prevent flipping into self-edges
  if (Fk_unique.second == Fl_unique.second)
  {
    return exit(-1);
  }
  // get ordered intrinsic edge lengths from neighborhood
  std::array<int, 5> edge_indices = order_quad_edge_indices(edge, Fk_unique.second, Fl_unique.second);
  std::array<int, 4> vertex_indices = order_quad_vertex_indices(edge, Fk_unique.second, Fl_unique.second);
  double l_ij = iSData.recovered_L(edge_indices[0]);
  double l_ik = iSData.recovered_L(edge_indices[1]);
  double l_jk = iSData.recovered_L(edge_indices[2]);
  double l_il = iSData.recovered_L(edge_indices[3]);
  double l_jl = iSData.recovered_L(edge_indices[4]);
  Quad2D unfolded = unfold(l_ij, l_ik, l_jk, l_il, l_jl);

  if (!is_edge_flippable(iSData, unfolded, edge))
  {
    std::cout << RED << "intrinsic edge-flip replay failed!!! (edge is not flippable)" << RESET << std::endl;
    exit(-1);
  }
  double l_flipped = flipped_edgelength(unfolded);
  // prevent flipping into degenerate triangles
  if (!satisfies_triangle_ineq(l_flipped, l_ik, l_il))
  {
    std::cout << RED << "intrinsic edge-flip triangle inequality violation in replay!!!" << std::endl;
    exit(-1);
  }
  if (!satisfies_triangle_ineq(l_flipped, l_jk, l_jl))
  {
    std::cout << RED << "intrinsic edge-flip triangle inequality violation in replay!!!" << std::endl;
    exit(-1);
  }
  bool couldFlip = iSData.recoveredMesh->flip(edge, false);
  if (!couldFlip)
  {
    std::cout << RED << "intrinsic edge-flip replay failed!!! (half-edge mesh edge flip failed)" << RESET << std::endl;
    exit(-1);
  }

  // update edge length
  iSData.recovered_L(edge.getIndex()) = l_flipped;
  return;
}

void replay_vertex_flattening(iSimpData &iSData, const int vertex_idx)
{
  gcs::Vertex v = iSData.recoveredMesh->vertex(vertex_idx);

  double u_i = compute_flattening_factor(iSData.recovered_L, v, false);
  if (u_i == -1.0)
  {
    exit(-1);
  }

  double eps = 0.000001;

  // check if length's satisfy the triangle inequality
  for (gcs::Face F : v.adjacentFaces())
  {
    // get edge lengths
    std::array<int, 3> edge_indices = order_triangle_edge_indices(F, vertex_idx);
    double l_ij = std::exp(u_i / 2.0) * iSData.recovered_L[edge_indices[0]];
    double l_ik = std::exp(u_i / 2.0) * iSData.recovered_L[edge_indices[1]];
    double l_jk = iSData.recovered_L[edge_indices[2]];

    if (!satisfies_triangle_ineq(l_ij, l_ik, l_jk))
    {
      // TODO: perform edge flips to try to mitigate it
      if (false)
      {
        std::cout << dye("Violated triangle inequality", RED) << std::endl;
      }
      exit(-1);
    }
  }

  for (gcs::Edge E : v.adjacentEdges())
  {
    iSData.recovered_L[E.getIndex()] = std::exp(u_i / 2.0) * iSData.L[E.getIndex()];
  }
  return;
}

void replay_vertex_flattening_with_heat(iSimpData &iSData, heatDiffData &hdData, const int vertex_idx)
{
  gcs::Vertex v = iSData.recoveredMesh->vertex(vertex_idx);
  double u_i = compute_flattening_factor(iSData.recovered_L, v, false);
  if (u_i == -1.0)
  {
    exit(-1);
  }

  double eps = 0.000001;

  // check if length's satisfy the triangle inequality
  for (gcs::Face F : v.adjacentFaces())
  {
    // get edge lengths
    std::array<int, 3> edge_indices = order_triangle_edge_indices(F, vertex_idx);
    double l_ij = std::exp(u_i / 2.0) * iSData.recovered_L[edge_indices[0]];
    double l_ik = std::exp(u_i / 2.0) * iSData.recovered_L[edge_indices[1]];
    double l_jk = iSData.recovered_L[edge_indices[2]];

    if (!satisfies_triangle_ineq(l_ij, l_ik, l_jk))
    {
      // TODO: perform edge flips to try to mitigate it
      if (false)
      {
        std::cout << dye("Violated triangle inequality", RED) << std::endl;
      }
      exit(-1);
    }
  }

  for (gcs::Edge E : v.adjacentEdges())
  {
    iSData.recovered_L[E.getIndex()] = std::exp(u_i / 2.0) * iSData.recovered_L[E.getIndex()];
  }
  return;
}

void replay_vertex_removal(iSimpData &iSData, const int vertex_idx)
{
  // check if vertex is intrinsically flat (we can only remove it, when it is)
  gcs::Vertex vertex = iSData.recoveredMesh->vertex(vertex_idx);
  double THETA_HAT = (!vertex.isBoundary()) ? 2.0 * M_PI : M_PI;
  if (vertex.isDead())
  {
    std::cout << RED << "replay_vertex_removal: Vertex has already been removed from mesh." << RESET << std::endl;
    exit(-1);
  }

  // NOTE: gcs::manifold_surface_mesh implementation does not allow boundary vertex removal
  if (vertex.isBoundary())
  {
    std::cout << RED << "replay_vertex_removal: Vertex is boundary vertex." << RESET << std::endl;
    exit(-1);
  }
  // check if vertex is degree 3
  if (vertex.degree() != 3)
  {
    exit(-1);
  }

  // remove vertex from intrinsicMesh
  iSData.recoveredMesh->removeVertex(vertex);
  return;
}

void replay_vertex_removal_with_heat(iSimpData &iSData, heatDiffData &hdData, const int vertex_idx)
{
  // check if vertex is intrinsically flat (we can only remove it, when it is)
  gcs::Vertex vertex = iSData.recoveredMesh->vertex(vertex_idx);
  double THETA_HAT = (!vertex.isBoundary()) ? 2.0 * M_PI : M_PI;
  if (vertex.isDead())
  {
    std::cout << RED << "replay_vertex_removal: Vertex has already been removed from mesh." << RESET << std::endl;
    exit(-1);
  }

  // NOTE: gcs::manifold_surface_mesh implementation does not allow boundary vertex removal
  if (vertex.isBoundary())
  {
    std::cout << RED << "replay_vertex_removal: Vertex is boundary vertex." << RESET << std::endl;
    exit(-1);
  }
  // check if vertex is degree 3
  if (vertex.degree() != 3)
  {
    exit(-1);
  }

  // heat redistribution
  int d = vertex.degree();
  double epsilon = 1e-3;
  heatRedist redist;
  redist.vertex_idx = vertex_idx;
  redist.neighbor_idcs = Eigen::VectorXi::Zero(d);
  redist.neighbor_mass = Eigen::VectorXd::Zero(d);
  redist.neighbor_mass_delta = Eigen::VectorXd::Zero(d);
  double heat = hdData.recovered_heat(vertex.getIndex(), 0);

  // Compute pre-removal circumcentric masses (dual area)
  int i = 0;
  double mass = circumcentric_dual_mass(iSData, vertex_idx);
  redist.vertex_mass = mass;
  for (gcs::Vertex V : vertex.adjacentVertices())
  {
    redist.neighbor_idcs(i) = V.getIndex();
    redist.neighbor_mass(i) = circumcentric_dual_mass(iSData, V.getIndex());
    i++;
  }

  // remove vertex from intrinsicMesh
  iSData.recoveredMesh->removeVertex(vertex);

  // Compute change in circumcentric masses post removal
  double total_change = 0.0;
  for (i=0; i<redist.neighbor_idcs.rows(); i++)
  {
    double mass_delta = circumcentric_dual_mass(iSData, redist.neighbor_idcs(i)) - redist.neighbor_mass(i);
    mass_delta = std::max(std::min(mass_delta, 1.0), 0.0);
    if (mass_delta >= -epsilon && mass_delta < 0) { mass_delta = 0.0; }
    redist.neighbor_mass_delta(i) = mass_delta;
    total_change += mass_delta;
    if (mass_delta < -epsilon)
    {
      std::cout << "encountered negative mass delta: " << mass_delta << " !!!!" << std::endl;
      std::cout << "previous mass: " << redist.neighbor_mass(i) << std::endl;
      std::cout << "post removal mass: " << (mass_delta + redist.neighbor_mass(i)) << std::endl;
      exit(-1);
    }
  }

  for (i=0; i<redist.neighbor_idcs.rows(); i++)
  {
    int neighbor_idx = redist.neighbor_idcs(i);
    double neighbor_heat = hdData.recovered_heat(neighbor_idx, 0);
    double weight = redist.neighbor_mass_delta(i) / (redist.neighbor_mass_delta(i) + redist.neighbor_mass(i));
    if (weight < 0.0)
    {
      std::cout << "encountered negative weight!!!!" << std::endl;
      exit(-1);
    }
    if (weight > 1.0)
    {
      std::cout << "encountered excessive weight > 1.0 !!!!" << std::endl;
      exit(-1);
    }
    hdData.recovered_heat(neighbor_idx, 0) = heat * weight + neighbor_heat * (1 - weight);
  }
  iSData.heat_mapping.push_back(redist);
  hdData.recovered_heat(vertex.getIndex(), 0) = 0.0;
  return;
}

// orders triangle edge indices of triangle ijk by a vertex i
// index 0: shared edge ij
// index 1: first edge ik of ijk-triangle
// index 2: second edge jk of ijk-triangle
std::array<int, 3> order_triangle_edge_indices(const gcs::Face &face, const int F_first_v)
{
  std::array<int, 3> indices = std::array<int, 3>();
  std::array<int, 3> v_ordered = order_triangle_vertex_indices(face, F_first_v);

  for (gcs::Edge E : face.adjacentEdges())
  {
    int e0 = E.firstVertex().getIndex();
    int e1 = E.secondVertex().getIndex();
    if (e0 == v_ordered[0] && e1 == v_ordered[1] || e0 == v_ordered[1] && e1 == v_ordered[0])
    {
      indices[0] = E.getIndex();
    }
    else if (e0 == v_ordered[0] && e1 == v_ordered[2] || e0 == v_ordered[2] && e1 == v_ordered[0])
    {
      indices[1] = E.getIndex();
    }
    else if (e0 == v_ordered[1] && e1 == v_ordered[2] || e0 == v_ordered[2] && e1 == v_ordered[1])
    {
      indices[2] = E.getIndex();
    }
  }

  return indices;
}

// index 0: shared edge ij
// index 1: first edge ik of ijk-triangle
// index 2: second edge jk of ijk-triangle
// index 3: first edge il of ijl-triangle
// index 4: second edge jl of ijl-triangle
// Note: this version assumes the first vertex of edge is i
std::array<int, 5> order_quad_edge_indices(const gcs::Edge &edge, const int Fk_unique_v, const int Fl_unique_v)
{
  std::array<int, 5> indices = std::array<int, 5>();
  indices[0] = edge.getIndex();

  int offset = 1;
  for (gcs::Face F : edge.adjacentFaces())
  {
    for (gcs::Edge E : F.adjacentEdges())
    {
      if ((E.firstVertex().getIndex() == edge.firstVertex().getIndex() && E.secondVertex().getIndex() == Fk_unique_v) || (E.secondVertex().getIndex() == edge.firstVertex().getIndex() && E.firstVertex().getIndex() == Fk_unique_v))
      {
        indices[1] = E.getIndex();
      }
      else if ((E.firstVertex().getIndex() == edge.secondVertex().getIndex() && E.secondVertex().getIndex() == Fk_unique_v) || (E.secondVertex().getIndex() == edge.secondVertex().getIndex() && E.firstVertex().getIndex() == Fk_unique_v))
      {
        indices[2] = E.getIndex();
      }
      else if ((E.firstVertex().getIndex() == edge.firstVertex().getIndex() && E.secondVertex().getIndex() == Fl_unique_v) || (E.secondVertex().getIndex() == edge.firstVertex().getIndex() && E.firstVertex().getIndex() == Fl_unique_v))
      {
        indices[3] = E.getIndex();
      }
      else if ((E.firstVertex().getIndex() == edge.secondVertex().getIndex() && E.secondVertex().getIndex() == Fl_unique_v) || (E.secondVertex().getIndex() == edge.secondVertex().getIndex() && E.firstVertex().getIndex() == Fl_unique_v))
      {
        indices[4] = E.getIndex();
      }
    }
  }
  return indices;
}

// index 0: shared edge ij
// index 1: first edge ik of ijk-triangle
// index 2: second edge jk of ijk-triangle
// index 3: first edge il of ijl-triangle
// index 4: second edge jl of ijl-triangle
std::array<int, 5> order_quad_edge_indices(const gcs::Edge &edge, const int Fk_unique_v, const int Fl_unique_v, const int vi_idx)
{
  std::array<int, 5> indices = std::array<int, 5>();
  indices[0] = edge.getIndex();

  int offset = 1;
  for (gcs::Face F : edge.adjacentFaces())
  {
    for (gcs::Edge E : F.adjacentEdges())
    {
      if ((E.firstVertex().getIndex() == vi_idx && E.secondVertex().getIndex() == Fk_unique_v) || (E.secondVertex().getIndex() == vi_idx && E.firstVertex().getIndex() == Fk_unique_v))
      {
        indices[1] = E.getIndex();
      }
      else if ((E.firstVertex().getIndex() != vi_idx && E.secondVertex().getIndex() == Fk_unique_v) || (E.secondVertex().getIndex() != vi_idx && E.firstVertex().getIndex() == Fk_unique_v))
      {
        indices[2] = E.getIndex();
      }
      else if ((E.firstVertex().getIndex() == vi_idx && E.secondVertex().getIndex() == Fl_unique_v) || (E.secondVertex().getIndex() == vi_idx && E.firstVertex().getIndex() == Fl_unique_v))
      {
        indices[3] = E.getIndex();
      }
      else if ((E.firstVertex().getIndex() != vi_idx && E.secondVertex().getIndex() == Fl_unique_v) || (E.secondVertex().getIndex() != vi_idx && E.firstVertex().getIndex() == Fl_unique_v))
      {
        indices[4] = E.getIndex();
      }
    }
  }
  return indices;
}

// index 0: left shared vertex (when orienting the quad such that the first triangle is ontop)
// index 1: right shared vertex
// index 2: first triangles unique vertex
// index 3: second triangles unique vertex
std::array<int, 4> order_quad_vertex_indices(const gcs::Edge &edge)
{
  // because in the ManifoldSurfaceMesh data structure the half-edges form counterclockwise cycles, this should work
  std::array<int, 4> indices = std::array<int, 4>();
  indices[0] = edge.firstVertex().getIndex();
  indices[1] = edge.secondVertex().getIndex();

  int offset = 2;
  for (gcs::Face F : edge.adjacentFaces())
  {
    for (gcs::Vertex V : F.adjacentVertices())
    {
      if (V.getIndex() != edge.firstVertex().getIndex() && V.getIndex() != edge.secondVertex().getIndex())
      {
        indices[offset] = V.getIndex();
        break;
      }
    }
    offset++;
  }
  return indices;
}

std::array<int, 4> order_quad_vertex_indices(const gcs::Edge &edge, const int Fk_unique_v, const int Fl_unique_v)
{
  // because in the ManifoldSurfaceMesh data structure the half-edges form counterclockwise cycles, this should work
  std::array<int, 4> indices = std::array<int, 4>();
  indices[0] = edge.firstVertex().getIndex();
  indices[1] = edge.secondVertex().getIndex();
  indices[2] = Fk_unique_v;
  indices[3] = Fl_unique_v;
  return indices;
}

std::array<int, 4> order_quad_vertex_indices(const gcs::Edge &edge, const int Fk_unique_v, const int Fl_unique_v, const int vi_idx)
{
  // because in the ManifoldSurfaceMesh data structure the half-edges form counterclockwise cycles, this should work
  std::array<int, 4> indices = std::array<int, 4>();
  indices[0] = vi_idx;
  if (edge.secondVertex().getIndex() != vi_idx)
  {
    indices[1] = edge.secondVertex().getIndex();
  }
  else
  {
    indices[1] = edge.firstVertex().getIndex();
  }
  indices[2] = Fk_unique_v;
  indices[3] = Fl_unique_v;
  return indices;
}

std::array<gcs::Face, 2> get_edge_faces(const gcs::Edge &edge)
{
  std::array<gcs::Face, 2> faces = std::array<gcs::Face, 2>();
  bool t = true;
  for (gcs::Face F : edge.adjacentFaces())
  {
    if (t)
    {
      faces[0] = F;
      t = false;
    }
    else
    {
      faces[1] = F;
      break;
    }
  }
  return faces;
}

std::pair<int, int> find_unique_vertex_index(const gcs::Face &F, const gcs::Edge &edge)
{
  int i = 0;
  for (gcs::Vertex V : F.adjacentVertices())
  {
    if (V.getIndex() != edge.firstVertex().getIndex() && V.getIndex() != edge.secondVertex().getIndex())
    {
      return std::pair<int, int>(i, V.getIndex());
    }
    i++;
  }
  return std::pair<int, int>(-1, -1);
}

int find_vertex_face_index(const gcs::Face &F, const int vertex)
{
  int i = 0;
  for (gcs::Vertex V : F.adjacentVertices())
  {
    if (V.getIndex() == vertex)
    {
      return i;
    }
    i++;
  }
  return -1;
}

double total_angle(const iSimpData &iSData, const gcs::Vertex &vertex)
{
  double total_angle = 0.0;
  int vertex_idx = vertex.getIndex();

  for (gcs::Face F : vertex.adjacentFaces())
  {
    // get edge lengths
    std::array<int, 3> edge_indices = order_triangle_edge_indices(F, vertex_idx);
    double l_ij = iSData.L[edge_indices[0]];
    double l_ik = iSData.L[edge_indices[1]];
    double l_jk = iSData.L[edge_indices[2]];
    if (l_jk < 0.0000000001)
    {
      std::cout << "Calling angle_i from total_angle with zero length" << std::endl;
    }

    total_angle += angle_i_from_lengths(l_ij, l_ik, l_jk, false, "total_angle");
  }
  // bool result = validate_intrinsic_edge_lengths(iSData);
  // if (!result) { std::cout << "occurred in total_angle: " << std::endl; }

  return total_angle;
}

// NOTE: Maybe the function name is a misnomer, bc it also calculates geodesic curvature for boundary vertices
double gaussian_curvature(const iSimpData &iSData, const gcs::Vertex &vertex)
{
  double THETA_HAT = (!vertex.isBoundary()) ? 2.0 * M_PI : M_PI;
  double curvature = THETA_HAT;
  int vertex_idx = vertex.getIndex();

  for (gcs::Face F : vertex.adjacentFaces())
  {
    // get edge lengths
    std::array<int, 3> edge_indices = order_triangle_edge_indices(F, vertex_idx);
    double l_ij = iSData.L[edge_indices[0]];
    double l_ik = iSData.L[edge_indices[1]];
    double l_jk = iSData.L[edge_indices[2]];
    // std::cout << "Calcing Curvature" << std::endl;

    if (l_jk < 0.0000000001)
    {
      std::cout << "gaussian_curvature: Calling angle_i with zero length" << std::endl;
    }
    curvature -= angle_i_from_lengths(l_ij, l_ik, l_jk, false, "gaussian_curvature");
  }
  // bool result = validate_intrinsic_edge_lengths(iSData);
  // if (!result) { std::cout << "occurred in gaussian_curvature: " << std::endl; }

  return curvature;
}

double find_tangent_space_angle(const iSimpData &iSData, const gcs::Edge &edge, const gcs::Vertex &vertex)
{
  // std::cout << "find_tangent_space_angle:" << std::endl;
  double angle = 0.0;
  // std::cout << "finding tangent space angle of edge: " << edge.getIndex() << " around vertex " << vertex.getIndex() << " with reference edge: " << iSData.intrinsicMesh->edge(iSData.tangent_space_reference_edges[vertex.getIndex()]).getIndex() << std::endl;

  // get reference edge
  int vertex_idx = vertex.getIndex();
  int ref_edge_idx = iSData.tangent_space_reference_edges[vertex_idx];
  gcs::Edge ref_edge = iSData.intrinsicMesh->edge(ref_edge_idx);
  if (edge.getIndex() == ref_edge_idx)
  {
    return angle;
  }
  std::array<int, 3> indices;

  bool found_sth = false;
  bool found_ref_edge_first = true;
  // iterate clockwise through the faces around vertex
  // NOTE: the edges around faces are ordered counter clockwise. Thus, we need to use indices[1] to get the "lagging behind" edge.
  for (gcs::Face F : vertex.adjacentFaces())
  {
    // printGCSFace(F);
    indices = order_triangle_edge_indices(F, vertex_idx);
    // std::cout << " determined indices: { " << indices[0] << ", " << indices[1] << ", " << indices[2] << "}" << std::endl;
    // for (int i = 0; i<3; i++)
    // {
    //   std::cout << "edge" << indices[i] << ": " << iSData.intrinsicMesh->edge(indices[i]).firstVertex() << " -> " << iSData.intrinsicMesh->edge(indices[i]).secondVertex() << std::endl;
    // }

    if (!found_sth)
    {
      // check for ref_edge
      if (indices[1] == ref_edge_idx)
      {
        found_ref_edge_first = true;
        found_sth = true;
        // std::cout << YELLOW << "found ref_edge " << ref_edge_idx << " first" << RESET << std::endl;
      }
      // check for edge
      else if (indices[1] == edge.getIndex())
      {
        found_ref_edge_first = false;
        found_sth = true;
        // std::cout << YELLOW << "found edge " << edge.getIndex() << " first" << RESET << std::endl;
      }
      else
      {
        continue;
      }
    }

    if (found_sth)
    {
      double l_ij = iSData.L[indices[0]];
      double l_ik = iSData.L[indices[1]];
      double l_jk = iSData.L[indices[2]];

      if (l_jk < 0.0000000001)
      {
        std::cout << "find_tangent_space_angle: Calling angle_i with approx. zero length (length: " << l_jk << " for edge: " << indices[0] << ")" << std::endl;
      }
      angle += angle_i_from_lengths(l_ij, l_ik, l_jk, false, "find_tangent_space_angle");
      // bool result = validate_intrinsic_edge_lengths(iSData);
      // if (!result) { std::cout << "occurred in find_tangent_space_angle: " << std::endl; }
      // std::cout << "Updated angle: " << angle << std::endl;

      // terminate, the second time, we find sth
      // check for ref_edge or edge
      if (indices[0] == ref_edge_idx || indices[0] == edge.getIndex())
      {
        // std::cout << "found sth the second time... Breaking." << std::endl;
        break;
      }
    }
  }

  double angle_sum = total_angle(iSData, vertex);
  if (angle_sum == 0.0)
  {
    std::cout << RED << "found 0.0 angle_sum!" << std::endl;
  }
  if (!found_ref_edge_first)
  {
    angle = angle_sum - angle;
  }
  // std::cout << BLUE << "tangent result: " << angle / angle_sum << RESET << std::endl;
  return angle / angle_sum;
}

PolarVector2D parallel_transport(const iSimpData &iSData, PolarVector2D tangent_vector, const gcs::Vertex origin, const gcs::Vertex target)
{
  gcs::Edge connection = iSData.intrinsicMesh->connectingEdge(origin, target);
  double origin_connection_angle = find_tangent_space_angle(iSData, connection, origin);
  double target_connection_angle = find_tangent_space_angle(iSData, connection, target);
  double target_tangent_space_angle = tangent_vector[0] - origin_connection_angle + target_connection_angle + 0.5;

  // constrain resulting angle to the half-open interval [0,1[
  while (target_tangent_space_angle >= 1.0)
  {
    target_tangent_space_angle -= 1.0;
  }
  while (target_tangent_space_angle < 1.0)
  {
    target_tangent_space_angle += 1.0;
  }
  tangent_vector[0] = target_tangent_space_angle;

  return tangent_vector;
}

bool validate_intrinsic_edge_lengths(const iSimpData &iSData)
{
  for (gcs::Face F : iSData.intrinsicMesh->faces())
  {
    int i = 0;
    int l_ij_idx, l_ik_idx, l_jk_idx;
    for (gcs::Edge E : F.adjacentEdges())
    {
      if (i == 0)
      {
        l_ij_idx = E.getIndex();
      }
      else if (i == 1)
      {
        l_ik_idx = E.getIndex();
      }
      else if (i == 2)
      {
        l_jk_idx = E.getIndex();
      }
      i++;
    }

    double eps = 1e-6;
    double l_ij = iSData.L[l_ij_idx];
    double l_ik = iSData.L[l_ik_idx];
    double l_jk = iSData.L[l_jk_idx];

    // Check for zero lengths
    if (l_ij < eps)
    {
      std::cout << dye("validate_intrinsic_edge_lengths: Found zero intrinsic length edge! edge " + std::to_string(l_ij_idx) + "(l_ij: " + std::to_string(l_ij) + " < 0.0 )", RED) << std::endl;
      return false;
    }
    if (l_ik < eps)
    {
      std::cout << dye("validate_intrinsic_edge_lengths: Found zero intrinsic length edge! edge " + std::to_string(l_ik_idx) + "(l_ik: " + std::to_string(l_ik) + " < 0.0 )", RED) << std::endl;
      return false;
    }
    if (l_jk < eps)
    {
      std::cout << dye("validate_intrinsic_edge_lengths: Found zero intrinsic length edge! edge " + std::to_string(l_jk_idx) + "(l_jk: " + std::to_string(l_jk) + " < 0.0 )", RED) << std::endl;
      return false;
    }

    // Check triangle inequality
    if (l_ij >= l_ik + l_jk)
    {
      std::cout << dye("validate_intrinsic_edge_lengths: Found triangle inequality violation! triangle " + std::to_string(F.getIndex()) + "(l_ij: " + std::to_string(l_ij) + " > l_ik: " + std::to_string(l_ik) + " + l_jk: " + std::to_string(l_jk) + ")", RED) << std::endl;
      return false;
    }
    if (l_ik >= l_ij + l_jk)
    {
      std::cout << dye("validate_intrinsic_edge_lengths: Found triangle inequality violation! triangle " + std::to_string(F.getIndex()) + "(l_ik: " + std::to_string(l_ik) + " > l_ij: " + std::to_string(l_ij) + " + l_jk: " + std::to_string(l_jk) + ")", RED) << std::endl;
      return false;
    }
    if (l_jk >= l_ik + l_ij)
    {
      std::cout << dye("validate_intrinsic_edge_lengths: Found triangle inequality violation! triangle " + std::to_string(F.getIndex()) + "(l_jk: " + std::to_string(l_jk) + " > l_ik: " + std::to_string(l_ik) + " + l_ij: " + std::to_string(l_ij) + ")", RED) << std::endl;
      return false;
    }
  }
  return true;
}