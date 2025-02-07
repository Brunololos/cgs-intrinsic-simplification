#include "isimp_utils.hpp"

void register_point(iSimpData& iSData, const BarycentricPoint point, const TexCoord texcoord, const int face_idx)
{
  int point_idx = iSData.tracked_points.size();
  iSData.tracked_by_triangle[face_idx].push_back(point_idx);
  iSData.tracked_points.push_back(point);
  iSData.tracked_texcoords.push_back(texcoord);
}

void map_registered(iSimpData& iSData)
{
  iSData.mapped_points = iSData.tracked_points;
  iSData.mapped_by_triangle = iSData.tracked_by_triangle;

  for(int i=0; i<iSData.mapping.size(); i++)
  {
    iSData.mapping[i]->map_to_coarse(iSData.mapped_points, iSData.mapped_by_triangle);
  }
}

void map_registered(iSimpData& iSData, const int n_vertices)
{
  iSData.mapped_points = iSData.tracked_points;
  iSData.mapped_by_triangle = iSData.tracked_by_triangle;
  int remove_n_vertices = iSData.inputMesh->nVertices() - n_vertices;

  int removal_count = 0;
  for(int i=0; i<iSData.mapping.size() && removal_count < remove_n_vertices; i++)
  {
    iSData.mapping[i]->map_to_coarse(iSData.mapped_points, iSData.mapped_by_triangle);
    if (iSData.mapping[i]->op_type == SIMP_OP::V_REMOVAL) { removal_count++; }
  }
}

// For an intrinsic edge to be (intrinsically) flippable two conditions have to hold:
//   1. the edge adjacent triangles have to form a convex quad
//   2. the edge endpoints have degree two or more
bool is_edge_flippable(const iSimpData& iSData, const Quad2D& quad, const gcs::Edge edge)
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
    // std::cout << "Quad is nonconvex" << std::endl;
    return false;
  }
  return true;
}

bool flip_intrinsic(iSimpData& iSData, const gcs::Edge edge)
{
  gcs::Face Fk;
  gcs::Face Fl;
  bool t = true;
  for (gcs::Face F : edge.adjacentFaces())
  {
    if (t)
    {
      Fk = F;
      t = false;
    }
    else {
      Fl = F;
    }
  }

  std::pair<int, int> Fk_unique = find_unique_vertex_index(Fk, edge);
  std::pair<int, int> Fl_unique = find_unique_vertex_index(Fl, edge);
  // get ordered intrinsic edge lengths from neighborhood
  std::array<int, 5> edge_indices = order_quad_edge_indices(edge, Fk_unique.second, Fl_unique.second);
  std::array<int, 4> vertex_indices = order_quad_vertex_indices(edge, Fk_unique.second, Fl_unique.second);
  double l_ij = iSData.L(edge_indices[0]);
  double l_ik = iSData.L(edge_indices[1]);
  double l_jk = iSData.L(edge_indices[2]);
  double l_il = iSData.L(edge_indices[3]);
  double l_jl = iSData.L(edge_indices[4]);
  Quad2D unfolded = unfold(l_ij, l_ik, l_jk, l_il, l_jl);

  if (!is_edge_flippable(iSData, unfolded, edge)) { return false; }
  // update tangent space reference edges that would be flipped
  if (iSData.tangent_space_reference_edges[edge.firstVertex().getIndex()] == edge.getIndex()) { rotate_reference_edge(iSData, edge.firstVertex()); }
  if (iSData.tangent_space_reference_edges[edge.secondVertex().getIndex()] == edge.getIndex()) { rotate_reference_edge(iSData, edge.secondVertex()); }
  bool couldFlip = iSData.intrinsicMesh->flip(edge, false);
  if (!couldFlip) { return false; }

  // update edge length
  iSData.L(edge.getIndex()) = flipped_edgelength(unfolded);
  // add intrinsic flip into mapping queue
  iSData.mapping.push_back(std::make_unique<Edge_Flip>(Fk.getIndex(), Fl.getIndex(), unfolded, Fk_unique.first, Fl_unique.first));

  return true;
}

bool flatten_vertex(iSimpData& iSData, const int vertex_idx)
{
  gcs::Vertex v = iSData.intrinsicMesh->vertex(vertex_idx);
  double threshold = 0.0001;

  // check if vertex is interior or boundary vertex
  double THETA_HAT = (!v.isBoundary()) ? 2.0 * M_PI : M_PI;
  double u_i = 0.0;

  // TODO: if vertex is boundary vertex, check if it is incident on a boundary self-edge, or is a vertex of a self-face. If it is, it is skipped and flattening fails.
  std::unordered_map<int, double> L_opt = std::unordered_map<int, double>();
  for (gcs::Edge E : v.adjacentEdges())
  {
    int idx = E.getIndex();
    L_opt.insert({idx, iSData.L[idx]});
  }

  double step = 1000.0;
  int i = 0;
  while (std::abs(step) > threshold)
  {
    // Calculate newton step
    double enumerator = THETA_HAT;
    double denominator = 0.0;
    for (gcs::Face F : v.adjacentFaces())
    {
      // get edge lengths
      std::array<int, 3> edge_indices = order_triangle_edge_indices(F, vertex_idx);
      double l_ij = L_opt[edge_indices[0]];
      double l_ik = L_opt[edge_indices[1]];
      double l_jk = iSData.L[edge_indices[2]];

      if (l_jk < 0.0000000001) { std::cout << "flatten_vertex: Calling angle_i with zero length" << std::endl; }
      enumerator -= angle_i_from_lengths(l_ij, l_ik, l_jk);
      if (l_ij < 0.0000000001) { std::cout << "flatten_vertex: Calling angle_i with zero length" << std::endl; }
      double theta_ij = angle_i_from_lengths(l_ik, l_jk, l_ij);
      if (l_ik < 0.0000000001) { std::cout << "flatten_vertex: Calling angle_i with zero length" << std::endl; }
      double theta_ik = angle_i_from_lengths(l_ij, l_jk, l_ik);
      denominator += 1.0/std::tan(theta_ij);
      denominator += 1.0/std::tan(theta_ik);
    }
    denominator /= 2.0;

    step = enumerator / denominator;
    u_i = u_i - step;

    // safety check for NaN
    if (u_i != u_i) { std::cout << dye("Encountered NaN!", RED) << std::endl; return false; }
    // update edge lengths in L_opt with u_i as a parameter
    for (const auto& pair : L_opt)
    {
      L_opt[pair.first] = std::exp(u_i / 2.0) * iSData.L[pair.first];
      gcs::Edge E = iSData.intrinsicMesh->edge(pair.first);
      // std::cout << "  (" << pair.first << ": (" << E.firstVertex().getIndex() << ", " << E.secondVertex().getIndex() << ")" << ", " << L_opt[pair.first] << ")" << std::endl;
    }
    i++;
  }

  double eps = 0.0000001;
  // if u_i is close to zero, the flattening is an identity operation (so we return a success, but don't insert a mapping operation into the mapping)
  if(u_i < eps) { return true; }

  // check if length's satisfy the triangle inequality
  for (gcs::Face F : v.adjacentFaces())
  {
    // get edge lengths
    std::array<int, 3> edge_indices = order_triangle_edge_indices(F, vertex_idx);
    double l_ij = L_opt[edge_indices[0]];
    double l_ik = L_opt[edge_indices[1]];
    double l_jk = iSData.L[edge_indices[2]];

    if ((l_ij > l_ik + l_jk + eps)
    ||  (l_ik > l_ij + l_jk + eps)
    ||  (l_jk > l_ij + l_ik + eps))
    {
      // TODO: perform edge flips to try to mitigate it
      std::cout << dye("Violated triangle inequality", RED) << std::endl;
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

  // add vertex flattening operation into mapping queue
  for (const auto& pair : L_opt)
  {
    iSData.L[pair.first] = pair.second;
  }
  std::vector<int> faces = std::vector<int>();
  std::vector<int> face_vertex_indices = std::vector<int>();
  for (gcs::Face F : v.adjacentFaces())
  {
    faces.push_back(F.getIndex());
    face_vertex_indices.push_back(find_vertex_face_index(F, vertex_idx));
  }
  iSData.mapping.push_back(std::make_unique<Vertex_Flattening>(faces, face_vertex_indices, u_i));

  // calculate weighting
  double total_curvature_change = 0.0;
  for (auto& pair : redistribution_weights)
  {
    gcs::Vertex neighbor = iSData.intrinsicMesh->vertex(pair.first);
    pair.second = std::abs(gaussian_curvature(iSData, neighbor) - pair.second);
    total_curvature_change += pair.second;
  }
  for (auto& pair : redistribution_weights)
  {
    // if already flat, distribute evenly
    if (total_curvature_change == 0.0)
    {
      pair.second = 1.0 / ((double) v.degree());
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

bool remove_vertex(iSimpData& iSData, const int vertex_idx)
{
      // check if vertex is intrinsically flat (we can only remove it, when it is)
      gcs::Vertex vertex = iSData.intrinsicMesh->vertex(vertex_idx);
      double THETA_HAT = (!vertex.isBoundary()) ? 2.0 * M_PI : M_PI;
      if (vertex.isDead())
      {
        std::cout << "Vertex has already been removed from mesh." << std::endl;
        return false;
      }

      // if (vertex.degree() < 3)
      // {
      //   std::cout << "Vertex has insufficient degree for removal." << std::endl;
      //   return false;
      // }

      // NOTE: gcs::manifold_surface_mesh implementation does not allow boundary vertex removal
      if (vertex.isBoundary())
      {
        std::cout << "Vertex is boundary vertex." << std::endl;
        return false;
      }

      // if (THETA_HAT - angle_sum != 0.0)
      // {
      //   std::cout << "Vertex has bad angle sum: " << angle_sum << std::endl;
      //   return false;
      // }

      // check if vertex is degree 3
      if (vertex.degree() != 3)
      {
        // TODO: try to perform edge flips to get to degree 3
        return false;
      }

      int F_jk = -1;
      int F_jl = -1;
      int F_kl = -1;
      int F_res;
      int F_res_permutation;
      for (gcs::Face F : vertex.adjacentFaces())
      {

        if(F_jk == -1) { F_jk = F.getIndex(); continue; }
        if(F_jl == -1) { F_jl = F.getIndex(); continue; }
        if(F_kl == -1) { F_kl = F.getIndex(); break; }
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
      // std::cout << "Determined F_jk: " << F_jk << ", F_jl: " << F_jl << std::endl;

      int i = 0;
      // printGCSFace(iSData.intrinsicMesh->face(F_jk));
      // printGCSFace(iSData.intrinsicMesh->face(F_jl));
      // std::cout << "Determined e_ij: " << e_ij.getIndex() << " (" << e_ij.firstVertex().getIndex() << "," << e_ij.secondVertex().getIndex() << ")" << std::endl;
      // std::cout << "AdjacentFaces: " << std::endl;
      // for (gcs::Face F : vertex.adjacentFaces())
      // {
      //   printGCSFace(F);
      // }
      i = 0;

      std::pair<int, int> Fk_unique = find_unique_vertex_index(iSData.intrinsicMesh->face(F_jk), e_ij);
      std::pair<int, int> Fl_unique = find_unique_vertex_index(iSData.intrinsicMesh->face(F_jl), e_ij);
      // std::cout << "Determined unique_v's" << std::endl;
      // std::cout << "Fk_unique: " << Fk_unique.second << std::endl;
      // std::cout << "Fl_unique: " << Fl_unique.second << std::endl;
      // get ordered intrinsic edge lengths from neighborhood
      std::array<int, 5> edge_indices = order_quad_edge_indices(e_ij, Fk_unique.second, Fl_unique.second, vertex_idx);
      double l_ij = iSData.L(edge_indices[0]);
      double l_ik = iSData.L(edge_indices[1]);
      double l_jk = iSData.L(edge_indices[2]);
      double l_il = iSData.L(edge_indices[3]);
      double l_jl = iSData.L(edge_indices[4]);
      Quad2D unfolded = unfold(l_ij, l_ik, l_jk, l_il, l_jl);
      // std::cout << "Determined unfolding" << std::endl;

      // determining at which position in each face the shared vertex resides
      std::vector<int> shared_face_vertex_indices = std::vector<int>();
      // std::cout << "Removing vertex: " << vertex_idx << " with faces: {";
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

      // std::cout << "predicted F_res: " << std::endl;
      // switch (F_res_permutation)
      // {
      //     default:
      //     case 0:
      //       std::cout << "vj, vk, vl: (" << ordered_verts[1] << ", " << ordered_verts[2] << ", " << ordered_verts[3] << ")" << std::endl;
      //       break;
      //     case 1:
      //       std::cout << "vl, vj, vk: (" << ordered_verts[3] << ", " << ordered_verts[1] << ", " << ordered_verts[2] << ")" << std::endl;
      //       break;
      //     case 2:
      //       std::cout << "vk, vl, vj: (" << ordered_verts[2] << ", " << ordered_verts[3] << ", " << ordered_verts[1] << ")" << std::endl;
      //       break;
      // }
      // std::cout << "Face after removal: " << std::endl;
      // printGCSFace(iSData.intrinsicMesh->face(F_res));
      // std::cout << "F_jk: ";
      // switch (F_jk_shared)
      // {
      //     default:
      //     case 0:
      //       std::cout << "vi, vj, vk: (" << ordered_verts[0] << ", " << ordered_verts[1] << ", " << ordered_verts[2] << ")" << std::endl;
      //       break;
      //     case 1:
      //       std::cout << "vk, vi, vj: (" << ordered_verts[2] << ", " << ordered_verts[0] << ", " << ordered_verts[1] << ")" << std::endl;
      //       break;
      //     case 2:
      //       std::cout << "vj, vk, vi: (" << ordered_verts[1] << ", " << ordered_verts[2] << ", " << ordered_verts[0] << ")" << std::endl;
      //       break;
      // }
      // std::cout << "F_jl: ";
      // switch (F_jl_shared)
      // {
      //     default:
      //     case 0:
      //       std::cout << "vi, vl, vj: (" << ordered_verts[0] << ", " << ordered_verts[3] << ", " << ordered_verts[1] << ")" << std::endl;
      //       break;
      //     case 1:
      //       std::cout << "vj, vi, vl: (" << ordered_verts[1] << ", " << ordered_verts[0] << ", " << ordered_verts[3] << ")" << std::endl;
      //       break;
      //     case 2:
      //       std::cout << "vl, vj, vi: (" << ordered_verts[3] << ", " << ordered_verts[1] << ", " << ordered_verts[0] << ")" << std::endl;
      //       break;
      // }
      // std::cout << "F_kl: ";
      // switch (F_kl_shared)
      // {
      //     default:
      //     case 0:
      //       std::cout << "vi, vk, vl: (" << ordered_verts[0] << ", " << ordered_verts[2] << ", " << ordered_verts[3] << ")" << std::endl;
      //       break;
      //     case 1:
      //       std::cout << "vl, vi, vk: (" << ordered_verts[3] << ", " << ordered_verts[0] << ", " << ordered_verts[2] << ")" << std::endl;
      //       break;
      //     case 2:
      //       std::cout << "vk, vl, vi: (" << ordered_verts[2] << ", " << ordered_verts[3] << ", " << ordered_verts[0] << ")" << std::endl;
      //       break;
      // }
      // NOTE: I tried to figure out how geometry central determines the permutation of vertex indices of the resulting face (e.g. a face ijk could be stored as ijk, jki or kij, and we need to be exact here because this affects our barycentric coordinates)
      // but I gave up. It mostly keeps the order of F_jk, but some unaccounted edge cases remain.
      // Instead I just look what order gcs has determined and pass that.

      // insert vertex removal into mapping queue
      iSData.mapping.push_back(std::make_unique<Vertex_Removal>(F_idxs, F_res, unfolded, shared_face_vertex_indices, F_res_permutation));
      return true;
}

// greedily flip edges, until delaunay
void flip_to_delaunay(iSimpData& iSData)
{
  for (gcs::Edge E : iSData.intrinsicMesh->edges())
  {
    if(E.isBoundary()) { continue; }
    gcs::Face Fk;
    gcs::Face Fl;
    bool t = true;
    for (gcs::Face F : E.adjacentFaces())
    {
      if (t) {
        Fk = F;
        t = false;
      } else {
        Fl = F;
      }
    }
    std::pair<int, int> Fk_unique = find_unique_vertex_index(Fk, E);
    std::pair<int, int> Fl_unique = find_unique_vertex_index(Fl, E);
    std::array<int, 5> edge_indices = order_quad_edge_indices(E, Fk_unique.second, Fl_unique.second);
    double l_ij = iSData.L(edge_indices[0]);
    double l_ik = iSData.L(edge_indices[1]);
    double l_jk = iSData.L(edge_indices[2]);
    double l_il = iSData.L(edge_indices[3]);
    double l_jl = iSData.L(edge_indices[4]);
    if (l_ij < 0.0000000001) { std::cout << "flip_to_delaunay: Calling angle_i with zero length" << std::endl; }
    double theta_k = angle_i_from_lengths(l_ik, l_jk, l_ij);
    if (l_ij < 0.0000000001) { std::cout << "flip_to_delaunay: Calling angle_i with zero length" << std::endl; }
    double theta_l = angle_i_from_lengths(l_il, l_jl, l_ij);
    if(theta_k + theta_l > M_PI)
    {
      flip_intrinsic(iSData, E);
    }
  }

  // bool found_flip = true;
  // while(found_flip)
  // {
  //   found_flip = false;
  //   for (gcs::Edge E : iSData.intrinsicMesh->edges())
  //   {
  //     if(E.isBoundary()) { continue; }
  //     gcs::Face Fk;
  //     gcs::Face Fl;
  //     bool t = true;
  //     for (gcs::Face F : E.adjacentFaces())
  //     {
  //       if (t) {
  //         Fk = F;
  //         t = false;
  //       } else {
  //         Fl = F;
  //       }
  //     }
  //     std::pair<int, int> Fk_unique = find_unique_vertex_index(Fk, E);
  //     std::pair<int, int> Fl_unique = find_unique_vertex_index(Fl, E);
  //     std::array<int, 5> edge_indices = order_quad_edge_indices(E, Fk_unique.second, Fl_unique.second);
  //     double l_ij = iSData.L(edge_indices[0]);
  //     double l_ik = iSData.L(edge_indices[1]);
  //     double l_jk = iSData.L(edge_indices[2]);
  //     double l_il = iSData.L(edge_indices[3]);
  //     double l_jl = iSData.L(edge_indices[4]);
  //     double theta_k = angle_i_from_lengths(l_ik, l_jk, l_ij);
  //     double theta_l = angle_i_from_lengths(l_il, l_jl, l_ij);
  //     if(theta_k + theta_l > M_PI)
  //     {
  //       found_flip = flip_intrinsic(iSData, E);
  //       if (found_flip) { break; }
  //     }
  //   }
  // }
}

bool flip_vertex_to_deg3(iSimpData& iSData, const int vertex_idx)
{
  // stack of performed edge flips in case we need to revert them
  std::vector<int> flips = std::vector<int>();

  gcs::Vertex V = iSData.intrinsicMesh->vertex(vertex_idx);
  if (V.degree() == 3) { return true; }

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
    }
    // found valid edge flip sequence
    if (V.degree() == 3)
    {
      return true;
    }
  }

  // TODO: the edge flip reversion currently leads to some issues
  // I think that a case (possibly due to numerics) can arise, where a first edge flip can be performed, but not the three edge flips to undo it.
  // Because in my reasoning a quad that can be flipped once, should be able to be flipped any number of times.
  // // revert edge-flips
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

void rotate_reference_edge(iSimpData& iSData, const gcs::Vertex& vertex)
{
  if (vertex.degree() < 2) { return; }
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
  if (iSData.T_minus(vertex_idx, 0) < 0.0) { iSData.T_minus(vertex_idx, 0) += 1.0; }
  if (iSData.T_plus(vertex_idx, 0) < 0.0) { iSData.T_plus(vertex_idx, 0) += 1.0; }
  iSData.tangent_space_reference_edges[vertex_idx] = next_ref_edge.getIndex();
}

PolarVector2D next_error_vector(const iSimpData& iSData, const double alpha, const gcs::Vertex vertex, const gcs::Vertex neighbor, bool for_T_minus)
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
  // printEigenVector2d(T);
  // printEigenVector2d(nT);
  // std::cout << "mi: " << mi << ", mj: " << mj << std::endl;

  gcs::Edge connection = iSData.intrinsicMesh->connectingEdge(vertex, neighbor);
  // TODO:
  // std::cout << connection << std::endl;
  // if (connection.isDead()) { std::cout << RED << "Chose dead connection!" << RESET << std::endl; }
  if (connection.isBoundary()) { std::cout << RED << "Chose dead connection!" << RESET << std::endl; }
  if (connection.firstVertex().isDead()) { std::cout << RED << "Chose dead connection!" << RESET << std::endl; return PolarVector2D(0.0, 0.0); }
  if (connection.secondVertex().isDead()) { std::cout << RED << "Chose dead connection!" << RESET << std::endl; return PolarVector2D(0.0, 0.0); }
  double edge_angle = find_tangent_space_angle(iSData, connection, neighbor);
  // TODO: 
  double edge_length = iSData.L[connection.getIndex()];
  Vector2D transport = to_cartesian({edge_angle, edge_length});
  Vector2D transported = to_cartesian(parallel_transport(iSData, T, vertex, neighbor));

  if (mi == 0.0 && mj == 0.0) { return PolarVector2D(0.0, 0.0); }
  PolarVector2D next_T = to_polar((alpha * mi * (transported + transport) + mj * nT) / (alpha * mi + mj));
  return next_T;
}

double intrinsic_curvature_error(const iSimpData& iSData, const gcs::Vertex vertex)
{
  if (vertex.isBoundary()) { return INFINITY; } // guard, because we chose to forgo boundary vertex simplification
  double ice = 0.0;
  double threshold = 0.0001;
  int vertex_idx = vertex.getIndex();

  // check if vertex is interior or boundary vertex
  double THETA_HAT = (!vertex.isBoundary()) ? 2.0 * M_PI : M_PI;
  double u_i = 0.0;

  // TODO: if vertex is boundary vertex, check if it is incident on a boundary self-edge, or is a vertex of a self-face. If it is, it is skipped and flattening fails.
  std::unordered_map<int, double> L_opt = std::unordered_map<int, double>();
  for (gcs::Edge E : vertex.adjacentEdges())
  {
    int idx = E.getIndex();
    L_opt.insert({idx, iSData.L[idx]});
  }

  double step = 1000.0;
  int i = 0;
  while (std::abs(step) > threshold)
  {
    // Calculate newton step
    double enumerator = THETA_HAT;
    double denominator = 0.0;
    for (gcs::Face F : vertex.adjacentFaces())
    {
      // get edge lengths
      std::array<int, 3> edge_indices = order_triangle_edge_indices(F, vertex_idx);
      double l_ij = L_opt[edge_indices[0]];
      double l_ik = L_opt[edge_indices[1]];
      double l_jk = iSData.L[edge_indices[2]];

      if (l_jk < 0.0000000001) { std::cout << "intrinsic_curvature_error: Calling angle_i with zero length" << std::endl; }
      enumerator -= angle_i_from_lengths(l_ij, l_ik, l_jk);
      if (l_ij < 0.0000000001) { std::cout << "intrinsic_curvature_error: Calling angle_i with zero length" << std::endl; }
      double theta_ij = angle_i_from_lengths(l_ik, l_jk, l_ij);
      if (l_ik < 0.0000000001) { std::cout << "intrinsic_curvature_error: Calling angle_i with zero length" << std::endl; }
      double theta_ik = angle_i_from_lengths(l_ij, l_jk, l_ik);
      denominator += 1.0/std::tan(theta_ij);
      denominator += 1.0/std::tan(theta_ik);
    }
    denominator /= 2.0;

    step = enumerator / denominator;
    u_i = u_i - step;

    // safety check for NaN
    if (u_i != u_i) { std::cout << dye("Encountered NaN!", RED) << std::endl; return INFINITY; }
    // update edge lengths in L_opt with u_i as a parameter
    for (const auto& pair : L_opt)
    {
      L_opt[pair.first] = std::exp(u_i / 2.0) * iSData.L[pair.first];
      gcs::Edge E = iSData.intrinsicMesh->edge(pair.first);
    }
    i++;
  }

  double eps = 0.0000001;

  // check if length's satisfy the triangle inequality
  for (gcs::Face F : vertex.adjacentFaces())
  {
    // get edge lengths
    std::array<int, 3> edge_indices = order_triangle_edge_indices(F, vertex_idx);
    double l_ij = L_opt[edge_indices[0]];
    double l_ik = L_opt[edge_indices[1]];
    double l_jk = iSData.L[edge_indices[2]];

    if ((l_ij > l_ik + l_jk + eps)
    ||  (l_ik > l_ij + l_jk + eps)
    ||  (l_jk > l_ij + l_ik + eps))
    {
      // TODO: perform edge flips to try to mitigate it
      std::cout << dye("Violated triangle inequality", RED) << std::endl;
      return INFINITY;
    }
  }

  // std::cout << "previous curvatures: " << std::endl;
  // save previous curvatures for weighted mass redistribution
  std::unordered_map<int, double> redistribution_weights = std::unordered_map<int, double>();
  for (gcs::Vertex neighbor : vertex.adjacentVertices())
  {
    int neighbor_idx = neighbor.getIndex();
    redistribution_weights.insert({neighbor_idx, gaussian_curvature(iSData, neighbor)});
    // std::cout << " - " << neighbor_idx << ": " << gaussian_curvature(iSData, neighbor) << std::endl;
  }

  // calculate weighting
  double total_curvature_change = 0.0;
  for (auto& pair : redistribution_weights)
  {
    gcs::Vertex neighbor = iSData.intrinsicMesh->vertex(pair.first);
    pair.second = std::abs(gaussian_curvature(iSData, neighbor) - pair.second);
    total_curvature_change += pair.second;
  }
  // std::cout << "total curvature change: " << total_curvature_change << std::endl;
  // std::cout << "alphas: " << std::endl;
  for (auto& pair : redistribution_weights)
  {
    // if already flat, distribute evenly
    if (total_curvature_change == 0.0)
    {
      pair.second = 1.0 / ((double) vertex.degree());
      // std::cout << " -> " << pair.second << std::endl;
      continue;
    }
    pair.second = pair.second / total_curvature_change;
    // std::cout << " -> " << pair.second << std::endl;
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
    // std::cout << "calculate ice term for vertex neighbor: " << neighbor_idx << std::endl;
    // std::cout << " -> alpha = " << alpha << std::endl;
    // printEigenVector2d(T_minus);
    // printEigenVector2d(T_plus);

    double m_minus = iSData.M(neighbor_idx, 0) + alpha * iSData.M(vertex_idx, 0);
    double m_plus = iSData.M(neighbor_idx, 1) + alpha * iSData.M(vertex_idx, 1);
    ice += m_minus * T_minus[1];
    ice += m_plus * T_plus[1];
    // std::cout << "summing up ice: " << ice << std::endl;
  }
  return ice;
}

// orders triangle vertex indices of triangle ijk by a vertex i
// index 0: vertex i
// index 1: vertex j of ijk-triangle
// index 2: vertex k of ijk-triangle
std::array<int, 3> order_triangle_vertex_indices(const gcs::Face& face, const int F_first_v)
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

// orders triangle edge indices of triangle ijk by a vertex i
// index 0: shared edge ij
// index 1: first edge ik of ijk-triangle
// index 2: second edge jk of ijk-triangle
std::array<int, 3> order_triangle_edge_indices(const gcs::Face& face, const int F_first_v)
{
  std::array<int, 3> indices = std::array<int, 3>();
  std::array<int, 3> v_ordered = order_triangle_vertex_indices(face, F_first_v);

  for (gcs::Edge E : face.adjacentEdges())
  {
    int e0 = E.firstVertex().getIndex();
    int e1 = E.secondVertex().getIndex();
    if (e0 == v_ordered[0] && e1 == v_ordered[1]
    ||  e0 == v_ordered[1] && e1 == v_ordered[0])
    {
      indices[0] = E.getIndex();
    }
    else if (e0 == v_ordered[0] && e1 == v_ordered[2]
         ||  e0 == v_ordered[2] && e1 == v_ordered[0])
    {
      indices[1] = E.getIndex();
    }
    else if (e0 == v_ordered[1] && e1 == v_ordered[2]
         ||  e0 == v_ordered[2] && e1 == v_ordered[1])
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
std::array<int, 5> order_quad_edge_indices(const gcs::Edge& edge, const int Fk_unique_v, const int Fl_unique_v)
{
  std::array<int, 5> indices = std::array<int, 5>();
  indices[0] = edge.getIndex();

  int offset = 1;
  for (gcs::Face F : edge.adjacentFaces())
  {
    for (gcs::Edge E : F.adjacentEdges())
    {
      if ((E.firstVertex().getIndex() == edge.firstVertex().getIndex()
      && E.secondVertex().getIndex() == Fk_unique_v)
      || (E.secondVertex().getIndex() == edge.firstVertex().getIndex()
      && E.firstVertex().getIndex() == Fk_unique_v))
      {
        indices[1] = E.getIndex();
      }
      else if ((E.firstVertex().getIndex() == edge.secondVertex().getIndex()
      && E.secondVertex().getIndex() == Fk_unique_v)
      || (E.secondVertex().getIndex() == edge.secondVertex().getIndex()
      && E.firstVertex().getIndex() == Fk_unique_v))
      {
        indices[2] = E.getIndex();
      }
      else if ((E.firstVertex().getIndex() == edge.firstVertex().getIndex()
      && E.secondVertex().getIndex() == Fl_unique_v)
      || (E.secondVertex().getIndex() == edge.firstVertex().getIndex()
      && E.firstVertex().getIndex() == Fl_unique_v))
      {
        indices[3] = E.getIndex();
      }
      else if ((E.firstVertex().getIndex() == edge.secondVertex().getIndex()
      && E.secondVertex().getIndex() == Fl_unique_v)
      || (E.secondVertex().getIndex() == edge.secondVertex().getIndex()
      && E.firstVertex().getIndex() == Fl_unique_v))
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
std::array<int, 5> order_quad_edge_indices(const gcs::Edge& edge, const int Fk_unique_v, const int Fl_unique_v, const int vi_idx)
{
  std::array<int, 5> indices = std::array<int, 5>();
  indices[0] = edge.getIndex();

  int offset = 1;
  for (gcs::Face F : edge.adjacentFaces())
  {
    for (gcs::Edge E : F.adjacentEdges())
    {
      if ((E.firstVertex().getIndex() == vi_idx
      && E.secondVertex().getIndex() == Fk_unique_v)
      || (E.secondVertex().getIndex() == vi_idx
      && E.firstVertex().getIndex() == Fk_unique_v))
      {
        indices[1] = E.getIndex();
      }
      else if ((E.firstVertex().getIndex() != vi_idx
      && E.secondVertex().getIndex() == Fk_unique_v)
      || (E.secondVertex().getIndex() != vi_idx
      && E.firstVertex().getIndex() == Fk_unique_v))
      {
        indices[2] = E.getIndex();
      }
      else if ((E.firstVertex().getIndex() == vi_idx
      && E.secondVertex().getIndex() == Fl_unique_v)
      || (E.secondVertex().getIndex() == vi_idx
      && E.firstVertex().getIndex() == Fl_unique_v))
      {
        indices[3] = E.getIndex();
      }
      else if ((E.firstVertex().getIndex() != vi_idx
      && E.secondVertex().getIndex() == Fl_unique_v)
      || (E.secondVertex().getIndex() != vi_idx
      && E.firstVertex().getIndex() == Fl_unique_v))
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
std::array<int, 4> order_quad_vertex_indices(const gcs::Edge& edge)
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
      if (V.getIndex() != edge.firstVertex().getIndex()
      && V.getIndex() != edge.secondVertex().getIndex())
      {
        indices[offset] = V.getIndex();
        break;
      }
    }
    offset++;
  }
  return indices;
}

std::array<int, 4> order_quad_vertex_indices(const gcs::Edge& edge, const int Fk_unique_v, const int Fl_unique_v)
{
  // because in the ManifoldSurfaceMesh data structure the half-edges form counterclockwise cycles, this should work
  std::array<int, 4> indices = std::array<int, 4>();
  indices[0] = edge.firstVertex().getIndex();
  indices[1] = edge.secondVertex().getIndex();
  indices[2] = Fk_unique_v;
  indices[3] = Fl_unique_v;
  return indices;
}

std::array<int, 4> order_quad_vertex_indices(const gcs::Edge& edge, const int Fk_unique_v, const int Fl_unique_v, const int vi_idx)
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

std::pair<int, int> find_unique_vertex_index(const gcs::Face& F, const gcs::Edge& edge)
{
  int i = 0;
  for (gcs::Vertex V : F.adjacentVertices())
  {
    if (V.getIndex() != edge.firstVertex().getIndex()
    && V.getIndex() != edge.secondVertex().getIndex())
    {
      return std::pair<int, int>(i, V.getIndex());
    }
    i++;
  }
  return std::pair<int, int>(-1, -1);
}

int find_vertex_face_index(const gcs::Face& F, const int vertex)
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

double total_angle(const iSimpData& iSData, const gcs::Vertex& vertex)
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
    if (l_jk < 0.0000000001) { std::cout << "Calling angle_i from total_angle with zero length" << std::endl; }

    total_angle += angle_i_from_lengths(l_ij, l_ik, l_jk);
  }

  return total_angle;
}


// NOTE: Maybe the function name is a misnomer, bc it also calculates geodesic curvature for boundary vertices
double gaussian_curvature(const iSimpData& iSData, const gcs::Vertex& vertex)
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

    if (l_jk < 0.0000000001) { std::cout << "gaussian_curvature: Calling angle_i with zero length" << std::endl; }
    curvature -= angle_i_from_lengths(l_ij, l_ik, l_jk);
  }

  return curvature;
}

double find_tangent_space_angle(const iSimpData& iSData, const gcs::Edge& edge, const gcs::Vertex& vertex)
{
  double angle = 0.0;
  // std::cout << "finding tangent space angle of edge: " << edge.getIndex() << " around vertex " << vertex.getIndex() << " with reference edge: " << iSData.intrinsicMesh->edge(iSData.tangent_space_reference_edges[vertex.getIndex()]).getIndex() << std::endl;

  // get reference edge
  int vertex_idx = vertex.getIndex();
  int ref_edge_idx = iSData.tangent_space_reference_edges[vertex_idx];
  gcs::Edge ref_edge = iSData.intrinsicMesh->edge(ref_edge_idx);
  if (edge.getIndex() == ref_edge_idx) { return angle; }
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
      else { continue; }
    }

    if (found_sth)
    {
      double l_ij = iSData.L[indices[0]];
      double l_ik = iSData.L[indices[1]];
      double l_jk = iSData.L[indices[2]];

      if (l_jk < 0.0000000001) { std::cout << "find_tangent_space_angle: Calling angle_i with approx. zero length (length: " << l_jk << " for edge: " << indices[0] << ")" << std::endl; }
      angle += angle_i_from_lengths(l_ij, l_ik, l_jk);
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
  if (angle_sum == 0.0) { std::cout << RED << "found 0.0 angle_sum!" << std::endl; }
  if (!found_ref_edge_first) { angle = angle_sum - angle; }
  // std::cout << BLUE << "tangent result: " << angle / angle_sum << RESET << std::endl;
  return angle / angle_sum;
}

PolarVector2D parallel_transport(const iSimpData& iSData, PolarVector2D tangent_vector, const gcs::Vertex origin, const gcs::Vertex target)
{
  gcs::Edge connection = iSData.intrinsicMesh->connectingEdge(origin, target);
  double origin_connection_angle = find_tangent_space_angle(iSData, connection, origin);
  double target_connection_angle = find_tangent_space_angle(iSData, connection, target);
  double target_tangent_space_angle = tangent_vector[0] - origin_connection_angle + target_connection_angle + 0.5;

  // constrain resulting angle to the half-open interval [0,1[
  while (target_tangent_space_angle >= 1.0) { target_tangent_space_angle -= 1.0; }
  while (target_tangent_space_angle < 1.0) { target_tangent_space_angle += 1.0; }
  tangent_vector[0] = target_tangent_space_angle;

  return tangent_vector;
}