#include "isimp_utils.hpp"

void register_point(iSimpData& iSData, const BarycentricPoint point, const int face_idx)
{
  int point_idx = iSData.tracked_points.size();
  iSData.tracked_by_triangle[face_idx].push_back(point_idx);
  iSData.tracked_points.push_back(point);
}

void map_registered(iSimpData& iSData)
{
  iSData.mapped_points = iSData.tracked_points;
  iSData.mapped_by_triangle = iSData.tracked_by_triangle;

  for(int i=0; i<iSData.mapping.size(); i++)
  {
    // Mapping_operation op = iSData.mapping[i];
    // Edge_Flip ef;
    // Vertex_Flattening vf;
    // Vertex_Removal vr;
    // switch(iSData.mapping[i].op_type)
    // {
    //   default:
    //   case SIMP_OP::E_FLIP:
    //     Edge_Flip ef = iSData.mapping[i];
    //     break;
    // }
    // op.map_to_coarse(iSData.mapped_points, iSData.mapped_by_triangle);
    iSData.mapping[i]->map_to_coarse(iSData.mapped_points, iSData.mapped_by_triangle);
  }
}

// bool flip_intrinsic(/* std::unique_ptr<gcs::ManifoldSurfaceMesh> intrinsicMesh */std::shared_ptr<gcs::ManifoldSurfaceMesh> intrinsicMesh, Eigen::Vector<double, -1> L, gcs::Edge edge)
bool flip_intrinsic(iSimpData& iSData, const gcs::Edge edge)
{
  // -------------------------------------------------------------v
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
  std::cout << "Flipping edge: " << edge.getIndex() << " with incident faces Fk: " << Fk.getIndex() << ", Fl: " << Fl.getIndex() << std::endl;
  printGCSFace("Fk", Fk);
  printGCSFace("Fl", Fl);
  // ------------------------------------------------------------------------^
  std::pair<int, int> Fk_unique = find_unique_vertex_index(Fk, edge);
  std::pair<int, int> Fl_unique = find_unique_vertex_index(Fl, edge);
  std::cout << "Fk_unique_idx: " << Fk_unique.first << std::endl;
  std::cout << "Fl_unique_idx: " << Fl_unique.first << std::endl;
  // get ordered intrinsic edge lengths from neighborhood
  std::array<int, 5> edge_indices = order_quad_edge_indices(edge, Fk_unique.second, Fl_unique.second);
  std::array<int, 4> vertex_indices = order_quad_vertex_indices(edge, Fk_unique.second, Fl_unique.second);
  std::cout << "ordered quad edge indices: " << edge_indices[0]  << edge_indices[1]  << edge_indices[2]  << edge_indices[3]  << edge_indices[4] << std::endl;
  double l_ij = iSData.L(edge_indices[0]);
  double l_ik = iSData.L(edge_indices[1]);
  double l_jk = iSData.L(edge_indices[2]);
  double l_il = iSData.L(edge_indices[3]);
  double l_jl = iSData.L(edge_indices[4]);
  std::cout << "lengths: l_ij: " << l_ij << ", l_ik: " << l_ik << ", l_jk: " << l_jk << ", l_il: " << l_il << ", l_jl: " << l_jl << std::endl;
  Quad2D unfolded = unfold(l_ij, l_ik, l_jk, l_il, l_jl);
  std::cout << "Planar unfolding of edge lengths:\n";
  printEigenVector2d("v0", unfolded.row(0));
  printEigenVector2d("v1", unfolded.row(1));
  printEigenVector2d("v2", unfolded.row(2));
  printEigenVector2d("v3", unfolded.row(3));
  std::cout << "ordered quad vertex indices: " << vertex_indices[0]  << vertex_indices[1]  << vertex_indices[2]  << vertex_indices[3] << std::endl;

  iSData.inputMesh->flip(edge, false); // TODO: remove this line when done testing
  bool couldFlip = iSData.intrinsicMesh->flip(edge, false);
  if (!couldFlip) { return false; }

  // update edge length
  iSData.L(edge.getIndex()) = flipped_edgelength(unfolded); // TODO: check if flipping messes with the edge index
  // std::cout << "(flip_intrinsic) Edge " << edge.getIndex() << " length: " << iSData.L[edge.getIndex()] << std::endl;
  // add intrinsic flip into mapping queue
  iSData.mapping.push_back(std::make_unique<Edge_Flip>(Fk.getIndex(), Fl.getIndex(), unfolded, Fk_unique.first, Fl_unique.first));
  // TODO: the mapping also needs to swap the barycentric coordinate indices during transformation accordingly

  // -------------------------------------------------------------v
  t = true;
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
  printGCSFace("New Fk", Fk);
  printGCSFace("New Fl", Fl);

  std::cout << "Predicted Fk = {" << vertex_indices[2] << ", " << vertex_indices[3] << ", " << vertex_indices[1] << "}" << std::endl;
  std::cout << "Predicted Fl = {" << vertex_indices[3] << ", " << vertex_indices[2] << ", " << vertex_indices[0] << "}" << std::endl;
  // ------------------------------------------------------------------------^

  return true;
}

// index 0: shared edge ij
// index 1: first edge ik of ijk-triangle
// index 2: second edge jk of ijk-triangle
// index 3: first edge il of ijl-triangle
// index 4: second edge jl of ijl-triangle
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
      // if ((E.firstVertex().getIndex() == edge.firstVertex().getIndex()
      // && E.secondVertex().getIndex() != edge.secondVertex().getIndex())
      // || (E.secondVertex().getIndex() == edge.firstVertex().getIndex()
      // && E.firstVertex().getIndex() != edge.secondVertex().getIndex()))
      // {
      //   indices[offset] = E.getIndex();
      // }
      // else if ((E.firstVertex().getIndex() == edge.secondVertex().getIndex()
      // && E.secondVertex().getIndex() != edge.firstVertex().getIndex())
      // || (E.secondVertex().getIndex() == edge.secondVertex().getIndex()
      // && E.firstVertex().getIndex() != edge.firstVertex().getIndex()))
      // {
      //   indices[offset+1] = E.getIndex();
      // }
    }
    // offset += 2;
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