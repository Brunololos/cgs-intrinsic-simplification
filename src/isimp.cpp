#include "isimp.hpp"

#include <iostream>
#include <fstream>

#include <igl/parallel_for.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/polar_svd3x3.h>

#include <Eigen/Sparse>

typedef Eigen::SparseVector<double>::InnerIterator SVIter;

void initMesh(const Eigen::Matrix<double, -1, 3>& V, const Eigen::Matrix<int, -1, 3>& F, iSimpData& data){
  data.V = V;
  data.F = F;
  data.M = Eigen::MatrixXd::Zero(V.rows(), 2);
  data.T_minus = Eigen::MatrixXd::Zero(V.rows(), 2);
  data.T_plus = Eigen::MatrixXd::Zero(V.rows(), 2);
  data.Q = std::set<std::shared_ptr<std::pair<int, double>>, customLess>();
  data.Q_elems = std::vector<std::pair<int, double>>(V.rows());

  data.inputMesh.reset(new gcs::ManifoldSurfaceMesh(F));
  data.inputGeometry.reset(new gcs::VertexPositionGeometry(*data.inputMesh, V));

  // initialize intrinsic connectivity
  data.intrinsicMesh.reset(new gcs::ManifoldSurfaceMesh(F));
  data.recoveredMesh.reset(new gcs::ManifoldSurfaceMesh(F));
  data.tangent_space_reference_edges = Eigen::MatrixXi::Zero(V.rows(), 1);

  // initialize intrinsic edge lengths
  int nEdges = data.intrinsicMesh->nEdges();
  data.inputL = Eigen::VectorXd::Zero(nEdges);
  int i = 0;
  for (gcs::Edge e : data.inputMesh->edges())
  {
    data.inputL(e.getIndex()) = data.inputGeometry->edgeLength(e);
    if(data.inputL(e.getIndex()) == 0.0) {
      std::cout << RED << "edge(" << e.firstVertex().getIndex() << ", " << e.secondVertex().getIndex() << ") L(" << i << ") = " << data.inputL(i) << RESET << std::endl;
    }
    i++;
  }
  data.L = data.inputL;
  data.recovered_L = data.inputL;
  bool result = validate_intrinsic_edge_lengths(data);
  if (!result) { std::cout << "occurred in iSData initialization " << std::endl; }

  // init masses & error vectors
  for (gcs::Vertex v : data.intrinsicMesh->vertices())
  {
    // calculate curvature
    double curvature = gaussian_curvature(data, v);
    double K_minus = - std::min(curvature, 0.0);
    double K_plus = std::max(curvature, 0.0);
    data.M.row(v.getIndex()) = Eigen::Vector2d(K_minus, K_plus);

    for (gcs::Edge e : v.adjacentEdges())
    {
      data.tangent_space_reference_edges(v.getIndex(), 0) = e.getIndex();
      break;
    }
  }

  // insert vertices with intrinsic curvature error into priority queue
  for (gcs::Vertex v : data.intrinsicMesh->vertices())
  {
    double ice = intrinsic_curvature_error(data, v);
    std::pair<int, double> v_elem = std::pair<int, double>(v.getIndex(), ice);
    data.Q_elems[v.getIndex()] = v_elem;
    data.Q.insert(std::make_shared<std::pair<int, double>>(v_elem));
  }

  data.mapping = std::vector<std::unique_ptr<Mapping_operation>>();
  data.tracked_points = std::vector<BarycentricPoint>();
  data.tracked_texcoords = std::vector<TexCoord>();
  data.tracked_by_triangle = std::vector<std::vector<int>>(data.F.rows());
  data.hasConverged = false;
}

bool iSimp_step(iSimpData& data, const bool verbose)
{
  std::vector<gcs::Vertex> temp_neighbors = std::vector<gcs::Vertex>();
  int vertex_idx;
  double ice;
  bool could_flip, could_flatten, could_remove, could_flip_to_deg3;
  bool result;

  if(data.Q.empty()) { data.hasConverged = true; return false; }
  std::pair<int, double>* current = data.Q.begin()->get();   // this is Q.top()
  data.Q.erase(data.Q.begin());                              // this is Q.pop()
  vertex_idx = current->first;
  ice = current->second;

  if(std::isinf(ice)) { if(verbose) { std::cout << "Converged!" << std::endl; } data.hasConverged = true; return false; }
  if (data.intrinsicMesh->vertex(vertex_idx).isDead()) { return false; }
  could_flatten = flatten_vertex(data, vertex_idx, verbose);
  if(could_flatten)
  {
    could_flip_to_deg3 = flip_vertex_to_deg3(data, vertex_idx, verbose);
    if(could_flip_to_deg3)
    {
      // save neighbors, because we cant iterate over vertex neighborhood after vertex removal
      temp_neighbors.clear();
      for (gcs::Vertex neighbor : data.intrinsicMesh->vertex(vertex_idx).adjacentVertices()) { temp_neighbors.push_back(neighbor); }
      could_remove = remove_vertex(data, vertex_idx, verbose);
      // try to enforce delaunay by iterating over all edges and flipping them, if necessary
      // NOTE: one could try to be smarter and only iterate over edges incident to the changed vertices/edges
      // NOTE2: I tried to be smarter, but not iterating over all edges yielded worse results
      flip_to_delaunay(data, verbose);
    }
    else if(verbose) { std::cout << "Flipping Vertex " + std::to_string(vertex_idx) + " to degree 3 failed!" << std::endl; }

  }
  else if(verbose) { std::cout << "Vertex Flattening of Vertex " + std::to_string(vertex_idx) + " failed!" << std::endl; }

  // on success: update values of neighbor ice's
  // on failure: reinsert vertex with infinite cost on failure
  if (could_flatten && could_flip_to_deg3 && could_remove)
  {
    for (gcs::Vertex neighbor : temp_neighbors)
    {
      // because we are working with delta complexes, we can have self edges and have to check not to reinsert already removed vertices
      if (neighbor.isDead()) { continue; }
      if (neighbor.getIndex() == vertex_idx) { continue; }
      double ice = intrinsic_curvature_error(data, neighbor);
      data.Q.erase(std::make_shared<std::pair<int, double>>(data.Q_elems[neighbor.getIndex()])); // erase the element to be sure the order gets updated
      data.Q_elems[neighbor.getIndex()].second = ice;
      data.Q.insert(std::make_shared<std::pair<int, double>>(data.Q_elems[neighbor.getIndex()]));
    }
    return true;
  } else {
    data.Q_elems[vertex_idx].second = INFINITY;
    data.Q.insert(std::make_shared<std::pair<int, double>>(data.Q_elems[vertex_idx]));
    return false;
  }
}

void recover_intrinsics_at(const int mapping_idx, iSimpData& data)
{
  data.recoveredMesh.reset(new gcs::ManifoldSurfaceMesh(data.F));
  data.recovered_L = data.inputL;

  std::cout << "Recovering intrinsics:        ";
  progressBar bar;
  start_progressBar(bar);

  // replay simplification operations
  for(int i=0; i<=mapping_idx; i++)
  {
    switch(data.mapping[i]->op_type)
    {
      default:
      case SIMP_OP::E_FLIP:
        replay_intrinsic_flip(data, data.mapping[i]->reduced_primitive_idx());
        break;
      case SIMP_OP::V_FLATTEN:
        replay_vertex_flattening(data, data.mapping[i]->reduced_primitive_idx());
        break;
      case SIMP_OP::V_REMOVAL:
        replay_vertex_removal(data, data.mapping[i]->reduced_primitive_idx());
        break;
    }
    double progress_percent = ((double) i + 1) / ((double) mapping_idx + 1);
    progress_progressBar(bar, progress_percent);
  }
  finish_progressBar(bar);
}

void recover_heat_and_intrinsics_at(const int mapping_idx, heatDiffData& hdData, iSimpData& data)
{
  data.recoveredMesh.reset(new gcs::ManifoldSurfaceMesh(data.F));
  data.recovered_L = data.inputL;
  hdData.recovered_heat = hdData.initial_heat;

  std::cout << "Recovering heat & intrinsics: ";
  progressBar bar;
  start_progressBar(bar);

  // replay simplification operations
  for(int i=0; i<=mapping_idx; i++)
  {
    switch(data.mapping[i]->op_type)
    {
      default:
      case SIMP_OP::E_FLIP:
        replay_intrinsic_flip(data, data.mapping[i]->reduced_primitive_idx());
        break;
      case SIMP_OP::V_FLATTEN:
        replay_vertex_flattening_with_heat(data, hdData, data.mapping[i]->reduced_primitive_idx()); // TODO: pass the scaling factor
        // replay_vertex_flattening_with_heat(data, hdData, data.mapping[i]->reduced_primitive_idx(), data.mapping[i]->get_u()); // TODO: pass the scaling factor
        break;
      case SIMP_OP::V_REMOVAL:
        replay_vertex_removal_with_heat(data, hdData, data.mapping[i]->reduced_primitive_idx());
        break;
    }
    double progress_percent = ((double) i + 1) / ((double) mapping_idx + 1);
    progress_progressBar(bar, progress_percent);
  }
  finish_progressBar(bar);
}