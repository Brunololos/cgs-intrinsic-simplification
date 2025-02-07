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
  data.tangent_space_reference_edges = Eigen::MatrixXi::Zero(V.rows(), 1);

  // initialize intrinsic edge lengths
  int nEdges = data.intrinsicMesh->nEdges();
  data.L = Eigen::VectorXd::Zero(nEdges);
  int i = 0;
  // std::cout << "Calculating L:\n" << std::endl; // TODO: remove
  for (gcs::Edge e : data.inputMesh->edges())
  {
    data.L(e.getIndex()) = data.inputGeometry->edgeLength(e);
    // /* TODO: remove */ std::cout << "edge(" << e.firstVertex().getIndex() << ", " << e.secondVertex().getIndex() << ") L(" << i << ") = " << data.L(i) << std::endl;
    // if(data.L(e.getIndex()) == 0.0) {
    //   std::cout << RED << "edge(" << e.firstVertex().getIndex() << ", " << e.secondVertex().getIndex() << ") L(" << i << ") = " << data.L(i) << RESET << std::endl;
    // }
    i++;
  }

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

void iSimp_step(iSimpData& data)
{
  std::vector<gcs::Vertex> temp_neighbors = std::vector<gcs::Vertex>();
  int vertex_idx;
  double ice;
  bool could_flip, could_flatten, could_remove, could_flip_to_deg3;
  // NOTE: Performing random intrinsic edge-flips
  // vertex_idx = rand() % data.inputMesh->nEdges();
  // vertex_idx = 1;
  // // for (gcs::Edge E : data.intrinsicMesh->edges())
  // TODO: just add this check into can_be_flipped() method (don't allow boundary flips)
  // while (true)
  // {
  //   vertex_idx = rand() % data.intrinsicMesh->nEdges();
  //   gcs::Edge E = data.intrinsicMesh->edge(r);
  //   int nk = 0;
  //   for (gcs::Face F : E.adjacentFaces())
  //   {
  //     nk++;
  //   }
  //   if (nk == 2)
  //   {
  //     vertex_idx = E.getIndex();
  //     break;
  //   }
  // }
  // std::cout << "flip edge: " << vertex_idx << " => (" << data.intrinsicMesh->edge(r).firstVertex().getIndex() << ", " << data.intrinsicMesh->edge(r).secondVertex().getIndex() << ")" << std::endl;
  // could_flip = flip_intrinsic(data, data.intrinsicMesh->edge(r));
  // if (could_flip)
  // {
  //   std::cout << "Edge flip successful!" << std::endl;
  // }
  // if (!could_flip)
  // {
  //   std::cout << "Edge flip failed!" << std::endl;
  // }

  // NOTE: Performing random intrinsic vertex flattenings
  // vertex_idx = rand() % data.intrinsicMesh->nVertices();
  // while (true)
  // {
  //   vertex_idx = rand() % data.intrinsicMesh->nVertices();
  //   gcs::Vertex V = data.intrinsicMesh->vertex(r);
  //   if (!V.isBoundary()) { break; }
  // }
  // // vertex_idx = 115;
  // std::cout << "flatten vertex: " << vertex_idx << std::endl;
  // could_flatten = flatten_vertex(data, vertex_idx);
  // if(could_flatten) { std::cout << "Vertex Flattening successful!" << std::endl; }
  // if(!could_flatten) { std::cout << "Vertex Flattening failed!" << std::endl; }

  // vertex_idx = rand() % V.rows(); // data.intrinsicMesh->nVertices();
  // i = 0;
  // while (i < 100)
  // {
  //   vertex_idx = rand() % V.rows(); // data.intrinsicMesh->nVertices();
  //   gcs::Vertex V = data.intrinsicMesh->vertex(r);
  //   // std::cout << "Checking vertex: " << vertex_idx << std::endl;
  //   if (!V.isDead() && !V.isBoundary())
  //   {
  //     // vertex_idx = 115;
  //     std::cout << "\n\n\nflatten vertex: " << vertex_idx << std::endl;
  //     could_flatten = flatten_vertex(data, vertex_idx);
  //     if(!could_flatten) { std::cout << "Vertex Flattening failed!" << std::endl; break; }
  //     std::cout << "\n\n\nflip vertex: " << vertex_idx << " to degree 3: " << std::endl;
  //     could_flip_to_deg3 = flip_vertex_to_deg3(data, r);
  //     if(!could_flip_to_deg3) { std::cout << dye("Vertex Flipping to degree 3 failed!", RED) << std::endl; break; }
  //     std::cout << "\n\n\nremove vertex: " << vertex_idx << std::endl;
  //     could_remove = remove_vertex(data, r);
  //     if(could_remove) { std::cout << "Vertex Removal successful!" << std::endl; }
  //     if(!could_remove) { std::cout << "Vertex Removal failed!" << std::endl; }
  //     break;
  //   }
  //   i++;
  // }

  if(data.Q.empty()) { data.hasConverged = true; return; }
  std::pair<int, double>* current = data.Q.begin()->get();   // this is Q.top()
  data.Q.erase(data.Q.begin());                              // this is Q.pop()
  vertex_idx = current->first;
  ice = current->second;

  if(std::isinf(ice)) { std::cout << "Converged!" << std::endl; data.hasConverged = true; return; }
  if (data.intrinsicMesh->vertex(vertex_idx).isDead()) { return; }
  // if (data.intrinsicMesh->vertex(r).isBoundary()) { return true; }
  // std::cout << "\n\n\nflatten vertex: " << vertex_idx << std::endl;
  could_flatten = flatten_vertex(data, vertex_idx);
  if(could_flatten)
  {
    // std::cout << "\n\n\nflip vertex: " << vertex_idx << " to degree 3: " << std::endl;
    could_flip_to_deg3 = flip_vertex_to_deg3(data, vertex_idx);
    if(could_flip_to_deg3)
    {
      // std::cout << "\n\n\nremove vertex: " << vertex_idx << std::endl;
      // save neighbors, because we cant iterate over vertex neighborhood after vertex removal
      temp_neighbors.clear();
      // std::cout << "cleared temp_neighbors: ";
      // for (gcs::Vertex neighbor : temp_neighbors) { std::cout << neighbor.getIndex() << ", "; }
      // std::cout << std::endl;
      for (gcs::Vertex neighbor : data.intrinsicMesh->vertex(vertex_idx).adjacentVertices()) { /* std::cout << "inserting temp_neighbor: " << neighbor.getIndex() << std::endl; */ temp_neighbors.push_back(neighbor); }
      could_remove = remove_vertex(data, vertex_idx);
      // if(could_remove) { std::cout << "Vertex Removal successful!" << std::endl; }
      // else { std::cout << "Vertex Removal failed!" << std::endl; }
      // try to enforce delaunay by iterating over all edges and flipping them, if necessary
      // NOTE: one could try to be smarter and only iterate over edges incident to the changed vertices/edges
      flip_to_delaunay(data);
    }
    else { std::cout << dye("Vertex Flipping to degree 3 failed!", RED) << std::endl; }

  }
  else { std::cout << "Vertex Flattening failed!" << std::endl; }

  // TODO: on success: update values of neighbor ice's
  // on failure: reinsert vertex with infinite value on failure
  if (could_flatten && could_flip_to_deg3 && could_remove)
  {
    for (gcs::Vertex neighbor : temp_neighbors)
    {
      // because we are working with delta complexes, we can have self edges and have to check not to reinsert already removed vertices
      if (neighbor.isDead()) { continue; }
      if (neighbor.getIndex() == vertex_idx) { continue; }
      double ice = intrinsic_curvature_error(data, neighbor);
      // data.Q.erase(data.Q.find(std::make_shared<std::pair<int, double>>(data.Q_elems[neighbor.getIndex()]))); // erase the element to be sure the order gets updated
      data.Q.erase(std::make_shared<std::pair<int, double>>(data.Q_elems[neighbor.getIndex()])); // erase the element to be sure the order gets updated
      data.Q_elems[neighbor.getIndex()].second = ice;
      data.Q.insert(std::make_shared<std::pair<int, double>>(data.Q_elems[neighbor.getIndex()]));
    }
  } else {
    data.Q_elems[vertex_idx].second = INFINITY;
    data.Q.insert(std::make_shared<std::pair<int, double>>(data.Q_elems[vertex_idx]));
  }
}