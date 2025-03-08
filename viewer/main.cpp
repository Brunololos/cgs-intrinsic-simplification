#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/surface_parameterization_quantity.h"

// #define MIN_QUAD_WITH_FIXED_CPP_DEBUG
#include <stdlib.h>
// #include <igl/png/readPNG.h>
#include <igl/read_triangle_mesh.h>
#include <stdlib.h>
#include <chrono>
#include <ctime>

#include "isimp.hpp"
#include "heat_diff.hpp"
#include "screenshot.hpp"
#include "lscm_wrapper.hpp"
#include "query_texture_barycentrics.hpp"

typedef Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> MatrixXuc;

// TODO: move this somewhere else
enum class ISIMP_MODE
{
  ON_INPUT,
  EACH_FRAME
};

int main(int argc, char *argv[])
{
  std::string filename = "../../meshes/spot.obj";
  int coarsen_to_n_vertices;
  int min_n_vertices = 0;
  bool using_texture = true;

  // int dot_pos = ((std::string) argv[1]).find_last_of('.');
  int dot_pos = ((std::string) argv[1]).size() - 4;
  if (argc > 1 && ((std::string) argv[1]).substr(dot_pos).compare(".obj"))
  {
    std::cout << "Please pass a .obj mesh as the first argument" << std::endl;
    exit(-1);
  }
  if (argc > 1) { filename = argv[1]; }
  // if (argc > 2) { coarsen_to_n_vertices = std::to_integer(string(argv[2])); }

  Eigen::MatrixXd V, UV, NV;
  Eigen::MatrixXi F, UF, NF;
  // store all colorings for all different simplification states (for each vertex removal) (This is pretty inefficient.)
  std::vector<Eigen::VectorXi> colorings;
  // points tracked by triangle for each simplification state (I did this just to make the polyscope slider movement smoother, but this takes a lot of storage.)
  std::vector<std::vector<std::vector<int>>> mapped_by_triangles;
  std::vector<std::vector<BarycentricPoint>> mapped_pointss;
  std::vector<int> snapshot_mapping_indices;
  std::vector<int> simp_step_mapping_indices;
  MatrixXuc R, G, B, A;
  int TEXTURE_WIDTH = 2000;
  int TEXTURE_HEIGHT = 2000;
  std::vector<unsigned char> polyscope_texture = std::vector<unsigned char>();
  polyscope::SurfaceCornerParameterizationQuantity* q = nullptr;
  polyscope::SurfaceMesh* psCoarseMesh;
  polyscope::SurfaceMesh* psMesh;

  igl::readOBJ(filename, V, UV, NV, F, UF, NF);

  if (UF.rows() <= 0)
  {
    std::cout << "Please pass a .obj mesh with texture coordinates included." << std::endl;
    exit(-1);
  }

  coarsen_to_n_vertices = V.rows();
  colorings = std::vector<Eigen::VectorXi>();
  mapped_by_triangles = std::vector<std::vector<std::vector<int>>>();
  mapped_pointss = std::vector<std::vector<BarycentricPoint>>();
  snapshot_mapping_indices = std::vector<int>();
  simp_step_mapping_indices = std::vector<int>();
  int snapshot_interval = 50;
  bool do_texture_update = false;

  // modified vertices & faces
  Eigen::MatrixXd U = V;
  Eigen::MatrixXi H = F;

  // settings
  int mode = 0;
  bool virgin = true;
  bool view_constraints = true;
  bool view_all_points = false;
  bool optimise = false;

  // ISIMP (intrinsic-simplification)
  ISIMP_MODE simpMode = ISIMP_MODE::ON_INPUT;

  iSimpData iSData;
  initMesh(V, F, iSData);

  // randomize seed
  srand(time(NULL));

  // TODO:
  std::ofstream logFile;                     // logfile
  bool doLogging = false;

  R = 255 * MatrixXuc::Ones(TEXTURE_WIDTH, TEXTURE_HEIGHT);
  G = 255 * MatrixXuc::Ones(TEXTURE_WIDTH, TEXTURE_HEIGHT);
  B = 255 * MatrixXuc::Ones(TEXTURE_WIDTH, TEXTURE_HEIGHT);
  A = 255 * MatrixXuc::Ones(TEXTURE_WIDTH, TEXTURE_HEIGHT);

  Eigen::MatrixXd CM = (Eigen::MatrixXd(16, 3) <<
                            // 0.8,0.8,0.8,                        // 0 gray
                            //  1,1,1,
                            //  1.0/255.0,31.0/255.0,255.0/255.0,   // 1 blue
                            //  0/255.0,201.0/255.0,195.0/255.0,    // 2 cyan
                            //  7.0/255.0,152.0/255.0,59.0/255.0,   // 3 green
                            //  190.0/255.0,255.0/255.0,0.0/255.0,  // 4 lime
                            //  255.0/255.0,243.0/255.0,0.0/255.0,  // 5 yellow
                            //  245.0/255.0,152.0/255.0,0.0/255.0,  // 6 orange
                            //  224.0/255.0,18.0/255.0,0.0/255.0,   // 7 red
                            //  220.0/255.0,13.0/255.0,255.0/255.0  // 8 magenta
                        1, 1, 1,                                    // 1 white
                        178.0 / 255.0, 31.0 / 255.0, 53.0 / 255.0,  // 1 dark red
                        216.0 / 255.0, 39.0 / 255.0, 53.0 / 255.0,  // 2 light red
                        255.0 / 255.0, 116.0 / 255.0, 53.0 / 255.0, // 3 orange
                        255.0 / 255.0, 161.0 / 255.0, 53.0 / 255.0, // 4 yellowish orange
                        255.0 / 255.0, 203.0 / 255.0, 53.0 / 255.0, // 5 redish yellow
                        255.0 / 255.0, 249.0 / 255.0, 53.0 / 255.0, // 6 yellow
                        0.0 / 255.0, 117.0 / 255.0, 58.0 / 255.0,   // 7 dark green
                        0.0 / 255.0, 158.0 / 255.0, 71.0 / 255.0,   // 8 green
                        22.0 / 255.0, 221.0 / 255.0, 53.0 / 255.0,  // 9 light green
                        0.0 / 255.0, 82.0 / 255.0, 165.0 / 255.0,   // 10 dark blue
                        0.0 / 255.0, 121.0 / 255.0, 231.0 / 255.0,  // 11 blue
                        0.0 / 255.0, 169.0 / 255.0, 252.0 / 255.0,  // 12 light blue
                        104.0 / 255.0, 30.0 / 255.0, 126.0 / 255.0, // 13 dark purple
                        125.0 / 255.0, 60.0 / 255.0, 181.0 / 255.0, // 14 purple
                        189.0 / 255.0, 122.0 / 255.0, 246.0 / 255.0 // 15 light purple
                        )
                           .finished();

  Eigen::MatrixXd CM2 = (Eigen::MatrixXd(12, 3) <<
    166.0 / 255.0, 206.0 / 255.0, 227.0 / 255.0,
    31.0 / 255.0, 120.0 / 255.0, 180.0 / 255.0,
    178.0 / 255.0, 223.0 / 255.0, 138.0 / 255.0,
    51.0 / 255.0, 160.0 / 255.0, 44.0 / 255.0,
    251.0 / 255.0, 154.0 / 255.0, 153.0 / 255.0,
    227.0 / 255.0, 26.0 / 255.0, 28.0 / 255.0,
    253.0 / 255.0, 191.0 / 255.0, 111.0 / 255.0,
    255.0 / 255.0, 127.0 / 255.0, 0.0 / 255.0,
    202.0 / 255.0, 178.0 / 255.0, 214.0 / 255.0,
    106.0 / 255.0, 61.0 / 255.0, 154.0 / 255.0,
    255.0 / 255.0, 255.0 / 255.0, 153.0 / 255.0,
    177.0 / 255.0, 89.0 / 255.0, 40.0 / 255.0
  ).finished();

  // coloring
  int numVorF = (true /* true: visualize over vertices, false: visualize over faces */) ? V.rows() : F.rows();
  // TODO: C can be removed, we now render colors differently
  Eigen::MatrixXd C = Eigen::MatrixXd(numVorF, 3);
  for (int i = 0; i < numVorF; i++)
  {
    C.row(i) = CM.row((i % (CM.rows() - 1)) + 1);
  }

  // NOTE: this approach assumes, that all triangles lie completely within the texture. (No triangles wrap around the outside.)
  // I note this here, because lscm sometimes yields uv-coordinates that exceed the [0, 1] bounds. (They can be negative and their absolute value can exceed 1.)
  // To fix this I translated & rescaled all uv's to be withing [0, 1]. Because we don't have a set texture and build it ourselves, this is fine.
  auto register_texels = [&]()
  {
    // std::cout << "Registering texels: " << std::endl;
    for (gcs::Face F : iSData.inputMesh->faces())
    {
      // std::cout << "-> for extrinsic face: " << F << std::endl;
      std::array<int, 3> verts = std::array<int, 3>();
      int i = 0;
      for (gcs::Vertex V : F.adjacentVertices())
      {
        verts[i] = V.getIndex();
        i++;
      }

      // mapping from uv- to texture space
      Eigen::Matrix2d TexExt;
      TexExt << TEXTURE_WIDTH, 0,
          0, TEXTURE_HEIGHT;

      // triangle corners in texture space
      Point2D t0 = UV.row(verts[0]) * TexExt;
      Point2D t1 = UV.row(verts[1]) * TexExt;
      Point2D t2 = UV.row(verts[2]) * TexExt;

      // triangle bounding box
      double minx = std::min(t0[0], std::min(t1[0], t2[0]));
      double maxx = std::max(t0[0], std::max(t1[0], t2[0]));
      double miny = std::min(t0[1], std::min(t1[1], t2[1]));
      double maxy = std::max(t0[1], std::max(t1[1], t2[1]));

      // iterate over triangle bounding box & register points contained in the triangle
      for (int i = std::floor(minx); i <= std::floor(maxx); i++)
      {
        for (int j = std::floor(miny); j <= std::floor(maxy); j++)
        {
          if (lies_inside_triangle(i, j, t0, t1, t2))
          {
            // std::cout << to_str(Point2D(i, j)) << std::endl;
            BarycentricPoint bp = to_barycentric(Point2D(i, j), t0, t1, t2);
            // std::cout << "Expressing texture coordinate: (" << i << ", " << j << ") barycentrically of: A=" << to_str(t0) << ", B=" << to_str(t1) << ", C=" << to_str(t2) << std::endl;
            // std::cout << "Created barycentric point: " << to_str(bp) << std::endl;
            // std::cout << "Converting back into explicit form 012, we get: " << to_str(to_explicit(bp, t0, t1, t2)) << std::endl;
            // std::cout << "Converting back into explicit form 021, we get: " << to_str(to_explicit(bp, t0, t2, t1)) << std::endl;
            // std::cout << "Converting back into explicit form 102, we get: " << to_str(to_explicit(bp, t1, t0, t2)) << std::endl;
            // std::cout << "Converting back into explicit form 120, we get: " << to_str(to_explicit(bp, t1, t2, t0)) << std::endl;
            // std::cout << "Converting back into explicit form 201, we get: " << to_str(to_explicit(bp, t2, t0, t1)) << std::endl;
            // std::cout << "Converting back into explicit form 210, we get: " << to_str(to_explicit(bp, t2, t1, t0)) << std::endl;
            TexCoord tc = TexCoord(i, j);
            register_point(iSData, bp, tc, F.getIndex());
          }
        }
      }
    }
  };

  auto init_viewer_data = [&]()
  {
    // polyscope_texture = std::vector<unsigned char>(4 * TEXTURE_WIDTH * TEXTURE_HEIGHT);
    polyscope_texture.resize(4 * TEXTURE_WIDTH * TEXTURE_HEIGHT);
  //   // prev face-based coloring: viewer.data().set_face_based(true);
  //   // viewer.data().set_texture(R,G,B,A);
  //   // viewer.data().set_colors(C);

  //   // Calculate texture mapping
  //   viewer.data().set_uv(UV);
  //   viewer.data().set_texture(R, G, B);
  //   viewer.data().show_texture = true;

  //   // viewer.data().use_matcap = true;
  //   viewer.data().point_size = 15;
  //   viewer.data().line_width = 1;
  //   viewer.data().show_lines = F.rows() < 20000;
  };

  // Selection
  bool translate = false;
  RAXIS rot_axis = RAXIS::YAW;
  bool only_visible = false;
  // igl::opengl::glfw::imgui::SelectionWidget selection;
  Eigen::Array<double, Eigen::Dynamic, 1> W = Eigen::Array<double, Eigen::Dynamic, 1>::Zero(V.rows());
  Eigen::Array<double, Eigen::Dynamic, 1> and_visible = Eigen::Array<double, Eigen::Dynamic, 1>::Zero(V.rows());

  Eigen::MatrixXd all_constraints = V;

  Eigen::Array<int, -1, 1> constraints_mask = Eigen::Array<int, -1, 1>::Zero(V.rows());
  Eigen::Array<int, -1, 1> selected_mask = Eigen::Array<int, -1, 1>::Zero(V.rows());
  Eigen::VectorXi selected_ids;

  Eigen::MatrixXd constraints;
  Eigen::VectorXi bi;
  Eigen::VectorXi indices = Eigen::VectorXd::LinSpaced(V.rows(), 0.0, V.rows() - 1).cast<int>();

  // igl::AABB<Eigen::MatrixXd, 3> tree;
  // tree.init(V, F);

  const auto update_edges = [&]()
  {
    // int UR = U.rows();
    // int UC = U.cols();

    // Eigen::MatrixXd grad_offset = U - gradients;
    // Eigen::MatrixXi E = Eigen::MatrixXi(UR, 2);
    // for (int i = 0; i<UR; i++)
    // {
    //   E(i, 0) = i;
    //   E(i, 1) = UR + i;
    // }
    // // vertical stacking
    // Eigen::MatrixXd P(U.rows()+grad_offset.rows(), U.cols());
    // P << U, grad_offset;
    // //viewer.data().line_width = 2.0; // NOTE: SUPPOSEDLY NOT SUPPORTED ON MAC & WINDOWS
    // viewer.data().set_edges(P, E, CM.row(4));
  };
  auto compute_coloring = [&]() {
    // update coloring
    int num_colors = CM.rows() - 1;
    Eigen::VectorXi Coloring = Eigen::VectorXi::Zero(iSData.F.rows());
    for (int i = 0; i < iSData.inputMesh->nFaces(); i++)
    {
      // Determine color
      int color_idx = (i % (num_colors - 1)) + 1;
      gcs::Face current_face = iSData.intrinsicMesh->face(i);
      if (current_face.isDead()) { continue; }
      int max_color_idx = 0;
      for(gcs::Face neighbor : current_face.adjacentFaces())
      {
        int neighbor_idx = neighbor.getIndex();
        // std::cout << "Checking neighbor face: " << neighbor_idx << std::endl;
        int neighbor_color_idx = Coloring(neighbor_idx);
        if (neighbor_color_idx > max_color_idx) { max_color_idx = neighbor_color_idx; }
        if (color_idx == neighbor_color_idx)
        {
          color_idx = max_color_idx + 1;
        }
      }
      if (color_idx >= num_colors) { color_idx = 0; }
      Coloring(i) = color_idx;
      // std::cout << "Setting coloring index: " << color_idx << std::endl;
    }
    colorings.push_back(Coloring);
  };


  const auto &update_texture = [&]()
  {
    int num_triangles = iSData.mapped_by_triangle.size();
    int num_colors = CM.rows() - 1;
    int n_removals = iSData.inputMesh->nVertices() - coarsen_to_n_vertices;

    // draw texture
    for (int i = 0; i < num_triangles; i++)
    {
      for (int j = 0; j < iSData.mapped_by_triangle[i].size(); j++)
      {
        int original_idx = iSData.mapped_by_triangle[i][j];
        // BarycentricPoint original_point = iSData.tracked_points[original_idx];
        // Point2D tex_coord = to_explicit(original_point, );
        TexCoord texcoord = iSData.tracked_texcoords[original_idx];
        // R(texcoord[0], texcoord[1]) = CM((i % num_colors) + 1, 0) * 255;
        // G(texcoord[0], texcoord[1]) = CM((i % num_colors) + 1, 1) * 255;
        // B(texcoord[0], texcoord[1]) = CM((i % num_colors) + 1, 2) * 255;
        polyscope_texture[4*(texcoord[0] + TEXTURE_WIDTH*texcoord[1]) + 0] = CM(colorings[n_removals](i), 0) * 255;
        polyscope_texture[4*(texcoord[0] + TEXTURE_WIDTH*texcoord[1]) + 1] = CM(colorings[n_removals](i), 1) * 255;
        polyscope_texture[4*(texcoord[0] + TEXTURE_WIDTH*texcoord[1]) + 2] = CM(colorings[n_removals](i), 2) * 255;
        // polyscope_texture[4*(texcoord[0] + TEXTURE_WIDTH*texcoord[1]) + 0] = CM(i % num_colors, 0) * 255;
        // polyscope_texture[4*(texcoord[0] + TEXTURE_WIDTH*texcoord[1]) + 1] = CM(i % num_colors, 1) * 255;
        // polyscope_texture[4*(texcoord[0] + TEXTURE_WIDTH*texcoord[1]) + 2] = CM(i % num_colors, 2) * 255;
        polyscope_texture[4*(texcoord[0] + TEXTURE_WIDTH*texcoord[1]) + 3] = 255;
        // std::cout << "Setting texcoord: (" << texcoord[0] << ", " << texcoord[1] << "): R=" << CM(i%9, 0) << ", G=" << CM(i%9, 1) << ", B=" << CM(i%9, 2) << std::endl;
      }
    }
    // viewer.data().set_texture(R, G, B);
    // std::cout << "write_success: " << igl::png::writePNG(R, G, B, A, "H:/GIT/cgs-intrinsic-simplification/textures/object_texture.png") << std::endl;
  };

  const auto update_mesh_points = [&](Eigen::VectorXi &selected_indices, Eigen::MatrixXd &selected_constraints)
  {
    for (int i = 0; i < selected_indices.size(); i++)
    {
      int index = selected_indices(i);
      float c0 = selected_constraints(i, 0);
      float c1 = selected_constraints(i, 1);
      float c2 = selected_constraints(i, 2);
      if (c0 != c0 || c1 != c1 || c2 != c2)
      {
        std::cout << "found nan value: " << constraints(index) << std::endl;
      }
      U(index, 0) = c0;
      U(index, 1) = c1;
      U(index, 2) = c2;
      // std::cout << "Set mesh point " << i << std::endl;
    }
  };

  // TODO: maybe make selected points translate with vertices during optimisation
  const auto update_points = [&]()
  {
    // constrained
    // viewer.data().clear_points();
    // if (view_constraints || view_all_points)
    // {
    //   viewer.data().set_points(igl::slice_mask(all_constraints, (selected_mask * constraints_mask).cast<bool>(), 1), CM.row(6));
    //   viewer.data().add_points(igl::slice_mask(/* all_constraints */ U, (selected_mask * (1 - constraints_mask)).cast<bool>(), 1), CM.row(5));
    //   viewer.data().add_points(igl::slice_mask(all_constraints, ((1 - selected_mask) * constraints_mask).cast<bool>(), 1), CM.row(7));
    // }
    // if (view_all_points)
    // {
    //   viewer.data().add_points(igl::slice_mask(/* all_constraints */ U, ((1 - selected_mask) * (1 - constraints_mask)).cast<bool>(), 1), CM.row(2));
    // }
  };

  // update entire mesh
  auto refresh_coarse_mesh = [&]()
  {
    if (psCoarseMesh != nullptr)
    {
      // polyscope::unre
    }
    // TODO: update connectivity
    // H = iSData.intrinsicMesh->faces();
    psCoarseMesh  = polyscope::registerSurfaceMesh("coarse mesh", U, H);
    // viewer.data().clear();
    // viewer.data().set_mesh(U, H);
    init_viewer_data();
    update_texture();
    // update_points();
    // viewer.draw();
  };

  // update mesh vertices
  auto refresh_mesh_vertices = [&]()
  {
    // viewer.data().set_vertices(U);
    init_viewer_data();
    update_texture();
    // update_points();
    // viewer.draw();
  };

  auto callback = [&]() {
    if(ImGui::SliderInt("remaining vertices", &coarsen_to_n_vertices, min_n_vertices, iSData.inputMesh->nVertices()))
    {
      do_texture_update = true;
    }
    if (ImGui::Button("Remove Vertex"))
    {
        coarsen_to_n_vertices = std::max(std::min(coarsen_to_n_vertices - 1, (int) iSData.inputMesh->nVertices()), min_n_vertices);
        do_texture_update = true;
    }
    if (ImGui::Button("Add Vertex"))
    {
        coarsen_to_n_vertices = std::max(std::min(coarsen_to_n_vertices + 1, (int) iSData.inputMesh->nVertices()), min_n_vertices);
        do_texture_update = true;
    }
    if (do_texture_update)
    {
      do_texture_update = false;
      int snapshot_idx = std::floor((iSData.inputMesh->nVertices() - coarsen_to_n_vertices) / snapshot_interval);
      int coarsen_from_n_vertices = iSData.inputMesh->nVertices() - snapshot_idx * snapshot_interval;
      // std::cout << "Coarsening from " << coarsen_from_n_vertices << std::endl;
      // std::cout << "Loading snapshots with index: " << snapshot_idx << std::endl;
      // std::cout << "mapped_by_triangles size: " << mapped_by_triangles.size() << std::endl;
      // std::cout << "mapped_pointss size: " << mapped_pointss.size() << std::endl;
      // TODO: commented snapshots
      iSData.mapped_by_triangle = mapped_by_triangles.at(snapshot_idx);
      iSData.mapped_points = mapped_pointss.at(snapshot_idx);
      int snapshot_mapping_idx = snapshot_mapping_indices.at(snapshot_idx);
      // std::cout << "Coarsening from: " << coarsen_from_n_vertices << " to " << coarsen_to_n_vertices << std::endl;
      // map_current_from_to(iSData, coarsen_from_n_vertices, coarsen_to_n_vertices);
      // if (snapshot_mapping_indices.size() > snapshot_idx + 1)
      if (simp_step_mapping_indices.size() > iSData.inputMesh->nVertices() - coarsen_to_n_vertices + 1)
      {
        // int next_snapshot_mapping_idx = snapshot_mapping_indices.at(snapshot_idx + 1);
        int target_removal_idx = iSData.inputMesh->nVertices() - coarsen_to_n_vertices;
        int simp_step_mapping_idx = simp_step_mapping_indices.at(iSData.inputMesh->nVertices() - coarsen_to_n_vertices + 1);
        map_current_from_to(iSData, snapshot_mapping_idx, simp_step_mapping_idx);
        // std::cout << "Mapping from " << snapshot_mapping_idx << " to " << simp_step_mapping_idx << std::endl;
        // std::cout << "target_removal_idx: " << target_removal_idx << std::endl;
      } else {
        map_current_from(iSData, snapshot_mapping_idx);
        // std::cout << "Mapping from " << snapshot_mapping_idx << " to the end." << std::endl;
      }
      // std::cout << "Mapping size: " << iSData.mapping.size() << std::endl;
      // std::cout << "Step mapping size: " << simp_step_mapping_indices.size() << std::endl;
      // std::cout << "Last step mapping idx: " << simp_step_mapping_indices.at(simp_step_mapping_indices.size() - 1) << std::endl;
      // map_registered(iSData, coarsen_to_n_vertices);
      // std::cout << iSData.mapped_points.size() << std::endl;
      update_texture();
      if (q != nullptr) {
        q->setTexture(TEXTURE_WIDTH, TEXTURE_HEIGHT, polyscope_texture, polyscope::TextureFormat::RGBA8);
        q->setEnabled(true);
        q->setStyle(polyscope::ParamVizStyle::TEXTURE);
        q->setCheckerSize(1);
      }
    }
  };

  // register_texels();
  query_texture_barycentrics(UV, UF, TEXTURE_WIDTH, TEXTURE_HEIGHT, iSData);
  // // These are inserted:
  map_registered(iSData);
  mapped_by_triangles.push_back(iSData.mapped_by_triangle);
  mapped_pointss.push_back(iSData.mapped_points);
  snapshot_mapping_indices.push_back(0);
  simp_step_mapping_indices.push_back(0);
  // update_texture();
  // std::cout << "write_success: " << igl::png::writePNG(R, G, B, A, "H:/GIT/cgs-intrinsic-simplification/textures/object_texture.png") << std::endl;
  // // end insertion
  init_viewer_data();
  compute_coloring();
  update_texture();
  // viewer.core().is_animating = true;
  // viewer.core().background_color.head(3) = CM.row(0).head(3).cast<float>();
  // // update_points()
  // viewer.launch();

  // TODO: this is neat, but for the intrinsic case, still the custom edge lengths need to be used.
  // iSData.inputGeometry->requireCotanLaplacian();
  // iSData.inputGeometry->requireVertexGalerkinMassMatrix();
  // Eigen::SparseMatrix<double> L = iSData.inputGeometry->cotanLaplacian;
  // Eigen::SparseMatrix<double> M = iSData.inputGeometry->vertexGalerkinMassMatrix;

  std::cout << "\n\nBeginning Simplification..." << std::endl;
  // while (!iSData.hasConverged && iSData.intrinsicMesh->nVertices() > coarsen_to_n_vertices)
  while (!iSData.hasConverged)
  {
    int step_mapping_idx = iSData.mapping.size();
    // make simplification step
    bool success = iSimp_step(iSData);
    // if case, because the iSimp step can fail & not introduce any new simplification operations
    // if(iSData.mapping.size() > step_mapping_idx) { simp_step_mapping_indices.push_back(step_mapping_idx); }
    if(success) { simp_step_mapping_indices.push_back(step_mapping_idx); }
    // TODO: commented snapshots
    // map_registered(iSData); // TODO: This should be done step-wise not remapping through the whole sequence
    // // map_current_from(iSData, iSData.intrinsicMesh->nVertices() + 1);
    // // std::cout << "Mapping from " << iSData.intrinsicMesh->nVertices() << std::endl;
    if ((iSData.inputMesh->nVertices() - iSData.intrinsicMesh->nVertices()) % snapshot_interval == 0) {
      map_current_from(iSData, snapshot_mapping_indices.at(snapshot_mapping_indices.size() - 1));
      std::cout << "Took snapshot at " << iSData.inputMesh->nVertices() - iSData.intrinsicMesh->nVertices() << " removed vertices" << std::endl;
      mapped_by_triangles.push_back(iSData.mapped_by_triangle);
      mapped_pointss.push_back(iSData.mapped_points);
      snapshot_mapping_indices.push_back(iSData.mapping.size());
    } // inefficiently store all simplification mappings
    compute_coloring();
    // TODO: remove after debugging
    // bool are_edge_lengths_valid = validate_intrinsic_edge_lengths(iSData);
    // if(!are_edge_lengths_valid) {
    //   std::cout << RED << "Found invalid edge-lengths after step_mapping_idx: " << step_mapping_idx << RESET << std::endl;
    //   iSData.hasConverged = true;
    // }
  }
  // refresh_coarse_mesh(); //TODO: update coarse mesh
  min_n_vertices = iSData.intrinsicMesh->nVertices();
  // min_n_vertices = 200; // TODO: remove this when done with visualizations

  // coarsen_to_n_vertices = 5;
  // map_registered(iSData, coarsen_to_n_vertices);
  // update_texture();

  polyscope::init();
  // psCoarseMesh  = polyscope::registerSurfaceMesh("coarse mesh", U, H); // TODO: visualize coarse mesh
  psMesh = polyscope::registerSurfaceMesh("input mesh", V, F);

  if (using_texture) {
    // convert parameterization to polyscope's desired input format
    Eigen::Matrix<glm::vec2, Eigen::Dynamic, 1> parameterization(3 * UF.rows());
    for (int iF = 0; iF < UF.rows(); iF++) {
      for (int iC = 0; iC < 3; iC++) {
        parameterization(3 * iF + iC) = glm::vec2{UV(UF(iF, iC), 0), UV(UF(iF, iC), 1)};
      }
    }
    q = psMesh->addParameterizationQuantity("intrinsic triangulation", parameterization);

    q->setTexture(TEXTURE_WIDTH, TEXTURE_HEIGHT, polyscope_texture, polyscope::TextureFormat::RGBA8);
    q->setEnabled(true);
    q->setStyle(polyscope::ParamVizStyle::TEXTURE);
    q->setCheckerSize(1);
  }

  polyscope::state::userCallback = callback;
  polyscope::show();
}