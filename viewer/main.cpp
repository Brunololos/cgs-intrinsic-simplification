#define igl_viewer_viewer_quiet
// #define MIN_QUAD_WITH_FIXED_CPP_DEBUG
#include <igl/opengl/glfw/Viewer.h>
#include <igl/AABB.h>
#include <igl/screen_space_selection.h>
#include <igl/opengl/glfw/imgui/ImGuiPlugin.h>
#include <igl/opengl/glfw/imgui/SelectionWidget.h>
#include <stdlib.h>
#include <igl/arap.h>
#include <igl/png/readPNG.h>
#include <stdlib.h>
#include <chrono>
#include <ctime>

#include "isimp.hpp"
#include "screenshot.hpp"
#include "lscm_wrapper.hpp"

typedef Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> MatrixXuc;

// TODO: move this somewhere else
enum class ISIMP_MODE
{
  ON_INPUT,
  EACH_FRAME
};

int main(int argc, char *argv[])
{
  if (argc < 2)
  {
    std::cout << "Please pass an input mesh as first argument" << std::endl;
    exit(1);
  }

  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  Eigen::MatrixXd UV;
  MatrixXuc R, G, B, A;
  int TEXTURE_WIDTH = 1000;
  int TEXTURE_HEIGHT = 1000;

  // int dot_pos = ((std::string) argv[1]).find_last_of('.');
  // std::cout << ((std::string) argv[1]).substr(dot_pos) << std::endl;
  // if (((std::string) argv[1]).substr(dot_pos) == ".obj")
  // {
  //   igl::readOBJ(argv[1], V, F);
  // }
  // else
  // {
  //   igl::read_triangle_mesh(argv[1], V, F);
  // }
  igl::read_triangle_mesh(argv[1], V, F);

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

  // TODO: remove
  int eid = 0;
  int vid = 0;
  // randomize seed
  srand(time(NULL));

  // TODO:
  std::pair<Eigen::MatrixXd, double> result; // result of the last simplification step
  std::ofstream logFile;                     // logfile
  bool doLogging = false;

  // Viewer
  igl::opengl::glfw::Viewer viewer;
  igl::opengl::glfw::imgui::ImGuiPlugin plugin;

  // Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> R,G,B,A;
  // igl::png::readPNG(argc > 2 ? argv[2] : "../textures/texture.png", R,G,B,A);
  // igl::png::readPNG("../textures/deformedCube (during optimisation).png", R,G,B,A);
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
                            1,
                        1, 1,                                       // 0 white
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
    // prev face-based coloring: viewer.data().set_face_based(true);
    // viewer.data().set_texture(R,G,B,A);
    // viewer.data().set_colors(C);

    // Calculate texture mapping
    viewer.data().set_uv(UV);
    viewer.data().set_texture(R, G, B);
    viewer.data().show_texture = true;

    // viewer.data().use_matcap = true;
    viewer.data().point_size = 15;
    viewer.data().line_width = 1;
    viewer.data().show_lines = F.rows() < 20000;
  };

  // Selection
  bool translate = false;
  RAXIS rot_axis = RAXIS::YAW;
  bool only_visible = false;
  igl::opengl::glfw::imgui::SelectionWidget selection;
  Eigen::Array<double, Eigen::Dynamic, 1> W = Eigen::Array<double, Eigen::Dynamic, 1>::Zero(V.rows());
  Eigen::Array<double, Eigen::Dynamic, 1> and_visible = Eigen::Array<double, Eigen::Dynamic, 1>::Zero(V.rows());

  Eigen::MatrixXd all_constraints = V;

  Eigen::Array<int, -1, 1> constraints_mask = Eigen::Array<int, -1, 1>::Zero(V.rows());
  Eigen::Array<int, -1, 1> selected_mask = Eigen::Array<int, -1, 1>::Zero(V.rows());
  Eigen::VectorXi selected_ids;

  Eigen::MatrixXd constraints;
  Eigen::VectorXi bi;
  Eigen::VectorXi indices = Eigen::VectorXd::LinSpaced(V.rows(), 0.0, V.rows() - 1).cast<int>();

  igl::AABB<Eigen::MatrixXd, 3> tree;
  tree.init(V, F);

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

  const auto &update_texture = [&]()
  {
    int num_triangles = iSData.mapped_by_triangle.size();
    int num_colors = CM.rows() - 1;
    for (int i = 0; i < num_triangles; i++)
    {
      // Color = colors(i%10)
      for (int j = 0; j < iSData.mapped_by_triangle[i].size(); j++)
      {
        int original_idx = iSData.mapped_by_triangle[i][j];
        // BarycentricPoint original_point = iSData.tracked_points[original_idx];
        // Point2D tex_coord = to_explicit(original_point, );
        TexCoord texcoord = iSData.tracked_texcoords[original_idx];
        R(texcoord[0], texcoord[1]) = CM((i % num_colors) + 1, 0) * 255;
        G(texcoord[0], texcoord[1]) = CM((i % num_colors) + 1, 1) * 255;
        B(texcoord[0], texcoord[1]) = CM((i % num_colors) + 1, 2) * 255;
        // std::cout << "Setting texcoord: (" << texcoord[0] << ", " << texcoord[1] << "): R=" << CM(i%9, 0) << ", G=" << CM(i%9, 1) << ", B=" << CM(i%9, 2) << std::endl;
      }
    }
    viewer.data().set_texture(R, G, B);
    std::cout << "write_success: " << igl::png::writePNG(R, G, B, A, "H:/GIT/cgs-intrinsic-simplification/textures/object_texture.png") << std::endl;
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
    viewer.data().clear_points();
    if (view_constraints || view_all_points)
    {
      viewer.data().set_points(igl::slice_mask(all_constraints, (selected_mask * constraints_mask).cast<bool>(), 1), CM.row(6));
      viewer.data().add_points(igl::slice_mask(/* all_constraints */ U, (selected_mask * (1 - constraints_mask)).cast<bool>(), 1), CM.row(5));
      viewer.data().add_points(igl::slice_mask(all_constraints, ((1 - selected_mask) * constraints_mask).cast<bool>(), 1), CM.row(7));
    }
    if (view_all_points)
    {
      viewer.data().add_points(igl::slice_mask(/* all_constraints */ U, ((1 - selected_mask) * (1 - constraints_mask)).cast<bool>(), 1), CM.row(2));
    }
  };

  // update entire mesh
  auto refresh_mesh = [&]()
  {
    viewer.data().clear();
    viewer.data().set_mesh(U, H);
    init_viewer_data();
    // update_texture();
    update_points();
    viewer.draw();
  };

  // update mesh vertices
  auto refresh_mesh_vertices = [&]()
  {
    viewer.data().set_vertices(U);
    init_viewer_data();
    // update_texture();
    update_points();
    viewer.draw();
  };

  selection.callback = [&]()
  {
    screen_space_selection(/* all_constraints */ U, F, tree, viewer.core().view, viewer.core().proj, viewer.core().viewport, selection.L, W, and_visible);
    if (only_visible)
    {
      W = (W * and_visible);
    }
    selected_mask = (W.abs() > 0.5).cast<int>();
    selected_ids = igl::slice_mask(indices, selected_mask.cast<bool>(), 1);
    selection.mode = selection.OFF;
    update_points();
    viewer.draw();
  };

  auto update_constraints = [&](bool draw)
  {
    if (virgin)
      return;
    bi = igl::slice_mask(indices, constraints_mask.cast<bool>(), 1);
    constraints = igl::slice(all_constraints, bi, 1);
    // does all precomputations at once
    igl::parallel_for(3, [&](const int i)
                      {
      // isimp_set_constraints(bi, constraints, iSData);
      update_mesh_points(bi, constraints); }, 1);
    update_points();
    if (draw)
    {
      viewer.draw();
    }
  };

  viewer.callback_pre_draw = [&](igl::opengl::glfw::Viewer &viewer) -> bool
  {
    // TODO: implement simplification
    // if (simpMode == ISIMP_MODE::EACH_FRAME && optimise && !iSData.hasConverged)
    // {
    //   // Simplification
    //   std::cout << "simplification iteration: " << iSData.iteration;
    //   // TODO: reduce
    //   // std::pair<double, double> optEnergies = asdap_energies(iSData, U);
    //   // // logging
    //   // if (doLogging)
    //   // {
    //   //   logFile << std::to_string(optEnergies.first + optEnergies.second) << ";" << std::to_string(optEnergies.first) << ";" << std::to_string(optEnergies.second) << "\n";
    //   // }
    //   // TODO:
    //   // std::cout << " => Energy: " << (optEnergies.first + optEnergies.second) << " (Scaling: " << optEnergies.first << ", Rotation: " << optEnergies.second << ")" << std::endl;
    //   // update viewer
    //   viewer.data().set_vertices(U);
    //   update_points();
    // }
    // if (optimise && iSData.hasConverged)
    // {
    //   if (doLogging) { (logFile).close(); doLogging = false; }
    //   std::cout << "optimisation converged!" << std::endl;
    //   optimise = false;
    // }
    return false;
  };

  viewer.callback_mouse_down = [&](igl::opengl::glfw::Viewer &viewer, int, int) -> bool
  {
    translate = false;
    selection.mode = selection.OFF;
    update_constraints(true);
    refresh_mesh_vertices();
    return false;
  };

  Eigen::RowVector3f last_mouse;
  viewer.callback_mouse_move = [&](igl::opengl::glfw::Viewer &viewer, int x, int y)
  {
    if (!translate)
      return false;
    Eigen::RowVector3f drag_mouse(viewer.current_mouse_x, viewer.core().viewport(3) - viewer.current_mouse_y, last_mouse(2));
    Eigen::RowVector3f drag_scene, last_scene;
    igl::unproject(drag_mouse, viewer.core().view, viewer.core().proj, viewer.core().viewport, drag_scene);
    igl::unproject(last_mouse, viewer.core().view, viewer.core().proj, viewer.core().viewport, last_scene);
    for (size_t i = 0; i < selected_ids.size(); i++)
    {
      all_constraints.row(selected_ids(i)) += (drag_scene - last_scene).cast<double>();
      // iSData.changedV[selected_ids(i)] = true; // TODO: implement constraint forwarding to ISIMP
    }
    constraints = igl::slice(all_constraints, bi, 1);
    last_mouse = drag_mouse;
    update_points();
    viewer.draw();
    return false;
  };

  viewer.callback_key_pressed = [&](decltype(viewer) &, unsigned int key, int mod)
  {
    // std::pair<double, double> optEnergies;
    std::pair<Eigen::MatrixXd, double> result;
    std::chrono::system_clock::time_point timestamp;
    std::time_t time;
    std::string filename;
    double theta = 0;
    double U22, U21;
    Eigen::MatrixXi F_new = Eigen::MatrixXi(F.rows(), F.cols());
    std::vector<gcs::Vertex> temp_neighbors = std::vector<gcs::Vertex>();
    std::pair<int, double>* current;
    int i = 0;
    int r;
    double ice;
    bool could_flip, could_flatten, could_remove, could_flip_to_deg3;
    // TODO: remove later
    gcs::Edge eddge = iSData.inputMesh->edge(12);
    std::vector<int> e_idxs = std::vector<int>();
    std::vector<int> v_idxs = std::vector<int>();
    // e_idxs.push_back(15);
    // e_idxs.push_back(11);
    // e_idxs.push_back(16);
    // e_idxs.push_back(-1);
    // // e_idxs.push_back(8);
    // // e_idxs.push_back(3);
    // // e_idxs.push_back(-1);
    // e_idxs.push_back(22);
    // e_idxs.push_back(20);
    // e_idxs.push_back(16);
    // e_idxs.push_back(19);
    // e_idxs.push_back(-1);
    // v_idxs.push_back(7);
    // // v_idxs.push_back(4);
    // v_idxs.push_back(10);
    e_idxs.push_back(40);
    e_idxs.push_back(14);
    e_idxs.push_back(20);
    e_idxs.push_back(-1);
    e_idxs.push_back(14);
    e_idxs.push_back(13);
    e_idxs.push_back(12);
    e_idxs.push_back(18);
    e_idxs.push_back(-1);
    // e_idxs.push_back(49);
    // e_idxs.push_back(45);
    // e_idxs.push_back(50);
    // e_idxs.push_back(-1);
    // e_idxs.push_back(43);
    // e_idxs.push_back(-1);
    e_idxs.push_back(23);
    e_idxs.push_back(-1);
    // e_idxs.push_back(33);
    // e_idxs.push_back(7);
    // e_idxs.push_back(-1);
    // e_idxs.push_back(37);
    // e_idxs.push_back(-1);
    // e_idxs.push_back(45);
    // e_idxs.push_back(43);
    // e_idxs.push_back(37);
    // e_idxs.push_back(33);
    e_idxs.push_back(18);
    // e_idxs.push_back(27);
    // e_idxs.push_back(-1);

    v_idxs.push_back(8);
    v_idxs.push_back(7);
    // v_idxs.push_back(21);
    // v_idxs.push_back(19);
    v_idxs.push_back(10);
    // v_idxs.push_back(5);
    // v_idxs.push_back(17);
    // v_idxs.push_back(11);

    // TODO: remove lateer
    BarycentricPoint bp;
    // bp << 0.25, 0.25, 0.5; // Don't forget that the barycentric coordinates have to sum to unity.
    // bp << 1.0/3.0, 1.0/3.0, 1.0/3.0; // Don't forget that the barycentric coordinates have to sum to unity.
    bp << 1.0, 0.0, 0.0; // Don't forget that the barycentric coordinates have to sum to unity.

    switch (key)
    {
    case ' ':
      // TODO: later change .inputMesh to .intrinsicMesh

      // NOTE: Performing random intrinsic edge-flips
      // r = rand() % iSData.inputMesh->nEdges();
      // r = 1;
      // // for (gcs::Edge E : iSData.intrinsicMesh->edges())
      // TODO: just add this check into can_be_flipped() method (don't allow boundary flips)
      // while (true)
      // {
      //   r = rand() % iSData.intrinsicMesh->nEdges();
      //   gcs::Edge E = iSData.intrinsicMesh->edge(r);
      //   int nk = 0;
      //   for (gcs::Face F : E.adjacentFaces())
      //   {
      //     nk++;
      //   }
      //   if (nk == 2)
      //   {
      //     r = E.getIndex();
      //     break;
      //   }
      // }
      // std::cout << "flip edge: " << r << " => (" << iSData.intrinsicMesh->edge(r).firstVertex().getIndex() << ", " << iSData.intrinsicMesh->edge(r).secondVertex().getIndex() << ")" << std::endl;
      // could_flip = flip_intrinsic(iSData, iSData.intrinsicMesh->edge(r));
      // if (could_flip)
      // {
      //   std::cout << "Edge flip successful!" << std::endl;
      // }
      // if (!could_flip)
      // {
      //   std::cout << "Edge flip failed!" << std::endl;
      // }

      // NOTE: Performing random intrinsic vertex flattenings
      // r = rand() % iSData.intrinsicMesh->nVertices();
      // while (true)
      // {
      //   r = rand() % iSData.intrinsicMesh->nVertices();
      //   gcs::Vertex V = iSData.intrinsicMesh->vertex(r);
      //   if (!V.isBoundary()) { break; }
      // }
      // // r = 115;
      // std::cout << "flatten vertex: " << r << std::endl;
      // could_flatten = flatten_vertex(iSData, r);
      // if(could_flatten) { std::cout << "Vertex Flattening successful!" << std::endl; }
      // if(!could_flatten) { std::cout << "Vertex Flattening failed!" << std::endl; }

      // r = rand() % V.rows(); // iSData.intrinsicMesh->nVertices();
      // i = 0;
      // while (i < 100)
      // {
      //   r = rand() % V.rows(); // iSData.intrinsicMesh->nVertices();
      //   gcs::Vertex V = iSData.intrinsicMesh->vertex(r);
      //   // std::cout << "Checking vertex: " << r << std::endl;
      //   if (!V.isDead() && !V.isBoundary())
      //   {
      //     // r = 115;
      //     std::cout << "\n\n\nflatten vertex: " << r << std::endl;
      //     could_flatten = flatten_vertex(iSData, r);
      //     if(!could_flatten) { std::cout << "Vertex Flattening failed!" << std::endl; break; }
      //     std::cout << "\n\n\nflip vertex: " << r << " to degree 3: " << std::endl;
      //     could_flip_to_deg3 = flip_vertex_to_deg3(iSData, r);
      //     if(!could_flip_to_deg3) { std::cout << dye("Vertex Flipping to degree 3 failed!", RED) << std::endl; break; }
      //     std::cout << "\n\n\nremove vertex: " << r << std::endl;
      //     could_remove = remove_vertex(iSData, r);
      //     if(could_remove) { std::cout << "Vertex Removal successful!" << std::endl; }
      //     if(!could_remove) { std::cout << "Vertex Removal failed!" << std::endl; }
      //     break;
      //   }
      //   i++;
      // }

      if(iSData.Q.empty()) { return true; } // TODO: could also later be replaced by a finished flag
      current = iSData.Q.begin()->get();   // this is Q.top()
      iSData.Q.erase(iSData.Q.begin());    // this is Q.pop()
      r = current->first;
      ice = current->second;

      if(std::isinf(ice)) { return true; }
      std::cout << "\n\n\nflatten vertex: " << r << std::endl;
      could_flatten = flatten_vertex(iSData, r);
      if(could_flatten)
      {
        std::cout << "\n\n\nflip vertex: " << r << " to degree 3: " << std::endl;
        could_flip_to_deg3 = flip_vertex_to_deg3(iSData, r);
        if(could_flip_to_deg3)
        {
          std::cout << "\n\n\nremove vertex: " << r << std::endl;
          // save neighbors, because we cant iterate over vertex neighborhood after vertex removal
          temp_neighbors.clear();
          std::cout << "cleared temp_neighbors: ";
          for (gcs::Vertex neighbor : temp_neighbors) { std::cout << neighbor.getIndex() << ", "; }
          std::cout << std::endl;
          for (gcs::Vertex neighbor : iSData.intrinsicMesh->vertex(r).adjacentVertices()) { std::cout << "inserting temp_neighbor: " << neighbor.getIndex() << std::endl; temp_neighbors.push_back(neighbor); }
          could_remove = remove_vertex(iSData, r);
          if(could_remove) { std::cout << "Vertex Removal successful!" << std::endl; }
          else { std::cout << "Vertex Removal failed!" << std::endl; }

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
          double ice = intrinsic_curvature_error(iSData, neighbor);
          iSData.Q.erase(iSData.Q.find(std::make_shared<std::pair<int, double>>(iSData.Q_elems[neighbor.getIndex()]))); // erase the element to be sure the order gets updated
          iSData.Q_elems[neighbor.getIndex()].second = ice;
          iSData.Q.insert(std::make_shared<std::pair<int, double>>(iSData.Q_elems[neighbor.getIndex()]));
        }
      } else {
        iSData.Q_elems[r].second = INFINITY;
        iSData.Q.insert(std::make_shared<std::pair<int, double>>(iSData.Q_elems[r]));
      }

      // if (e_idxs[eid] == -1)
      // {
      //   std::cout << "removing vertex " << v_idxs[vid] << std::endl;
      //   could_remove = remove_vertex(iSData, v_idxs[vid]);
      //   if(!could_remove) { std::cout << "Vertex Removal failed!" << std::endl; }
      //   vid++;
      // }
      // else
      // {
      //   std::cout << "flipping edge " << e_idxs[eid] << "(" << iSData.intrinsicMesh->edge(e_idxs[eid]).firstVertex().getIndex() << ", " << iSData.intrinsicMesh->edge(e_idxs[eid]).secondVertex().getIndex() << ")" << std::endl;
      //   could_flip = flip_intrinsic(iSData, iSData.intrinsicMesh->edge(e_idxs[eid]));
      //   if(!could_flip) { std::cout << "Edge Flip failed!" << std::endl; }
      // }
      // eid++;

      i = 0;
      // I think this was updating the displayed faces for debugging purposes
      for (gcs::Face F : iSData.intrinsicMesh->faces())
      {
        Eigen::Vector3i ffface = Eigen::Vector3i();
        int j = 0;
        for (gcs::Vertex VV : F.adjacentVertices())
        {
          ffface(j) = VV.getIndex();
          j++;
        }
        F_new.row(i) = ffface;
        i++;
      }
      H = F_new;
      // TODO: figure out a way to store the mapping with varying mapping objects in a vector
      // ((std::unique_ptr<Edge_Flip>) iSData.mapping[0])->map_to_coarse(iSData.mapped_points, iSData.mapped_by_triangle);
      // iSData.mapping[0]->map_to_coarse(iSData.mapped_points, iSData.mapped_by_triangle);
      // map_registered(iSData);
      // refresh_mesh();
      map_registered(iSData);
      update_texture();
      // if (iSData.hasConverged) return true; // TODO: replace
      // step manually
      if (simpMode == ISIMP_MODE::ON_INPUT)
      {
        // make simplification step
        // TODO: implement simplification step

        // Logging
        // TODO: implement logging
        // optEnergies = asdap_energies(iSData, U);
        // if(doLogging)
        // {
        //   logFile << std::to_string(optEnergies.first + optEnergies.second) << ";" << std::to_string(optEnergies.first) << ";" << std::to_string(optEnergies.second) << "\n";
        // }
        // std::cout << " => Energy: " << (optEnergies.first + optEnergies.second) << " (Scaling: " << optEnergies.first << ", Rotation: " << optEnergies.second << ")" << std::endl;
        // if (iSData.hasConverged)
        // {
        //   if (doLogging) { (logFile).close(); doLogging = false; }
        //   std::cout << "optimisation converged!" << std::endl;
        // }
      }
      else
      {
        optimise = !optimise;
      }
      viewer.draw();
      return true;
    case 'H':
    case 'h':
      only_visible = !only_visible;
      update_points();
      viewer.draw();
      return true;
    case 'N':
    case 'n':
      // TODO: find new use
      printEigenMatrixXd("L", iSData.L);
      return true;
    case 'A':
    case 'a':
      if (doLogging)
      {
        (logFile).close();
      }
      else
      {
        // TODO: implement logging
        // doLogging = true;
        // timestamp = std::chrono::system_clock::now();
        // time = std::chrono::system_clock::to_time_t(timestamp);
        // filename = "../logs/EnergyCurve" + std::to_string(time) + ".csv";
        // logFile.open(filename);
        // logFile << std::setprecision(15);
      }
      return true;
    case 'W':
    case 'w':
      // // update texture
      // std::cout << "Updating texture: " << std::endl;
      // for (gcs::Face F : iSData.intrinsicMesh->faces())
      // {
      //   std::cout << "-> texturing intrinsic face: " << F << std::endl;
      //   for (gcs::Vertex V : F.adjacentVertices())
      //   {
      //     std::cout << "---> texturing intrinsic vertex: " << V << std::endl;
      //     int v_idx = V.getIndex();
      //     Eigen::Vector2d uv = UV.row(v_idx);
      //     std::cout << "---> that maps to u: " << uv[0] << ", v: " << uv[1] << std::endl;
      //     // int x = TEXTURE_WIDTH - (int) ((uv[0] + 1.0) / 2.0 * TEXTURE_WIDTH);
      //     // int y = TEXTURE_HEIGHT - (int) ((uv[1] + 1.0) / 2.0 * TEXTURE_HEIGHT);
      //     int x = (int) ((uv[0]) * TEXTURE_WIDTH);
      //     int y = (int) ((uv[1]) * TEXTURE_HEIGHT);
      //     // NOTE: negative checks not needed anymore, because we constrain uv's to be within the range [0, 1]
      //     // if (x < 0) { x += TEXTURE_WIDTH; }
      //     // if (y < 0) { y += TEXTURE_HEIGHT; }
      //     std::cout << "---> Calculated texture position of uv. x: " << x << ", y: " << y << std::endl;
      //     G(x, y) = 0;
      //   }
      // }
      // // G(22, 88) = 0;
      // // G(22, 12) = 0;
      // // G(43, 23) = 0;
      // for (int i=0; i<TEXTURE_WIDTH; i++)
      // {
      //   for (int j=0; j<TEXTURE_HEIGHT; j++)
      //   {
      //     R(i, j) = (int) ((((double) i) * 255.0) / ((double) TEXTURE_WIDTH));
      //     B(i, j) = (int) ((((double) j) * 255.0) / ((double) TEXTURE_HEIGHT));
      //   }
      // }
      refresh_mesh();
      // igl::png::writePNG(R, G, B, A, "textures/object_texture.png");
      // std::cout << "write_success: " << igl::png::writePNG(R, G, B, A, "H:/GIT/cgs-intrinsic-simplification/textures/object_texture.png") << std::endl;
      return true;
    case 'E':
    case 'e':
      // Map registered points
      map_registered(iSData);
      std::cout << "Mapped points" << std::endl;
      // std::cout << "Mapped points:\n";
      // for(int i=0; i<iSData.mapped_by_triangle.size(); i++)
      // {
      //   for(int j=0; j<iSData.mapped_by_triangle[i].size(); j++)
      //   {
      //     std::cout << to_str(iSData.tracked_points[j]) << " => " << to_str(iSData.mapped_points[j]) << std::endl;
      //   }
      //   if(iSData.mapped_by_triangle[i].size() > 0)
      //   {
      //     std::cout << " to triangle: " << i << std::endl;
      //   }
      // }
      update_texture();
      return true;
    case 'J':
    case 'j':
      simpMode = (simpMode == ISIMP_MODE::ON_INPUT) ? ISIMP_MODE::EACH_FRAME : ISIMP_MODE::ON_INPUT;
      return true;
    case 'Y':
    case 'y':
      (logFile).close();
      return true;
    case 'C':
    case 'c':
      // register_point(iSData, bp, 0);
      // std::cout << "Registered point barycentric point: (" << bp[0] << ", " << bp[1] << ", " << bp[2] << ") at ";
      // printGCSFace(iSData.inputMesh->face(0));
      // TODO: we don't need constraints
      //   virgin = false;
      //   iSData.hasConverged = false;
      //   constraints_mask = ((constraints_mask+selected_mask)>0).cast<int>();
      //   update_constraints(true);
      //   refresh_mesh();
      //   // std::cout << constraints << "\n\n\n";
      //   // std::cout << bi << "\n\n\n";
      //   ---------------------------------------------------------------
      //   std::array<int, 5> fp_indices = order_quad_edge_indices(eddge);
      //   std::cout << fp_indices[0] << ", " << fp_indices[1] << ", " << fp_indices[2] << ", " << fp_indices[3] << ", " << fp_indices[4] << std::endl;
      //   for (int i=0; i<5; i++)
      //   {
      //     std::cout << "Edge " << fp_indices[i] << ": (" << iSData.inputMesh->edge(fp_indices[i]).firstVertex().getIndex() << ", " << iSData.inputMesh->edge(fp_indices[i]).secondVertex().getIndex() << ")" << std::endl;
      //   }
      //   // flip_intrinsic(iSData.inputMesh, iSData.L, eddge);
      //   flip_intrinsic(iSData, eddge);
      //   for (gcs::Face F : iSData.inputMesh->faces())
      //   {
      //     Eigen::Vector3i ffface = Eigen::Vector3i();
      //     int j = 0;
      //     for (gcs::Vertex VV : F.adjacentVertices())
      //     {
      //       ffface(j) = VV.getIndex();
      //       j++;
      //     }
      //     F_new.row(i) = ffface;
      //     i++;
      //   }
      //   H = F_new;
      //   refresh_mesh();
      //   std::cout << "Edge " << fp_indices[0] << " length: " << iSData.L[fp_indices[0]] << std::endl;
      return true;
    case 'U':
    case 'u':
      constraints_mask = ((constraints_mask - selected_mask) == 1).cast<int>();
      if (constraints_mask.sum() == 0)
      {
        virgin = true;
        return true;
      }
      update_constraints(true);
      return true;
    case 'G':
    case 'g':
      if (virgin)
        return true;
      translate = !translate;
      {
        Eigen::MatrixXd CP;
        Eigen::Vector3d mean = constraints.colwise().mean();
        igl::project(mean, viewer.core().view, viewer.core().proj, viewer.core().viewport, CP);
        last_mouse = Eigen::RowVector3f(viewer.current_mouse_x, viewer.core().viewport(3) - viewer.current_mouse_y, CP(0, 2));
        update_constraints(true);
      }
      return true;
    case 'D':
    case 'd':
      if (virgin)
        return true;
      selected_mask = (W.abs() > 0.5).cast<int>();
      selected_ids = igl::slice_mask(indices, selected_mask.cast<bool>(), 1);
      U = rotate_selected(U, selected_ids, rot_axis, 10.0);
      all_constraints = U;
      update_constraints(true);
      refresh_mesh_vertices();
      break;
    case 'K':
    case 'k':
      selection.mode = selection.RECTANGULAR_MARQUEE;
      selection.clear();
      return true;
    case 'S':
    case 's':
    {
      auto t = std::time(nullptr);
      auto tm = *std::localtime(&t);
      std::ostringstream oss;
      oss << std::put_time(&tm, "%d-%m-%Y %H-%M-%S");
      screenshot(viewer, oss.str() + ".png", 1);
    }
      return true;
    case 'T':
    case 't':
      // if (mode == 0)
      // {
      //   viewer.data().set_texture(R,G,B,A);
      // } else {
      //   Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> R,G,B,A;
      //   igl::png::readPNG(argc > 3 ? argv[3] : "../textures/texture.png", R,G,B,A);
      //   viewer.data().set_texture(R,G,B,A);
      // }
      return true;
    case 'R':
    case 'r':
      // reset mesh
      U = V;
      refresh_mesh();
      return true;
    case 'F':
    case 'f':
      // filter
      selected_mask = selected_mask * constraints_mask;
      selected_ids = igl::slice_mask(indices, selected_mask.cast<bool>(), 1);
      update_points();
      viewer.draw();
      return true;
    case 'P':
    case 'p':
      view_all_points = !view_all_points;
      update_points();
      viewer.draw();
      return true;
    case '.':
      view_constraints = !view_constraints;
      if (!view_constraints)
      {
        view_all_points = false;
      }
      update_points();
      viewer.draw();
      return true;
    case '-':
      switch (rot_axis)
      {
      default:
      case RAXIS::ROLL:
        rot_axis = RAXIS::YAW;
        break;
      case RAXIS::YAW:
        rot_axis = RAXIS::PITCH;
        break;
      case RAXIS::PITCH:
        rot_axis = RAXIS::ROLL;
        break;
      }
      return true;
    case 'X':
    case 'x':
      if (virgin)
        return true;
      // reset constrainted positions of selected
      for (size_t i = 0; i < selected_mask.size(); i++)
      {
        if (selected_mask(i))
          all_constraints.row(i) = V.row(i);
      }
      update_constraints(true);
      return true;
    case 'Q':
    case 'q':
      std::vector<Eigen::VectorXd> points;
      std::vector<std::array<int, 2>> edges;
      Eigen::MatrixXd P1, P2;
      P1.resize(0, 3);
      P2.resize(0, 3);
      viewer.data().set_edges(P1, P2.cast<int>(), CM.row(4));
      viewer.data().add_edges(P1, P2, CM.row(4));
      return true;
    }
    return false;
  };

  viewer.plugins.push_back(&plugin);
  plugin.widgets.push_back((igl::opengl::glfw::imgui::ImGuiWidget *)&selection);
  selection.mode = selection.OFF;

  std::cout << R"( igl::opengl::glfw::Viewer usage:
  I,i     Toggle invert normals
  L       Toggle wireframe
  O,o     Toggle orthographic/perspective projection
  Z       Snap to canonical view
  [,]     Toggle between rotation control types (trackball, two-axis
          valuator with fixed up, 2D mode with no rotation))
  <,>     Toggle between models
  ;       Toggle vertex labels
  :       Toggle face labels/

ASDAP Usage:
  [space]  Run one ISIMP Step (hold for animation)
  J,j      Toggle single Step/continual simplification (changes [space] behaviour)
  N,n      Toggle gradient visibility // TODO: not needed
  A,a      Toggle logging
  R,r      Reset mesh position

  H,h      Toggle whether to take visibility into account for selection
  F,f      Remove non constrained points from selection

  C,c      Add selection to constraints // TODO: not needed
  U,u      Remove selection from constraints

  X,x      Reset constrain position of selected vertices
  G,g      Move selection
  D,d      Rotate selection
  P,p      View all point positions
  .        View constrained point positions
  -        Toggle rotation axis (Yaw -> Pitch -> Roll)

  S,s      Make screenshot
)";

  viewer.data().clear();
  viewer.data().set_mesh(V, F);
  lscm(V, F, UV);
  // scale factor, such that UV's range from [-1, 1]
  double absmax = std::max(UV.maxCoeff(), -UV.minCoeff());
  // transforming UV's to range from [0, 1]
  UV = (UV / absmax + Eigen::MatrixXd::Ones(UV.rows(), UV.cols())) / 2.0;
  register_texels();
  init_viewer_data();
  viewer.core().is_animating = true;
  viewer.core().background_color.head(3) = CM.row(0).head(3).cast<float>();
  // update_points()
  viewer.launch();
}
