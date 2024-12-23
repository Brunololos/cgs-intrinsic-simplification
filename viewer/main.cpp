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
enum class ISIMP_MODE {
  ON_INPUT,
  EACH_FRAME
};

int main(int argc, char *argv[])
{
  if (argc<2) {
    std::cout << "Please pass an input mesh as first argument" << std::endl;
    exit(1);
  }

  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  Eigen::MatrixXd UV;
  MatrixXuc R, G, B, A;
  int TEXTURE_WIDTH = 1000;
  int TEXTURE_HEIGHT = 1000;
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

  // TODO:
  std::pair<Eigen::MatrixXd, double> result; // result of the last simplification step
  std::ofstream logFile; // logfile
  bool doLogging = false;

  // Viewer
  igl::opengl::glfw::Viewer viewer;
  igl::opengl::glfw::imgui::ImGuiPlugin plugin;

  // Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> R,G,B,A;
  // igl::png::readPNG(argc > 2 ? argv[2] : "../textures/texture.png", R,G,B,A);
  // igl::png::readPNG("../textures/deformedCube (during optimisation).png", R,G,B,A);
  R = 255 * MatrixXuc::Ones(TEXTURE_HEIGHT, TEXTURE_WIDTH);
  G = 255 * MatrixXuc::Ones(TEXTURE_HEIGHT, TEXTURE_WIDTH);
  B = 255 * MatrixXuc::Ones(TEXTURE_HEIGHT, TEXTURE_WIDTH);
  A = 255 * MatrixXuc::Ones(TEXTURE_HEIGHT, TEXTURE_WIDTH);

  Eigen::MatrixXd CM = (Eigen::MatrixXd(9,3)<<
      //0.8,0.8,0.8,                        // 0 gray
      1,1,1,
      1.0/255.0,31.0/255.0,255.0/255.0,   // 1 blue
      0/255.0,201.0/255.0,195.0/255.0,    // 2 cyan
      7.0/255.0,152.0/255.0,59.0/255.0,   // 3 green
      190.0/255.0,255.0/255.0,0.0/255.0,  // 4 lime
      255.0/255.0,243.0/255.0,0.0/255.0,  // 5 yellow
      245.0/255.0,152.0/255.0,0.0/255.0,  // 6 orange
      224.0/255.0,18.0/255.0,0.0/255.0,   // 7 red
      220.0/255.0,13.0/255.0,255.0/255.0  // 8 magenta
    ).finished();

  // coloring
  int numVorF = (true /* true: visualize over vertices, false: visualize over faces */) ? V.rows() : F.rows();
  Eigen::MatrixXd C = Eigen::MatrixXd(numVorF, 3);
  for ( int i=0; i<numVorF; i++)
  {
    C.row(i) = CM.row(i % 9);
  }

  auto register_texels = [&](){
        std::cout << "Registering texels: " << std::endl;
        for (gcs::Face F : iSData.inputMesh->faces())
        {
          std::cout << "-> for extrinsic face: " << F << std::endl;
          for (gcs::Vertex V : F.adjacentVertices())
          {
            std::cout << "---> texturing intrinsic vertex: " << V << std::endl;
            int v_idx = V.getIndex();
            Eigen::Vector2d uv = UV.row(v_idx);
            std::cout << "---> that maps to u: " << uv[0] << ", v: " << uv[1] << std::endl;
            // int x = TEXTURE_WIDTH - (int) ((uv[0] + 1.0) / 2.0 * TEXTURE_WIDTH);
            // int y = TEXTURE_HEIGHT - (int) ((uv[1] + 1.0) / 2.0 * TEXTURE_HEIGHT);
            int x = (int) ((uv[0]) * TEXTURE_WIDTH);
            int y = (int) ((uv[1]) * TEXTURE_HEIGHT);
            if (x < 0) { x += TEXTURE_WIDTH; }
            if (y < 0) { y += TEXTURE_HEIGHT; }
          }
          // NOTE: this approach assumes, that all triangles lie completely within the texture. (No triangles wrap around the outside.)
          // I note this here, because lscm sometimes yields uv-coordinates that exceed the [0, 1] bounds. (They can be negative and their absolute value can exceed 1.)
          // To fix this I translated & rescaled all uv's to be withing [0, 1]. Because we don't have a set texture and build it ourselves, this is fine.

          // TODO: find triangle spanned by the uv-coordinates of the three corner vertices
          // TODO: find bounds on texels that might lie within the triangle
          // TODO: iterate through the texel-candidates and check for each if it lies within the triangle.
          // TODO: ---> IF it lies within the triangle, reexpress it barycentrically vie the uv-coordinates and register it for the mapping.
        }
  };

  auto init_viewer_data = [&](){
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
  Eigen::Array<double,Eigen::Dynamic,1> W = Eigen::Array<double,Eigen::Dynamic,1>::Zero(V.rows());
  Eigen::Array<double,Eigen::Dynamic,1> and_visible = Eigen::Array<double,Eigen::Dynamic,1>::Zero(V.rows());

  Eigen::MatrixXd all_constraints = V;

  Eigen::Array<int,-1,1> constraints_mask = Eigen::Array<int,-1,1>::Zero(V.rows());
  Eigen::Array<int,-1,1> selected_mask = Eigen::Array<int,-1,1>::Zero(V.rows());
  Eigen::VectorXi selected_ids;

  Eigen::MatrixXd constraints;
  Eigen::VectorXi bi;
  Eigen::VectorXi indices = Eigen::VectorXd::LinSpaced(V.rows(), 0.0, V.rows() - 1).cast<int>();


  igl::AABB<Eigen::MatrixXd, 3> tree;
  tree.init(V,F);

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
  
  const auto & update_texture = [&]()
  {
    // Eigen::MatrixXuc M;
    // TODO: set R, G, B depending on mapping to coarse mesh
    // YImage yimg;
    // yimg.load(argv[2]);
    // R.resize(yimg.width(),yimg.height());
    // B.resize(R.rows(),R.cols());
    // G.resize(R.rows(),R.cols());
    // for(int i = 0;i<yimg.width();i++)
    // {
    //   for(int j = 0;j<yimg.height();j++)
    //   {
    //     const auto & p = yimg.at(i,j);
    //     R(i,yimg.height()-1-j) = p.r;
    //     G(i,yimg.height()-1-j) = p.g;
    //     B(i,yimg.height()-1-j) = p.b;
    //   }
    // }
    // viewer.data().set_texture(R, G, B);
  };

  const auto update_mesh_points = [&](Eigen::VectorXi& selected_indices, Eigen::MatrixXd& selected_constraints)
  {
    for (int i=0; i<selected_indices.size(); i++)
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
      //std::cout << "Set mesh point " << i << std::endl;
    }
  };

  // TODO: maybe make selected points translate with vertices during optimisation
  const auto update_points = [&]()
  {
    // constrained
    viewer.data().clear_points();
    if (view_constraints || view_all_points) {
      viewer.data().set_points(igl::slice_mask(all_constraints, (selected_mask * constraints_mask).cast<bool>(), 1), CM.row(6));
      viewer.data().add_points(igl::slice_mask(/* all_constraints */ U, (selected_mask * (1-constraints_mask)).cast<bool>(), 1), CM.row(5));
      viewer.data().add_points(igl::slice_mask(all_constraints, ((1-selected_mask) * constraints_mask).cast<bool>(), 1), CM.row(7));
    }
    if (view_all_points) {
      viewer.data().add_points(igl::slice_mask(/* all_constraints */ U, ((1-selected_mask) * (1-constraints_mask)).cast<bool>(), 1), CM.row(2));
    }
  };

  // update entire mesh
  auto refresh_mesh = [&](){
    viewer.data().clear();
    viewer.data().set_mesh(U,H);
    init_viewer_data();
    update_texture();
    update_points();
    viewer.draw();
  };

  // update mesh vertices
  auto refresh_mesh_vertices = [&](){
    viewer.data().set_vertices(U);
    init_viewer_data();
    update_texture();
    update_points();
    viewer.draw();
  };

  selection.callback = [&]()
  {
    screen_space_selection(/* all_constraints */ U, F, tree, viewer.core().view, viewer.core().proj, viewer.core().viewport, selection.L, W, and_visible);
    if (only_visible){
      W = (W * and_visible);
    }
    selected_mask = (W.abs()>0.5).cast<int>();
    selected_ids = igl::slice_mask(indices, selected_mask.cast<bool>(), 1);
    selection.mode=selection.OFF;
    update_points();
    viewer.draw();
  };

  auto update_constraints = [&](bool draw){
    if (virgin) return;
    bi = igl::slice_mask(indices, constraints_mask.cast<bool>(), 1);
    constraints = igl::slice(all_constraints, bi, 1);
    // does all precomputations at once
    igl::parallel_for(3, [&](const int i) {
      // isimp_set_constraints(bi, constraints, iSData);
      update_mesh_points(bi, constraints);
    }, 1);
    update_points();
    if(draw) { viewer.draw(); }
  };

  viewer.callback_pre_draw = [&](igl::opengl::glfw::Viewer& viewer)->bool
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

  viewer.callback_mouse_down = [&](igl::opengl::glfw::Viewer& viewer, int, int)->bool
  {
    translate = false;
    selection.mode=selection.OFF;
    update_constraints(true);
    refresh_mesh_vertices();
    return false;
  };

  Eigen::RowVector3f last_mouse;
  viewer.callback_mouse_move = [&](igl::opengl::glfw::Viewer& viewer, int x, int y)
  {
    if (!translate) return false;
    Eigen::RowVector3f drag_mouse(viewer.current_mouse_x, viewer.core().viewport(3) - viewer.current_mouse_y, last_mouse(2));
    Eigen::RowVector3f drag_scene,last_scene;
    igl::unproject(drag_mouse, viewer.core().view, viewer.core().proj, viewer.core().viewport, drag_scene);
    igl::unproject(last_mouse, viewer.core().view, viewer.core().proj, viewer.core().viewport, last_scene);
    for (size_t i = 0; i < selected_ids.size(); i++) {
      all_constraints.row(selected_ids(i)) += (drag_scene-last_scene).cast<double>();
      // iSData.changedV[selected_ids(i)] = true; // TODO: implement constraint forwarding to ISIMP
    }
    constraints = igl::slice(all_constraints, bi, 1);
    last_mouse = drag_mouse;
    update_points();
    viewer.draw();
    return false;
  };



  viewer.callback_key_pressed = [&](decltype(viewer) &,unsigned int key, int mod)
  {
    // std::pair<double, double> optEnergies;
    std::pair<Eigen::MatrixXd, double> result;
    std::chrono::system_clock::time_point timestamp;
    std::time_t time;
    std::string filename;
    double theta = 0;
    double U22, U21;
    Eigen::MatrixXi F_new = Eigen::MatrixXi(F.rows(), F.cols());
    int i = 0;
    int r;
    // TODO: remove later
    gcs::Edge eddge = iSData.inputMesh->edge(12);

    // TODO: remove lateer
    BarycentricPoint bp;
    // bp << 0.25, 0.25, 0.5; // Don't forget that the barycentric coordinates have to sum to unity.
    // bp << 1.0/3.0, 1.0/3.0, 1.0/3.0; // Don't forget that the barycentric coordinates have to sum to unity.
    bp << 1.0, 0.0, 0.0; // Don't forget that the barycentric coordinates have to sum to unity.

    switch(key)
    {
      case ' ':
        // TODO: later change .inputMesh to .intrinsicMesh
        r = rand() % iSData.inputMesh->nEdges();
        r = 1;
        std::cout << "flip edge: " << r << std::endl;
        flip_intrinsic(iSData, iSData.inputMesh->edge(r));
        for (gcs::Face F : iSData.inputMesh->faces())
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
        map_registered(iSData);
        refresh_mesh();
        if (virgin) return true;
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
        } else {
          optimise = !optimise;
        }
        viewer.draw();
        return true;
      case 'H': case 'h': only_visible = !only_visible; update_points(); viewer.draw(); return true;
      case 'N': case 'n':
        // TODO: find new use
        return true;
      case 'A': case 'a':
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
      case 'W': case 'w':
        // update texture
        std::cout << "Updating texture: " << std::endl;
        for (gcs::Face F : iSData.intrinsicMesh->faces())
        {
          std::cout << "-> texturing intrinsic face: " << F << std::endl;
          for (gcs::Vertex V : F.adjacentVertices())
          {
            std::cout << "---> texturing intrinsic vertex: " << V << std::endl;
            int v_idx = V.getIndex();
            Eigen::Vector2d uv = UV.row(v_idx);
            std::cout << "---> that maps to u: " << uv[0] << ", v: " << uv[1] << std::endl;
            // int x = TEXTURE_WIDTH - (int) ((uv[0] + 1.0) / 2.0 * TEXTURE_WIDTH);
            // int y = TEXTURE_HEIGHT - (int) ((uv[1] + 1.0) / 2.0 * TEXTURE_HEIGHT);
            int x = (int) ((uv[0]) * TEXTURE_WIDTH);
            int y = (int) ((uv[1]) * TEXTURE_HEIGHT);
            if (x < 0) { x += TEXTURE_WIDTH; }
            if (y < 0) { y += TEXTURE_HEIGHT; }
            std::cout << "---> Calculated texture position of uv. x: " << x << ", y: " << y << std::endl;
            G(x, y) = 0;
          }
        }
        // G(22, 88) = 0;
        // G(22, 12) = 0;
        // G(43, 23) = 0;
        for (int i=0; i<TEXTURE_WIDTH; i++)
        {
          for (int j=0; j<TEXTURE_HEIGHT; j++)
          {
            R(i, j) = (int) ((((double) i) * 255.0) / ((double) TEXTURE_WIDTH));
            B(i, j) = (int) ((((double) j) * 255.0) / ((double) TEXTURE_HEIGHT));
          }
        }
        refresh_mesh();
        // igl::png::writePNG(R, G, B, A, "textures/object_texture.png");
        std::cout << "write_success: " << igl::png::writePNG(R, G, B, A, "H:/GIT/cgs-intrinsic-simplification/textures/object_texture.png") << std::endl;
        return true;
      case 'E': case 'e':
        // Map registered points
        map_registered(iSData);
        std::cout << "Mapped points:\n";
        for(int i=0; i<iSData.mapped_by_triangle.size(); i++)
        {
          for(int j=0; j<iSData.mapped_by_triangle[i].size(); j++)
          {
            std::cout << to_str(iSData.tracked_points[j]) << " => " << to_str(iSData.mapped_points[j]) << std::endl;
          }
          if(iSData.mapped_by_triangle[i].size() > 0)
          {
            std::cout << " to triangle: " << i << std::endl;
          }
        }
        return true;
      case 'J': case 'j':
        simpMode = (simpMode == ISIMP_MODE::ON_INPUT) ? ISIMP_MODE::EACH_FRAME : ISIMP_MODE::ON_INPUT;
        return true;
      case 'Y': case 'y':
        (logFile).close();
        return true;
      case 'C': case 'c':
        register_point(iSData, bp, 0);
        std::cout << "Registered point barycentric point: (" << bp[0] << ", " << bp[1] << ", " << bp[2] << ") at ";
        printGCSFace(iSData.inputMesh->face(0));
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
      case 'U': case 'u':
        constraints_mask = ((constraints_mask-selected_mask)==1).cast<int>();
        if (constraints_mask.sum() == 0){
          virgin = true;
          return true;
        }
        update_constraints(true);
        return true;
      case 'G': case 'g':
        if (virgin) return true;
        translate = !translate;
        {
          Eigen::MatrixXd CP;
          Eigen::Vector3d mean = constraints.colwise().mean();
          igl::project(mean, viewer.core().view, viewer.core().proj, viewer.core().viewport, CP);
          last_mouse = Eigen::RowVector3f(viewer.current_mouse_x, viewer.core().viewport(3) - viewer.current_mouse_y, CP(0,2));
          update_constraints(true);
        }
        return true;
      case 'D': case 'd':
        if (virgin) return true;
        selected_mask = (W.abs()>0.5).cast<int>();
        selected_ids = igl::slice_mask(indices, selected_mask.cast<bool>(), 1);
        U = rotate_selected(U, selected_ids, rot_axis, 10.0);
        all_constraints = U;
        update_constraints(true);
        refresh_mesh_vertices();
        break;
      case 'K': case 'k':
        selection.mode = selection.RECTANGULAR_MARQUEE;
        selection.clear();
        return true;
      case 'S': case 's':
        {
          auto t = std::time(nullptr);
          auto tm = *std::localtime(&t);
          std::ostringstream oss;
          oss << std::put_time(&tm, "%d-%m-%Y %H-%M-%S");
          screenshot(viewer, oss.str()+".png", 1);
        }
        return true;
      case 'T': case 't':
        // if (mode == 0)
        // {
        //   viewer.data().set_texture(R,G,B,A);
        // } else {
        //   Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> R,G,B,A;
        //   igl::png::readPNG(argc > 3 ? argv[3] : "../textures/texture.png", R,G,B,A);
        //   viewer.data().set_texture(R,G,B,A);
        // }
        return true;
      case 'R': case 'r':
        // reset mesh
        U = V; refresh_mesh(); return true;
      case 'F': case 'f':
        // filter
        selected_mask = selected_mask*constraints_mask;
        selected_ids = igl::slice_mask(indices, selected_mask.cast<bool>(), 1);
        update_points();
        viewer.draw();
        return true;
      case 'P': case 'p':
        view_all_points = !view_all_points;
        update_points();
        viewer.draw();
        return true;
      case '.':
        view_constraints = !view_constraints;
        if (!view_constraints) { view_all_points = false; }
        update_points();
        viewer.draw();
        return true;
      case '-':
        switch(rot_axis)
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
      case 'X': case 'x':
        if (virgin) return true;
        // reset constrainted positions of selected
        for (size_t i = 0; i < selected_mask.size(); i++) {
          if (selected_mask(i))
            all_constraints.row(i) = V.row(i);
        }
        update_constraints(true);
        return true;
      case 'Q': case 'q':
        std::vector<Eigen::VectorXd> points;
        std::vector<std::array<int, 2>> edges;
        Eigen::MatrixXd P1, P2;
        P1.resize(0,3);
        P2.resize(0,3);
        viewer.data().set_edges(P1, P2.cast<int>(), CM.row(4));
        viewer.data().add_edges(P1, P2, CM.row(4));
        return true;
    }
    return false;
  };

  viewer.plugins.push_back(&plugin);
  plugin.widgets.push_back((igl::opengl::glfw::imgui::ImGuiWidget*) &selection);
  selection.mode = selection.OFF;

  std::cout<<R"( igl::opengl::glfw::Viewer usage:
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
  viewer.data().set_mesh(V,F);
  lscm(V, F, UV);
  // scale factor, such that UV's range from [-1, 1]
  double absmax = std::max(UV.maxCoeff(), -UV.minCoeff());
  // transforming UV's to range from [0, 1]
  UV = (UV / absmax + Eigen::MatrixXd::Ones(UV.rows(), UV.cols()))/ 2.0;
  register_texels();
  init_viewer_data();
  viewer.core().is_animating = true;
  viewer.core().background_color.head(3) = CM.row(0).head(3).cast<float>();
  // update_points()
  viewer.launch();
}
