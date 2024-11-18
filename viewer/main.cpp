#define IGL_VIEWER_VIEWER_QUIET
// #define MIN_QUAD_WITH_FIXED_CPP_DEBUG
#include <igl/opengl/glfw/Viewer.h>
#include <igl/AABB.h>
#include <igl/screen_space_selection.h>
#include <igl/opengl/glfw/imgui/ImGuiPlugin.h>
#include <igl/opengl/glfw/imgui/SelectionWidget.h>
#include <igl/arap.h>
#include <igl/png/readPNG.h>
#include <stdlib.h>
#include <chrono>
#include <ctime>

#include "isimp.hpp"
#include "screenshot.hpp"

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
  igl::read_triangle_mesh(argv[1], V, F);

  Eigen::MatrixXd U = V;

  // settings
  int mode = 0;
  bool virgin = true;
  bool view_constraints = true;
  bool view_all_points = false;
  bool optimise = false;

  // ISIMP (intrinsic-simplification)
  ISIMP_MODE simpMode = ISIMP_MODE::ON_INPUT;

  iSimpData Adata;
  initMesh(V, F, Adata);

  // TODO:
  std::pair<Eigen::MatrixXd, double> result; // result of the last simplification step
  std::ofstream logFile; // logfile
  bool doLogging = false;

  // Viewer
  igl::opengl::glfw::Viewer viewer;
  igl::opengl::glfw::imgui::ImGuiPlugin plugin;

  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> R,G,B,A;
  igl::png::readPNG(argc > 2 ? argv[2] : "../textures/texture.png", R,G,B,A);

  Eigen::MatrixXd CM = (Eigen::MatrixXd(9,3)<<
      //0.8,0.8,0.8,                        // 0 gray
      1,1,1,
      1.0/255.0,31.0/255.0,255.0/255.0,   // 1 blue
      0/255.0,201.0/255.0,195.0/255.0,    // 2 cyan
      7.0/255.0,152.0/255.0,59.0/255.0,   // 3 green
      190.0/255.0,255.0/255.0,0.0/255.0,  // 4 lime
      255.0/255.0,243.0/255.0,0.0/255.0,  // 5 yellow
      245.0/255.0,152.0/255.0,0.0/255.0,  // 6 orange
      245.0/255.0,152.0/255.0,0.0/255.0,  // 7 orange
      // 224.0/255.0,18.0/255.0,0.0/255.0,   // 7 red
      220.0/255.0,13.0/255.0,255.0/255.0  // 8 magenta
    ).finished();

  auto init_viewer_data = [&](){
    viewer.data().set_face_based(true);
    viewer.data().set_texture(R,G,B,A);
    viewer.data().show_texture = true;
    viewer.data().use_matcap = true;
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
    viewer.data().set_mesh(U,F);
    init_viewer_data();
    update_points();
    viewer.draw();
  };

  // update mesh vertices
  auto refresh_mesh_vertices = [&](){
    viewer.data().set_vertices(U);
    init_viewer_data();
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
      // isimp_set_constraints(bi, constraints, Adata);
      update_mesh_points(bi, constraints);
    }, 1);
    update_points();
    if(draw) { viewer.draw(); }
  };

  viewer.callback_pre_draw = [&](igl::opengl::glfw::Viewer& viewer)->bool
  {
    // TODO: implement simplification
    // if (simpMode == ISIMP_MODE::EACH_FRAME && optimise && !Adata.hasConverged)
    // {
    //   // Simplification
    //   std::cout << "simplification iteration: " << Adata.iteration;
    //   // TODO: reduce
    //   // std::pair<double, double> optEnergies = asdap_energies(Adata, U);
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
    // if (optimise && Adata.hasConverged)
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
      // Adata.changedV[selected_ids(i)] = true; // TODO: implement constraint forwarding to ISIMP
    }
    constraints = igl::slice(all_constraints, bi, 1);
    last_mouse = drag_mouse;
    update_points();
    viewer.draw();
    return false;
  };



  viewer.callback_key_pressed = [&](decltype(viewer) &,unsigned int key, int mod)
  {
    std::pair<double, double> optEnergies;
    std::pair<Eigen::MatrixXd, double> result;
    std::chrono::system_clock::time_point timestamp;
    std::time_t time;
    std::string filename;
    double theta = 0;
    double U22, U21;

    switch(key)
    {
      case ' ':
        if (virgin) return true;
        // if (Adata.hasConverged) return true; // TODO: replace
        // step manually
        if (simpMode == ISIMP_MODE::ON_INPUT)
        {
          // make simplification step
          // TODO: implement simplification step

          // Logging
          // TODO: implement logging
          // optEnergies = asdap_energies(Adata, U);
          // if(doLogging)
          // {
          //   logFile << std::to_string(optEnergies.first + optEnergies.second) << ";" << std::to_string(optEnergies.first) << ";" << std::to_string(optEnergies.second) << "\n";
          // }
          // std::cout << " => Energy: " << (optEnergies.first + optEnergies.second) << " (Scaling: " << optEnergies.first << ", Rotation: " << optEnergies.second << ")" << std::endl;
          // if (Adata.hasConverged)
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
        // TODO: reserved for test examples
        return true;
      case 'E': case 'e':
        // TODO: reserved for test examples
        return true;
      case 'J': case 'j':
        simpMode = (simpMode == ISIMP_MODE::ON_INPUT) ? ISIMP_MODE::EACH_FRAME : ISIMP_MODE::ON_INPUT;
        return true;
      case 'Y': case 'y':
        (logFile).close();
        return true;
      case 'C': case 'c':
      // TODO: we don't need constraints
      //   virgin = false;
      //   Adata.hasConverged = false;
      //   constraints_mask = ((constraints_mask+selected_mask)>0).cast<int>();
      //   update_constraints(true);
      //   refresh_mesh();
      //   // std::cout << constraints << "\n\n\n";
      //   // std::cout << bi << "\n\n\n";
      //   return true;
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
        if (mode == 0)
        {
          viewer.data().set_texture(R,G,B,A);
        } else {
          Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> R,G,B,A;
          igl::png::readPNG(argc > 3 ? argv[3] : "../textures/texture.png", R,G,B,A);
          viewer.data().set_texture(R,G,B,A);
        }
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
  init_viewer_data();
  viewer.core().is_animating = true;
  viewer.core().background_color.head(3) = CM.row(0).head(3).cast<float>();
  // update_points()
  viewer.launch();
}
