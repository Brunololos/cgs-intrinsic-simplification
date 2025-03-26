#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/surface_parameterization_quantity.h"

#include <stdlib.h>
#include <igl/read_triangle_mesh.h>
#include <stdlib.h>
#include <chrono>
#include <ctime>

#include "isimp.hpp"
#include "heat_diff.hpp"
#include "screenshot.hpp"
#include "lscm_wrapper.hpp"
#include "query_texture_barycentrics.hpp"

// TODO: move this somewhere else
enum class ISIMP_MODE
{
  ON_INPUT,
  EACH_FRAME
};

int main(int argc, char *argv[])
{
  // PROGRAM SETTINGS
  std::string filename = "../../meshes/spot.obj";
  int TEXTURE_WIDTH = 2000;
  int TEXTURE_HEIGHT = 2000;
  int snapshot_interval = 50;
  bool verbose_simplification = false;

  // UI SETTINGS
  int coarsen_to_n_vertices;
  int diffuse_over_n_steps;
  double diffusion_stepsize;
  double min_diff_stepsize = 0.0001;
  double max_diff_stepsize = 0.01;

  // UI TRIGGERS/STATE
  bool do_intrinsics_recovery = true;
  bool do_diffuse = false;
  bool do_triang_texture_update = false;
  bool do_heat_texture_update = false;

  bool is_triang_map_up_to_date = true;
  bool is_heat_map_up_to_date = false;
  bool is_heat_map_initialized = false;

  // POLYSCOPE
  polyscope::SurfaceMesh* psMesh;
  polyscope::SurfaceMesh* psCoarseMesh;
  polyscope::SurfaceVertexScalarQuantity* g = nullptr;
  polyscope::SurfaceCornerParameterizationQuantity* triang_map = nullptr;
  polyscope::SurfaceCornerParameterizationQuantity* heat_map = nullptr;


  // TEXTURING
  std::vector<unsigned char> ps_triang_texture = std::vector<unsigned char>();
  std::vector<unsigned char> ps_heat_texture = std::vector<unsigned char>();

  // store all colorings for all different simplification states (for each vertex removal).
  std::vector<Eigen::VectorXi> colorings;

  // SIMPLIFICATION
  Eigen::MatrixXd V, UV, NV;
  Eigen::MatrixXi F, UF, NF;
  int min_n_vertices = 0;

  // SIMPLIFICATION - SNAPSHOTS
  // points tracked by triangle for each simplification state (I did this just to make the polyscope slider movement smoother, but this takes a lot of storage.)
  std::vector<std::vector<std::vector<int>>> mapped_by_triangle_snapshots;
  std::vector<std::vector<BarycentricPoint>> mapped_points_snapshots;
  std::vector<int> snapshot_mapping_indices;
  std::vector<int> simp_step_mapping_indices;
  // Eigen::VectorXd current_heat;

  // RUDIMENTARY ARGUMENT PARSING
  int arg_idx = 1;
  std::string mesh_arg = "";
  std::string tex_arg = "";
  std::string snap_arg = "";
  std::string verb_arg = "";

  // allow custom order args via "[arg_name]=[arg]"
  while (arg_idx < argc)
  {
    std::string arg = std::string(argv[arg_idx]);
    int eq_idx = arg.find_first_of("=");

    // argument not specified
    if (eq_idx == -1)
    {
      // std::cout << "no arg specifier provided for argument " << arg_idx << ":" << std::endl;
      switch(arg_idx)
      {
        default:
        case 1:
          // std::cout << "-> assuming " << arg << " is the mesh." << std::endl;
          mesh_arg = arg;
          break;
        case 2:
          // std::cout << "-> assuming " << arg << " is the texture size." << std::endl;
          tex_arg = arg;
          break;
        case 3:
          // std::cout << "-> assuming " << arg << " is the snapshot interval." << std::endl;
          snap_arg = arg;
          break;
        case 4:
          // std::cout << "-> assuming " << arg << " is the verbosity." << std::endl;
          verb_arg = arg;
          break;
      }
    } else {
      std::string arg_name = arg.substr(0, eq_idx);
      if (arg_name == "mesh") {
        mesh_arg = arg.substr(eq_idx+1);
      } else if (arg_name == "tex") {
        tex_arg = arg.substr(eq_idx+1);
      } else if (arg_name == "snap") {
        snap_arg = arg.substr(eq_idx+1);
      } else if (arg_name == "verb") {
        verb_arg = arg.substr(eq_idx+1);
      } else {
        std::cout << "Failed to parse argument!" << std::endl;
      }
    }
    arg_idx++;
  }

  // setting args
  int dot_pos = mesh_arg.size() - 4;
  if (mesh_arg != "" && (mesh_arg.substr(dot_pos).compare(".obj")))
  {
    std::cout << "Please pass a .obj mesh as the first argument" << std::endl;
    exit(-1);
  }
  if (mesh_arg != "") { filename = mesh_arg; }
  if (tex_arg != "") {
    TEXTURE_WIDTH = std::stoi(tex_arg);
    TEXTURE_HEIGHT = TEXTURE_WIDTH;
  }
  if (snap_arg != "") { snapshot_interval = std::stoi(snap_arg); }
  if (verb_arg != "") { verbose_simplification = (verb_arg == "verbose"); }


  igl::readOBJ(filename, V, UV, NV, F, UF, NF);

  if (UF.rows() <= 0)
  {
    std::cout << "Please pass a .obj mesh with texture coordinates included." << std::endl;
    exit(-1);
  }

  coarsen_to_n_vertices = V.rows();
  diffuse_over_n_steps = 0;
  diffusion_stepsize = 0.005;
  colorings = std::vector<Eigen::VectorXi>();
  mapped_by_triangle_snapshots = std::vector<std::vector<std::vector<int>>>();
  mapped_points_snapshots = std::vector<std::vector<BarycentricPoint>>();
  snapshot_mapping_indices = std::vector<int>();
  simp_step_mapping_indices = std::vector<int>();

  // modified vertices & faces
  Eigen::MatrixXd U = V;
  Eigen::MatrixXi H = F;

  // settings
  int mode = 0;
  bool virgin = true;
  bool view_constraints = true;
  bool view_all_points = false;
  int selected_cmap = 9;
  const char* cmaps[]{"coolwarm", "reds", "blues", "viridis", "pink-green", "phase", "spectral", "rainbow", "jet", "turbo"};
  bool optimise = false;

  // ISIMP (intrinsic-simplification)
  ISIMP_MODE simpMode = ISIMP_MODE::ON_INPUT;

  iSimpData iSData;
  initMesh(V, F, iSData);

  // HEAT DIFFUSION
  heatDiffData hdData;
  initHeatDiff(hdData, iSData);

  // randomize seed
  srand(time(NULL));

  // COLORING/VISUALIZATION
  Eigen::MatrixXd CM = (Eigen::MatrixXd(16, 3) <<
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
  ).finished();

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

  polyscope::init();
  std::cout << RESET << R"(  Intrinsic Simplification Settings:
    Mesh:                       )" << filename << R"(
    texture size:               )" << TEXTURE_WIDTH << "x" << TEXTURE_HEIGHT << R"(px
    snapshot interval:          every )" << snapshot_interval << R"( vertex removals
    verbose simplification:     )" << (verbose_simplification ? "true" : "false") << R"(
  )" << std::endl;

  // NOTE: This is only for the polyscope VertexScalarQuantity. Not needed for the heat texture mapping.
  // auto update_heat = [&]()
  // {
  //   current_heat = Eigen::VectorXd::Zero(iSData.inputMesh->nVertices());
  //   for (int i=0; i<hdData.heat_sample_vertices.size(); i++)
  //   {
  //     int vertex_idx = hdData.heat_sample_vertices[i];
  //     current_heat(vertex_idx, 0) = get_heat(hdData, vertex_idx);
  //   }
  // };

  // NOTE: this approach assumes, that all triangles lie completely within the texture. (No triangles wrap around the outside.)
  // I note this here, because lscm sometimes yields uv-coordinates that exceed the [0, 1] bounds. (They can be negative and their absolute value can exceed 1.)
  // To fix this uv's can be translated & rescaled to be withing [0, 1]. Because we don't have a set texture and build it ourselves, this is fine.
  auto register_texels = [&]()
  {
    for (gcs::Face F : iSData.inputMesh->faces())
    {
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
            BarycentricPoint bp = to_barycentric(Point2D(i, j), t0, t1, t2);
            TexCoord tc = TexCoord(i, j);
            register_point(iSData, bp, tc, F.getIndex());
          }
        }
      }
    }
  };

  auto init_viewer_data = [&]()
  {
    ps_triang_texture.resize(4 * TEXTURE_WIDTH * TEXTURE_HEIGHT);
    ps_heat_texture.resize(4 * TEXTURE_WIDTH * TEXTURE_HEIGHT);
  };

  // Selection
  bool translate = false;
  RAXIS rot_axis = RAXIS::YAW;
  bool only_visible = false;
  Eigen::Array<double, Eigen::Dynamic, 1> W = Eigen::Array<double, Eigen::Dynamic, 1>::Zero(V.rows());
  Eigen::Array<double, Eigen::Dynamic, 1> and_visible = Eigen::Array<double, Eigen::Dynamic, 1>::Zero(V.rows());

  Eigen::MatrixXd all_constraints = V;

  Eigen::Array<int, -1, 1> constraints_mask = Eigen::Array<int, -1, 1>::Zero(V.rows());
  Eigen::Array<int, -1, 1> selected_mask = Eigen::Array<int, -1, 1>::Zero(V.rows());
  Eigen::VectorXi selected_ids;

  Eigen::MatrixXd constraints;
  Eigen::VectorXi bi;
  Eigen::VectorXi indices = Eigen::VectorXd::LinSpaced(V.rows(), 0.0, V.rows() - 1).cast<int>();

  const auto build_F_snapshot = [&](const iSimpData& iSData)
  {
    std::vector<std::array<int, 3>> snapshot = std::vector<std::array<int, 3>>();

    int n = iSData.inputMesh->nFaces();
    for (int i=0; i<n; i++)
    {
      gcs::Face F = iSData.intrinsicMesh->face(i);
      if(!F.isDead())
      {
        int v1_idx, v2_idx, v3_idx;
        int k = 0;
        for (gcs::Vertex V : F.adjacentVertices())
        {
          if (k==0) { v1_idx = V.getIndex(); }
          else if (k==1) { v2_idx = V.getIndex(); }
          else if (k==2) { v3_idx = V.getIndex(); }
          k++;
        }
        std::array<int, 3> idcs = std::array<int, 3>();
        idcs[0] = v1_idx;
        idcs[1] = v2_idx;
        idcs[2] = v3_idx;
        snapshot.push_back(idcs);
      } else {
        snapshot.push_back(std::array<int, 3>({-1, -1, -1}));
      }
    }

    return snapshot;
  };

  // TODO: If I wanted I could implement a visualization of the tangent space reference edges and/or the error vectors
  const auto update_edges = [&]()
  {

  };

  auto compute_coloring = [&]() {
    int num_colors = CM.rows() - 1;
    Eigen::VectorXi Coloring = Eigen::VectorXi::Zero(iSData.F.rows());
    for (int i = 0; i < iSData.inputMesh->nFaces(); i++)
    {
      // Determine color
      int seeded_color_idx = (i % (num_colors - 1));
      gcs::Face current_face = iSData.intrinsicMesh->face(i);
      if (current_face.isDead()) { continue; }

      int color_idx = 0;
      for (int j=0; j<num_colors-1; j++)
      {
        bool valid = true;
        int current_color_idx = ((seeded_color_idx + j) % (num_colors - 1)) + 1;
        for(gcs::Face neighbor : current_face.adjacentFaces())
        {
          int neighbor_idx = neighbor.getIndex();
          int neighbor_color_idx = Coloring(neighbor_idx);
          if (current_color_idx == neighbor_color_idx)
          {
            valid = false;
            break;
          }
        }
        if (valid) { color_idx = current_color_idx; break; }
      }

      Coloring(i) = color_idx;
    }
    colorings.push_back(Coloring);
  };

  const auto &reset_triang_texture = [&]()
  {
    for(int i=0; i<TEXTURE_WIDTH; i++)
    {
      for(int j=0; j<TEXTURE_HEIGHT; j++)
      {
        ps_triang_texture[4*(i + TEXTURE_WIDTH*j) + 0] = 0;
        ps_triang_texture[4*(i + TEXTURE_WIDTH*j) + 1] = 0;
        ps_triang_texture[4*(i + TEXTURE_WIDTH*j) + 2] = 0;
        ps_triang_texture[4*(i + TEXTURE_WIDTH*j) + 3] = 255;
      }
    }
  };

  const auto &update_triang_texture = [&]()
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
        TexCoord texcoord = iSData.tracked_texcoords[original_idx];
        ps_triang_texture[4*(texcoord[0] + TEXTURE_WIDTH*texcoord[1]) + 0] = CM(colorings[n_removals](i), 0) * 255;
        ps_triang_texture[4*(texcoord[0] + TEXTURE_WIDTH*texcoord[1]) + 1] = CM(colorings[n_removals](i), 1) * 255;
        ps_triang_texture[4*(texcoord[0] + TEXTURE_WIDTH*texcoord[1]) + 2] = CM(colorings[n_removals](i), 2) * 255;
        ps_triang_texture[4*(texcoord[0] + TEXTURE_WIDTH*texcoord[1]) + 3] = 255;
      }
    }
  };

  const auto &update_heat_texture = [&]()
  {
    int num_triangles = iSData.mapped_by_triangle.size();

    // draw texture
    for (int i = 0; i < num_triangles; i++)
    {
      if (iSData.recoveredMesh->face(i).isDead()) { continue; }
      int v1_idx, v2_idx, v3_idx;
      int k = 0;
      for (gcs::Vertex V : iSData.recoveredMesh->face(i).adjacentVertices())
      {
        if (k==0) { v1_idx = V.getIndex(); }
        else if (k==1) { v2_idx = V.getIndex(); }
        else if (k==2) { v3_idx = V.getIndex(); }
        k++;
      }
      double heat1 = get_heat(hdData, v1_idx);
      double heat2 = get_heat(hdData, v2_idx);
      double heat3 = get_heat(hdData, v3_idx);
      for (int j = 0; j < iSData.mapped_by_triangle[i].size(); j++)
      {
        int original_idx = iSData.mapped_by_triangle[i][j];
        TexCoord texcoord = iSData.tracked_texcoords[original_idx];

        // calculate heat value via interpolation
        BarycentricPoint p = iSData.mapped_points[original_idx];
        double psum = p[0] + p[1] + p[2];
        if (psum > 1.0) { p = p / psum; }

        double heat = heat1*p[0] + heat2*p[1] + heat3*p[2];
        heat = std::max(std::min(heat, 1.0), 0.0);
        if (p[0] + p[1] + p[2] > 1.0) { std::cout << "encountered excessive barycentric coordinates: (" << p[0] << ", " << p[1] << ", " << p[2] << ") => " << p[0] + p[1] + p[2] << std::endl; }
        glm::vec3 C;
        polyscope::render::ValueColorMap cmap = polyscope::render::engine->getColorMap(cmaps[selected_cmap]);
        C = cmap.getValue(heat);

        ps_heat_texture[4*(texcoord[0] + TEXTURE_WIDTH*texcoord[1]) + 0] = (C.r) * 255;
        ps_heat_texture[4*(texcoord[0] + TEXTURE_WIDTH*texcoord[1]) + 1] = (C.g) * 255;
        ps_heat_texture[4*(texcoord[0] + TEXTURE_WIDTH*texcoord[1]) + 2] = (C.b) * 255;
        ps_heat_texture[4*(texcoord[0] + TEXTURE_WIDTH*texcoord[1]) + 3] = 255;
      }
    }
  };

  const auto &update_textures = [&]()
  {
    update_triang_texture();
    update_heat_texture();
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
    }
  };

  // update entire mesh
  auto refresh_coarse_mesh = [&]()
  {
    // if (psCoarseMesh != nullptr)
    // {
    //   // polyscope::unregister...
    // }
    // TODO: update connectivity
    psCoarseMesh  = polyscope::registerSurfaceMesh("coarse mesh", U, H);
    init_viewer_data();
    update_textures();
  };

  // update mesh vertices
  auto refresh_mesh_vertices = [&]()
  {
    init_viewer_data();
    update_textures();
  };

  auto callback = [&]() {
    ImGui::Text("Intrinsic Simplification");
    ImGui::Separator();
    if(ImGui::SliderInt("remaining vertices", &coarsen_to_n_vertices, min_n_vertices, iSData.inputMesh->nVertices()))
    {
      do_intrinsics_recovery = true;
      if (triang_map->isEnabled()) {
        is_heat_map_up_to_date = false;
        do_triang_texture_update = true;
      } else if(heat_map->isEnabled()) {
        is_triang_map_up_to_date = false;
        do_diffuse = true;
        do_heat_texture_update = true;
      }
    }
    if (ImGui::Button("Add Vertex"))
    {
      coarsen_to_n_vertices = std::max(std::min(coarsen_to_n_vertices + 1, (int) iSData.inputMesh->nVertices()), min_n_vertices);
      do_intrinsics_recovery = true;
      if (triang_map->isEnabled()) {
        is_heat_map_up_to_date = false;
        do_triang_texture_update = true;
      } else if(heat_map->isEnabled()) {
        is_triang_map_up_to_date = false;
        do_diffuse = true;
        do_heat_texture_update = true;
      }
    }
    ImGui::SameLine();
    if (ImGui::Button("Remove Vertex"))
    {
      coarsen_to_n_vertices = std::max(std::min(coarsen_to_n_vertices - 1, (int) iSData.inputMesh->nVertices()), min_n_vertices);
      do_intrinsics_recovery = true;
      if (triang_map->isEnabled()) {
        is_heat_map_up_to_date = false;
        do_triang_texture_update = true;
      } else if(heat_map->isEnabled()) {
        is_triang_map_up_to_date = false;
        do_diffuse = true;
        do_heat_texture_update = true;
      }
    }
    ImGui::Dummy(ImVec2(0.0f, 10.0f));
    ImGui::Text("Heat Diffusion");
    ImGui::Separator();
    if(ImGui::SliderScalar("Diffusion Stepsize", ImGuiDataType_Double, &diffusion_stepsize, &min_diff_stepsize, &max_diff_stepsize))
    {
      do_diffuse = true;
      do_heat_texture_update = true;
    }
    if(ImGui::SliderInt("Diffusion Steps", &diffuse_over_n_steps, 0, 200))
    {
      do_diffuse = true;
      do_heat_texture_update = true;
    }
    if (ImGui::Button("Add Diffusion Step"))
    {
      diffuse_over_n_steps++;
      do_diffuse = true;
      do_heat_texture_update = true;
    }
    ImGui::SameLine();
    if (ImGui::Button("Remove Diffusion Step"))
    {
        diffuse_over_n_steps = std::max(diffuse_over_n_steps - 1, 0);
        do_diffuse = true;
        do_heat_texture_update = true;
    }
    if(ImGui::Combo("Heat color map", &selected_cmap, cmaps, 10))
    {
      do_heat_texture_update = true;
    }
    ImGui::Dummy(ImVec2(0.0f, 10.0f));
    ImGui::Text("Debug");
    ImGui::Separator();
    if (ImGui::Button("Print intrinsic triangulation"))
    {
      std::cout << "Remaining intrinsic triangles:" << std::endl;
        // assuming recoveredMesh ist recovered and consistent with mapped points
      for (gcs::Face F : iSData.recoveredMesh->faces())
      {
        std::cout << "> Triangle " << F.getIndex() << ": {";
        bool first = true;
        for (gcs::Vertex V : F.adjacentVertices())
        {
          if (!first) { std::cout << ", "; }
          std::cout << V.getIndex();
          first = false;
        }
        std::cout << "},      N: {";
        first = true;
        for (gcs::Face NF : F.adjacentFaces())
        {
          if (!first) { std::cout << ", "; }
          std::cout << NF.getIndex();
          if (first) { first = false; }
        }
        std::cout << "}" << std::endl;
      }
    }
    if (ImGui::Button("Print tracked points"))
    {
      std::cout << "Tracked points by triangle:" << std::endl;
      for (int i=0; i<iSData.tracked_by_triangle.size(); i++)
      {
        if (iSData.tracked_by_triangle[i].size() <= 0) { continue; }
        // assuming recoveredMesh ist recovered and consistent with mapped points
        gcs::Face F = iSData.inputMesh->face(i);
        std::cout << "> Triangle " << i << ": {";
        bool first = true;
        for (gcs::Vertex V : F.adjacentVertices())
        {
          if (!first) { std::cout << ", "; }
          std::cout << V.getIndex();
          first = false;
        }
        std::cout << "}:" << std::endl;

        for(int j=0; j<iSData.tracked_by_triangle[i].size(); j++)
        {
          BarycentricPoint bp = iSData.tracked_points[iSData.tracked_by_triangle[i][j]];
          std::cout << "  > Barycentric Point: {" << bp[0] << ", " << bp[1] << ", " << bp[2] << "}" << std::endl;
        }
      }
    }
    ImGui::SameLine();
    if (ImGui::Button("Print mapped points"))
    {
      std::cout << "Mapped points by triangle:" << std::endl;
      for (int i=0; i<iSData.mapped_by_triangle.size(); i++)
      {
        if (iSData.mapped_by_triangle[i].size() <= 0) { continue; }
        // assuming recoveredMesh ist recovered and consistent with mapped points
        gcs::Face F = iSData.recoveredMesh->face(i);
        std::cout << "> Triangle " << i << ": {";
        bool first = true;
        for (gcs::Vertex V : F.adjacentVertices())
        {
          if (!first) { std::cout << ", "; }
          std::cout << V.getIndex();
          first = false;
        }
        std::cout << "}:" << std::endl;

        for(int j=0; j<iSData.mapped_by_triangle[i].size(); j++)
        {
          BarycentricPoint bp = iSData.mapped_points[iSData.mapped_by_triangle[i][j]];
          std::cout << "  > Barycentric Point: {" << bp[0] << ", " << bp[1] << ", " << bp[2] << "}" << std::endl;
        }
      }
    }
    if (ImGui::Button("Print initial heat"))
    {
      std::cout << "Initial heat by vertex:" << std::endl;
      for (int i=0; i<hdData.final_heat.size(); i++)
      {
        // assuming recoveredMesh ist recovered and consistent with mapped points
        gcs::Face F = iSData.recoveredMesh->face(i);
        int vertex_idx = hdData.heat_sample_vertices[i];
        std::cout << "> Vertex " << vertex_idx << " heat: " << hdData.initial_heat[vertex_idx] << std::endl;
      }
    }
    ImGui::SameLine();
    if (ImGui::Button("Print heat"))
    {
      std::cout << "Current heat by vertex:" << std::endl;
      for (int i=0; i<hdData.final_heat.size(); i++)
      {
        // assuming recoveredMesh ist recovered and consistent with mapped points
        gcs::Face F = iSData.recoveredMesh->face(i);
        int vertex_idx = hdData.heat_sample_vertices[i];
        std::cout << "> Vertex " << vertex_idx << " heat: " << hdData.final_heat[i] << std::endl;
      }
    }
    if (!is_heat_map_initialized) { is_heat_map_initialized = true; do_heat_texture_update = true; }
    if (!is_triang_map_up_to_date && triang_map->isEnabled()) { do_triang_texture_update = true; }
    if (!is_heat_map_up_to_date && heat_map->isEnabled()) { do_diffuse = true; do_heat_texture_update = true; }
    if (do_intrinsics_recovery)
    {
      // load latest snapshot
      do_intrinsics_recovery = false;
      std::cout << "\nUpdate for " << coarsen_to_n_vertices << " remaining vertices:" << std::endl;
      int snapshot_idx = std::floor((iSData.inputMesh->nVertices() - coarsen_to_n_vertices) / snapshot_interval);
      int coarsen_from_n_vertices = iSData.inputMesh->nVertices() - snapshot_idx * snapshot_interval;
      iSData.mapped_by_triangle = mapped_by_triangle_snapshots.at(snapshot_idx);
      iSData.mapped_points = mapped_points_snapshots.at(snapshot_idx);

      // map from latest snapshot to desired simplification state & recover necessary information like connectivity & heat (those are replayed from the initial state as opposed to the latest snapshot)
      int snapshot_mapping_idx = snapshot_mapping_indices.at(snapshot_idx);
      if (simp_step_mapping_indices.size() > iSData.inputMesh->nVertices() - coarsen_to_n_vertices + 1)
      {
        int target_removal_idx = iSData.inputMesh->nVertices() - coarsen_to_n_vertices;
        int simp_step_mapping_idx = simp_step_mapping_indices.at(iSData.inputMesh->nVertices() - coarsen_to_n_vertices + 1);
        map_current_from_to(iSData, snapshot_mapping_idx, simp_step_mapping_idx);
        // get intrinsic connectivity/edge-lengths for heat diffusion
        recover_heat_and_intrinsics_at(simp_step_mapping_idx, hdData, iSData);
      } else {
        map_current_from(iSData, snapshot_mapping_idx);
        // get intrinsic connectivity/edge-lengths for heat diffusion
        recover_heat_and_intrinsics_at(iSData.mapping.size()-1, hdData, iSData);
      }
      build_heat_sample_indices(hdData, iSData);
      build_cotan_mass(hdData, iSData);
      build_cotan_laplacian(hdData, iSData);
      // if (g != nullptr) {
      //   g = psMesh->addVertexScalarQuantity("ps heat map", current_heat);
      // }
    }
    if (do_diffuse)
    {
      do_diffuse = false;
      diffuse_heat(hdData, diffusion_stepsize, diffuse_over_n_steps);
      // update_heat()
    }
    if (do_triang_texture_update)
    {
      do_triang_texture_update = false;
      update_triang_texture();
      if (triang_map != nullptr) {
        triang_map->setTexture(TEXTURE_WIDTH, TEXTURE_HEIGHT, ps_triang_texture, polyscope::TextureFormat::RGBA8);
        triang_map->setStyle(polyscope::ParamVizStyle::TEXTURE);
        triang_map->setCheckerSize(1);
      }
      is_triang_map_up_to_date = true;
    }
    if (do_heat_texture_update)
    {
      do_heat_texture_update = false;
      update_heat_texture();
      if (heat_map != nullptr) {
        heat_map->setTexture(TEXTURE_WIDTH, TEXTURE_HEIGHT, ps_heat_texture, polyscope::TextureFormat::RGBA8);
        heat_map->setStyle(polyscope::ParamVizStyle::TEXTURE);
        heat_map->setCheckerSize(1);
      }
      is_heat_map_up_to_date = true;
    }
  };

  // register initial snapshots & build textures for no simplification
  query_texture_barycentrics(UV, UF, TEXTURE_WIDTH, TEXTURE_HEIGHT, iSData);
  map_registered(iSData);
  mapped_by_triangle_snapshots.push_back(iSData.mapped_by_triangle);
  mapped_points_snapshots.push_back(iSData.mapped_points);
  snapshot_mapping_indices.push_back(0);
  simp_step_mapping_indices.push_back(0);

  init_viewer_data();
  compute_coloring();

  update_triang_texture();
  // NOTE: This is only for the polyscope VertexScalarQuantity. Not needed for the heat texture mapping.
  // update_heat();

  progressBar bar;
  if(!verbose_simplification) {
    std::cout << "\n\nPerforming Simplification:    ";
    start_progressBar(bar);
  } else {
    std::cout << "\n\nBeginning Simplification..." << std::endl;
  }
  while (!iSData.hasConverged)
  {
    int step_mapping_idx = iSData.mapping.size();
    // make simplification step
    bool success = iSimp_step(iSData, verbose_simplification);
    // compute coloring & save simplification state index for replays
    if(success) { simp_step_mapping_indices.push_back(step_mapping_idx); compute_coloring(); }

    // take snapshot
    if ((iSData.inputMesh->nVertices() - iSData.intrinsicMesh->nVertices()) % snapshot_interval == 0) {
      if(verbose_simplification) {
        std::cout << "Took snapshot at " << iSData.inputMesh->nVertices() - iSData.intrinsicMesh->nVertices() << " removed vertices" << std::endl;
      }

      // save texel mapping state for intrinsic triangulation texture mapping
      map_current_from(iSData, snapshot_mapping_indices.at(snapshot_mapping_indices.size() - 1));
      mapped_by_triangle_snapshots.push_back(iSData.mapped_by_triangle);
      mapped_points_snapshots.push_back(iSData.mapped_points);
      snapshot_mapping_indices.push_back(iSData.mapping.size());
    }
    if(!verbose_simplification)
    {
      double progress_percent = ((double) iSData.inputMesh->nVertices() - iSData.intrinsicMesh->nVertices()) / ((double) iSData.inputMesh->nVertices());
      progress_progressBar(bar, progress_percent);
    }
  }
  if(!verbose_simplification) { finish_progressBar(bar); }
  // refresh_coarse_mesh(); //TODO: update coarse mesh
  min_n_vertices = iSData.intrinsicMesh->nVertices();

  // psCoarseMesh  = polyscope::registerSurfaceMesh("coarse mesh", U, H); // TODO: visualize coarse mesh
  psMesh = polyscope::registerSurfaceMesh("input mesh", V, F);

  // convert parameterization to polyscope's desired input format
  Eigen::Matrix<glm::vec2, Eigen::Dynamic, 1> parameterization(3 * UF.rows());
  for (int iF = 0; iF < UF.rows(); iF++) {
    for (int iC = 0; iC < 3; iC++) {
      parameterization(3 * iF + iC) = glm::vec2{UV(UF(iF, iC), 0), UV(UF(iF, iC), 1)};
    }
  }
  triang_map = psMesh->addParameterizationQuantity("intrinsic triangulation", parameterization);

  triang_map->setTexture(TEXTURE_WIDTH, TEXTURE_HEIGHT, ps_triang_texture, polyscope::TextureFormat::RGBA8);
  triang_map->setEnabled(true);
  triang_map->setStyle(polyscope::ParamVizStyle::TEXTURE);
  triang_map->setCheckerSize(1);

  heat_map = psMesh->addParameterizationQuantity("heat map", parameterization);
  heat_map->setTexture(TEXTURE_WIDTH, TEXTURE_HEIGHT, ps_heat_texture, polyscope::TextureFormat::RGBA8);
  heat_map->setStyle(polyscope::ParamVizStyle::TEXTURE);
  heat_map->setCheckerSize(1);

  // g = psMesh->addVertexScalarQuantity("ps heat map", current_heat);

  polyscope::state::userCallback = callback;
  polyscope::show();
}