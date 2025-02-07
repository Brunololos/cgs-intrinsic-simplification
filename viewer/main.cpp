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
  // if (argc < 2)
  // {
  //   std::cout << "Please pass an input mesh as first argument" << std::endl;
  //   exit(1);
  // }
  std::string filename = "../../meshes/spot.obj";
  int coarsen_to_n_vertices;
  int min_n_vertices = 0;
  bool using_texture = true;

  if (argc > 1) { filename = argv[1]; }
  // if (argc > 2) { coarsen_to_n_vertices = std::to_integer(string(argv[2])); }

  Eigen::MatrixXd V, UV, NV;
  Eigen::MatrixXi F, UF, NF;
  MatrixXuc R, G, B, A;
  int TEXTURE_WIDTH = 1000;
  int TEXTURE_HEIGHT = 1000;
  std::vector<unsigned char> polyscope_texture = std::vector<unsigned char>();
  polyscope::SurfaceCornerParameterizationQuantity* q = nullptr;
  polyscope::SurfaceMesh* psCoarseMesh;
  polyscope::SurfaceMesh* psMesh;

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
  igl::readOBJ(argv[1], V, UV, NV, F, UF, NF);
  coarsen_to_n_vertices = V.rows();

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
        // R(texcoord[0], texcoord[1]) = CM((i % num_colors) + 1, 0) * 255;
        // G(texcoord[0], texcoord[1]) = CM((i % num_colors) + 1, 1) * 255;
        // B(texcoord[0], texcoord[1]) = CM((i % num_colors) + 1, 2) * 255;
        polyscope_texture[4*(texcoord[0] + TEXTURE_WIDTH*texcoord[1]) + 0] = CM((i % num_colors) + 1, 0) * 255;
        polyscope_texture[4*(texcoord[0] + TEXTURE_WIDTH*texcoord[1]) + 1] = CM((i % num_colors) + 1, 1) * 255;
        polyscope_texture[4*(texcoord[0] + TEXTURE_WIDTH*texcoord[1]) + 2] = CM((i % num_colors) + 1, 2) * 255;
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
      map_registered(iSData, coarsen_to_n_vertices);
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
  // // THese are inserted:
  map_registered(iSData);
  // update_texture();
  // std::cout << "write_success: " << igl::png::writePNG(R, G, B, A, "H:/GIT/cgs-intrinsic-simplification/textures/object_texture.png") << std::endl;
  // // end insertion
  init_viewer_data();
  update_texture();
  // viewer.core().is_animating = true;
  // viewer.core().background_color.head(3) = CM.row(0).head(3).cast<float>();
  // // update_points()
  // viewer.launch();

  // while (!iSData.hasConverged && iSData.intrinsicMesh->nVertices() > coarsen_to_n_vertices)
  while (!iSData.hasConverged)
  {
    // make simplification step
    iSimp_step(iSData);
  }
  // refresh_coarse_mesh(); //TODO: update coarse mesh
  min_n_vertices = iSData.intrinsicMesh->nVertices();

  // map_registered(iSData, iSData.inputMesh->nVertices());
  // update_texture();

  polyscope::init();
  // psCoarseMesh  = polyscope::registerSurfaceMesh("coarse mesh", U, H); // TODO: visualize coarse mesh
  psMesh  = polyscope::registerSurfaceMesh("input mesh", V, F);

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