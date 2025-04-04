#include <igl/AABB.h>
#include <igl/readOBJ.h>
#include <igl/writePLY.h>

#include <Eigen/Eigenvalues>
#include "isimp.hpp"

int main(int argc, char *argv[])
{
  // load mesh
  
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  Eigen::MatrixXd U;
  std::cout << "Read mesh" << std::endl;
  igl::readOBJ("../meshes/horse.obj", V, F);
  U = V;

  int max_iter = 10000;

  //add constraints
  std::cout << "Define constraints" << std::endl;
  Eigen::VectorXi constraints_indices(59);
  constraints_indices << 45, 103, 177, 204, 283, 312, 343, 408, 465, 507, 601, 631, 827, 868, 1494, 1803, 1896, 2067, 2902, 3429, 3492, 3567, 3570, 3612, 3661, 4127, 4394, 4442, 4465, 4515, 4544, 4870, 4886, 4900, 4901, 4924, 4928, 4931, 4932, 4938, 4939, 4941, 4945, 4952, 4957, 4959, 4962, 4964, 4969, 4974, 4977, 4978, 4980, 4982, 4983, 4984, 4985, 4987, 4988;
  Eigen::MatrixXd constraints(59,3);
  constraints << 0.241461, -0.277445,  -0.42548, 
                -0.0820156, -0.392467, -0.472385, 
                0.226665, -0.281277,  -0.37685, 
                -0.0597339, -0.392607, -0.476943, 
                -0.044731, -0.759223,  0.805673, 
                0.264775, -0.275534,  -0.42396, 
                0.246686, -0.278313,  -0.40283, 
                -0.0450468, -0.391733, -0.463577, 
                0.320483, -0.746594,  0.618826, 
                -0.0533566, -0.758245,  0.840553, 
                0.347788, -0.746674,  0.631189, 
                0.28076,  0.174832,  -0.54724, 
                -0.048603, -0.757855,  0.853794, 
                -0.0501507, -0.759285,  0.805814, 
                0.22085, -0.279699, -0.417232, 
                -0.0327781, -0.758129,  0.857802, 
                0.220514, -0.280898, -0.389586, 
                -0.0456697, -0.392223, -0.479773, 
                0.224928, -0.278192, -0.447134, 
                -0.0926333, -0.392652, -0.494451, 
                0.340997, -0.746671,  0.654024, 
                0.302125, -0.746488,  0.635107, 
                -0.0698434, -0.393035, -0.520588, 
                0.269187,  0.148073, -0.546521, 
                0.266767, -0.276607, -0.403145, 
                -0.101492, -0.389092,  -0.52791, 
                0.249272,  0.152527, -0.546577, 
                0.328425, -0.746918, 0.63969, 
                0.24503,  0.171966, -0.552968, 
                -0.0598982, 0.0429398, -0.541973, 
                -0.0788974, 0.0262442, -0.542602, 
                0.341203, -0.746608,  0.612959, 
                0.325149,  -0.74425,  0.680691, 
                0.323408, -0.746312,  0.590353, 
                0.282271, -0.745395,  0.626948, 
                0.301812, -0.741624,  0.673186, 
                0.310509, -0.739648,  0.572693, 
                0.287972, -0.743229,  0.656179, 
                0.348575, -0.740295,  0.571123, 
                0.245492, -0.275553, -0.459654, 
                0.284876, -0.745289,  0.602403, 
                -0.1132, -0.389669, -0.456985, 
                0.198556, -0.276806, -0.408464, 
                0.240477, -0.280049, -0.377989, 
                -0.0170756, -0.758117,  0.831104, 
                -0.119904, -0.388288, -0.477635, 
                0.204259, -0.278105, -0.434468, 
                -0.0368408, -0.758565,  0.84483, 
                -0.115442, -0.388172,  -0.51278, 
                -0.0488478,  -0.38663, -0.538671, 
                -0.0464777, -0.759305,  0.812379, 
                -0.0310131, -0.758796,  0.815747, 
                -0.0861562, -0.386498, -0.541731, 
                -0.0859028, -0.391868,  -0.53696, 
                -0.0346011,  -0.38839, -0.523012, 
                -0.0557857, -0.392328, -0.501448, 
                -0.0657594, -0.392396, -0.453397, 
                -0.0639448, -0.390876,  -0.54407, 
                -0.0865047, -0.392429, -0.456581;

  // TODO: implement iSIMP (intrinsic-simplification)
  /* // precomputation
  iSimpData data;
  std::cout << "Precompute matrix decomposition" << std::endl;
  asdap_precompute(constraints_indices, data);
  Eigen::MatrixXd U;
  U = V;
  std::cout << "Do 50 solve steps" << std::endl;
  // do 50 solving steps
  for (size_t i = 0; i < 50; i++) {
    asdap_solve(constraints, data, U);
  }
  std::cout << "Write result to output.ply" << std::endl;
  igl::writePLY("output.ply", U, F); */

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

  // precomputation
  iSimpData Adata;
  initMesh(V, F, Adata);
  igl::parallel_for(3, [&](const int i) {
    // TODO: asdap_set_constraints(constraints_indices, constraints, Adata);
    update_mesh_points(constraints_indices, constraints);
  }, 1);

  // TODO: implement simplification
  std::cout << "Writing result to output.ply" << std::endl;
  igl::writePLY("output.ply", U, F);
}