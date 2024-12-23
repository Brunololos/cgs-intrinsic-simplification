/* // Used to test edge flips
Eigen::MatrixXi F_new = Eigen::MatrixXi(F.rows(), F.cols());
int i = 0;
gcs::Edge eddge = iSData.inputMesh->edge(12);
std::array<int, 5> fp_indices = order_face_pair_indices(eddge);
std::cout << fp_indices[0] << ", " << fp_indices[1] << ", " << fp_indices[2] << ", " << fp_indices[3] << ", " << fp_indices[4] << std::endl;
for (int i=0; i<5; i++)
{
  std::cout << "Edge " << fp_indices[i] << ": (" << iSData.inputMesh->edge(fp_indices[i]).firstVertex().getIndex() << ", " << iSData.inputMesh->edge(fp_indices[i]).secondVertex().getIndex() << ")" << std::endl;
}
flip_intrinsic(iSData.inputMesh, iSData.L, eddge);
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
refresh_mesh();
std::cout << "Edge " << fp_indices[0] << " length: " << iSData.L[fp_indices[0]] << std::endl;
*/