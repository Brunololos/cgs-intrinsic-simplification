#include "query_texture_barycentrics.hpp"

void query_texture_barycentrics(
    const Eigen::MatrixXd& UV,
    const Eigen::MatrixXi& UF,
    const int& TEXTURE_WIDTH,
    const int& TEXTURE_HEIGHT,
    iSimpData& iSData)
{
    igl::AABB<Eigen::MatrixXd, 2> tree;
    tree.init(UV, UF);

    int TEXTURE_SIZE = TEXTURE_WIDTH * TEXTURE_HEIGHT;

    // TODO: setup fields

    for (int j = 0; j < TEXTURE_HEIGHT; j++) {
        for (int i = 0; i < TEXTURE_WIDTH; i++) {
            int coordinate = j*TEXTURE_WIDTH + i;

            double u = ((double) i + 0.5) / (double) TEXTURE_WIDTH;
            double v = ((double) j + 0.5) / (double) TEXTURE_HEIGHT;
            Eigen::RowVector2d uv = Eigen::RowVector2d(u, v);

            int f;
            Eigen::RowVector2d C;
            double distance = tree.squared_distance(UV, UF, uv, f, C);
            if (distance < 1e-4) {
                Eigen::Vector3d bc = to_barycentric(C, UV.row(UF(f, 0)), UV.row(UF(f, 1)), UV.row(UF(f, 2)));

                TexCoord tc = TexCoord(i, j);
                register_point(iSData, bc, tc, f);
            }
        }

    }
}