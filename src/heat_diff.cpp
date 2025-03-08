#include "heat_diff.hpp"

#include <random>

void initHeat(heatDiffData& hdData, const iSimpData& iSData)
{
    int n = iSData.inputMesh->nVertices();
    hdData.initial_heat = Eigen::MatrixXd::Zero(n);

    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    for (int i=0; i<n; i++)
    {
        hdData.initial_heat(i, 0) = distribution(generator);
    }
}

void diffuse_heat(heatDiffData& hdData, const double stepsize, const int steps)
{
    // TODO: check if Laplacian matrices already have been calculated in heatDiffData
    Eigen::MatrixXd A, B;
    A = hdData.M - stepsize * hdData.L;

    hdData.final_heat = hdData.initial_heat;
    for (int i=0; i<steps; i++)
    {
        // TODO: in the case for a simplified mesh, we need to only consider a subset of the initial heat, i.e. select those entries into a new matrix
        B = hdData.M * hdData.final_heat;
        hdData.final_heat = A.ldlt().solve(B);
    }
}

void build_cotan_mass(heatDiffData& hdData, iSimpData& iSData)
{
    // TODO: generalize this for multiple different vertexcounts
    int n = hdData.initial_heat.rows();
    hdData.M = Eigen::MatrixXd::Zero(n, n);
    for (int i=0; i < n; i++)
    {
        // TODO: replace with cotan mass
        hdData.M(i, i) = barycentric_lumped_mass(iSData, i);
    }
}

void build_cotan_laplacian(heatDiffData& hdData, iSimpData& iSData)
{

}

double barycentric_lumped_mass(iSimpData& iSData, const int vertex_idx)
{
    if (iSData.intrinsicMesh->nVertices() > vertex_idx) { std::cout << RED << "passed out of bounds vertex idx to barycentric_lumped_mass!" << RESET << std::endl; exit(-1); }
    double totalArea = 0.0;
    for (gcs::Face F : iSData.intrinsicMesh->vertex(vertex_idx).adjacentFaces())
    {
        int i = 0;
        double l_ij, l_ik, l_jk;
        for (gcs::Edge E : F.adjacentEdges())
        {
            if (i == 0) { l_ij = iSData.L(E.getIndex(), 0); }
            else if (i == 1) { l_ik = iSData.L(E.getIndex(), 0); }
            else if (i == 2) { l_jk = iSData.L(E.getIndex(), 0); }
        }
        totalArea += triangle_area_from_lengths(l_ij, l_ik, l_jk);
    }
    return totalArea / 3.f;
}

double cotan_weight(const iSimpData& iSData, const int edge_idx)
{
    gcs::Edge edge = iSData.intrinsicMesh->edge(edge_idx);
    double deg_epsilon = 0.000001*M_PI/180.0;
    double weight = 0.0;

    gcs::Face Fk;
    gcs::Face Fl;
    int faces = 0;
    for (gcs::Face F : edge.adjacentFaces())
    {
      if (faces == 0)
      {
        Fk = F;
        faces++;
      }
      else {
        Fl = F;
        faces++;
      }
    }
    if (faces < 2) { return 1.0; }

    std::pair<int, int> Fk_unique = find_unique_vertex_index(Fk, edge);
    std::pair<int, int> Fl_unique = find_unique_vertex_index(Fl, edge);
    std::array<int, 5> edge_indices = order_quad_edge_indices(edge, Fk_unique.second, Fl_unique.second);
    int l_ij = iSData.L(edge_indices[0]);

    if (faces >= 1)
    {
        int l_ik = iSData.L(edge_indices[1]);
        int l_jk = iSData.L(edge_indices[2]);
        double alpha = angle_i_from_lengths(l_ij, l_ik, l_jk);
        if (alpha > 2*M_PI) { std::cout << RED << "cotan_weight encountered angle larger than 2PI or equivalently 180 degrees!" << RESET << std::endl; exit(-1); }

        if (alpha <= M_PI/2.0)
        {
            weight += (alpha > deg_epsilon) ? (cos(alpha)/sin(alpha)) : l_ik;
        }
        if (alpha > M_PI/2.0)
        {
            weight += (alpha < M_PI - deg_epsilon) ? -(cos(alpha)/sin(alpha)) : -l_ik;
        }
    }

    if (faces >= 2)
    {
        int l_il = iSData.L(edge_indices[3]);
        int l_jl = iSData.L(edge_indices[4]);
        double beta = angle_i_from_lengths(l_ij, l_il, l_jl);
        if (beta > 2*M_PI) { std::cout << RED << "cotan_weight encountered angle larger than 2PI or equivalently 180 degrees!" << RESET << std::endl; exit(-1); }

        if (beta <= M_PI/2.0)
        {
            weight += (beta > deg_epsilon) ? (cos(beta)/sin(beta)) : l_il;
        }
        if (beta > M_PI/2.0)
        {
            weight += (beta < M_PI - deg_epsilon) ? -(cos(beta)/sin(beta)) : -l_il;
        }
    }

    return weight / 2.0;
}