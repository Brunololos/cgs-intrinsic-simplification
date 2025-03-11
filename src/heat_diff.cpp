#include "heat_diff.hpp"

#include <random>

void initHeatDiff(heatDiffData& hdData, const iSimpData& iSData)
{
    initHeat(hdData, iSData);
    build_cotan_mass(hdData, iSData);
    build_cotan_laplacian(hdData, iSData);
}

void initHeat(heatDiffData& hdData, const iSimpData& iSData)
{
    int n = iSData.inputMesh->nVertices();
    hdData.initial_heat = Eigen::MatrixXd::Zero(n, 1);

    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    for (int i=0; i<n; i++)
    {
        hdData.initial_heat(i, 0) = distribution(generator);
    }
}

void diffuse_heat(heatDiffData& hdData, const double stepsize, const int steps)
{
    std::cout << "diffusing heat..." << std::endl;
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
    std::cout << "Successfully diffused heat!" << std::endl;
}

void build_cotan_mass(heatDiffData& hdData, const iSimpData& iSData)
{
    // TODO: generalize this for multiple different vertexcounts
    int n = hdData.initial_heat.rows();
    hdData.M = Eigen::MatrixXd::Zero(n, n);
    for (int i=0; i < n; i++)   // TODO: this iteration does not work for reduced meshes
    {
        hdData.M(i, i) = circumcentric_dual_mass(iSData, i);
        // std::cout << "vertex " << i << " area mass: " << hdData.M(i, i) << std::endl;
    }
    // std::cout << "computed total area mass: " << total_mass(hdData) << std::endl;
}

void build_cotan_laplacian(heatDiffData& hdData, const iSimpData& iSData)
{
    int n = iSData.intrinsicMesh->nVertices();  // TODO: this doesnt work for reduced meshes, because you can't iterate through
    hdData.L = Eigen::MatrixXd::Zero(n, n);

    for (int i=0; i<n; i++) // TODO: This iteration doesnt work for reduced meshes
    {
        double vertex_weight = 0.0;
        gcs::Vertex vi = iSData.intrinsicMesh->vertex(i);
        for (gcs::Vertex vj : vi.adjacentVertices())
        {
            int edge_idx = iSData.intrinsicMesh->connectingEdge(vi, vj).getIndex();
            double edge_weight = cotan_weight(iSData, edge_idx);
            int j = vj.getIndex();

            // std::cout << "determined connecting edge: " << edge_idx << ": (" << iSData.intrinsicMesh->edge(edge_idx).firstVertex() << ", " << iSData.intrinsicMesh->edge(edge_idx).secondVertex() << ") between the vertices: " << i << ", " << j << std::endl;
            if (i == j) { std::cout << RED << "build_cotan_laplacian: Encountered self-edge!" << RESET << std::endl; exit(-1); }
            hdData.L(i, j) = edge_weight;
            vertex_weight += edge_weight;
        }
        hdData.L(i, i) = -vertex_weight;
    }
}

double barycentric_lumped_mass(const iSimpData& iSData, const int vertex_idx)
{
    // if (iSData.intrinsicMesh->nVertices() <= vertex_idx) { std::cout << RED << "passed out of bounds vertex idx to barycentric_lumped_mass!" << RESET << std::endl; exit(-1); }
    double total_area = 0.0;
    for (gcs::Face F : iSData.intrinsicMesh->vertex(vertex_idx).adjacentFaces())
    {
        int i = 0;
        double l_ij, l_ik, l_jk;
        for (gcs::Edge E : F.adjacentEdges())
        {
            if (i == 0) { l_ij = iSData.L(E.getIndex(), 0); }
            else if (i == 1) { l_ik = iSData.L(E.getIndex(), 0); }
            else if (i == 2) { l_jk = iSData.L(E.getIndex(), 0); }
            i++;
        }
        total_area += triangle_area_from_lengths(l_ij, l_ik, l_jk);
    }
    return total_area / 3.f;
}

double circumcentric_dual_mass(const iSimpData& iSData, const int vertex_idx)
{
    double total_area = 0.0;
    gcs::Vertex V = iSData.intrinsicMesh->vertex(vertex_idx);
    double epsilon = std::numeric_limits<double>().epsilon();

    for (gcs::Face F : V.adjacentFaces())
    {
        CircumcentricShape calculation_case = CircumcentricShape::ZERO;
        double area_contribution;
        std::array<int, 3> face_edges = order_triangle_edge_indices(F, vertex_idx);

        int e_ij = face_edges[0];
        int e_ik = face_edges[1];
        int e_jk = face_edges[2];
        double l_ij = iSData.L(e_ij);
        double l_ik = iSData.L(e_ik);
        double l_jk = iSData.L(e_jk);
        double midpoint_l_ij = l_ij/2.0;
        double midpoint_l_ik = l_ik/2.0;

        // Swap so the first edge is always the shorter one
        if (midpoint_l_ij > midpoint_l_ik)
        {
            double tmp_len = l_ij;
            l_ij = l_ik;
            l_ik = tmp_len;

            tmp_len = midpoint_l_ij;
            midpoint_l_ij = midpoint_l_ik;
            midpoint_l_ik = tmp_len;

            int tmp_e = e_ij;
            e_ij = e_ik;
            e_ik = tmp_e;
        }

        // The angle beta is only is reasonable, when considering a right triangle, i.e. in the TRIANGLE/CHOPPED TRIANGLE CASES
        double alpha = angle_i_from_lengths(l_ij, l_ik, l_jk);
        double beta = M_PI/2.0 - alpha;

        // calc triangle side lengths
        double b = midpoint_l_ij;
        double a = tan(alpha)*b;
        double c = sqrt(a*a + b*b);

        // NOTE: If the angle is 0 or 180 degrees, the circumcentric area of p0 is 0.
        // NOTE: If the angle is greater than 0 degrees and less than 90 degrees, the circumcentric area of p0 can be calculated based on the area of the triangle between the midpoint1 and midpoint2 edges and the line containing midpoint1 with the direction of the normal of the midpoint1 edge.
        // NOTE: If the angle is 90 degrees, the circumcentric area of p0 can be calculated as the square of the length of the edge between p0 and midpoin1.
        // NOTE: If the angle is greater than 90 degrees and less than 180 degrees, the circumcentric area of p0 can be calculated based on the quadrilateral between the midpoint1 and midpoint2 edges and the two lines each containing either midpoint1 or midpoint2 with either the direction of the normal of the midpoint1 or midpoint2 edge.
        if (alpha < M_PI/2 - epsilon && c <= midpoint_l_ik)                     { calculation_case = CircumcentricShape::TRIANGLE; }
        if (alpha < M_PI/2 - epsilon && c > midpoint_l_ik)                      { calculation_case = CircumcentricShape::CHOPPED_TRIANGLE; }
        if (alpha >= M_PI/2 - epsilon && alpha <= M_PI/2 + epsilon)             { calculation_case = CircumcentricShape::RECTANGLE; }
        if (alpha > M_PI/2 + epsilon && alpha < M_PI - epsilon)                 { calculation_case = CircumcentricShape::CHOPPED_QUADRILATERAL; }

        switch(calculation_case)
        {
            default:
            case CircumcentricShape::ZERO:
                area_contribution = 0;
                break;
            case CircumcentricShape::TRIANGLE:
                area_contribution = a*b/2.0;
                break;
            case CircumcentricShape::CHOPPED_TRIANGLE:
            {
                // calc chopped off triangle side lengths
                double sa = c - midpoint_l_ik;
                double sb = tan(beta)*sa;

                area_contribution = (a*b - sa*sb)/2.0;
                break;
            }
            case CircumcentricShape::RECTANGLE:
                area_contribution = midpoint_l_ij*midpoint_l_ik;
                break;
            case CircumcentricShape::CHOPPED_QUADRILATERAL:
            {
                double triangle_area = 2.0*midpoint_l_ij*midpoint_l_ik*sin(alpha);

                double beta1 = angle_i_from_lengths(l_ij, l_jk, l_ik);
                double beta2 = angle_i_from_lengths(l_ik, l_jk, l_ij);
                double alpha1 = M_PI/2.0 - beta1;
                double alpha2 = M_PI/2.0 - beta2;
                double b1 = midpoint_l_ij / tan(alpha1);
                double b2 = midpoint_l_ik / tan(alpha2);
                double sta1 = abs(midpoint_l_ij*b1/2.0);
                double sta2 = abs(midpoint_l_ik*b2/2.0);

                area_contribution = triangle_area - sta1 - sta2;
                break;
            }
        }

        total_area += area_contribution;
    }

    return total_area;
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
    // std::cout << "e_idxs: " << edge_indices[0] << ", " << edge_indices[1] << ", " << edge_indices[2] << ", " << edge_indices[3] << ", " << edge_indices[4] << std::endl;
    double l_ij = iSData.L(edge_indices[0]);

    if (faces >= 1)
    {
        double l_ik = iSData.L(edge_indices[1]);
        double l_jk = iSData.L(edge_indices[2]);
        // std::cout << "Calling angle_i_from_lengths for alpha with: l_ij=" << l_ij << ", l_ik=" << l_ik << ", l_jk=" << l_jk << std::endl;
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
        double l_il = iSData.L(edge_indices[3]);
        double l_jl = iSData.L(edge_indices[4]);
        // std::cout << "Calling angle_i_from_lengths for beta with: l_ij=" << l_ij << ", l_il=" << l_il << ", l_jl=" << l_jl << std::endl;
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

double total_mass(const heatDiffData& hdData)
{
    double total_area = 0.0;
    for (int i=0; i<hdData.M.rows(); i++)
    {
        total_area += hdData.M(i, i);
    }
    return total_area;
}