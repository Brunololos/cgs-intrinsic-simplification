#include "heat_diff.hpp"

#include <random>

void initHeatDiff(heatDiffData& hdData, const iSimpData& iSData)
{
    initHeat(hdData, iSData);
    build_heat_sample_indices(hdData, iSData);
    build_cotan_mass(hdData, iSData);
    build_cotan_laplacian(hdData, iSData);
}

// generate random heat values for each vertex
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

void build_heat_sample_indices(heatDiffData& hdData, const iSimpData& iSData)
{
    hdData.heat_sample_vertices = std::vector<int>();
    hdData.vertex_idcs_to_heat_idcs = std::unordered_map<int, int>();
    int i=0;
    for (gcs::Vertex V : iSData.recoveredMesh->vertices())
    {
        // if (V.isDead()) { continue; }
        hdData.heat_sample_vertices.push_back(V.getIndex());
        hdData.vertex_idcs_to_heat_idcs.insert({V.getIndex(), i});
        i++;
    }
}

void diffuse_heat(heatDiffData& hdData, const double stepsize, const int steps)
{
    // TODO: guard check that Laplacian matrices already have been calculated in heatDiffData
    Eigen::MatrixXd A, B;
    A = hdData.M - stepsize * hdData.L;

    // choose only relevant heat values from initial heat
    int n = hdData.heat_sample_vertices.size();
    hdData.final_heat = Eigen::Matrix<double, -1, 1>::Zero(n, 0);
    for (int i=0; i<n; i++)
    {
        int vertex_idx = hdData.heat_sample_vertices[i];
        hdData.final_heat(i, 0) = hdData.initial_heat(vertex_idx);
        // std::cout << "iteration " << i << " progress_percent: " << progress_percent << ", next_bar_percent: " << next_bar_percent << std::endl;
    }
    // printEigenMatrixXd("initial heat", hdData.final_heat);

    // implicit Euler steps
    // std::cout << "vertex 383 laplacian idx: " << hdData.vertex_idcs_to_heat_idcs.at(383) << std::endl;
    // auto Lrow = hdData.L.row(hdData.vertex_idcs_to_heat_idcs.at(383));
    // for (int i=0; i<Lrow.cols(); i++)
    // {
    //     if (Lrow(0, i) != 0.0) {
    //         std::cout << "Vertex 383 adjacency entry: " << i << " for vertex " << hdData.heat_sample_vertices[i] << " has weight: " << Lrow(0, i) << std::endl;
    //     }
    // }
    // std::cout << "Vertex 383 mass: " << hdData.M(hdData.vertex_idcs_to_heat_idcs.at(383), hdData.vertex_idcs_to_heat_idcs.at(383)) << std::endl;
    // printEigenMatrixXd("vertex 383 laplacian row", hdData.L.row(hdData.vertex_idcs_to_heat_idcs.at(383)));
    // std::cout << "updated vertex 383 heat: " << get_heat(hdData, 383) << std::endl;
    std::cout << "diffusing heat: [" << std::flush;
    int total_bars = 24;
    int current_bars = 0;
    if(steps == 0) { for (int i=0; i<total_bars; i++) { std::cout << "#"; } }

    for (int i=0; i<steps; i++)
    {
        B = hdData.M * hdData.final_heat;
        hdData.final_heat = A.ldlt().solve(B);

        // progress bar
        float progress_percent = ((float) i+1)/((float) steps);
        float next_bar_percent = ((float) current_bars+1)/((float) total_bars);
        // std::cout << "iteration " << i << " progress_percent: " << progress_percent << ", next_bar_percent: " << next_bar_percent << std::endl;
        for (; progress_percent >= next_bar_percent;)
        {
            current_bars++;
            next_bar_percent = ((float) current_bars+1)/((float) total_bars);
            std::cout << "#" << std::flush;
        }
    }
    // printEigenMatrixXd("final heat", hdData.final_heat);
    // std::cout << "Successfully diffused heat!" << std::endl;
    std::cout << "]" << std::endl;
}

double get_heat(heatDiffData& hdData, const int vertex_idx)
{
    return hdData.final_heat(hdData.vertex_idcs_to_heat_idcs.at(vertex_idx), 0);
}

void build_cotan_mass(heatDiffData& hdData, const iSimpData& iSData)
{
    int n = hdData.heat_sample_vertices.size();
    hdData.M = Eigen::MatrixXd::Zero(n, n);
    for (int i=0; i < n; i++)
    {
        int vertex_idx = hdData.heat_sample_vertices[i];
        hdData.M(i, i) = circumcentric_dual_mass(iSData, vertex_idx);
        // std::cout << "vertex " << i << " area mass: " << hdData.M(i, i) << std::endl;
    }
    // std::cout << "computed total area mass: " << total_mass(hdData) << std::endl;
    // std::cout << "Built cotan mass" << std::endl;
}

void build_cotan_laplacian(heatDiffData& hdData, const iSimpData& iSData)
{
    int n = hdData.heat_sample_vertices.size();
    hdData.L = Eigen::MatrixXd::Zero(n, n);

    for (int i=0; i<n; i++)
    {
        double vertex_weight = 0.0;
        int vertex_idx = hdData.heat_sample_vertices[i];
        gcs::Vertex vi = iSData.recoveredMesh->vertex(vertex_idx);
        // std::cout << "Got vertex" << std::endl;
        for (gcs::Vertex vj : vi.adjacentVertices())
        {
            int neighbor_idx = vj.getIndex();
            int j = hdData.vertex_idcs_to_heat_idcs.at(neighbor_idx);
            if (vertex_idx == neighbor_idx) { std::cout << RED << "build_cotan_laplacian: Encountered self-edge!" << RESET << std::endl; exit(-1); }

            // std::cout << "Getting connecting edge idx... " << std::endl;
            int edge_idx = iSData.recoveredMesh->connectingEdge(vi, vj).getIndex();
            // std::cout << "Got connecting edge idx: " << edge_idx << std::endl;
            double edge_weight = cotan_weight(iSData, edge_idx);
            // std::cout << "Calculated cotan weight" << std::endl;

            // std::cout << "determined connecting edge: " << edge_idx << ": (" << iSData.intrinsicMesh->edge(edge_idx).firstVertex() << ", " << iSData.intrinsicMesh->edge(edge_idx).secondVertex() << ") between the vertices: " << i << ", " << j << std::endl;
            hdData.L(i, j) = edge_weight;
            vertex_weight += edge_weight;
        }
        hdData.L(i, i) = -vertex_weight;
        // std::cout << "Calculated vertex: " << i << std::endl;
    }
    // std::cout << "Built cotan laplacian" << std::endl;
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
    gcs::Vertex V = iSData.recoveredMesh->vertex(vertex_idx);
    // double epsilon = std::numeric_limits<double>().epsilon();
    double epsilon = 0.000001;

    for (gcs::Face F : V.adjacentFaces())
    {
        CircumcentricShape calculation_case = CircumcentricShape::ZERO;
        double area_contribution;
        std::array<int, 3> face_edges = order_triangle_edge_indices(F, vertex_idx);

        int e_ij = face_edges[0];
        int e_ik = face_edges[1];
        int e_jk = face_edges[2];
        double l_ij = iSData.recovered_L(e_ij);
        double l_ik = iSData.recovered_L(e_ik);
        double l_jk = iSData.recovered_L(e_jk);
        double midpoint_l_ij = l_ij/2.0;
        double midpoint_l_ik = l_ik/2.0;
        // if (vertex_idx == 383) { std::cout << "Got edge-lengths: l_ij: " << l_ij << ", l_ik: " << l_ik << ", l_jk: " << l_jk << std::endl; }

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
            // if (vertex_idx == 383) { std::cout << "swapping sides" << std::endl; }
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
        // NOTE: If the angle is 90 degrees, the circumcentric area of p0 can be calculated as the square of the length of the edge between p0 and midpoint1.
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

                // if (vertex_idx == 383)
                // {
                //     std::cout << "alpha = " << alpha << std::endl;
                //     std::cout << "beta = " << beta << std::endl;
                //     std::cout << "a = " << a << std::endl;
                //     std::cout << "b = " << b << std::endl;
                //     std::cout << "sa = " << sa << std::endl;
                //     std::cout << "sb = " << sb << std::endl;
                // }

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
        // if (vertex_idx == 383)
        // {
        //     std::cout << "Vertex 383 mass contribution of: " << area_contribution << " for calculation case: ";
        //     switch(calculation_case)
        //     {
        //         default:
        //         case CircumcentricShape::ZERO:
        //             std::cout << "ZERO!" << std::endl;
        //             break;
        //         case CircumcentricShape::TRIANGLE:
        //             std::cout << "TRIANGLE" << std::endl;
        //             break;
        //         case CircumcentricShape::CHOPPED_TRIANGLE:
        //         {
        //             std::cout << "CHOPPED_TRIANGLE" << std::endl;
        //             break;
        //         }
        //         case CircumcentricShape::RECTANGLE:
        //             std::cout << "RECTANGLE" << std::endl;
        //             break;
        //         case CircumcentricShape::CHOPPED_QUADRILATERAL:
        //         {
        //             std::cout << "CHOPPED_QUADRILATERAL" << std::endl;
        //             break;
        //         }
        //     }
        // }

        total_area += area_contribution;
    }

    return total_area;
}

double cotan_weight(const iSimpData& iSData, const int edge_idx)
{
    gcs::Edge edge = iSData.recoveredMesh->edge(edge_idx);
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
    double l_ij = iSData.recovered_L(edge_indices[0]);

    if (faces >= 1)
    {
        double l_ik = iSData.recovered_L(edge_indices[1]);
        double l_jk = iSData.recovered_L(edge_indices[2]);
        // std::cout << "Calling angle_i_from_lengths for alpha with: l_ij=" << l_ij << ", l_ik=" << l_ik << ", l_jk=" << l_jk << std::endl;
        // double alpha = angle_i_from_lengths(l_ij, l_ik, l_jk);
        double alpha = angle_i_from_lengths(l_ik, l_jk, l_ij);
        // if (edge.firstVertex().getIndex() == 383 || edge.secondVertex().getIndex() == 383)
        // {
        //     std::cout << "Angle alpha: " << (alpha*180)/M_PI << " for edge: (" << edge.firstVertex().getIndex() << ", " << edge.secondVertex().getIndex() << ")" << std::endl;
        // }
        if (alpha > 2*M_PI) { std::cout << RED << "cotan_weight encountered angle larger than 2PI or equivalently 180 degrees!" << RESET << std::endl; exit(-1); }

        if (alpha <= M_PI/2.0)
        {
            weight += (alpha > deg_epsilon) ? (cos(alpha)/sin(alpha)) : l_ik;
        }
        if (alpha > M_PI/2.0)
        {
            // weight += (alpha < M_PI - deg_epsilon) ? -(cos(alpha)/sin(alpha)) : -l_ik;
            weight += (alpha < M_PI - deg_epsilon) ? (cos(alpha)/sin(alpha)) : -l_ik;
        }

        // if (edge.firstVertex().getIndex() == 383 || edge.secondVertex().getIndex() == 383)
        // {
        //     std::cout << "alpha weight: " << weight << std::endl;
        // }
    }
    // TODO: remove print temp thing
    double prev_weight = weight;

    if (faces >= 2)
    {
        double l_il = iSData.recovered_L(edge_indices[3]);
        double l_jl = iSData.recovered_L(edge_indices[4]);
        // std::cout << "Calling angle_i_from_lengths for beta with: l_ij=" << l_ij << ", l_il=" << l_il << ", l_jl=" << l_jl << std::endl;
        // double beta = angle_i_from_lengths(l_ij, l_il, l_jl);
        double beta = angle_i_from_lengths(l_il, l_jl, l_ij);
        // if (edge.firstVertex().getIndex() == 383 || edge.secondVertex().getIndex() == 383)
        // {
        //     std::cout << "Angle beta: " << (beta*180)/M_PI << " for edge: (" << edge.firstVertex().getIndex() << ", " << edge.secondVertex().getIndex() << ")" << std::endl;
        // }
        if (beta > 2*M_PI) { std::cout << RED << "cotan_weight encountered angle larger than 2PI or equivalently 180 degrees!" << RESET << std::endl; exit(-1); }

        if (beta <= M_PI/2.0)
        {
            weight += (beta > deg_epsilon) ? (cos(beta)/sin(beta)) : l_il;
        }
        if (beta > M_PI/2.0)
        {
            // weight += (beta < M_PI - deg_epsilon) ? -(cos(beta)/sin(beta)) : -l_il;
            weight += (beta < M_PI - deg_epsilon) ? (cos(beta)/sin(beta)) : -l_il;
        }
        // if (edge.firstVertex().getIndex() == 383 || edge.secondVertex().getIndex() == 383)
        // {
        //     std::cout << "beta weight: " << weight - prev_weight << std::endl;
        // }
    }

    // if (edge.firstVertex().getIndex() == 383 || edge.secondVertex().getIndex() == 383)
    // {
    //     std::cout << "With the total weight: " << weight/2.0 << std::endl;
    // }

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