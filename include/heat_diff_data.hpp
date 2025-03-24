#ifndef HEAT_DIFF_DATA_H
#define HEAT_DIFF_DATA_H

#include <Eigen/Core>
#include <Eigen/Dense>
#include <unordered_map>

struct heatDiffData {
    Eigen::Matrix<double, -1, 1> initial_heat;
    Eigen::Matrix<double, -1, 1> recovered_heat;
    Eigen::Matrix<double, -1, 1> final_heat;

    // ordered list of vertex indices, selected as a subset from the vertices of the original mesh
    // NOTE: this is to ensure a consistent ordering & allows for calculating heat diff on a simplified mesh,
    // while being able to map the results back to the original mesh via the original vertex indices.
    std::vector<int> heat_sample_vertices;
    std::unordered_map<int, int> vertex_idcs_to_heat_idcs;

    // Laplacian // TODO: store laplacian for multiple different simplification states
    Eigen::MatrixXd M;
    Eigen::MatrixXd L;
};

#endif