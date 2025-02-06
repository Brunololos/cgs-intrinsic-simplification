#ifndef QUERY_TEXTURE_BARYCENTRICS
#define QUERY_TEXTURE_BARYCENTRICS

#include <Eigen/Core>
#include <Eigen/Dense>

#include <assert.h>
#include <iostream>
#include <vector>
#include <queue>
#include <cmath>
#include <tuple>

#include <igl/AABB.h>
#include "isimp_utils.hpp"

// Code roughly adopted from https://github.com/HTDerekLiu/intrinsic-simplification/src/query_texture_barycentric.cpp
void query_texture_barycentrics(
                                const Eigen::MatrixXd& UV,
                                const Eigen::MatrixXi& UF,
                                const int& TEXTURE_WIDTH,
                                const int& TEXTURE_HEIGHT,
                                iSimpData& iSData);

#endif