#ifndef GEOM_HELPERS_H
#define GEOM_HELPERS_H

#include <Eigen/Core>
#include <Eigen/Dense>
#include <vector>
#include <math.h>
#include <iomanip>
#include <cassert>

// TODO: remove later
#include "print.hpp"

#define _USE_MATH_DEFINES
#define assertm(exp, msg) assert(((void)msg, exp))

// semantic typedefs
typedef double Bias;
typedef Eigen::Vector2i TexCoord;
typedef Eigen::Vector2d Point2D, Normal2D, Vector2D, PolarVector2D;
typedef Eigen::Vector3d Point3D, BarycentricPoint;
typedef Eigen::Matrix<double, 4, 2> Quad2D; // first two positions are the common vertices and the last two are the unique ones

enum class RAXIS {
    YAW = 0,
    PITCH = 1,
    ROLL = 2,
    NONE = 3
};

Eigen::MatrixXd rotate_selected(Eigen::MatrixXd U, const Eigen::VectorXi selected_ids, const RAXIS rot_axis, double angle);

// unfolds intrinsic triangles ijk, ijl into an extrinsic planar representation
Quad2D unfold(const double l_ij, const double l_ik, const double l_jk, const double l_il, const double l_jl);
double flipped_edgelength(const Quad2D& trianglepair);
bool is_convex(const Quad2D& trianglepair);

bool lies_inside_triangle(const double px, const double py, const Point2D v0, const Point2D v1, const Point2D v2);

Point2D to_explicit(const BarycentricPoint& bary_point, const Point2D& A, const Point2D& B, const Point2D& C);
BarycentricPoint to_barycentric(const Point2D& point, const Point2D& A, const Point2D& B, const Point2D& C);

// calculate interior angle at corner i of a triangle ijk with edge lengths l_ij, l_ik, l_jk (via law of cosines)
double angle_i_from_lengths(const double l_ij, const double l_ik, const double l_jk);

double scalar_cross(const Vector2D& v, const Vector2D& w);

PolarVector2D to_polar(const Vector2D& v);
Vector2D to_cartesian(const PolarVector2D& v);

double to_degrees(const double angle_in_radians);
double to_radians(const double angle_in_degrees);

#endif