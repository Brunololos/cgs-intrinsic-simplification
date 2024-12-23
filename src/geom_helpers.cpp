#include "geom_helpers.hpp"

Eigen::MatrixXd rotate_selected(Eigen::MatrixXd U, const Eigen::VectorXi selected_ids, const RAXIS rot_axis, double angle)
{
  angle = (0.5 / 360.0) * 2.0 * 3.14159265;
  Eigen::Vector3d mean = Eigen::Vector3d::Zero();
  for (int i=0; i<selected_ids.size(); i++)
  {
    double id = selected_ids(i);
    mean += U.row(id);
  }
  mean = mean / selected_ids.size();

  double u0, u1, u2, p0, p1, p2;
  for (int i=0; i<selected_ids.size(); i++)
  {
    double id = selected_ids(i);
    u0 = U(id, 0) - mean(0);
    u1 = U(id, 1) - mean(1);
    u2 = U(id, 2) - mean(2);
    switch (rot_axis)
    {
      case RAXIS::ROLL:
        p0 = u0;
        p1 = cos(angle)*u1 - sin(angle)*u2;
        p2 = sin(angle)*u1 + cos(angle)*u2;
        break;
      default:
      case RAXIS::YAW:
        p0 = cos(angle)*u0 + sin(angle)*u2;
        p1 = u1;
        p2 = - sin(angle)*u0 + cos(angle)*u2;
        break;
      case RAXIS::PITCH:
        p0 = cos(angle)*u0 - sin(angle)*u1;
        p1 = sin(angle)*u0 + cos(angle)*u1;
        p2 = u2;
        break;
    }
    U(id, 0) = p0 + mean(0);
    U(id, 1) = p1 + mean(1);
    U(id, 2) = p2 + mean(2);
  }
  return U;
}

Quad2D unfold(const double l_ij, const double l_ik, const double l_jk, const double l_il, const double l_jl)
{
  double theta_i_k = angle_i_from_lengths(l_ij, l_ik, l_jk);
  double theta_i_l = angle_i_from_lengths(l_ij, l_il, l_jl);
  std::cout << "calced theta_i_k: " << theta_i_k << std::endl;
  std::cout << "calced theta_i_l: " << theta_i_l << std::endl;
  Quad2D unfolded = Quad2D();
  unfolded.row(0) = Point2D(0.0, 0.0);
  unfolded.row(1) = Point2D(l_ij, 0.0);
  unfolded.row(2) = Point2D(l_ik * cos(theta_i_k), l_ik * sin(theta_i_k));
  unfolded.row(3) = Point2D(l_il * cos(theta_i_l), -l_il * sin(theta_i_l));
  return unfolded;
}

double flipped_edgelength(const Quad2D& trianglepair)
{
  return (trianglepair.row(2) - trianglepair.row(3)).norm();
}

Point2D to_explicit(const BarycentricPoint& bary_point, const Point2D& A, const Point2D& B, const Point2D& C)
{
  std::cout << "making explicit, barycentric point " << to_str(bary_point) << " as convex combination of A: " << to_str(A) << ", B: " << to_str(B) << " and C: " << to_str(C) << std::endl;
  return A*bary_point(0) + B*bary_point(1) + C*bary_point(2);
}

// TODO: one could also store the values that do not depend on point in the simplification operation
// Barycentric point calculation adopted from: https://gamedev.stackexchange.com/questions/23743/whats-the-most-efficient-way-to-find-barycentric-coordinates
BarycentricPoint to_barycentric(const Point2D& point, const Point2D& A, const Point2D& B, const Point2D& C)
{
  std::cout << "making barycentric, explicit point " << to_str(point) << " as convex combination of A: " << to_str(A) << ", B: " << to_str(B) << " and C: " << to_str(C) << std::endl;
  Vector2D AB = B - A;
  Vector2D AC = C - A;
  Vector2D AP = point - A;

  double D00 = AB.dot(AB);
  double D01 = AB.dot(AC);
  double D11 = AC.dot(AC);
  double D20 = AB.dot(AP);
  double D21 = AC.dot(AP);
  double denominator = D00 * D11 - D01 * D01;
  double x2 = (D11 * D20 - D01 * D21) / denominator;
  double x3 = (D00 * D21 - D01 * D20) / denominator;
  double x1 = (1.0f - x2 - x3) / denominator;
  return BarycentricPoint(x1, x2, x3);
}

double angle_i_from_lengths(const double l_ij, const double l_ik, const double l_jk)
{
  // TODO: It might be that the denominator becomes zero because we work on Delta complexes where distances can become zero.
  double enumerator = l_ij*l_ij + l_ik*l_ik - l_jk*l_jk;
  double denominator = 2 * l_ij * l_ik;
  std::cout << "Calcing theta_i_... enumerator=" << enumerator << ", denominator=" << denominator << std::endl; // TODO: remove
  assertm(denominator != 0, "angle_i_from_lengths: Denominator became zero!");
  assertm(enumerator / denominator < 0, "angle_i_from_lengths: cos(theta) became less than zero! We need to implement a case for obtuse triangles.");
  return acos(enumerator / denominator);
}

double to_degrees(const double angle_in_radians) { return (angle_in_radians * 180.0) / M_PI; }
double to_radians(const double angle_in_degrees) { return (angle_in_degrees * M_PI) / 180.0; }