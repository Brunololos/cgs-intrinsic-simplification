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
  // TODO: remove prints
  // std::cout << "calced theta_i_k: " << theta_i_k << std::endl;
  // std::cout << "calced theta_i_l: " << theta_i_l << std::endl;
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

// assuming quad ijkl of triangles ijk, ijl
bool is_convex(const Quad2D& trianglepair)
{
  // calculate cycle around quad boundary
  Vector2D ik = trianglepair.row(2) - trianglepair.row(0);
  Vector2D kj = trianglepair.row(1) - trianglepair.row(2);
  Vector2D jl = trianglepair.row(3) - trianglepair.row(1);
  Vector2D li = trianglepair.row(0) - trianglepair.row(3);

  Vector2D ij = trianglepair.row(1) - trianglepair.row(0);
  Vector2D kl = trianglepair.row(3) - trianglepair.row(2);
  Vector2D ji = -ij;
  Vector2D lk = -kl;

  double prev = 0.0;
  double curr = 0.0;
  for (int i=0; i<4; i++)
  {
    // if (i == 0) { curr = ik.cross(ij).sum(); }
    // if (i == 1) { curr = kj.cross(kl).sum(); }
    // if (i == 2) { curr = jl.cross(ji).sum(); }
    // if (i == 3) { curr = li.cross(lk).sum(); }

    if (i == 0) { curr = scalar_cross(ik, ij); }
    if (i == 1) { curr = scalar_cross(kj, kl); }
    if (i == 2) { curr = scalar_cross(jl, ji); }
    if (i == 3) { curr = scalar_cross(li, lk); }

    if (curr != 0.0)
    {
      if (curr * prev < 0.0)
      {
        return false;
      }
      else
      {
        prev = curr;
      }
    }
  }
}

// Point in triangle check logic adopted from: https://stackoverflow.com/questions/2049582/how-to-determine-if-a-point-is-in-a-2d-triangle
double sign_helper(const double px, const double py, const Point2D p1, const Point2D p2)
{
    return (px - p2[0]) * (p1[1] - p2[1]) - (p1[0] - p2[0]) * (py - p2[1]);
}

bool lies_inside_triangle(const double px, const double py, const Point2D v0, const Point2D v1, const Point2D v2)
{
  double d0, d1, d2;
  bool found_neg, found_pos;

  d0 = sign_helper(px, py, v0, v1);
  d1 = sign_helper(px, py, v1, v2);
  d2 = sign_helper(px, py, v2, v0);

  found_neg = (d0 < 0) || (d1 < 0) || (d2 < 0);
  found_pos = (d0 > 0) || (d1 > 0) || (d2 > 0);

  return !(found_neg && found_pos);
}

Point2D to_explicit(const BarycentricPoint& bary_point, const Point2D& A, const Point2D& B, const Point2D& C)
{
  // TODO: remove print
  // std::cout << "making explicit, barycentric point " << to_str(bary_point) << " as convex combination of A: " << to_str(A) << ", B: " << to_str(B) << " and C: " << to_str(C) << std::endl;
  return A*bary_point(0) + B*bary_point(1) + C*bary_point(2);
}

// TODO: one could also store the values that do not depend on point in the simplification operation
// Barycentric point calculation adopted from: https://gamedev.stackexchange.com/questions/23743/whats-the-most-efficient-way-to-find-barycentric-coordinates
BarycentricPoint to_barycentric(const Point2D& point, const Point2D& A, const Point2D& B, const Point2D& C)
{
  // TODO: remove print
  // std::cout << "making barycentric, explicit point " << to_str(point) << " as convex combination of A: " << to_str(A) << ", B: " << to_str(B) << " and C: " << to_str(C) << std::endl;
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
  double x1 = (1.0f - x2 - x3);

  return BarycentricPoint(x1, x2, x3);
}

double angle_i_from_lengths(const double l_ij, const double l_ik, const double l_jk)
{
  double epsilon = 0.000001;
  // TODO: It might be that the denominator becomes zero because we work on Delta complexes where distances can become zero.
  double enumerator = l_ij*l_ij + l_ik*l_ik - l_jk*l_jk;
  double denominator = 2 * l_ij * l_ik;
  // TODO: remove print
  // std::cout << "Calcing theta_i_... enumerator=" << enumerator << ", denominator=" << denominator << std::endl; // TODO: remove
  // assertm(denominator != 0, "angle_i_from_lengths: Denominator became zero!");
  // assertm(enumerator / denominator < 0, "angle_i_from_lengths: cos(theta) became less than zero! We need to implement a case for obtuse triangles.");
  if (denominator == 0) { std::cout << dye("angle_i_from_lengths: Denominator became zero!", RED) << std::endl; }
  // if (enumerator / denominator < 0) { std::cout << dye("angle_i_from_lengths: cos(theta) became less than zero! We need to implement a case for obtuse triangles.", RED) << std::endl; }
  // TODO: Here I tried to introduce an epsilon for precisions slackness to remedy numerical issues (wasn't successful)
  if (enumerator / denominator >= 1.0 && enumerator / denominator <= 1.0 + epsilon) { return 0.0; }
  if (enumerator / denominator <= -1.0 && enumerator / denominator >= -1.0 - epsilon) { return M_PI; }
  if (enumerator / denominator < -1 || enumerator / denominator > 1) { std::cout << dye("angle_i_from_lengths: cos(theta) = " + std::to_string(enumerator / denominator) + " of the lengths l_ij: " + std::to_string(l_ij) + ", l_ik: " + std::to_string(l_ik) + ", l_jk: " + std::to_string(l_jk) + " is outside the interval [-1, 1]! The arccos would yield NaN. => The passed edge length's were invalid and likely fail to satisfy the triangle inequality!", RED) << std::endl; }
  return acos(enumerator / denominator);
}

double scalar_cross(const Vector2D& v, const Vector2D& w)
{
  return v[0] * w[1] - v[1] * w[0];
}

PolarVector2D to_polar(const Vector2D& v)
{
  double magnitude = v.norm();
  double angle = std::atan2(v[1], v[0]);
  return PolarVector2D({angle, magnitude});
}

Vector2D to_cartesian(const PolarVector2D& v)
{
  double x = v[1] * std::cos(v[0]);
  double y = v[1] * std::sin(v[0]);
  return Vector2D({x, y});
}

double to_degrees(const double angle_in_radians) { return (angle_in_radians * 180.0) / M_PI; }
double to_radians(const double angle_in_degrees) { return (angle_in_degrees * M_PI) / 180.0; }