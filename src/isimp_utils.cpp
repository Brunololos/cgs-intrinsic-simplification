#include "isimp_utils.hpp"

Eigen::MatrixXd rotate_selected(Eigen::MatrixXd U, Eigen::VectorXi selected_ids, RAXIS rot_axis, double angle)
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