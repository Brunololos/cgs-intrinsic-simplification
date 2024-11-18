#ifndef ISIMP_UTILS_H
#define ISIMP_UTILS_H

#define _USE_MATH_DEFINES

#include <iomanip>
#include <math.h>
#include "tinyad_defs.hpp"
#include "print.hpp"

enum class RAXIS {
    YAW = 0,
    PITCH = 1,
    ROLL = 2,
    NONE = 3
};

Eigen::MatrixXd rotate_selected(Eigen::MatrixXd U, Eigen::VectorXi selected_ids, RAXIS rot_axis, double angle);

#endif