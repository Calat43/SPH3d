#pragma once

#include <vector>

#include "MathUtils.h"

class Grid;

Grid squared_ball_init_state(double radius, double step);

namespace Sod_tube_1d
{
    Grid init();
}

void fill_initial_sod_coord(std::vector<double> & coord, std::vector<double> & image_coord, uint real_left_p_num,
    uint real_right_p_num, uint image_left_p_num, uint image_right_p_num);