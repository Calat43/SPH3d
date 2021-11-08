#pragma once

#include "Particle.h"

class Cell;
class Grid;

namespace Dusty_shock_1d
{
    Grid init();

    Grid do_time_step(Grid & old_grid, int step_num);
}