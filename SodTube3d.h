#pragma once

#include "Solver.h"
#include "InitStates.h"
#include "CompareFiles.h"

namespace Sod_tube_3d
{
    Grid init();

    Point find_new_velocity(Particle * particle, Cell * cell);

    double find_new_energy(Particle * particle, Cell * cell);

    Grid do_time_step(Grid & old_grid, int step_num);
}