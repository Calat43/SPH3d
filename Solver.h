#pragma once

#include <iostream>

#include "Params.h"
#include "Common.h"
#include "Cell.h"

Point find_new_coordinates(Particle const & particle);

double find_density_no_sort(Particle const & particle, Grid const & grid);

void recalc_density(Grid & grid);

Point find_new_velocity(Particle * particle, Cell * cell);

Point find_new_velocity_no_sort(Particle * particle, Grid * grid);

// 1D, x axis
namespace Sod_tube
{
    Point find_new_velocity(Particle * particle, Cell * cell);

    double find_new_pressure(Particle * particle);

    double find_new_energy(Particle * particle, Cell * cell);

    Grid do_time_step(Grid & old_grid, int step_num);

    void recalc_pressure(Grid & grid);
}