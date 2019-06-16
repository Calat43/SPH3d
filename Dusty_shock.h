#pragma once

#include "SodTube3d.h"
#include "InitStates.h"
#include "Viscosity.h"
#include "Common.h"

namespace Dusty_shock_3d
{
    Grid init();

    Grid do_time_step(Grid & old_grid, int step_num, bool isIDIC);
}

namespace monaghan
{
    Point find_gas_velocity(Particle * particle, Cell * cell);
    Point find_dust_velocity(Particle * particle, Cell * cell);
}

namespace idic
{
    Point vel_asterisk(Cell * cell, Particle::Kind kind);

    double eps_asterisk(Cell * cell);

    double density_asterisk(Cell * cell, Particle::Kind kind);

    double t_stop_asterisk(Cell * cell);

    Point pressure_term_asterisk(Cell * cell);

    Point x_through_vel(Cell * cell);
    Point y_through_vel(Cell * cell);

    Point find_x(Cell * cell);
    Point find_y(Cell * cell);

    Point find_gas_vel_asterisk(Cell * cell);
    Point find_dust_vel_asterisk(Cell * cell);

    Point find_gas_velocity(Particle * particle, Cell * cell);
    Point find_dust_velocity(Particle * particle, Cell * cell);
}
