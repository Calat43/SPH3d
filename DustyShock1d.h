#pragma once

#include "Solver.h"
#include "InitStates.h"
#include "Cell.h"
#include "Grid.h"
#include "Particle.h"
#include "NonLinear.h"
#include "Viscosity.h"

namespace Dusty_shock_1d
{
    double pressure_term(Particle * particle, Cell * cell);

    Grid init();

    Grid do_time_step(Grid & old_grid, int step_num);
}

namespace idic_1d
{
    double pressure_term_asterisk(Cell * cell);

    double vel_asterisk(Cell * cell, Particle::Kind kind);

    double eps_asterisk(Cell * cell);

    double x_through_vel(Cell * cell);

    double y_through_vel(Cell * cell);

    double find_x(Cell * cell);

    double find_y(Cell * cell);

    double find_gas_vel_asterisk(Cell * cell);

    double find_dust_vel_asterisk(Cell * cell);

    double find_gas_velocity(Particle * particle, Cell * cell);

    double find_dust_velocity(Particle * particle, Cell * cell);
}