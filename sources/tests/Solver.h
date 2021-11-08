#pragma once

class Cell;
class Grid;

Point find_new_coordinates(Particle const & particle);

//double find_density(Particle * particle, Cell * cell);

//void recalc_density(Grid & grid, Particle::Kind kind);

double theta(Particle * particle, Grid & grid);

namespace idic_1d
{
    //double pressure_term_asterisk(Cell * cell);

    double vel_asterisk(Cell * cell, Particle::Kind kind);

    double eps_asterisk(Cell * cell);

    double t_stop_asterisk(Cell *cell);

    double x_through_vel(Cell * cell);

    double y_through_vel(Cell * cell);

    double find_x(Cell * cell);

    double find_y(Cell * cell);

    double find_gas_vel_asterisk(Cell * cell);

    double find_dust_vel_asterisk(Cell * cell);

    //double pressure_term(Particle * particle, Cell * cell);

    double find_gas_velocity(Particle * particle, Cell * cell);

    double find_dust_velocity(Particle * particle, Cell * cell);
}

double gas_velocity_MK(Particle * particle, Cell * cell);

double dust_velocity_MK(Particle * particle, Cell * cell);

double gas_velocity_shock_MK(Particle * particle, Cell * cell);

double dust_velocity_shock_MK(Particle * particle, Cell * cell);