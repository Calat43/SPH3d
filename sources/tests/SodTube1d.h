#pragma once

class Cell;
class Grid;
class Point;

double find_density(Particle * particle, Cell * cell);

double find_density_no_sort(Particle const & particle, Grid const & grid);

void recalc_density(Grid & grid, Particle::Kind kind);

Point find_new_velocity(Particle * particle, Cell * cell);

Point find_new_velocity_no_sort(Particle * particle, Grid * grid);

// 1D, x axis
namespace Sod_tube_1d
{
    Point find_new_velocity(Particle * particle, Cell * cell);

    double find_new_pressure(Particle * particle);

    double find_new_energy(Particle * particle, Cell * cell);

    Grid do_time_step(Grid & old_grid, int step_num);

    void recalc_pressure(Grid & grid, Particle::Kind kind);
}
