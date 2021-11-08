#pragma once

class Grid;
class Particle;
class Cell;

namespace Sod_tube_3d
{
    Point find_new_velocity(Particle * particle, Cell * cell, int step_num);

    double find_new_energy(Particle * particle, Cell * cell);

    Grid init_with_boundaries();
}

