#pragma once

class Grid;
class Point;
class Particle;


class SPHSolver
{
public:
    SPHSolver();
    virtual ~SPHSolver();

    int generate_initial_distribution();
    int do_time_step();

private:
    int generate_uniform_distribution();
    int generate_ball_distribution();

    void swap_grids_and_clear_next();
    
    Grid* current_grid = nullptr;
    Grid* next_grid = nullptr;
};

Point find_new_coordinates_(Particle const & particle);
