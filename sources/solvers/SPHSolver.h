#pragma once

class Grid;


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
    
    Grid* current_grid = 0;
    Grid* next_grid = 0;
};