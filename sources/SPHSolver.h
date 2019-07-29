#pragma once

class Grid;


class SPHSolver
{
public:
    SPHSolver();
    virtual ~SPHSolver();

    int generate_initial_conditions();
    int do_time_step();

private:
    Grid* grid;

};