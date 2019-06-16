#pragma once

#include <cassert>
#include <iostream>
#include <fstream>

#include "Params.h"
#include "Common.h"
#include "Cell.h"
#include "Grid.h"
#include "SodTube3d.h"
#include "Solver.h"

namespace ball_analytic
{
    //t(R)
    double time_through_R(double R);

    double R_bisection(double exact_time, double step, double defect);

    double density(double exact_time, double step, double defect);

    double velocity(double r, double f);

    double pressure(double r, double R, double ddot_R);

    void print_solution(double exact_time, double step_r, double step_R, double defect);
}

namespace ball_in_vacuum
{
    double find_pressure(Particle * particle, Cell * cell);

    Point init_velocity(Point center, Point particle_coordinates);

    Grid init(double step);

    Grid do_time_step(Grid & old_grid, int step_num);
}