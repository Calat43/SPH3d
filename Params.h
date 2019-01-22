#pragma once

#include <cmath>
#include "Point.h"

// Singleton
class Params
{
public:
    static Params & get_instance()
    {
        static Params p; // lazy initialized
        return p;
    }

    double t = 0.01;
    double c_s = 1;
    double K = NAN;
    double kx = NAN;
    double ky = NAN;
    double kz = NAN;
    int n_gas = 10000;
    int n_dust = 0;
    int dimensions = 3;
    double h = 0.02;
    double smooth_radius = 2. * h;
    double tau = 0.001;

    double middle_gas_dens = NAN;
    double d2g = NAN;
    double delta = NAN;

    // boundary of area
    Point border1 = Point(0, 0, 0);
    Point border2 = Point(1, 1, 1);

    double grid_step_x = 0.01;
    double grid_step_y = 0.01;
    double grid_step_z = 0.01;

    static int DBG_ID;

    Params(Params const &) = delete;
    Params & operator=(Params const &) = delete;

private:
    Params() = default;

    ~Params() = default;
};
