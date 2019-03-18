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

    double t = 0.1;
    double c_s = 1;
    double gamma = 7. / 5.;
    double K = NAN;
    double kx = NAN;
    double ky = NAN;
    double kz = NAN;
    int n_gas = 10000;
    int n_dust = 0;
    int dimensions = 3;
    double h = 0.01;
    double smooth_radius = 2. * h;
    double tau = 0.001;

    double middle_gas_dens = NAN;
    double d2g = NAN;
    double delta = NAN;

    //Sod tube
    uint real_particles = 200;//304;
    uint image_particles = 20;//32;
    double membrane = 0.5;
    double left = 0.2;
    double right = 0.8;

    double dens_right = 1.;//0.125;
    double press_right = 0.8;//0.1;
    double vel_right = 0.;
    double ener_right = 2.;

    double dens_left = 1.;//3. * dens_right; //0.375
    double press_left = 1.;//0.375;
    double vel_left = 0.;
    double ener_left = 2.5;


    //viscosity
    double alpha = 1.;
    double beta = 2.;
    double nu = 0.1 * h;
    bool have_viscosity = true;

    // boundary of area
    Point border1 = Point(0, -0.1, -0.1);
    Point border2 = Point(1, 0.2, 0.2);

    double grid_step_x = 0.04;
    double grid_step_y = 0.04;
    double grid_step_z = 0.04;

    static int DBG_ID;

    Params(Params const &) = delete;
    Params & operator=(Params const &) = delete;

private:
    Params() = default;

    ~Params() = default;
};
