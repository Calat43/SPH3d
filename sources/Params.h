#pragma once

#include <cmath>
#include "Point.h"
#include "Common.h"

const std::string OUTPUT_PATH = "./output/";

const bool PRINT_DENSITY = false;

const bool PRINT_STUFF = false;

const bool PRINT_FILE = true;

const bool IS_IDIC = true;

// Singleton
class Params
{
public:
    static Params & get_instance()
    {
        static Params p; // lazy initialized
        return p;
    }

    double t = 0.2;
    double c_s = 1;
    double gamma = 4. / 3.; //7. / 5.;
    double K = 0;
    double kx = NAN;
    double ky = NAN;
    double kz = NAN;
    static const int n_gas = 10000;
    int n_dust = 0;
    int dimensions = 1;
    double h = 0.01;
    double smooth_radius = 2. * h;
    double tau = 0.001;

    double middle_gas_dens = NAN;
    double d2g = NAN;
    double delta = NAN;

    //non-linear drag
    double a = 5. * pow(10, -3);
    double particle_density = 2.3;

    //ball in vacuum
    double A = 1.;
    double ball_mass = 1.;
    double ball_init_radius = 0.1;
    double ball_init_velocity = 0.1;

    //viscosity
    double alpha = 1.;
    double beta = 2.;
    double nu = 0.1 * h;
    bool have_viscosity = true;

    //monaghan
    double sigma = 1. / 3.;
    double eta_squared = 0.001 * h * h;

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

class Dusty_shock_params
{
public:
    static Dusty_shock_params & get_instance()
    {
        static Dusty_shock_params p; // lazy initialized
        return p;
    }

    double membrane = 0.5;
    double x_left = 0.;
    double x_right = 1.;
    double yz_left = 0.;
    double yz_right = 0.1;

    //gas
    uint gas_real_particles = 70;
    uint gas_image_particles = 16;
    uint gas_yz_particles = 16;

    double gas_dens_left = 1.;
    double gas_press_left = 1.;
    double gas_vel_left = 0.;
    double gas_ener_left = 2.5;

    double gas_dens_right = 1.;
    double gas_press_right = 0.8;
    double gas_vel_right = 0.;
    double gas_ener_right = 2.;

    //dust
    uint dust_real_particles = 70;
    uint dust_image_particles = 16;
    uint dust_yz_particles = 16;

    double dust_dens_left = 1.;
    double dust_vel_left = 0.;

    double dust_dens_right = 1.;
    double dust_vel_right = 0.;

    Dusty_shock_params(Params const &) = delete;
    Dusty_shock_params & operator=(Dusty_shock_params const &) = delete;

private:
    Dusty_shock_params() = default;

    ~Dusty_shock_params() = default;
};
