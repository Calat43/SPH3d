#pragma once

#include <string>

#include "Point.h"
#include "MathUtils.h"

const std::string OUTPUT_PATH = "/home/calat/documents/output/";

const bool PRINT_STUFF = false;

const bool PRINT_FILE = true;

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
    double gamma = 4. / 3.;
    int dimensions = 1;
    double h = 0.1;
    double smooth_radius = 2. * h;
    double tau = 0.01;

    double sound_speed = 1.;
    double a = 5. * pow(10, -6);
    double particle_density = 0.1;

    //ball in vacuum
    double A = 1.;
    double ball_mass = 1.;
    double ball_init_radius = 1.;
    double ball_init_velocity = 0.1;
    //TODO remove
    int n_gas = 1000;
    int n_dust = 1000;

    //viscosity
    double alpha = 1.;
    double beta = 2.;
    double nu = 0.1 * h;
    bool have_viscosity = false;

    // boundary of area
    Point border1 = Point(-2, -2, -2);
    Point border2 = Point(2, 2, 2);
    //Point border1 = Point(1, 1, 1);
    //Point border2 = Point(5, 5, 5);

    double grid_step_x = 2. * h;
    double grid_step_y = 2. * h;
    double grid_step_z = 2. * h;

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