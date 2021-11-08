#pragma once

#include <string>
#include <cmath>

#include "Point.h"
#include "MathUtils.h"

const std::string OUTPUT_PATH = "/home/calat/documents/SPH3D-reworked/output/";

const bool PRINT_DENSITY = false; // TODO restored from "unmerged..."

const bool PRINT_STUFF = false;

const bool PRINT_FILE = true;

const bool IS_IDIC = true; // TODO restored from "unmerged..."

// Singleton
class Params
{
public:
    static Params & get_instance()
    {
        static Params p; // lazy initialized
        return p;
    }

    ///////////// TODO restored from "unmerged..."
    double t = 0.2; // total time of numerical experiment
    double t_diagnostics = 0.05; // time interval to save intermediate results

    double c_s = 1;
    double gamma = 4. / 3.;
    static const int n_gas = 10000;
    int n_dust = 0;
    int amount = 600;
    double delta = pow(10, -2);
    ////////////

    int dimensions = 1;
    double h = 0.01;
    double smooth_radius = 2. * h;
    double tau = 0.001;

    double gas_sound_speed = 1.;
    double dust_sound_speed = 0.;
    double a = 5. * pow(10, -6);
    double K = 0.01;
    //TODO delete from t_stop_asterisk
    double t_stop = 1;
    double d2g = 0.01;
    double middle_rho_gas = 1;
    double theta_const = 0.001;
    double particle_density = middle_rho_gas * d2g / theta_const;

    double dens_coeff1 = 0.01;
    double dens_coeff2 = 0;
    double dens_coeff3 = 0.0000122;
    double dens_coeff4 = -0.000014;

    double gas_vel_coeff1 = 0.0099995;
    double gas_vel_coeff2 = 0.0000063;
    double dust_vel_coeff1 = 0.0012238;
    double dust_vel_coeff2 = -0.0013977;

    //ball in vacuum
    double A = 1.;
    double ball_mass = 100.;
    double ball_init_radius = 0.5;
    double ball_init_velocity = 0.1;
    //TODO remove
//    int n_gas = 1000;
//    int n_dust = 1000;

    //viscosity
    double alpha = 1.;
    double beta = 2.;
    double nu = 0.1 * h;
    bool have_viscosity = false;

    //monaghan // TODO restored from "unmerged..."
    double sigma = 1. / dimensions;
    double eta_squared = 0.001 * h * h;

    // boundary of area
    Point border1 = Point(-1, 0, 0);
    Point border2 = Point(2, 0, 0);
    //Point border1 = Point(1, 1, 1);
    //Point border2 = Point(5, 5, 5);

    double grid_step_x = h / 2;
    double grid_step_y = 0;
    double grid_step_z = 0;

    static int DBG_ID;

    Params(Params const &) = delete;
    Params & operator=(Params const &) = delete;

    enum DistributionType // TODO restored from "unmerged..."
    {
        dtUniform = 0,
        dtBall = 1,
    };
    DistributionType dt = dtUniform;

    enum BoundaryConditions // TODO restored from "unmerged..."
    {
        bcIsolated = 0, // particles do not leave global computational domain and reflect from boundaries
        bcPeriodic = 1, // particles move as in periodic global computational domain
    };
    BoundaryConditions bc = bcIsolated;

    size_t sz_gas_particles_total; // total number of gas particles for all processors
    size_t sz_dust_particles_total; // total number of dust particles for all processors

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
    uint gas_real_particles = 45000;
    uint gas_image_particles = 450;
    uint gas_yz_particles = 0;

    double gas_dens_left = 1.;
    double gas_press_left = 1.;
    double gas_vel_left = 0.;
    double gas_ener_left = 2.5;

    double gas_dens_right = 0.125;
    double gas_press_right = 0.1;
    double gas_vel_right = 0.;
    double gas_ener_right = 2.;

    //dust
    uint dust_real_particles = 45000;
    uint dust_image_particles = 450;
    uint dust_yz_particles = 0;

    double dust_dens_left = 1.46;
    double dust_vel_left = 0.;

    double dust_dens_right = 0.125 * 1.46;
    double dust_vel_right = 0.;

    Dusty_shock_params(Params const &) = delete;
    Dusty_shock_params & operator=(Dusty_shock_params const &) = delete;

private:
    Dusty_shock_params() = default;

    ~Dusty_shock_params() = default;
};