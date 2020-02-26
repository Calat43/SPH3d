#pragma once

#include <string>

#include "Point.h"
#include "MathUtils.h"

const std::string OUTPUT_PATH = "/home/calat/documents/SPH3D-debug/output/";

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

    double t = 0.3;
    double gamma = 4. / 3.;
    int dimensions = 3;
    double h = 0.1;
    double smooth_radius = 2. * h;
    double tau = 0.01;

    double sound_speed = 1.;

    //ball in vacuum
    double A = 1.;
    double ball_mass = 1.;
    double ball_init_radius = 1.;
    double ball_init_velocity = 0.1;

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