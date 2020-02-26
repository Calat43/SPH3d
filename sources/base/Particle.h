#pragma once

#include <string>
#include <cmath>
#include "Numbered.h"

class Particle : public Numbered
{
public:
    enum Kind {
        Gas,
        Dust
    };

    int dbg_state = 0;

    Kind kind;

    double vx = NAN;
    double vy = NAN;
    double vz = NAN;

    double mass = NAN;

    double pressure = NAN;
    double density = NAN;
    double energy = NAN;

    ////TODO remove after debug
    double right_sum_x = 1024;
    double right_sum_y = 1024;
    double right_sum_z = 1024;

    double right_an_x = 2048;
    double right_an_y = 2048;
    double right_an_z = 2048;
    /////

    double x = NAN;
    double y = NAN;
    double z = NAN;
    
    Particle(Kind kind, double x, double y, double z);

    Particle(Particle const & that);

    Particle & operator=(Particle const & that);

    void set_coordinates(double new_x, double new_y, double new_z); // TODO maybe remove
    void set_velocities(double new_vx, double new_vy, double new_vz); // TODO maybe remove

    friend std::ostream & operator<<(std::ostream & stream, Particle const & p);

};
