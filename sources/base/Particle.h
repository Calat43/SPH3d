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
