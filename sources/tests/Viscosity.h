#pragma once

class Particle;
class Point;

namespace viscosity_3d
{
    double find_viscosity(Particle * p1, Particle * p2);
}

namespace viscosity_1d
{
    //double find_mu(double vel_ab, double coord_ab); // no function definition
    double find_viscosity(Particle * p1, Particle * p2);
}

namespace viscosity
{
    double find_sound_speed(Particle * particle);

    double find_mu(Point vel_ab, Point coord_ab);
    double find_mu(double vel_ab, double coord_ab);

    //double find_viscosity(Particle * p1, Particle * p2);
}