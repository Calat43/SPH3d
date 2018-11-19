#include <iostream>
#include "Particle.h"

Particle::Particle(Particle::Kind kind, double x, double y, double z) : Numbered()
{
    this->kind = kind;
    this->x = x;
    this->y = y;
    this->z = z;
}

Particle::Particle(Particle const & that) : Numbered(that)
{
    kind = that.kind;
    x = that.x;
    y = that.y;
    z = that.z;
    mass = that.mass;
    pressure = that.pressure;
    density = that.density;
    energy = that.energy;
    vx = that.vx;
    vy = that.vy;
    vz = that.vz;
}

Particle & Particle::operator=(Particle const & that)
{
    Numbered::operator=(that);
    kind = that.kind;
    x = that.x;
    y = that.y;
    z = that.z;
    mass = that.mass;
    pressure = that.pressure;
    density = that.density;
    energy = that.energy;
    vx = that.vx;
    vy = that.vy;
    vz = that.vz;
    return *this;
}

void Particle::set_coordinates(double new_x, double new_y, double new_z) { // TODO maybe remove
    x = new_x;
    y = new_y;
    z = new_z;
}

std::ostream & operator<<(std::ostream & stream, Particle const & p)
{
    stream << "Particle#" << p.get_id() << "(" << p.x << ", " << p.y << ", " << p.z << ")";
    return stream;
}
