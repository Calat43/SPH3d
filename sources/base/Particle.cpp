
#include <iostream>

#include "Particle.h"

Particle::Particle(Particle::Kind kind, double x, double y, double z) : Numbered()
{
    this->kind = kind;
    this->x = x;
    this->y = y;
    this->z = z;
}

Particle::Particle(Particle const & that) : Numbered(that) //WHY Numbered(that)??
{
    kind = that.kind;
    is_border = that.is_border;
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
    dbg_state = that.dbg_state;
}

Particle & Particle::operator=(Particle const & that) = default;

void Particle::set_coordinates(double new_x, double new_y, double new_z)
{ // TODO maybe remove
    x = new_x;
    y = new_y;
    z = new_z;
}

void Particle::set_velocities(double new_vx, double new_vy, double new_vz)
{ // TODO maybe remove
    vx = new_vx;
    vy = new_vy;
    vz = new_vz;
}

std::ostream & operator<<(std::ostream & stream, Particle const & p)
{
    stream << "Particle#" << p.get_id() << "(" << p.x << ", " << p.y << ", " << p.z << ")";
    return stream;
}
