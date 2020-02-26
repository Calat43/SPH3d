
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

    ////TODO remove after debug
    right_sum_x = that.right_sum_x;
    right_sum_y = that.right_sum_y;
    right_sum_z = that.right_sum_z;

    right_an_x = that.right_an_x;
    right_an_y = that.right_an_y;
    right_an_z = that.right_an_z;
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
    dbg_state = that.dbg_state;

    ////TODO remove after debug
    right_sum_x = that.right_sum_x;
    right_sum_y = that.right_sum_y;
    right_sum_z = that.right_sum_z;

    right_an_x = that.right_an_x;
    right_an_y = that.right_an_y;
    right_an_z = that.right_an_z;
    ////

    return *this;
}

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
