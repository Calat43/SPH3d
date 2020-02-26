
#include <fstream>
#include <iostream>

#include "Params.h"
#include "Cell.h"
#include "Grid.h"
#include "Viscosity.h"
#include "Solver.h"
#include "SodTube3d.h"

Point Sod_tube_3d::find_new_velocity(Particle * particle, Cell * cell, int step_num)
{
    Point sum = {0, 0, 0};
    Point vel = {0, 0, 0};
    Point prev_vel = {particle->vx, particle->vy, particle->vz};
    Point kernel_grad = {0, 0, 0};
    Point kernel_grad_ba = {0, 0, 0};

    Params & params = Params::get_instance();

    assert(particle->density != 0);
    assert(!__isnan(particle->density));

    int nears = 0;

    //полуаналитика для p = 1 - 100 * r * r;
/*
    double center_x = (params.border2.x + params.border1.x) / 2;
    double center_y = (params.border2.y + params.border1.y) / 2;
    double center_z = (params.border2.z + params.border1.z) / 2;

    double coeff = 2. * params.tau * 100. / particle->density;
    vel.x = coeff * (particle->x - center_x) + particle->vx;
    vel.y = coeff * (particle->y - center_y) + particle->vy;
    vel.z = coeff * (particle->z - center_z) + particle->vz;
*/

    assert(cell->grid->find_cell(*particle) == cell);

    for (Cell * neighbour : cell->get_neighbours())
    {
        for (Particle & p : particle->kind == Particle::Kind::Gas ? neighbour->gas_particles : neighbour->dust_particles)
        { // TODO remove direct access

            kernel_grad = kernel_gradient(*(particle), p, params.dimensions);

            assert(!__isnan(p.density));
            assert(!__isnan(kernel_grad.x));

            double term1 = (particle->pressure / pow(particle->density, 2) + p.pressure / pow(p.density, 2));

            //sum = sum + kernel_grad * term1;
            sum.x += kernel_grad.x * term1;
            sum.y += kernel_grad.y * term1;
            sum.z += kernel_grad.z * term1;
        }
    }


    double term2 = params.tau * particle->mass;

    vel.x = prev_vel.x - term2 * sum.x;
    vel.y = prev_vel.y - term2 * sum.y;
    vel.z = prev_vel.z - term2 * sum.z;

    //vel = prev_vel - sum * term2;
    return vel;
}

double Sod_tube_3d::find_new_energy(Particle * particle, Cell * cell)
{
    double pres_member = 0;
    double visc_member = 0;
    double result = 0;

    for (Cell * neighbour : cell->get_neighbours())
    {
        for (Particle & p : particle->kind == Particle::Kind::Gas ? neighbour->gas_particles : neighbour->dust_particles)
        {
            //TODO calculate one time

            double vel_dot_gr_k = vel_dot_grad_kernel(*particle, p);

            pres_member += vel_dot_gr_k;
            visc_member += vel_dot_gr_k * viscosity_3d::find_viscosity(particle, &p);
        }
    }

    result = particle->pressure /  pow(particle->density, 2) * particle->mass * Params::get_instance().tau * pres_member +
             Params::get_instance().tau * particle->mass / 2. * visc_member + particle->energy;
    assert(visc_member >= 0);
    assert(result >= 0);

    return result;
}