#include <fstream>
#include <iostream>

#include "Params.h"

#include "Grid.h"
#include "Cell.h"
#include "Viscosity.h"
#include "SPHSolver.h"

#include "SodTube1d.h"

double find_density(Particle * particle, Cell * cell)
{
    double density = 0;
    int nears = 0;

    //Params & params = Params::get_instance();

    for (Cell * neighbour : cell->get_neighbours())
    {
        // std::cout << *neighbour << std::endl;
        for (Particle & p : particle->kind == Particle::Kind::Gas ? neighbour->gas_particles : neighbour->dust_particles)
        { // TODO remove direct access
            density += p.mass * kernel(*particle, p, Params::get_instance().dimensions);

            if(PRINT_STUFF)
            {
                if(kernel(*particle, p, Params::get_instance().dimensions) != 0)
                {
                    nears += 1;
                }
            }
        }
    }

    if(PRINT_STUFF)
    {
        std::cout << particle->get_id() << " " << nears << std::endl;
    }
    assert(!__isnan(density));

    return density;
}

double find_density_no_sort(Particle const & particle, Grid const & grid)
{
    double density = 0;
    int nears = 0;

    for(int i = 0; i < grid.x_size; ++i)
    {
        for(int j = 0; j < grid.y_size; ++j)
        {
            for(int k = 0; k < grid.z_size; ++k)
            {
                for(Particle const & p : particle.kind == Particle::Kind::Gas ? grid.cells[i][j][k].gas_particles
                                                                      : grid.cells[i][j][k].dust_particles)
                {
                    density += p.mass * kernel(particle, p, Params::get_instance().dimensions);
                    if(kernel(particle, p, Params::get_instance().dimensions) != 0)
                    {
                        nears += 1;
                    }
                }
            }
        }
    }

    return density;
}

void recalc_density(Grid & grid, Particle::Kind kind)
{
    for(int i = 0; i < grid.x_size; ++i)
    {
        for(int j = 0; j < grid.y_size; ++j)
        {
            for(int k = 0; k < grid.z_size; ++k)
            {
                Cell & cell = grid.cells.at(i).at(j).at(k);
                for(Particle & particle : cell.particles_of_kind(kind))
                {
                    particle.density = find_density(&particle, &(cell));
                    //particle.density = find_density_no_sort(particle, grid);
                    assert(!__isnan(particle.density));
                }
            }
        }
    }
}

Point find_new_velocity(Particle * particle, Cell * cell)
{
    Point sum = {0, 0, 0};
    Point vel = {0, 0, 0};
    Point prev_vel = {particle->vx, particle->vy, particle->vz};
    Point kernel_grad = {0, 0, 0};

    Params & params = Params::get_instance();

    assert(particle->density != 0);

    for (Cell * neighbour : cell->get_neighbours())
    {
        for (Particle & p : particle->kind == Particle::Kind::Gas ? neighbour->gas_particles : neighbour->dust_particles)
        { // TODO remove direct access
            assert(!__isnan(p.density));
            assert(p.density != 42);

            double term1 = (1. / p.density + 1. / particle->density);

            kernel_grad = kernel_gradient(*(particle), p, params.dimensions);

            sum = sum + kernel_grad * term1;
        }
    }

    double term2 = - params.tau * particle->mass * params.gas_sound_speed * params.gas_sound_speed;

    vel = sum * term2 + prev_vel;

    assert(!__isnan(vel.x));
    assert(!__isnan(vel.y));
    assert(!__isnan(vel.z));

    return vel;
}

Point find_new_velocity_no_sort(Particle * particle, Grid * grid)
{
    Point sum = {0, 0, 0};
    Point vel = {0, 0, 0};
    Point prev_vel = {particle->vx, particle->vy, particle->vz};
    Point kernel_grad = {0, 0, 0};

    Params & params = Params::get_instance();

    for(int i = 0; i < grid->x_size; ++i)
    {
        for(int j = 0; j < grid->y_size; ++j)
        {
            for(int k = 0; k < grid->z_size; ++k)
            {
                for(Particle & p : particle->kind == Particle::Kind::Gas ? grid->cells[i][j][k].gas_particles
                                                                      : grid->cells[i][j][k].dust_particles)
                {
                    double term1 = (1. / p.density + 1. / particle->density);

                    kernel_grad = kernel_gradient(*(particle), p, params.dimensions);

                    sum = sum + kernel_grad * term1;
                }
            }
        }
    }

    double term2 = - params.tau * particle->mass * params.gas_sound_speed * params.gas_sound_speed;

    vel = sum * term2 + prev_vel;

    assert(!__isnan(vel.x));

    return vel;
}

//particle from prev step
Point Sod_tube_1d::find_new_velocity(Particle * particle, Cell * cell)
{
    double sum_x = 0;

    double vx = 0;

    Params & params = Params::get_instance();

    assert(particle->density != 0);

    for (Cell * neighbour : cell->get_neighbours())
    {
        for (Particle & p : particle->kind == Particle::Kind::Gas ? neighbour->gas_particles : neighbour->dust_particles)
        { // TODO remove direct access

            double kernel_grad = kernel_gradient(*(particle), p, params.dimensions).x;
            assert(!__isnan(p.density));
            assert(!__isnan(kernel_grad));
            sum_x += (particle->pressure / pow(particle->density, 2) + p.pressure / pow(p.density, 2)
                    + viscosity_1d::find_viscosity(particle, &p))
                    * kernel_grad;
        }
    }

    vx = particle->vx - params.tau * particle->mass * sum_x;

    return {vx, 0, 0};
}

double Sod_tube_1d::find_new_pressure(Particle * particle)
{
    double result = particle->density * particle->energy * (Params::get_instance().gamma - 1);
    assert(result >= 0);
    return result;
}

void Sod_tube_1d::recalc_pressure(Grid & grid, Particle::Kind kind)
{
    for(int i = 0; i < grid.x_size; ++i)
    {
        for(int j = 0; j < grid.y_size; ++j)
        {
            for(int k = 0; k < grid.z_size; ++k)
            {
                Cell & cell = grid.cells[i][j][k];
                for(Particle & particle : cell.particles_of_kind(kind))
                {
                    particle.pressure = Sod_tube_1d::find_new_pressure(&particle);
                    assert(!__isnan(particle.pressure));
                }
            }
        }
    }
}

//particle from prev step
double Sod_tube_1d::find_new_energy(Particle * particle, Cell * cell)
{
    double pres_member = 0;
    double visc_member = 0;
    double result = 0;
    double kern_grad_x = 0;

    for (Cell * neighbour : cell->get_neighbours())
    {
        for (Particle & p : particle->kind == Particle::Kind::Gas ? neighbour->gas_particles : neighbour->dust_particles)
        {
            kern_grad_x = kernel_gradient(*(particle), p, Params::get_instance().dimensions).x;
            pres_member += (particle->vx - p.vx) * kern_grad_x;
            visc_member += (particle->vx - p.vx) * kern_grad_x *
                            viscosity_1d::find_viscosity(particle, &p);
        }
    }

    result = particle->pressure /  pow(particle->density, 2) * particle->mass * Params::get_instance().tau * pres_member +
            Params::get_instance().tau * particle->mass / 2. * visc_member + particle->energy;
    assert(visc_member >= 0);
    assert(result >= 0);

    return result;
}

Grid Sod_tube_1d::do_time_step(Grid & old_grid, int step_num)
{
    std::ofstream part_file;
    if(PRINT_FILE)
    {
        part_file.open(OUTPUT_PATH + "part_" + std::to_string(step_num) + ".dat");
    }

    Params & params = Params::get_instance();

    Grid next_grid(params.grid_step_x, params.grid_step_y, params.grid_step_z, params.border1, params.border2);

    for (int i = 0; i < old_grid.x_size; ++i)
    { // TODO foreach
        for(int j = 0; j < old_grid.y_size; ++j)
        {
            for(int k = 0; k < old_grid.z_size; ++k)
            {
                Cell & cell = old_grid.cells[i][j][k];
                for(Particle * particle : cell.get_all_particles())
                {
                    if(PRINT_FILE)
                    {
                        part_file << particle->x << " " << particle->vx << " " << particle->density << " "
                                  << particle->pressure << " " << particle->energy
                                  << std::endl;
                    }

                    Particle new_particle(*particle);
                    Point new_coords = find_new_coordinates_(*particle);

                    Point new_vel = Sod_tube_1d::find_new_velocity(particle, &cell);
                    new_particle.density = NAN;

                    new_particle.set_coordinates(new_coords.x, new_coords.y, new_coords.z);
                    new_particle.set_velocities(new_vel.x, new_vel.y, new_vel.z);

                    new_particle.energy = find_new_energy(particle, &cell);

                    next_grid.with_copy_of(new_particle);
                }
            }
        }
    }

    recalc_density(next_grid, Particle::Kind::Gas);
    recalc_pressure(next_grid, Particle::Kind::Gas);

    return next_grid;
}