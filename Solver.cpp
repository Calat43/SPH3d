#include "Solver.h"

Point find_new_coordinates(Particle const & particle)
{
    Params & params = Params::get_instance();
    double x = particle.x + params.tau * particle.vx;
    double y = particle.y + params.tau * particle.vy;
    double z = particle.z + params.tau * particle.vz;
    return {x, y, z};
}

double find_density(Particle * particle, Cell * cell)
{
    double density = 0;
    int nears = 0;

    Params & params = Params::get_instance();
    double center_x = (params.border2.x - params.border1.x) / 2;
    double center_y = (params.border2.y - params.border1.y) / 2;
    double center_z = (params.border2.z - params.border1.z) / 2;

    // std::cout << "--------" << std::endl;
    // std::cout << *particle  << " in " << *cell << std::endl;
    // std::cout << "--------" << std::endl;
    for (Cell * neighbour : cell->get_neighbours())
    {
        // std::cout << *neighbour << std::endl;
        for (Particle & p : particle->kind == Particle::Kind::Gas ? neighbour->gas_particles : neighbour->dust_particles)
        { // TODO remove direct access
            density += p.mass * kernel(*particle, p, Params::get_instance().dimensions);

            if(kernel(*particle, p, Params::get_instance().dimensions) != 0)
            {
                nears += 1;
            }
        }
    }

    if((pow((particle->x - center_x), 2) + pow((particle->y - center_y), 2) + pow((particle->z - center_z), 2)) <= 0.0025)
    {
        std::cout << particle->get_id() << " " << nears << std::endl;
    }

    return density;
}

double find_density_no_sort(Particle const & particle, Grid const & grid)
{
    double density = 0;

    int nears = 0;

    Params & params = Params::get_instance();
    double center_x = (params.border2.x - params.border1.x) / 2;
    double center_y = (params.border2.y - params.border1.y) / 2;
    double center_z = (params.border2.z - params.border1.z) / 2;

    for(int i = 0; i < grid.x_size; ++i)
    {
        for(int j = 0; j < grid.y_size; ++j)
        {
            for(int k = 0; k < grid.z_size; ++k)
            {
                for(Particle p : particle.kind == Particle::Kind::Gas ? grid.cells[i][j][k].gas_particles
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

    if((pow((particle.x - center_x), 2) + pow((particle.y - center_y), 2) + pow((particle.z - center_z), 2)) <= 0.0025)
    {
        std::cout << particle.get_id() << " " << nears << std::endl;
    }


    return density;
}

void recalc_density(ParticlesState & state, double radius)
{

    Params & params = Params::get_instance();
    double center_x = (params.border2.x - params.border1.x) / 2;
    double center_y = (params.border2.y - params.border1.y) / 2;
    double center_z = (params.border2.z - params.border1.z) / 2;

    Grid & grid = state.grid;
    for(int i = 0; i < grid.x_size; ++i)
    {
        for(int j = 0; j < grid.y_size; ++j)
        {
            for(int k = 0; k < grid.z_size; ++k)
            {
                Cell & cell = grid.cells[i][j][k];
                for(Particle * particle : cell.get_all_particles())
                {
                    if((pow((particle->x - center_x), 2) + pow((particle->y - center_y), 2) +
                        pow((particle->z - center_z), 2)) <= radius * radius)
                    {
                        particle->density = find_density(particle, &(cell));
                        //particle->density = find_density_no_sort(*particle, state.grid);
                    }
                    assert(!__isnan(particle->density));
                }
            }
        }
    }
}


Point find_new_velocity(Particle * particle, Cell * cell)
{
    double sum_x = 0;
    double sum_y = 0;
    double sum_z = 0;

    double vx = 0;
    double vy = 0;
    double vz = 0;

    Params & params = Params::get_instance();

    assert(particle->density != 0);

    for (Cell * neighbour : cell->get_neighbours())
    {
        for (Particle & p : particle->kind == Particle::Kind::Gas ? neighbour->gas_particles : neighbour->dust_particles)
        { // TODO remove direct access
            assert(!__isnan(p.density));
            assert(p.density != 42);
            assert(!__isnan(kernel_gradient_x(*(particle), p, params.dimensions)));
            assert(!__isnan(kernel_gradient_y(*(particle), p, params.dimensions)));
            assert(!__isnan(kernel_gradient_z(*(particle), p, params.dimensions)));
            sum_x += (1. / p.density + 1. / particle->density) * kernel_gradient_x(*(particle), p, params.dimensions);
            sum_y += (1. / p.density + 1. / particle->density) * kernel_gradient_y(*(particle), p, params.dimensions);
            sum_z += (1. / p.density + 1. / particle->density) * kernel_gradient_z(*(particle), p, params.dimensions);
        }
    }

    vx = - params.tau * particle->mass * params.c_s * params.c_s * sum_x + particle->vx;
    vy = - params.tau * particle->mass * params.c_s * params.c_s * sum_y + particle->vy;
    vz = - params.tau * particle->mass * params.c_s * params.c_s * sum_z + particle->vz;

    assert(!__isnan(vx));
    assert(!__isnan(vy));
    assert(!__isnan(vz));

    return {vx, vy, vz};
}

Point find_new_velocity_no_sort(Particle * particle, Grid * grid)
{
    double sum_x = 0;
    double sum_y = 0;
    double sum_z = 0;

    double vx = 0;
    double vy = 0;
    double vz = 0;

    Params & params = Params::get_instance();

    for(int i = 0; i < grid->x_size; ++i)
    {
        for(int j = 0; j < grid->y_size; ++j)
        {
            for(int k = 0; k < grid->z_size; ++k)
            {
                for(Particle p : particle->kind == Particle::Kind::Gas ? grid->cells[i][j][k].gas_particles
                                                                      : grid->cells[i][j][k].dust_particles)
                {
                    sum_x += (1. / p.density + 1. / particle->density) * kernel_gradient_x(*(particle), p, params.dimensions);
                    sum_y += (1. / p.density + 1. / particle->density) * kernel_gradient_y(*(particle), p, params.dimensions);
                    sum_z += (1. / p.density + 1. / particle->density) * kernel_gradient_z(*(particle), p, params.dimensions);
                }
            }
        }
    }

    vx = - params.tau * particle->mass * params.c_s * params.c_s * sum_x + particle->vx;
    vy = - params.tau * particle->mass * params.c_s * params.c_s * sum_y + particle->vy;
    vz = - params.tau * particle->mass * params.c_s * params.c_s * sum_z + particle->vz;

    assert(!__isnan(vx));

    return {vx, vy, vz};
}


