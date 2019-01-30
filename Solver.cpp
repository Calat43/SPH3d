#include "Solver.h"

namespace viscosity_1d
{

    double find_sound_speed(Particle * particle)
    {
        assert(!__isnan(sqrt(Params::get_instance().gamma * particle->pressure / particle->density)));
        return sqrt(Params::get_instance().gamma * particle->pressure / particle->density);
    }

    double find_mu(double vel_ab, double coord_ab)
    {
        return (Params::get_instance().h * vel_ab * coord_ab) /
               (pow(coord_ab, 2) + pow(Params::get_instance().nu, 2));
    }

    double find_viscosity(Particle * p1, Particle * p2)
    {
        Params & params = Params::get_instance();

        if(params.have_viscosity)
        {
            double vel_ab = p1->vx - p2->vx;
            double coord_ab = p1->x - p2->x;

            if(vel_ab * coord_ab >= 0)
            {
                return 0;
            }
            if(vel_ab * coord_ab < 0)
            {
                double c_a = find_sound_speed(p1);
                double c_b = find_sound_speed(p2);

                double dens_ab = 1. / 2 * (p1->density + p2->density);
                double c_ab = 1. / 2 * (c_a + c_b);

                double mu_ab = find_mu(vel_ab, coord_ab);

                assert(!__isnan(c_a));
                assert(!__isnan(c_b));
                assert(!__isnan(mu_ab));

                return (- params.alpha * c_ab * mu_ab + params.beta * pow(mu_ab, 2)) / dens_ab;
            }
        }
        return 0;
    }

}

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

//    if((pow((particle->x - center_x), 2) + pow((particle->y - center_y), 2) + pow((particle->z - center_z), 2)) <= 0.0025)
//    {
//        std::cout << particle->get_id() << " " << nears << std::endl;
//    }

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

void recalc_density(Grid & grid)
{
    for(int i = 0; i < grid.x_size; ++i)
    {
        for(int j = 0; j < grid.y_size; ++j)
        {
            for(int k = 0; k < grid.z_size; ++k)
            {
                Cell & cell = grid.cells[i][j][k];
                for(Particle * particle : cell.get_all_particles())
                {
                    particle->density = find_density(particle, &(cell));
                    //particle->density = find_density_no_sort(*particle, state.grid);
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

//particle from prev step
Point Sod_tube::find_new_velocity(Particle * particle, Cell * cell)
{
    double sum_x = 0;

    double vx = 0;

    Params & params = Params::get_instance();

    assert(particle->density != 0);

    for (Cell * neighbour : cell->get_neighbours())
    {
        for (Particle & p : particle->kind == Particle::Kind::Gas ? neighbour->gas_particles : neighbour->dust_particles)
        { // TODO remove direct access
            assert(!__isnan(p.density));
            assert(!__isnan(kernel_gradient_x(*(particle), p, params.dimensions)));
            sum_x += (particle->pressure / pow(particle->density, 2) + p.pressure / pow(p.density, 2)
                    + viscosity_1d::find_viscosity(particle, &p))
                    * kernel_gradient_x(*(particle), p, params.dimensions);
        }
    }

    vx = particle->vx - params.tau * particle->mass * sum_x;

    return {vx, 0, 0};
}

double Sod_tube::find_new_pressure(Particle * particle)
{
    return particle->density * particle->energy * (Params::get_instance().gamma - 1);
}

void Sod_tube::recalc_pressure(Grid & grid)
{
    for(int i = 0; i < grid.x_size; ++i)
    {
        for(int j = 0; j < grid.y_size; ++j)
        {
            for(int k = 0; k < grid.z_size; ++k)
            {
                Cell & cell = grid.cells[i][j][k];
                for(Particle * particle : cell.get_all_particles())
                {
                    particle->pressure = Sod_tube::find_new_pressure(particle);
                    //particle->density = find_density_no_sort(*particle, state.grid);
                    assert(!__isnan(particle->pressure));
                }
            }
        }
    }
}

//particle from prev step
double Sod_tube::find_new_energy(Particle * particle, Cell * cell)
{
    double pres_member = 0;
    double visc_member = 0;
    double result = 0;

    for (Cell * neighbour : cell->get_neighbours())
    {
        for (Particle & p : particle->kind == Particle::Kind::Gas ? neighbour->gas_particles : neighbour->dust_particles)
        {
            pres_member += (particle->vx - p.vx) * kernel_gradient_x(*(particle), p, Params::get_instance().dimensions);
            visc_member += (particle->vx - p.vx) * kernel_gradient_x(*(particle), p, Params::get_instance().dimensions) *
                            viscosity_1d::find_viscosity(particle, &p);
        }
    }

    result = particle->pressure /  pow(particle->density, 2) * particle->mass * Params::get_instance().tau * pres_member +
            Params::get_instance().tau * particle->mass / 2. * visc_member + particle->energy;
    assert(visc_member >= 0);
    assert(result >= 0);

    return result;
}

Grid Sod_tube::do_time_step(Grid & old_grid, int step_num)
{
    char filename[512];
    sprintf(filename, "/home/calat/tmp/part_%0d.dat", step_num);
    FILE * f = fopen(filename, "w");

    Params & params = Params::get_instance();

    Grid next_grid(params.grid_step_x, params.grid_step_y, params.grid_step_z, params.border1, params.border2);

    for (int i = 0; i < old_grid.x_size; ++i)
    { // TODO foreach
        for(int j = 0; j < old_grid.y_size; ++j)
        {
            for(int k = 0; k < old_grid.z_size; ++k)
            {
                Cell cell = old_grid.cells[i][j][k];
                for(Particle * particle : cell.get_all_particles())
                {
                    fprintf(f, "%lf %lf %lf %lf %lf\n", particle->x, particle->vx, particle->density,
                            particle->pressure, particle->energy);

                    Particle new_particle(*particle);
                    Point new_coords = find_new_coordinates(*particle);

                    Point new_vel = Sod_tube::find_new_velocity(particle, &cell);
                    new_particle.density = NAN;

                    new_particle.set_coordinates(new_coords.x, new_coords.y, new_coords.z);
                    new_particle.set_velocities(new_vel.x, new_vel.y, new_vel.z);

                    new_particle.energy = find_new_energy(particle, &cell);

                    next_grid.with_copy_of(new_particle);
                }
            }
        }
    }

    recalc_density(next_grid);
    recalc_pressure(next_grid);

    return next_grid;
}