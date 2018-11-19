#include <iostream>
#include <string>
#include "Particle.h"
#include "ParticlesState.h"
#include "Common.h"
#include "Grid.h"
#include "Cell.h"
#include "Params.h"

Particle new_random_particle(Particle::Kind kind, Point border1, Point border2) {

    double radius = (border2.x - border1.x) / 50;

    double center_x = (border2.x - border1.x) / 2;
    double center_y = (border2.y - border1.y) / 2;
    double center_z = (border2.z - border1.z) / 2;

    double x = center_x + random_double(-radius, radius);
    double y = center_y + random_double(-radius, radius);
    double z = center_z + random_double(-radius, radius);

    double factor = random_double(0, 0.2);
    double vx = (x - center_x) * factor / 100.;
    double vy = (y - center_y) * factor / 100.;
    double vz = (z - center_z) * factor / 100.;

    Particle particle(kind, x, y, z);

    particle.vx = vx;
    particle.vy = vy;
    particle.vz = vz;

    return particle;
}

ParticlesState random_init_state()
{
    ParticlesState init_state;

    Params & params = Params::get_instance();

    for (int i = 0; i < params.n_gas; ++i)
    {
        Particle particle = new_random_particle(Particle::Kind::Gas, init_state.grid.border1, init_state.grid.border2);
        init_state.with_copy_of(particle);
    }
    for (int i = 0; i < params.n_dust; ++i)
    {
        Particle particle = new_random_particle(Particle::Kind::Dust, init_state.grid.border1, init_state.grid.border2);
        init_state.with_copy_of(particle);
        
    }
    return init_state;
}

double find_density(Particle * particle, Cell * cell)
{
    double density = 0;
    // std::cout << "--------" << std::endl;
    // std::cout << *particle  << " in " << *cell << std::endl;
    // std::cout << "--------" << std::endl;
    for (Cell * neighbour : cell->get_neighbours())
    {
        // std::cout << *neighbour << std::endl;
        for (Particle & p : particle->kind == Particle::Kind::Gas ? neighbour->gas_particles : neighbour->dust_particles)
        { // TODO remove direct access
            // std::cout << p  << " in " << *neighbour << std::endl;
            double tmp = kernel(*particle, p, Params::get_instance().dimensions);
            density += p.mass * kernel(*particle, p, Params::get_instance().dimensions);
        }
    }

    return density;
}

void recalc_density(ParticlesState & state)
{
    Grid & grid = state.grid;
    for (int i = 0; i < grid.x_size; ++i)
    {
        for (int j = 0; j < grid.y_size; ++j)
        {
            for (int k = 0; k < grid.z_size; ++k)
            {
                Cell & cell = grid.cells[i][j][k];
                for (Particle * particle : cell.get_all_particles())
                {
                    particle->density = find_density(particle, &(cell));
                }
            }
        }
    }
}

ParticlesState regular_init_state()
{
    ParticlesState init_state;

    Params & params = Params::get_instance();

    double radius = 0.2;
    double abs_velo = radius / 2;

    double center_x = (params.border2.x - params.border1.x) / 2;
    double center_y = (params.border2.y - params.border1.y) / 2;
    double center_z = (params.border2.z - params.border1.z) / 2;

    for(int i = 0; i < params.n_gas; i += 1)
    {
        double theta = acos(1 - 2 * (i + 0.5) / params.n_gas);
        double phi = PI * (1 + pow(5, 0.5)) * (i + 0.5);

        double rel_x = radius * cos(phi) * sin(theta);
        double rel_y = radius * sin(phi) * sin(theta);
        double rel_z = radius * cos(theta);

        Particle particle(Particle::Kind::Gas, rel_x + center_x, rel_y + center_y, rel_z + center_z);

        particle.vx = rel_x / radius * abs_velo;
        particle.vy = rel_y / radius * abs_velo;
        particle.vz = rel_z / radius * abs_velo;

        particle.mass = 1. / params.n_gas;

        init_state.with_copy_of(particle);
    }

    recalc_density(init_state);

    return init_state;
}

Point find_new_coordinates(Particle const & particle)
{
    Params & params = Params::get_instance();
    double x = particle.x + params.tau * particle.vx;
    double y = particle.y + params.tau * particle.vy;
    double z = particle.z + params.tau * particle.vz;
    return {x, y, z};
}

ParticlesState do_time_step(ParticlesState & old, int step_num)
{
    char filename[512];
    sprintf(filename, "/home/calat/tmp/part_%0d.dat", step_num);
    FILE * f = fopen(filename, "w");

    ParticlesState nextState;

    Grid old_grid = old.grid;
    for (int i = 0; i < old_grid.x_size; ++i)
    { // TODO foreach
        for (int j = 0; j < old_grid.y_size; ++j)
        {
            for (int k = 0; k < old_grid.z_size; ++k)
            {
                Cell cell = old_grid.cells[i][j][k];
                for (Particle * particle : cell.get_all_particles())
                {
                    fprintf(f, "%lf %lf %lf %lf\n", particle->x, particle->y, particle->z, particle->density);

                    Particle new_particle(*particle);
                    Point new_coords = find_new_coordinates(new_particle);
                    new_particle.set_coordinates(new_coords.x, new_coords.y, new_coords.z);

                    nextState.with_copy_of(new_particle);
                }
            }
        }
    }

    recalc_density(old); // TODO next
    fclose(f);

    return nextState;
}

int main()
{
    clock_t startTime = clock();
    ParticlesState state = regular_init_state();

    for (int i = 0; i < 10; ++i)
    {
        clock_t iter_start = clock();
        state = do_time_step(state, i);
        clock_t iter_fin = clock();

        double iter_time = (double)(iter_fin - iter_start) / CLOCKS_PER_SEC;
        std::cout << i << " " << iter_time << " s" << std::endl;
    }
    std::cout << "Done!" << std::endl;
    clock_t finishTime = clock();

    double executionTime = (double)(finishTime - startTime) / CLOCKS_PER_SEC;
    printf("Finished in %lf seconds.\n", executionTime);

    return 0;
}
