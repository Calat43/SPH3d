#include <iostream>
#include <string>
#include "Particle.h"
#include "ParticlesState.h"
#include "Common.h"
#include "Grid.h"
#include "Cell.h"
#include "Params.h"
#include "CompareFiles.h"
#include "Solver.h"


Particle new_random_particle(Particle::Kind kind, Point border1, Point border2)
{
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

ParticlesState sphere_init_state()
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

    recalc_density(init_state, 0.1);

    return init_state;
}

ParticlesState ball_init_state(double radius)
{
    ParticlesState init_state;

    Params & params = Params::get_instance();

    double center_x = (params.border2.x - params.border1.x) / 2;
    double center_y = (params.border2.y - params.border1.y) / 2;
    double center_z = (params.border2.z - params.border1.z) / 2;

    int rad_quantity = 0;
    int part_quant = 1;

    while(part_quant <= params.n_gas)
    {
        part_quant = part_quant + part_quant * 3;
        rad_quantity += 1;
    }

    part_quant = 1;
    double rad_step = radius / (double) rad_quantity;

    for(int i = 0; i <= rad_quantity; ++i)
    {
        double this_rad = rad_step * i;

        for(int k = 0; k < part_quant; ++k)
        {
            double theta = acos(1 - 2 * (k + 0.5) / part_quant);
            double phi = PI * (1 + pow(5, 0.5)) * (k + 0.5);

            double rel_x = this_rad * cos(phi) * sin(theta);
            double rel_y = this_rad * sin(phi) * sin(theta);
            double rel_z = this_rad * cos(theta);

            Particle particle(Particle::Kind::Gas, rel_x + center_x, rel_y + center_y, rel_z + center_z);

            particle.vx = 0;
            particle.vy = 0;
            particle.vz = 0;

            particle.mass = 1. / params.n_gas;

            init_state.with_copy_of(particle);
        }
        part_quant *= 3;

    }

    recalc_density(init_state, radius);

    return init_state;
}

ParticlesState ball_rand_init_state(double radius)
{
    ParticlesState init_state;

    Params & params = Params::get_instance();

    double center_x = (params.border2.x - params.border1.x) / 2;
    double center_y = (params.border2.y - params.border1.y) / 2;
    double center_z = (params.border2.z - params.border1.z) / 2;

    int i = 0;
    while(i < params.n_gas)
    {
        double x = center_x + random_double(-radius, radius);
        double y = center_y + random_double(-radius, radius);
        double z = center_z + random_double(-radius, radius);

        if(sqrt(pow((x - center_x), 2) + pow((y - center_y), 2) + pow((z - center_z), 2)) <= radius)
        {
            Particle particle(Particle::Kind::Gas, x, y, z);

            particle.vx = 0;
            particle.vy = 0;
            particle.vz = 0;

            particle.mass = 1. / params.n_gas;

            init_state.with_copy_of(particle);

            i += 1;
        }
    }

    recalc_density(init_state, radius);

    return init_state;
}

ParticlesState squared_ball_init_state(double radius, double step)
{
    Params & params = Params::get_instance();

    ParticlesState init_state;

    double center_x = (params.border2.x - params.border1.x) / 2;
    double center_y = (params.border2.y - params.border1.y) / 2;
    double center_z = (params.border2.z - params.border1.z) / 2;

    int count = 0;

    for(double x = params.border1.x; x < params.border2.x; x += step)
    {
        for(double y = params.border1.y; y < params.border2.y; y += step)
        {
            for(double z = params.border1.z; z < params.border2.z; z += step)
            {
                if (distance(Point(x, y, z), Point(center_x, center_y, center_z)) <= radius)
                {
                    Particle particle(Particle::Kind::Gas, x, y, z);

                    particle.vx = 0;
                    particle.vy = 0;
                    particle.vz = 0;

                    // particle.mass = 1. / params.n_gas;

                    init_state.with_copy_of(particle);
                    ++count;
                }
            }
        }
    }


    for (int i = 0; i < init_state.grid.x_size; ++i)
    {
        for(int j = 0; j < init_state.grid.y_size; ++j)
        {
            for(int k = 0; k < init_state.grid.z_size; ++k)
            {
                Cell * cell = &(init_state.grid.cells[i][j][k]);
                for (Particle * particle : cell->get_all_particles())
                {
                    particle->mass = 1. / count;

                }
            }
        }
    }
    std::cout << "------particles:" << count << std::endl;

    recalc_density(init_state, radius);

    return init_state;
}

ParticlesState do_time_step(ParticlesState & old, int step_num)
{
    char filename[512];
    sprintf(filename, "/home/calat/tmp/part_%0d.dat", step_num);
    FILE * f = fopen(filename, "w");

    ParticlesState nextState;

    Params & params = Params::get_instance();
    double center_x = (params.border2.x - params.border1.x) / 2;
    double center_y = (params.border2.y - params.border1.y) / 2;
    double center_z = (params.border2.z - params.border1.z) / 2;

    Grid & old_grid = old.grid;
    for (int i = 0; i < old_grid.x_size; ++i)
    { // TODO foreach
        for (int j = 0; j < old_grid.y_size; ++j)
        {
            for (int k = 0; k < old_grid.z_size; ++k)
            {
                Cell cell = old_grid.cells[i][j][k];
                for (Particle * particle : cell.get_all_particles())
                {
                    if(PRINT_DENSITY)
                    {
                        fprintf(f, "%lf %lf %lf %lf\n", particle->x, particle->y, particle->z, particle->density);
                    }
                    else
                    {
                        fprintf(f, "%lf %lf %lf %lf %lf %lf %lf\n", particle->x, particle->y, particle->z, particle->vx,
                                particle->vy, particle->vz, particle->density);
                    }


                    Particle new_particle(*particle);
                    Point new_coords = find_new_coordinates(*particle);
                    //Cell * new_particle_old_cell = old.grid.find_cell(new_particle);
                    //assert(new_particle_old_cell != nullptr);

                    Point new_vel = find_new_velocity(particle, &cell);
                    //Point new_vel = find_new_velocity_no_sort(&new_particle, &nextState.grid);
                    new_particle.density = NAN;
                    //new_particle.density = 42;
                    new_particle.set_coordinates(new_coords.x, new_coords.y, new_coords.z);
                    new_particle.set_velocities(new_vel.x, new_vel.y, new_vel.z);

                    nextState.with_copy_of(new_particle);
                }
            }
        }
    }
    recalc_density(nextState, 1);

    fclose(f);

    return nextState;
}

int main()
{
    clock_t startTime = clock();

    Params & params = Params::get_instance();
    //ParticlesState state = ball_rand_init_state(0.1);
    ParticlesState state = squared_ball_init_state(0.1, 0.01); // FIXME squared_ball does not take n_gas into account !!

    for (int frameId = 0; frameId < floor(params.t / params.tau); ++frameId)
    {
        clock_t step_start = clock();
        state = do_time_step(state, frameId);
        clock_t step_fin = clock();

        double step_time = (double)(step_fin - step_start) / CLOCKS_PER_SEC;
        std::cout << frameId << " " << step_time << " s" << std::endl;

    }

    std::cout << "Done!" << std::endl;
    clock_t finishTime = clock();

    double executionTime = (double)(finishTime - startTime) / CLOCKS_PER_SEC;
    printf("Finished in %lf seconds.\n", executionTime);

    return 0;
}
