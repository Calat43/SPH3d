#include <iostream>
#include <string>
#include "Particle.h"
#include "Common.h"
#include "Grid.h"
#include "Cell.h"
#include "Params.h"
#include "CompareFiles.h"
#include "Solver.h"
#include "InitStates.h"
#include "SodTube3d.h"
#include "Dusty_shock.h"
#include "Point.h"


//TODO empty grid constructor
Grid do_time_step(Grid & old_grid, int step_num)
{
    char filename[512];
    sprintf(filename, "/home/calat/tmp/part_%0d.dat", step_num);
    FILE * f = fopen(filename, "w");

    Params & params = Params::get_instance();

    Grid next_grid(params.grid_step_x, params.grid_step_y, params.grid_step_z, params.border1, params.border2);

    double center_x = (params.border2.x - params.border1.x) / 2;
    double center_y = (params.border2.y - params.border1.y) / 2;
    double center_z = (params.border2.z - params.border1.z) / 2;

    for (int i = 0; i < old_grid.x_size; ++i)
    { // TODO foreach
        for (int j = 0; j < old_grid.y_size; ++j)
        {
            for (int k = 0; k < old_grid.z_size; ++k)
            {
                Cell & cell = old_grid.cells[i][j][k];
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

                    next_grid.with_copy_of(new_particle);
                }
            }
        }
    }
    recalc_density(next_grid, Particle::Kind::Gas);

    fclose(f);

    return next_grid;
}

int main()
{
    clock_t startTime = clock();

    Params & params = Params::get_instance();
    //ParticlesState state = ball_rand_init_state(0.1);
    //Grid grid = squared_ball_init_state(0.1, 0.01); // FIXME squared_ball does not take n_gas into account !!
    clock_t init_start = clock();
   // Grid grid = Sod_tube_3d::init_with_boundaries();
    Grid grid = Dusty_shock_3d::init();
    clock_t init_fin = clock();
    std::cout << "Init: " << (double)(init_fin - init_start) / CLOCKS_PER_SEC << std::endl;

    for (int frameId = 0; frameId < floor(params.t / params.tau); ++frameId)
    {
        clock_t step_start = clock();
        //TODO PRINT
        //grid = Sod_tube_3d::do_time_step_with_boundaries(grid, frameId);
        grid = Dusty_shock_3d::do_time_step(grid, frameId, IS_IDIC);
        clock_t step_fin = clock();

        double step_time = (double)(step_fin - step_start) / CLOCKS_PER_SEC;
        std::cout << frameId << " " << step_time << " s" << std::endl;
    }

    std::cout << "Done!" << std::endl;
    clock_t finishTime = clock();

    double executionTime = (double)(finishTime - startTime) / CLOCKS_PER_SEC;
    printf("Finished in %lf seconds.\n", executionTime);

    //Grid grid = Dusty_shock_3d::init();

    //centering();

    return 0;
}
