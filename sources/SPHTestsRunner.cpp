#include <iostream>
#include <string>
#include <fstream>
#include <time.h>

#include "DustyWave.h"
#include "Particle.h"
#include "MathUtils.h"
#include "Grid.h"
#include "Cell.h"
#include "Params.h"
#include "CompareFiles.h"
#include "SodTube1d.h"
#include "InitStates.h"
#include "SodTube3d.h"
#include "DustyShock.h"
#include "Point.h"
#include "DustyShock1d.h"
#include "BallInVacuum.h"
#include "SPHSolver.h"


//TODO empty grid constructor
Grid do_time_step(Grid & old_grid, int step_num)
{
    char filename[512];
    sprintf(filename, (OUTPUT_PATH + "part_%0d.dat").c_str(), step_num);
    FILE * f = fopen(filename, "w");

    Params & params = Params::get_instance();

    Grid next_grid(params.grid_step_x, params.grid_step_y, params.grid_step_z, params.border1, params.border2);

    // unreferenced variables
    //double center_x = (params.border2.x - params.border1.x) / 2;
    //double center_y = (params.border2.y - params.border1.y) / 2;
    //double center_z = (params.border2.z - params.border1.z) / 2;

    for (int i = 0; i < old_grid.x_size; ++i)
    { // TODO foreach
        for (int j = 0; j < old_grid.y_size; ++j)
        {
            for (int k = 0; k < old_grid.z_size; ++k)
            {
                Cell & cell = old_grid.cells[i][j][k];
                for (Particle * particle : cell.get_all_particles())
                {
                    if(PRINT_FILE)
                    {
                        fprintf(f, "%lf %lf %lf %lf %lf %lf %lf\n", particle->x, particle->y, particle->z, particle->vx,
                                particle->vy, particle->vz, particle->density);
                    }

                    Particle new_particle(*particle);
                    Point new_coords = find_new_coordinates_(*particle);
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


void check_output_dir()
{
    std::ofstream test_file(OUTPUT_PATH + "test");
    if (!test_file.is_open())
    {
        std::cerr << "Cant write to the output directory" << std::endl;
        std::exit(1);
    }
    else
    {
        test_file.close();
    }
}

Particle scan_dust_particle(std::ifstream & fin)
{
    double x, y, z, vx, vy, vz, density;
    fin >> x >> y >> z >> vx >> vy >> vz >> density;
    Particle part(Particle::Kind::Dust, x, y, z);
    part.vx = vx;
    part.vy = vy;
    part.vz = vz;
    part.density = density;
    return part;
}

Particle scan_gas_particle(std::ifstream & fin)
{
    double x, y, z, vx, vy, vz, density, pressure, energy;
    fin >> x >> y >> z >> vx >> vy >> vz >> density >> pressure >> energy;
    Particle part(Particle::Kind::Gas, x, y, z);
    part.vx = vx;
    part.vy = vy;
    part.vz = vz;
    part.density = density;
    part.pressure = pressure;
    part.energy = energy;
    return part;
}

Grid restore_saved_grid()
{
    Params & params = Params::get_instance();
    Grid grid(params.grid_step_x, params.grid_step_y, params.grid_step_z, params.border1, params.border2);

    std::ifstream dust_parts("/home/haitaka/tania/SPH3d/output/dust_part_34.dat");
    std::ifstream gas_parts("/home/haitaka/tania/SPH3d/output/gas_part_34.dat");

    while (!dust_parts.eof())
    {
        grid.with_copy_of(scan_dust_particle(dust_parts));
    }

    while (!gas_parts.eof())
    {
        grid.with_copy_of(scan_gas_particle(gas_parts));
    }

    return grid;
}


int main()
{
    check_output_dir();

    Params & params = Params::get_instance();

/*
    std::cout << ball_analytic::time_through_R(0.9) << std::endl;
    std::cout << ball_analytic::time_through_R(0.87714) << std::endl;

    std::cout << ball_analytic::R_bisection(0.2, 0.00001, 0.000001) << std::endl;
    std::cout << ball_analytic::time_through_R(ball_analytic::R_bisection(0.2, 0.00001, 0.000001)) << std::endl;

 */
    //ball_analytic::print_solution(0.2, 0.01, 0.0001, 0.0001);


    //Grid grid = ball_in_vacuum::init(0.05);
    Grid grid = DustyWave::init();
    //Grid grid = Dusty_shock_1d::init();



    std::cerr << "Init done" << std::endl;
    clock_t startTime = clock();
    for (int frameId = 0; frameId < floor(params.t / params.tau); ++frameId)
    {
        clock_t step_start = clock();

        //grid = ball_in_vacuum::do_time_step(grid, frameId);
        grid = DustyWave::do_time_step(grid, frameId);
        //grid = Dusty_shock_1d::do_time_step(grid, frameId);
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
