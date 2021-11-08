#include <fstream>
#include <sstream>

#include "Params.h"
#include "Grid.h"
#include "Cell.h"
#include "Viscosity.h"
#include "NonLinear.h"

#include "SodTube1d.h"
#include "InitStates.h"
#include "DustyShock1d.h"
#include "SPHSolver.h"
#include "Solver.h"

Grid Dusty_shock_1d::init()
{
    Params & params = Params::get_instance();
    Dusty_shock_params & shock_params = Dusty_shock_params::get_instance();

    Grid init(params.grid_step_x, params.grid_step_y, params.grid_step_z, params.border1, params.border2);

    init = Sod_tube_1d::init();

    uint real_particles = shock_params.dust_real_particles;
    uint all_particles = shock_params.dust_real_particles + shock_params.dust_image_particles;

    double dust_l2r = shock_params.dust_dens_left / shock_params.dust_dens_right;

    // считаем количества реальных частиц
    uint dust_real_right_p_num = (uint) round((double) shock_params.dust_real_particles / (dust_l2r + 1));
    uint dust_real_left_p_num = shock_params.dust_real_particles - dust_real_right_p_num;

    // аналогично для виртуальных
    uint dust_image_right_p_num = (uint) round((double) shock_params.dust_image_particles / (dust_l2r + 1));
    uint dust_image_left_p_num = shock_params.dust_image_particles - dust_image_right_p_num;

    std::vector<double> coord(real_particles);
    std::vector<double> image_coord(all_particles);

    fill_initial_sod_coord(coord, image_coord, dust_real_left_p_num, dust_real_right_p_num, dust_image_left_p_num,
                           dust_image_right_p_num);

    double dust_image_left_lenght = image_coord.at(dust_image_left_p_num + dust_real_left_p_num) - image_coord.at(0);
    double mass = dust_image_left_lenght * shock_params.dust_dens_left / (dust_image_left_p_num + dust_real_left_p_num);

    uint i = 0;

    for(i = 0; i < dust_image_left_p_num; ++i)
    {
        Particle particle(Particle::Kind::Dust, image_coord.at(i), 0, 0);

        particle.mass = mass;

        particle.set_velocities(shock_params.dust_vel_left, 0, 0);
        particle.density = shock_params.dust_dens_left;
        particle.pressure = NAN;
        particle.energy = NAN;

        particle.is_border = true;

        init.with_copy_of(particle);
    }

    for(i = dust_image_left_p_num; i < real_particles + dust_image_left_p_num; ++i)
    {
        Particle particle(Particle::Kind::Dust, image_coord.at(i), 0, 0);

        particle.mass = mass;
        if(particle.x <= shock_params.membrane)
        {
            particle.set_velocities(shock_params.dust_vel_left, 0, 0);
            particle.density = shock_params.dust_dens_left;
            particle.pressure = NAN;
            particle.energy = NAN;
        }
        else
        {
            particle.set_velocities(shock_params.dust_vel_right, 0, 0);
            particle.density = shock_params.dust_dens_right;
            particle.pressure = NAN;
            particle.energy = NAN;
        }

        particle.is_border = false;

        init.with_copy_of(particle);
    }

    for(i = real_particles + dust_image_left_p_num; i < all_particles; ++i)
    {
        Particle particle(Particle::Kind::Dust, image_coord.at(i), 0, 0);

        particle.mass = mass;

        particle.set_velocities(shock_params.dust_vel_right, 0, 0);
        particle.density = shock_params.dust_dens_right;
        particle.pressure = NAN;
        particle.energy = NAN;

        particle.is_border = true;

        init.with_copy_of(particle);
    }
/*
    for(uint i = 0; i < all_particles; ++i)
    {
        Particle particle(Particle::Kind::Dust, image_coord.at(i), 0, 0);

        particle.mass = mass;
        if(particle.x <= shock_params.membrane)
        {
            particle.set_velocities(shock_params.dust_vel_left, 0, 0);
            particle.density = shock_params.dust_dens_left;
            particle.pressure = NAN;
            particle.energy = NAN;
        }
        else
        {
            particle.set_velocities(shock_params.dust_vel_right, 0, 0);
            particle.density = shock_params.dust_dens_right;
            particle.pressure = NAN;
            particle.energy = NAN;
        }

        init.with_copy_of(particle);
    }
*/
    recalc_density(init, Particle::Kind::Dust);
    //Sod_tube_1d::recalc_pressure(init, Particle::Kind::Dust);

    std::ofstream gas_file;
    std::ofstream dust_file;

    if(PRINT_FILE)
    {
        gas_file.open(OUTPUT_PATH + "gas_init" + ".dat");
        dust_file.open(OUTPUT_PATH + "dust_init" + ".dat");
    }

    for(int i = 0; i < init.x_size; ++i)
    { // TODO foreach
        for(int j = 0; j < init.y_size; ++j)
        {
            for(int k = 0; k < init.z_size; ++k)
            {
                Cell & cell = init.cells[i][j][k];

                for(Particle & particle : cell.gas_particles)
                {
                    if(PRINT_FILE)
                    {
                        gas_file.precision(10);
                        gas_file << particle.x << " " << particle.vx << " " << particle.density << " "
                                 << particle.energy << " " << particle.pressure << std::endl;
                        //gas_file << std::fixed << particle.density << std::endl;
                    }
                }

                for(Particle & particle : cell.dust_particles)
                {
                    if(PRINT_FILE)
                    {
                        dust_file.precision(10);
                        dust_file << particle.x << " " << particle.vx << " " << particle.density << std::endl;
                        //dust_file << std::fixed << particle.density  << std::endl;
                    }
                }
            }
        }
    }

    return init;
}

Grid Dusty_shock_1d::do_time_step(Grid & old_grid, int step_num)
{
    std::ofstream f_dust;
    std::ofstream f_gas;

    Params & params = Params::get_instance();
    Dusty_shock_params & shock_params = Dusty_shock_params::get_instance();
    int final_step = params.t / params.tau - 1;

    if (PRINT_FILE && (step_num == final_step || step_num == 0))
    {
        //char filename_gas[512];
        //sprintf(filename_gas, (OUTPUT_PATH + "dust_h%lg_cfl%lg_n%d_hc%lg_part_%0d.dat").c_str(), params.h, params.tau / params.h, shock_params.dust_real_particles, params.grid_step_x, step_num);
        //f_dust.open(filename_dust);

        std::stringstream filename_dust;
        filename_dust << OUTPUT_PATH
                     << "dust"
                     << "_h" << params.h
                     << "_tau" << params.tau
                     << "_part_" << step_num
                     << ".dat";
        f_dust.open(filename_dust.str());

        std::stringstream filename_gas;
        filename_gas << OUTPUT_PATH
                     << "gas"
                     << "_h" << params.h
                     << "_tau" << params.tau
                     << "_part_" << step_num
                     << ".dat";
        f_gas.open(filename_gas.str());
    }

    Grid next_grid(params.grid_step_x, params.grid_step_y, params.grid_step_z, params.border1, params.border2);

    for (int i = 0; i < old_grid.x_size; ++i)
    { // TODO foreach
        for(int j = 0; j < old_grid.y_size; ++j)
        {
            for(int k = 0; k < old_grid.z_size; ++k)
            {
                Cell & cell = old_grid.cells[i][j][k];
                for(Particle & particle : cell.gas_particles)
                {
                    if (PRINT_FILE && (step_num == final_step || step_num == 0))
                    {
                        f_gas << particle.x << " " << particle.vx << " " << particle.density << " "
                              << particle.pressure << " " << particle.energy << std::endl;
                    }

                    if(particle.is_border)
                    {
                        next_grid.with_copy_of(particle);
                    }
                    else {
                        Particle new_particle(particle);
                        Point new_coords = find_new_coordinates(particle);

                        double new_vel = gas_velocity_shock_MK(&particle, &cell);
                        new_particle.density = NAN;

                        new_particle.set_coordinates(new_coords.x, new_coords.y, new_coords.z);
                        new_particle.set_velocities(new_vel, particle.vy, particle.vz);

                        new_particle.energy = Sod_tube_1d::find_new_energy(&particle, &cell);

                        next_grid.with_copy_of(new_particle);
                    }
                }
                for(Particle & particle : cell.dust_particles)
                {
                    if (PRINT_FILE && (step_num == final_step || step_num == 0)) {
                        f_dust << particle.x << " " << particle.vx << " " << particle.density << std::endl;
                    }

                    if (particle.is_border) {
                        next_grid.with_copy_of(particle);
                    }
                    else {
                        Particle new_particle(particle);
                        Point new_coords = find_new_coordinates(particle);

                        double new_vel = dust_velocity_shock_MK(&particle, &cell);
                        new_particle.density = NAN;

                        new_particle.set_coordinates(new_coords.x, new_coords.y, new_coords.z);
                        new_particle.set_velocities(new_vel, particle.vy, particle.vz);

                        next_grid.with_copy_of(new_particle);
                    }
                }
            }
        }
    }

    recalc_density(next_grid, Particle::Kind::Gas);
    recalc_density(next_grid, Particle::Kind::Dust);
    Sod_tube_1d::recalc_pressure(next_grid, Particle::Kind::Gas);

    return next_grid;
}