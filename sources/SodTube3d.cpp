
#include <fstream>
#include <iostream>

#include "Params.h"
#include "Cell.h"
#include "Grid.h"
#include "Viscosity.h"

#include "InitStates.h"
#include "CompareFiles.h"
#include "SodTube1d.h"
#include "SodTube3d.h"


Grid Sod_tube_3d::init()
{
    Params & params = Params::get_instance();
    Dusty_shock_params & shock_params = Dusty_shock_params::get_instance();

    Grid init(params.grid_step_x, params.grid_step_y, params.grid_step_z, params.border1, params.border2);

    uint real_1d_particles = shock_params.gas_real_particles;
    uint all_1d_particles = shock_params.gas_real_particles + shock_params.gas_image_particles;

    double gas_l2r = shock_params.gas_dens_left / shock_params.gas_dens_right;

    // считаем количества реальных частиц
    uint real_right_1d_p_num = (uint) round((double) shock_params.gas_real_particles / (gas_l2r + 1));
    uint real_left_1d_p_num = shock_params.gas_real_particles - real_right_1d_p_num;

    // аналогично для виртуальных
    uint image_right_1d_p_num = (uint) round((double) shock_params.gas_image_particles / (gas_l2r + 1));
    uint image_left_1d_p_num = shock_params.gas_image_particles - image_right_1d_p_num;

    std::vector<double> coord(real_1d_particles);
    std::vector<double> image_coord(all_1d_particles);

    fill_initial_sod_coord(coord, image_coord, real_left_1d_p_num, real_right_1d_p_num, image_left_1d_p_num,
                           image_right_1d_p_num);

    //По осям z и у ставим по yz_particles частиц
    int yz_particles = shock_params.gas_yz_particles;
    double yz_length = shock_params.yz_right - shock_params.yz_left;
    double step_yz = yz_length / (double)yz_particles;

    double image_left_lenght = image_coord.at(image_left_1d_p_num + real_left_1d_p_num) - image_coord.at(0);
    double mass = image_left_lenght * yz_length * yz_length * shock_params.gas_dens_left /
                  (image_left_1d_p_num + real_left_1d_p_num) / (double)yz_particles / (double)yz_particles;

    for(int i = 0; i < (int) all_1d_particles; ++i)
    {
        for(int j = 0; j < yz_particles; ++j)
        {
            for(int k = 0; k < yz_particles; ++k)
            {
                Particle particle_main(Particle::Kind::Gas, image_coord.at(i), step_yz / 2. + j * step_yz,
                                       step_yz / 2. + k * step_yz);

                particle_main.mass = mass;
                particle_main.set_velocities(0, 0, 0);
                //делим объем плоскостью x = membrane
                if(particle_main.x <= shock_params.membrane)
                {
                    particle_main.density = shock_params.gas_dens_left;
                    particle_main.pressure = shock_params.gas_press_left;
                    particle_main.energy = shock_params.gas_ener_left;
                }
                else
                {
                    particle_main.density = shock_params.gas_dens_right;
                    particle_main.pressure = shock_params.gas_press_right;
                    particle_main.energy = shock_params.gas_ener_right;
                }
                init.with_copy_of(particle_main);
            }
        }
    }
    recalc_density(init, Particle::Kind::Gas);
    Sod_tube_1d::recalc_pressure(init, Particle::Kind::Gas);

    if(PRINT_STUFF)
    {
        char init_dens[] = "init_dens";
        print_grid_dens(init, init_dens);
    }

    return init;
}

Point Sod_tube_3d::find_new_velocity(Particle * particle, Cell * cell)
{
    Point sum = {0, 0, 0};
    Point vel = {0, 0, 0};
    Point prev_vel = {particle->vx, particle->vy, particle->vz};
    Point kernel_grad = {0, 0, 0};

    Params & params = Params::get_instance();

    assert(particle->density != 0);
    if(PRINT_STUFF)
    {
        if(__isnan(particle->density))
        {
            std::cout << "NAN particle: " << particle->dbg_state << " " << cell->get_i() << " " << cell->get_j() << " "
                      << cell->get_k() << std::endl;
        }
    }
    assert(!__isnan(particle->density));

    for (Cell * neighbour : cell->get_neighbours())
    {
        for (Particle & p : particle->kind == Particle::Kind::Gas ? neighbour->gas_particles : neighbour->dust_particles)
        { // TODO remove direct access

            kernel_grad = kernel_gradient(*(particle), p, params.dimensions);

            assert(!__isnan(p.density));
            assert(!__isnan(kernel_grad.x));
            double viscosity = viscosity_3d::find_viscosity(particle, &p);

            double term1 = (particle->pressure / pow(particle->density, 2) + p.pressure / pow(p.density, 2) + viscosity);

            sum = sum + kernel_grad * term1;
        }
    }

    double term2 = params.tau * particle->mass;

    vel = prev_vel - sum * term2;

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

Grid Sod_tube_3d::do_time_step(Grid & old_grid, int step_num)
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
                        part_file << particle->x << " " << particle->y << " " << particle->z << " "
                                  << particle->vx << " " << particle->vy << " " << particle->vz << " "
                                  << particle->density << " " << particle->pressure << " " << particle->energy
                                  << std::endl;
                    }

                    Particle new_particle(*particle);
                    Point new_coords = find_new_coordinates(*particle);

                    Point new_vel = Sod_tube_3d::find_new_velocity(particle, &cell);
                    new_particle.density = NAN;
                    new_particle.dbg_state = 2;

                    new_particle.set_coordinates(new_coords.x, new_coords.y, new_coords.z);
                    new_particle.set_velocities(new_vel.x, new_vel.y, new_vel.z);

                    new_particle.energy = Sod_tube_3d::find_new_energy(particle, &cell);

                    next_grid.with_copy_of(new_particle);
                }
            }
        }
    }

    recalc_density(next_grid, Particle::Kind::Gas);
    Sod_tube_1d::recalc_pressure(next_grid, Particle::Kind::Gas);

    return next_grid;
}

Grid Sod_tube_3d::init_with_boundaries()
{
    Params & params = Params::get_instance();
    Dusty_shock_params & shock_params = Dusty_shock_params::get_instance();

    Grid init(params.grid_step_x, params.grid_step_y, params.grid_step_z, params.border1, params.border2);

    uint real_1d_particles = shock_params.gas_real_particles;
    uint all_1d_particles = shock_params.gas_real_particles + shock_params.gas_image_particles;

    double gas_l2r = shock_params.gas_dens_left / shock_params.gas_dens_right;

    // считаем количества реальных частиц
    uint real_right_1d_p_num = (uint) round((double) shock_params.gas_real_particles / (gas_l2r + 1));
    uint real_left_1d_p_num = shock_params.gas_real_particles - real_right_1d_p_num;

    // аналогично для виртуальных
    uint image_right_1d_p_num = (uint) round((double) shock_params.gas_image_particles / (gas_l2r + 1));
    uint image_left_1d_p_num = shock_params.gas_image_particles - image_right_1d_p_num;

    std::vector<double> coord(real_1d_particles);
    std::vector<double> image_coord(all_1d_particles);

    fill_initial_sod_coord(coord, image_coord, real_left_1d_p_num, real_right_1d_p_num, image_left_1d_p_num,
                           image_right_1d_p_num);

    //По осям z и у ставим по yz_particles частиц на длину yz_length
    double yz_length = shock_params.yz_right - shock_params.yz_left;
    double step_yz = yz_length / (double)shock_params.gas_yz_particles;

    double image_left_lenght = image_coord.at(image_left_1d_p_num + real_left_1d_p_num) - image_coord.at(0);
    double mass = image_left_lenght * yz_length * yz_length * shock_params.gas_dens_left /
                  (image_left_1d_p_num + real_left_1d_p_num) / (double)shock_params.gas_yz_particles
                                                             / (double)shock_params.gas_yz_particles;

    for(int i = 0; i < (int) all_1d_particles; ++i)
    {
        for(int j = 0; j < (int) shock_params.gas_yz_particles; ++j)
        {
            for(int k = 0; k < (int) shock_params.gas_yz_particles; ++k)
            {
                Particle particle_main(Particle::Kind::Gas, image_coord.at(i), step_yz / 2. + j * step_yz,
                                       step_yz / 2. + k * step_yz);

                particle_main.mass = mass;
                particle_main.set_velocities(0, 0, 0);
                //делим объем плоскостью x = membrane
                if(particle_main.x <= shock_params.membrane)
                {
                    particle_main.density = shock_params.gas_dens_left;
                    particle_main.pressure = shock_params.gas_press_left;
                    particle_main.energy = shock_params.gas_ener_left;
                }
                else
                {
                    particle_main.density = shock_params.gas_dens_right;
                    particle_main.pressure = shock_params.gas_press_right;
                    particle_main.energy = shock_params.gas_ener_right;
                }
                init.with_copy_of(particle_main);

                Particle particle_zp(particle_main);
                particle_zp.set_coordinates(particle_main.x, particle_main.y, particle_main.z + yz_length);
                init.with_copy_of(particle_zp);

                Particle particle_zm(particle_main);
                particle_zm.set_coordinates(particle_main.x, particle_main.y, particle_main.z - yz_length);
                init.with_copy_of(particle_zm);

                Particle particle_yp(particle_main);
                particle_yp.set_coordinates(particle_main.x, particle_main.y + yz_length, particle_main.z);
                init.with_copy_of(particle_yp);

                Particle particle_ym(particle_main);
                particle_ym.set_coordinates(particle_main.x, particle_main.y - yz_length, particle_main.z);
                init.with_copy_of(particle_ym);

                Particle particle_yp_zp(particle_main);
                particle_yp_zp.set_coordinates(particle_main.x, particle_main.y + yz_length, particle_main.z + yz_length);
                init.with_copy_of(particle_yp_zp);

                Particle particle_ym_zp(particle_main);
                particle_ym_zp.set_coordinates(particle_main.x, particle_main.y - yz_length, particle_main.z + yz_length);
                init.with_copy_of(particle_ym_zp);

                Particle particle_yp_zm(particle_main);
                particle_yp_zm.set_coordinates(particle_main.x, particle_main.y + yz_length, particle_main.z - yz_length);
                init.with_copy_of(particle_yp_zm);

                Particle particle_ym_zm(particle_main);
                particle_ym_zm.set_coordinates(particle_main.x, particle_main.y - yz_length, particle_main.z - yz_length);
                init.with_copy_of(particle_ym_zm);
            }
        }
    }

    int x_left_border = (int) floor((shock_params.x_left - params.border1.x) / params.grid_step_x);
    int x_right_border = (int) ceil((shock_params.x_right - shock_params.x_left) / params.grid_step_x + x_left_border);

    int yz_left_border = (int) floor((shock_params.yz_left - params.border1.y) / params.grid_step_y);
    int yz_right_border = (int) ceil((shock_params.yz_right - shock_params.yz_left) / params.grid_step_y + yz_left_border);

    //частицы в расчетной области
    for(int i = x_left_border; i < x_right_border; ++i)
    {
        for(int j = yz_left_border; j < yz_right_border; ++j)
        {
            for(int k = yz_left_border; k < yz_right_border; ++k)
            {
                Cell & cell = init.cells[i][j][k];

                for(Particle & particle : cell.particles_of_kind(Particle::Kind::Gas))
                {
                    particle.density = find_density(&particle, &cell);
                    particle.pressure = Sod_tube_1d::find_new_pressure(&particle);
                }
            }
        }
    }


    for(int i = x_left_border; i < x_right_border; ++i)
    {
        for(int j = 0; j < init.y_size; ++j)
        {
            for(int k = yz_right_border; k < init.z_size; ++k)
            {
                Cell & cell = init.cells[i][j][k];

                cell.gas_particles.erase(cell.gas_particles.begin(), cell.gas_particles.end());
            }
        }
    }
    for(int i = x_left_border; i < x_right_border; ++i)
    {
        for(int j = 0; j < init.y_size; ++j)
        {
            for(int k = 0; k < yz_left_border; ++k)
            {
                Cell & cell = init.cells[i][j][k];

                cell.gas_particles.erase(cell.gas_particles.begin(), cell.gas_particles.end());
            }
        }
    }
    for(int i = x_left_border; i < x_right_border; ++i)
    {
        for(int j = 0; j < yz_left_border; ++j)
        {
            for(int k = yz_left_border; k < yz_right_border; ++k)
            {
                Cell & cell = init.cells[i][j][k];

                cell.gas_particles.erase(cell.gas_particles.begin(), cell.gas_particles.end());
            }
        }
    }
    for(int i = x_left_border; i < x_right_border; ++i)
    {
        for(int j = yz_right_border; j < init.y_size; ++j)
        {
            for(int k = yz_left_border; k < yz_right_border; ++k)
            {
                Cell & cell = init.cells[i][j][k];

                cell.gas_particles.erase(cell.gas_particles.begin(), cell.gas_particles.end());
            }
        }
    }

    for(int i = x_left_border; i < x_right_border; ++i)
    {
        for(int j = yz_left_border; j < yz_right_border; ++j)
        {
            for(int k = yz_left_border; k < yz_right_border; ++k)
            {
                Cell & cell = init.cells[i][j][k];

                for(Particle * particle : cell.get_all_particles())
                {
                    Particle ym_particle(*particle);
                    ym_particle.set_coordinates(particle->x, particle->y - yz_length, particle->z);
                    init.with_copy_of(ym_particle);

                    Particle yp_particle(*particle);
                    yp_particle.set_coordinates(particle->x, particle->y + yz_length, particle->z);
                    init.with_copy_of(yp_particle);

                    Particle zm_particle(*particle);
                    zm_particle.set_coordinates(particle->x, particle->y, particle->z - yz_length);
                    init.with_copy_of(zm_particle);

                    Particle zp_particle(*particle);
                    zp_particle.set_coordinates(particle->x, particle->y, particle->z + yz_length);
                    init.with_copy_of(zp_particle);

                    Particle particle_yp_zp(*particle);
                    particle_yp_zp.set_coordinates(particle->x, particle->y + yz_length, particle->z + yz_length);
                    init.with_copy_of(particle_yp_zp);

                    Particle particle_ym_zp(*particle);
                    particle_ym_zp.set_coordinates(particle->x, particle->y - yz_length, particle->z + yz_length);
                    init.with_copy_of(particle_ym_zp);

                    Particle particle_yp_zm(*particle);
                    particle_yp_zm.set_coordinates(particle->x, particle->y + yz_length, particle->z - yz_length);
                    init.with_copy_of(particle_yp_zm);

                    Particle particle_ym_zm(*particle);
                    particle_ym_zm.set_coordinates(particle->x, particle->y - yz_length, particle->z - yz_length);
                    init.with_copy_of(particle_ym_zm);
                }
            }
        }
    }

    if(PRINT_STUFF)
    {
        char init_dens[] = "init_dens";
        print_grid_dens(init, init_dens);
    }

    return init;
}

Grid Sod_tube_3d::do_time_step_with_boundaries(Grid & old_grid, int step_num)
{
    std::ofstream part_file;
    if (PRINT_FILE)
    {
        part_file.open(OUTPUT_PATH + "part_" + std::to_string(step_num) + ".dat");
    }

    Params & params = Params::get_instance();
    Dusty_shock_params & shock_params = Dusty_shock_params::get_instance();

    Grid next_grid(params.grid_step_x, params.grid_step_y, params.grid_step_z, params.border1, params.border2);

    int x_left_border = (int) floor((shock_params.x_left - params.border1.x) / params.grid_step_x);
    int x_right_border = (int) ceil((shock_params.x_right - shock_params.x_left) / params.grid_step_x + x_left_border);

    int yz_left_border = (int) floor((shock_params.yz_left - params.border1.y) / params.grid_step_y);
    int yz_right_border = (int) ceil((shock_params.yz_right - shock_params.yz_left) / params.grid_step_y + yz_left_border);

    double yz_length = shock_params.yz_right - shock_params.yz_left;

    for(int i = x_left_border; i < x_right_border; ++i)
    {
        for(int j = yz_left_border; j < yz_right_border; ++j)
        {
            for(int k = yz_left_border; k < yz_right_border; ++k)
            {
                Cell & cell = old_grid.cells[i][j][k];

                for(Particle * particle : cell.get_all_particles())
                {
                    if(step_num == 199 || step_num == 99)
                    {
                        if (PRINT_FILE)
                        {
                            part_file << particle->x << " " << particle->y << " " << particle->z << " "
                                      << particle->vx << " " << particle->vy << " " << particle->vz << " "
                                      << particle->density << " " << particle->pressure << " " << particle->energy
                                      << std::endl;
                        }
                    }

                    Particle new_particle(*particle);
                    Point new_coords = find_new_coordinates(*particle);

                    Point new_vel = Sod_tube_3d::find_new_velocity(particle, &cell);
                    new_particle.density = NAN;
                    new_particle.dbg_state = 2;

                    new_particle.set_coordinates(new_coords.x, new_coords.y, new_coords.z);
                    new_particle.set_velocities(new_vel.x, new_vel.y, new_vel.z);

                    new_particle.energy = Sod_tube_3d::find_new_energy(particle, &cell);

                    next_grid.with_copy_of(new_particle);

                    Particle ym_particle(*particle);
                    ym_particle.set_coordinates(particle->x, particle->y - yz_length, particle->z);
                    next_grid.with_copy_of(ym_particle);

                    Particle yp_particle(*particle);
                    yp_particle.set_coordinates(particle->x, particle->y + yz_length, particle->z);
                    next_grid.with_copy_of(yp_particle);

                    Particle zm_particle(*particle);
                    zm_particle.set_coordinates(particle->x, particle->y, particle->z - yz_length);
                    next_grid.with_copy_of(zm_particle);

                    Particle zp_particle(*particle);
                    zp_particle.set_coordinates(particle->x, particle->y, particle->z + yz_length);
                    next_grid.with_copy_of(zp_particle);

                    Particle particle_yp_zp(*particle);
                    particle_yp_zp.set_coordinates(particle->x, particle->y + yz_length, particle->z + yz_length);
                    next_grid.with_copy_of(particle_yp_zp);

                    Particle particle_ym_zp(*particle);
                    particle_ym_zp.set_coordinates(particle->x, particle->y - yz_length, particle->z + yz_length);
                    next_grid.with_copy_of(particle_ym_zp);

                    Particle particle_yp_zm(*particle);
                    particle_yp_zm.set_coordinates(particle->x, particle->y + yz_length, particle->z - yz_length);
                    next_grid.with_copy_of(particle_yp_zm);

                    Particle particle_ym_zm(*particle);
                    particle_ym_zm.set_coordinates(particle->x, particle->y - yz_length, particle->z - yz_length);
                    next_grid.with_copy_of(particle_ym_zm);
                }
            }
        }
    }

    //Частицы слева и справа (по х) от области не трогаем
    for(int i = 0; i < x_left_border; ++i)
    {
        for(int j = 0; j < old_grid.y_size; ++j)
        {
            for(int k = 0; k < old_grid.z_size; ++k)
            {
                Cell & cell = old_grid.cells[i][j][k];

                for(Particle * particle : cell.get_all_particles())
                {
                    next_grid.with_copy_of(*particle);
                }
            }
        }
    }
    for(int i = x_right_border; i < old_grid.x_size; ++i)
    {
        for(int j = 0; j < old_grid.y_size; ++j)
        {
            for(int k = 0; k < old_grid.z_size; ++k)
            {
                Cell & cell = old_grid.cells[i][j][k];

                for(Particle * particle : cell.get_all_particles())
                {
                    next_grid.with_copy_of(*particle);
                }
            }
        }
    }
    //частицы в расчетной области
    for(int i = x_left_border; i < x_right_border; ++i)
    {
        for(int j = yz_left_border; j < yz_right_border; ++j)
        {
            for(int k = yz_left_border; k < yz_right_border; ++k)
            {
                Cell & cell = next_grid.cells[i][j][k];

                for(Particle * particle : cell.get_all_particles())
                {
                    particle->density = find_density(particle, &cell);
                    particle->pressure = Sod_tube_1d::find_new_pressure(particle);
                    particle->dbg_state = 1;
                }
            }
        }
    }

    //PRINT DENSITIES
    if(PRINT_STUFF)
    {
        char before_copy[] = "before_copy";
        print_grid_dens(next_grid, before_copy);
    }

    for(int i = x_left_border; i < x_right_border; ++i)
    {
        for(int j = 0; j < next_grid.y_size; ++j)
        {
            for(int k = yz_right_border; k < next_grid.z_size; ++k)
            {
                Cell & cell = next_grid.cells[i][j][k];

                cell.gas_particles.erase(cell.gas_particles.begin(), cell.gas_particles.end());
            }
        }
    }
    for(int i = x_left_border; i < x_right_border; ++i)
    {
        for(int j = 0; j < next_grid.y_size; ++j)
        {
            for(int k = 0; k < yz_left_border; ++k)
            {
                Cell & cell = next_grid.cells[i][j][k];

                cell.gas_particles.erase(cell.gas_particles.begin(), cell.gas_particles.end());
            }
        }
    }
    for(int i = x_left_border; i < x_right_border; ++i)
    {
        for(int j = 0; j < yz_left_border; ++j)
        {
            for(int k = yz_left_border; k < yz_right_border; ++k)
            {
                Cell & cell = next_grid.cells[i][j][k];

                cell.gas_particles.erase(cell.gas_particles.begin(), cell.gas_particles.end());
            }
        }
    }
    for(int i = x_left_border; i < x_right_border; ++i)
    {
        for(int j = yz_right_border; j < next_grid.y_size; ++j)
        {
            for(int k = yz_left_border; k < yz_right_border; ++k)
            {
                Cell & cell = next_grid.cells[i][j][k];

                cell.gas_particles.erase(cell.gas_particles.begin(), cell.gas_particles.end());
            }
        }
    }

    for(int i = x_left_border; i < x_right_border; ++i)
    {
        for(int j = yz_left_border; j < yz_right_border; ++j)
        {
            for(int k = yz_left_border; k < yz_right_border; ++k)
            {
                Cell & cell = next_grid.cells[i][j][k];

                for(Particle * particle : cell.get_all_particles())
                {
                    Particle ym_particle(*particle);
                    ym_particle.set_coordinates(particle->x, particle->y - yz_length, particle->z);
                    next_grid.with_copy_of(ym_particle);

                    Particle yp_particle(*particle);
                    yp_particle.set_coordinates(particle->x, particle->y + yz_length, particle->z);
                    next_grid.with_copy_of(yp_particle);

                    Particle zm_particle(*particle);
                    zm_particle.set_coordinates(particle->x, particle->y, particle->z - yz_length);
                    next_grid.with_copy_of(zm_particle);

                    Particle zp_particle(*particle);
                    zp_particle.set_coordinates(particle->x, particle->y, particle->z + yz_length);
                    next_grid.with_copy_of(zp_particle);

                    Particle particle_yp_zp(*particle);
                    particle_yp_zp.set_coordinates(particle->x, particle->y + yz_length, particle->z + yz_length);
                    next_grid.with_copy_of(particle_yp_zp);

                    Particle particle_ym_zp(*particle);
                    particle_ym_zp.set_coordinates(particle->x, particle->y - yz_length, particle->z + yz_length);
                    next_grid.with_copy_of(particle_ym_zp);

                    Particle particle_yp_zm(*particle);
                    particle_yp_zm.set_coordinates(particle->x, particle->y + yz_length, particle->z - yz_length);
                    next_grid.with_copy_of(particle_yp_zm);

                    Particle particle_ym_zm(*particle);
                    particle_ym_zm.set_coordinates(particle->x, particle->y - yz_length, particle->z - yz_length);
                    next_grid.with_copy_of(particle_ym_zm);
                }
            }
        }
    }


    //PRINT DENSITIES
    if(PRINT_STUFF)
    {
        char after_copy[] = "after_copy";
        print_grid_dens(next_grid, after_copy);
    }

    return next_grid;
}