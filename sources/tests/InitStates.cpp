#include <vector>
#include <iostream>

#include "MathUtils.h"
#include "Params.h"
#include "Grid.h"
#include "Cell.h"
#include "SodTube1d.h"

#include "InitStates.h"


void fill_initial_sod_coord(std::vector<double> & coord, std::vector<double> & image_coord, uint real_left_p_num,
                        uint real_right_p_num, uint image_left_p_num, uint image_right_p_num)
{
    //Params & params = Params::get_instance();
    Dusty_shock_params & shock_params = Dusty_shock_params::get_instance();

    // step = length / number
    double step_left = (shock_params.membrane - shock_params.x_left) / real_left_p_num;
    double step_right = (shock_params.x_right - shock_params.membrane) / real_right_p_num;

    // write left "real" into both arrays
    coord.at(0) = shock_params.x_left + (step_left / 2);
    image_coord.at(image_left_p_num) = coord.at(0);
    for (uint i = 1; i < real_left_p_num; i++) {
        coord.at(i) = coord.at(i - 1) + step_left;
        image_coord.at(image_left_p_num + i) = coord.at(i);
    }
    // calculate separately left "imaginary" into image array
    for (int i = image_left_p_num - 1; i >= 0; i--) {
        image_coord.at(i) = image_coord.at(i + 1) - step_left;
    }

    // write right real into both arrays
    coord.at(real_left_p_num) = shock_params.membrane + (step_right / 2);
    image_coord.at(image_left_p_num + real_left_p_num) = coord.at(real_left_p_num);
    for (uint i = 1; i < real_right_p_num; i++) {
        coord.at(real_left_p_num + i) = coord.at(real_left_p_num + i - 1) + step_right;
        image_coord.at(image_left_p_num + real_left_p_num + i) = coord.at(real_left_p_num + i);
    }
    // calculate separately left "imaginary"
    for (uint i = 0; i < image_right_p_num; i++) {
        int r_image_id = image_left_p_num + real_left_p_num + real_right_p_num + i;
        image_coord.at(r_image_id) = image_coord.at(r_image_id - 1) + step_right;
    }
}

Grid Sod_tube_1d::init()
{
    Params & params = Params::get_instance();
    Dusty_shock_params & shock_params = Dusty_shock_params::get_instance();

    Grid init(params.grid_step_x, params.grid_step_y, params.grid_step_z, params.border1, params.border2);

    uint real_particles = shock_params.gas_real_particles;
    uint all_particles = shock_params.gas_real_particles + shock_params.gas_image_particles;

    double gas_l2r = shock_params.gas_dens_left / shock_params.gas_dens_right;

    // считаем количества реальных частиц
    uint gas_real_right_p_num = (uint) round((double) shock_params.gas_real_particles / (gas_l2r + 1));
    uint gas_real_left_p_num = shock_params.gas_real_particles - gas_real_right_p_num;

    // аналогично для виртуальных
    uint gas_image_right_p_num = (uint) round((double) shock_params.gas_image_particles / (gas_l2r + 1));
    uint gas_image_left_p_num = shock_params.gas_image_particles - gas_image_right_p_num;

    std::vector<double> coord(real_particles);
    std::vector<double> image_coord(all_particles);

    fill_initial_sod_coord(coord, image_coord, gas_real_left_p_num, gas_real_right_p_num, gas_image_left_p_num,
                       gas_image_right_p_num);

    double gas_image_left_lenght = image_coord.at(gas_image_left_p_num + gas_real_left_p_num) - image_coord.at(0);
    double mass = gas_image_left_lenght * shock_params.gas_dens_left / (gas_image_left_p_num + gas_real_left_p_num);

    for(uint i = 0; i < all_particles; ++i)
    {
        Particle particle(Particle::Kind::Gas, image_coord.at(i), 0, 0);

        particle.mass = mass;
        if(particle.x <= shock_params.membrane)
        {
            particle.set_velocities(shock_params.gas_vel_left, 0, 0);
            particle.density = shock_params.gas_dens_left;
            particle.pressure = shock_params.gas_press_left;
            particle.energy = shock_params.gas_ener_left;
        }
        else
        {
            particle.set_velocities(shock_params.gas_vel_right, 0, 0);
            particle.density = shock_params.gas_dens_right;
            particle.pressure = shock_params.gas_press_right;
            particle.energy = shock_params.gas_ener_right;
        }

        init.with_copy_of(particle);
    }

    recalc_density(init, Particle::Kind::Gas);
    Sod_tube_1d::recalc_pressure(init, Particle::Kind::Gas);

    return init;
}

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

Grid random_init_state()
{
    Params & params = Params::get_instance();

    Grid init(params.grid_step_x, params.grid_step_y, params.grid_step_z, params.border1, params.border2);

    for (int i = 0; i < params.n_gas; ++i)
    {
        Particle particle = new_random_particle(Particle::Kind::Gas, init.border1, init.border2);
        init.with_copy_of(particle);
    }
    for (int i = 0; i < params.n_dust; ++i)
    {
        Particle particle = new_random_particle(Particle::Kind::Dust, init.border1, init.border2);
        init.with_copy_of(particle);

    }
    return init;
}

Grid sphere_init_state()
{
    Params & params = Params::get_instance();

    Grid init(params.grid_step_x, params.grid_step_y, params.grid_step_z, params.border1, params.border2);

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

        init.with_copy_of(particle);
    }

    recalc_density(init, Particle::Kind::Gas);

    return init;
}

Grid ball_init_state(double radius)
{
    Params & params = Params::get_instance();

    Grid init(params.grid_step_x, params.grid_step_y, params.grid_step_z, params.border1, params.border2);

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

            init.with_copy_of(particle);
        }
        part_quant *= 3;

    }

    recalc_density(init, Particle::Kind::Gas);

    return init;
}

Grid ball_rand_init_state(double radius)
{
    Params & params = Params::get_instance();

    Grid init(params.grid_step_x, params.grid_step_y, params.grid_step_z, params.border1, params.border2);

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

            init.with_copy_of(particle);

            i += 1;
        }
    }

    recalc_density(init, Particle::Kind::Gas);

    return init;
}

Grid squared_ball_init_state(double radius, double step)
{
    Params & params = Params::get_instance();

    Grid init(params.grid_step_x, params.grid_step_y, params.grid_step_z, params.border1, params.border2);

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

                    init.with_copy_of(particle);
                    ++count;
                }
            }
        }
    }


    for (int i = 0; i < init.x_size; ++i)
    {
        for(int j = 0; j < init.y_size; ++j)
        {
            for(int k = 0; k < init.z_size; ++k)
            {
                Cell * cell = &(init.cells[i][j][k]);
                for (Particle * particle : cell->get_all_particles())
                {
                    particle->mass = 1. / count;

                }
            }
        }
    }
    std::cout << "------particles:" << count << std::endl;

    recalc_density(init, Particle::Kind::Gas);

    return init;
}