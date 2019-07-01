#include <fstream>
#include <sstream>
#include "DustyShock1d.h"

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

    recalc_density(init, Particle::Kind::Dust);

    return init;
}

double Dusty_shock_1d::pressure_term(Particle * particle, Cell * cell)
{
    Params & params = Params::get_instance();

    double a_p_x = 0;

    for (Cell * neighbour : cell->get_neighbours())
    {
        for (Particle & p : neighbour->gas_particles)
        { // TODO remove direct access
            assert(!__isnan(p.density));
            double viscosity = viscosity_1d::find_viscosity(particle, &p);

            a_p_x += (particle->pressure / pow(particle->density, 2) + p.pressure / pow(p.density, 2) + viscosity)
                     * kernel_gradient(*(particle), p, params.dimensions).x;

        }
    }
    a_p_x *= - particle->mass;

    return a_p_x;
}

double idic_1d::pressure_term_asterisk(Cell * cell)
{
    double gas_size = cell->gas_particles.size();

    double a_p_x_astr = 0;

    if(gas_size == 0)
    {
        return 0;
    }
    else
    {
        for (Particle & p : cell->gas_particles)
        {
            a_p_x_astr += Dusty_shock_1d::pressure_term(&p, cell);
        }
    }

    return a_p_x_astr / gas_size;
}

double idic_1d::vel_asterisk(Cell * cell, Particle::Kind kind)
{
    int size = (int)cell->particles_of_kind(kind).size();

    double vel_x = 0;

    if(size == 0)
    {
        return vel_x;
    }

    for(Particle & particle : cell->particles_of_kind(kind))
    {
        vel_x += particle.vx;
    }

    return vel_x / (double)size;
}

double idic_1d::eps_asterisk(Cell * cell)
{
    //TODO int or not?
    double dust_size = cell->dust_particles.size();
    double gas_size = cell->gas_particles.size();

    double dust_mass = 0;

    if(gas_size == 0)
    {
        return 0;
    }
    else
    {
        if(dust_size != 0)
        {
            dust_mass = cell->dust_particles.at(0).mass;
        }

        double gas_mass = cell->gas_particles.at(0).mass;

        return dust_mass * dust_size / (gas_mass * gas_size);
    }
}

double idic_1d::x_through_vel(Cell * cell)
{
    return idic_1d::vel_asterisk(cell, Particle::Kind::Gas) - idic_1d::vel_asterisk(cell, Particle::Kind::Dust);
}

double idic_1d::y_through_vel(Cell * cell)
{
    return idic_1d::vel_asterisk(cell, Particle::Kind::Gas) + idic_1d::vel_asterisk(cell, Particle::Kind::Dust)
            * idic_1d::eps_asterisk(cell);
}

double idic_1d::find_x(Cell * cell)
{
    double tau = Params::get_instance().tau;

    double x_prev = idic_1d::x_through_vel(cell);
    double a_p_asrt = idic_1d::pressure_term_asterisk(cell);

    double front_next_x = 1. / tau + (idic_1d::eps_asterisk(cell) + 1.) / nl_idic::t_stop_asterisk(cell);

    return (x_prev / tau + a_p_asrt) / front_next_x;
}

double idic_1d::find_y(Cell * cell)
{
    double tau = Params::get_instance().tau;

    double y_prev = idic_1d::y_through_vel(cell);
    double a_p_asrt = idic_1d::pressure_term_asterisk(cell);

    return y_prev + a_p_asrt * tau;
}

double idic_1d::find_gas_vel_asterisk(Cell * cell)
{
    return (idic_1d::find_y(cell) + idic_1d::find_x(cell) * idic_1d::eps_asterisk(cell)) / (1. + idic_1d::eps_asterisk(cell));
}

double idic_1d::find_dust_vel_asterisk(Cell * cell)
{
    return (idic_1d::find_y(cell) - idic_1d::find_x(cell)) / (1. + idic_1d::eps_asterisk(cell));
}

double idic_1d::find_gas_velocity(Particle * particle, Cell * cell)
{
    double tau = Params::get_instance().tau;

    double prev_vel = particle->vx;

    double result = 0;

    double front_next_vel = 0;

    if(nl_idic::t_stop_asterisk(cell) == 0)
    {
        front_next_vel = 1. / tau;

        result = (prev_vel / tau + Dusty_shock_1d::pressure_term(particle, cell)) / front_next_vel;
    }

    if(nl_idic::t_stop_asterisk(cell) != 0)
    {
        front_next_vel = 1. / tau + idic_1d::eps_asterisk(cell) / nl_idic::t_stop_asterisk(cell);

        result = (prev_vel / tau + idic_1d::find_dust_vel_asterisk(cell) * (idic_1d::eps_asterisk(cell)
                  / nl_idic::t_stop_asterisk(cell))
                  + Dusty_shock_1d::pressure_term(particle, cell)) / front_next_vel;
    }

    assert(!__isnan(result));

    return result;
}

//TODO assert Kind
double idic_1d::find_dust_velocity(Particle * particle, Cell * cell)
{
    double tau = Params::get_instance().tau;

    double front_next_vel = 1. / tau + 1. / nl_idic::t_stop_asterisk(cell);

    double prev_vel = particle->vx;
    double result = 0;

    if(nl_idic::t_stop_asterisk(cell) != 0)
    {
        result = (prev_vel / tau + idic_1d::find_gas_vel_asterisk(cell) / nl_idic::t_stop_asterisk(cell)) / front_next_vel;
    }

    assert(!__isnan(result));

    return result;
}

Grid Dusty_shock_1d::do_time_step(Grid & old_grid, int step_num)
{
    std::ofstream f_dust;
    std::ofstream f_gas;

    Params & params = Params::get_instance();
    Dusty_shock_params & shock_params = Dusty_shock_params::get_instance();

    if (PRINT_FILE)
    {
        //char filename_gas[512];
        //sprintf(filename_gas, (OUTPUT_PATH + "dust_h%lg_cfl%lg_n%d_hc%lg_part_%0d.dat").c_str(), params.h, params.tau / params.h, shock_params.dust_real_particles, params.grid_step_x, step_num);
        //f_dust.open(filename_dust);

        std::stringstream filename_dust;
        filename_dust << OUTPUT_PATH
                     << "dust_h" << params.h
                     << "_cfl" << (params.tau / params.h)
                     << "_n" << shock_params.dust_real_particles
                     << "_hc" << params.grid_step_x
                     << "_part_" << step_num;
        f_dust.open(filename_dust.str());

        std::stringstream filename_gas;
        filename_gas << OUTPUT_PATH
                     << "gas_h" << params.h
                     << "_cfl" << (params.tau / params.h)
                     << "_n" << shock_params.gas_real_particles
                     << "_hc" << params.grid_step_x
                     << "_part_" << step_num;
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
                    if (PRINT_FILE)
                    {
                        f_gas << particle.x << " " << particle.vx << " " << particle.density << " "
                              << particle.pressure << " " << particle.energy << std::endl;
                    }

                    Particle new_particle(particle);
                    Point new_coords = find_new_coordinates(particle);

                    double new_vel = idic_1d::find_gas_velocity(&particle, &cell);
                    new_particle.density = NAN;

                    new_particle.set_coordinates(new_coords.x, new_coords.y, new_coords.z);
                    new_particle.set_velocities(new_vel, particle.vy, particle.vz);

                    new_particle.energy = Sod_tube_1d::find_new_energy(&particle, &cell);

                    next_grid.with_copy_of(new_particle);
                }
                for(Particle & particle : cell.dust_particles)
                {
                    if (PRINT_FILE)
                    {
                        f_dust << particle.x << " " << particle.vx << " " << particle.density << std::endl;
                    }

                    Particle new_particle(particle);
                    Point new_coords = find_new_coordinates(particle);

                    double new_vel = idic_1d::find_dust_velocity(&particle, &cell);
                    new_particle.density = NAN;

                    new_particle.set_coordinates(new_coords.x, new_coords.y, new_coords.z);
                    new_particle.set_velocities(new_vel, particle.vy, particle.vz);

                    next_grid.with_copy_of(new_particle);
                }
            }
        }
    }

    recalc_density(next_grid, Particle::Kind::Gas);
    recalc_density(next_grid, Particle::Kind::Dust);
    Sod_tube_1d::recalc_pressure(next_grid, Particle::Kind::Gas);

    return next_grid;
}