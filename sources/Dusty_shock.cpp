#include <fstream>
#include "Dusty_shock.h"

Point pressure_term(Particle * particle, Cell * cell)
{
    Params & params = Params::get_instance();

    Point a_p = {0, 0, 0};
    Point kernel_grad = {0, 0, 0};

    for (Cell * neighbour : cell->get_neighbours())
    {
        for (Particle & p : neighbour->gas_particles)
        { // TODO remove direct access

            kernel_grad = kernel_gradient(*(particle), p, params.dimensions);

            assert(!__isnan(p.density));
            assert(!__isnan(kernel_grad.x));
            double viscosity = viscosity_3d::find_viscosity(particle, &p);

            double term = particle->pressure / pow(particle->density, 2) + p.pressure / pow(p.density, 2) + viscosity;

            a_p = a_p + kernel_grad * term;
        }
    }

    a_p = a_p * (- particle->mass);

    return a_p;
}

Point idic::pressure_term_asterisk(Cell * cell)
{
    double gas_size = cell->gas_particles.size();

    double a_p_x_astr = 0;
    double a_p_y_astr = 0;
    double a_p_z_astr = 0;

    if(gas_size == 0)
    {
        return {0, 0, 0};
    }
    else
    {
        for (Particle & p : cell->gas_particles)
        {
            a_p_x_astr += pressure_term(&p, cell).x;
            a_p_y_astr += pressure_term(&p, cell).y;
            a_p_z_astr += pressure_term(&p, cell).z;
        }
    }

    return {a_p_x_astr / gas_size, a_p_y_astr / gas_size, a_p_z_astr / gas_size};
}

Point idic::vel_asterisk(Cell * cell, Particle::Kind kind)
{
    int size = (int)cell->particles_of_kind(kind).size();

    double vel_x = 0;
    double vel_y = 0;
    double vel_z = 0;

    if(size == 0)
    {
        return {vel_x, vel_y, vel_z};
    }

    for(Particle & particle : cell->particles_of_kind(kind))
    {
        vel_x += particle.vx;
        vel_y += particle.vy;
        vel_z += particle.vz;
    }

    return {vel_x / (double)size, vel_y / (double)size, vel_z / (double)size};
}

double idic::eps_asterisk(Cell * cell)
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

double idic::density_asterisk(Cell * cell, Particle::Kind kind)
{
    int size = (int)cell->particles_of_kind(kind).size();

    double result = 0;

    if(size == 0)
    {
        return result;
    }

    for(Particle & particle : cell->particles_of_kind(kind))
    {
        result += particle.density;
    }

    return result / (double)size;
}

double idic::t_stop_asterisk(Cell * cell)
{
    if(Params::get_instance().K != 0)
    {
        return density_asterisk(cell, Particle::Kind::Dust) / Params::get_instance().K;
    }
    return 0;
}

Point idic::x_through_vel(Cell * cell)
{
    return vel_asterisk(cell, Particle::Kind::Gas) - vel_asterisk(cell, Particle::Kind::Dust);
}

Point idic::y_through_vel(Cell * cell)
{
    return vel_asterisk(cell, Particle::Kind::Gas) + vel_asterisk(cell, Particle::Kind::Dust) * eps_asterisk(cell);
}

Point idic::find_x(Cell * cell)
{
    double tau = Params::get_instance().tau;

    Point x_prev = x_through_vel(cell);
    Point a_p_asrt = pressure_term_asterisk(cell);

    double front_next_x = 1. / tau + (eps_asterisk(cell) + 1.) / t_stop_asterisk(cell);

    return (x_prev / tau + a_p_asrt) / front_next_x;
}

Point idic::find_y(Cell * cell)
{
    double tau = Params::get_instance().tau;

    Point y_prev = y_through_vel(cell);
    Point a_p_asrt = pressure_term_asterisk(cell);

    return y_prev + a_p_asrt * tau;
}

Point idic::find_gas_vel_asterisk(Cell * cell)
{
    return (find_y(cell) + find_x(cell) * eps_asterisk(cell)) / (1. + eps_asterisk(cell));
}

Point idic::find_dust_vel_asterisk(Cell * cell)
{
    return (find_y(cell) - find_x(cell)) / (1. + eps_asterisk(cell));
}

//TODO assert Kind
Point idic::find_gas_velocity(Particle * particle, Cell * cell)
{
    double tau = Params::get_instance().tau;

    Point prev_vel = {particle->vx, particle->vy, particle->vz};

    Point result = {0, 0, 0};

    double front_next_vel = 0;

    if(t_stop_asterisk(cell) == 0)
    {
        front_next_vel = 1. / tau;

        result = (prev_vel / tau + pressure_term(particle, cell)) / front_next_vel;
    }

    if(t_stop_asterisk(cell) != 0)
    {
        front_next_vel = 1. / tau + eps_asterisk(cell) / t_stop_asterisk(cell);

        result = (prev_vel / tau + find_dust_vel_asterisk(cell) * (eps_asterisk(cell) / t_stop_asterisk(cell))
                  + pressure_term(particle, cell)) / front_next_vel;
    }

    assert(!__isnan(result.x));
    assert(!__isnan(result.y));
    assert(!__isnan(result.z));

    return result;
}

//TODO assert Kind
Point idic::find_dust_velocity(Particle * particle, Cell * cell)
{
    double tau = Params::get_instance().tau;

    double front_next_vel = 1. / tau + 1. / t_stop_asterisk(cell);

    Point prev_vel = {particle->vx, particle->vy, particle->vz};
    Point result = {0, 0, 0};

    if(t_stop_asterisk(cell) != 0)
    {
        result = (prev_vel / tau + find_gas_vel_asterisk(cell) / t_stop_asterisk(cell)) / front_next_vel;
    }

    assert(!__isnan(result.x));
    assert(!__isnan(result.y));
    assert(!__isnan(result.z));

    return result;
}

Point monaghan::find_gas_velocity(Particle * particle, Cell * cell)
{
    Params & params = Params::get_instance();

    Point a_p = pressure_term(particle, cell);

    double drag_x = 0;
    double drag_y = 0;
    double drag_z = 0;

    double dust_mass = 0;
    if((int)cell->dust_particles.size() != 0)
    {
        dust_mass = cell->dust_particles.at(0).mass;
    }

    Point vel = {0, 0, 0};
    Point coord = {0, 0, 0};

    for(Cell * neighbour : cell->get_neighbours())
    {
        for(Particle & p : neighbour->dust_particles)
        {
            vel = {particle->vx - p.vx, particle->vy - p.vy, particle->vz - p.vz};
            coord = {particle->x - p.x, particle->y - p.y, particle->z - p.z};

            double v_dot_r = dot_product(vel, coord);
            double r_magnitude = magnitude(coord);

            drag_x += (1. / particle->density / p.density) * v_dot_r / (r_magnitude + params.eta_squared)
                    * coord.x * kernel(*particle, p, params.dimensions);

            drag_y += (1. / particle->density / p.density) * v_dot_r / (r_magnitude + params.eta_squared)
                      * coord.y * kernel(*particle, p, params.dimensions);

            drag_z += (1. / particle->density / p.density) * v_dot_r / (r_magnitude + params.eta_squared)
                      * coord.z * kernel(*particle, p, params.dimensions);
        }
    }

    drag_x *= params.sigma * dust_mass * params.K;
    drag_y *= params.sigma * dust_mass * params.K;
    drag_z *= params.sigma * dust_mass * params.K;

    return {params.tau * (a_p.x - drag_x) + particle->vx, params.tau * (a_p.y - drag_y) + particle->vy,
            params.tau * (a_p.z - drag_z) + particle->vz};
}

Point monaghan::find_dust_velocity(Particle * particle, Cell * cell)
{
    Params & params = Params::get_instance();

    double drag_x = 0;
    double drag_y = 0;
    double drag_z = 0;

    double gas_mass = 0;
    if((int)cell->gas_particles.size() != 0)
    {
        gas_mass = cell->gas_particles.at(0).mass;
    }

    Point vel = {0, 0, 0};
    Point coord = {0, 0, 0};

    for(Cell * neighbour : cell->get_neighbours())
    {
        for (Particle & p : neighbour->gas_particles)
        {
            vel = {particle->vx - p.vx, particle->vy - p.vy, particle->vz - p.vz};
            coord = {particle->x - p.x, particle->y - p.y, particle->z - p.z};

            double v_dot_r = dot_product(vel, coord);
            double r_magnitude = magnitude(coord);

            drag_x += (1. / particle->density / p.density) * v_dot_r / (r_magnitude + params.eta_squared)
                      * coord.x * kernel(*particle, p, params.dimensions);

            drag_y += (1. / particle->density / p.density) * v_dot_r / (r_magnitude + params.eta_squared)
                      * coord.y * kernel(*particle, p, params.dimensions);

            drag_z += (1. / particle->density / p.density) * v_dot_r / (r_magnitude + params.eta_squared)
                      * coord.z * kernel(*particle, p, params.dimensions);
        }
    }

    drag_x *= params.sigma * gas_mass * params.K;
    drag_y *= params.sigma * gas_mass * params.K;
    drag_z *= params.sigma * gas_mass * params.K;

    return {params.tau * drag_x + particle->vx, params.tau * drag_y + particle->vy, params.tau * drag_z + particle->vz};
}

Grid Dusty_shock_3d::init()
{
    Params & params = Params::get_instance();
    Dusty_shock_params & shock_params = Dusty_shock_params::get_instance();

    Grid init(params.grid_step_x, params.grid_step_y, params.grid_step_z, params.border1, params.border2);

    init = Sod_tube_3d::init_with_boundaries();

    uint real_1d_particles = shock_params.dust_real_particles;
    uint all_1d_particles = shock_params.dust_real_particles + shock_params.dust_image_particles;

    double dust_l2r = shock_params.dust_dens_left / shock_params.dust_dens_right;

    // считаем количества реальных частиц
    uint real_right_1d_p_num = (uint) round((double) shock_params.dust_real_particles / (dust_l2r + 1));
    uint real_left_1d_p_num = shock_params.dust_real_particles - real_right_1d_p_num;

    // аналогично для виртуальных
    uint image_right_1d_p_num = (uint) round((double) shock_params.dust_image_particles / (dust_l2r + 1));
    uint image_left_1d_p_num = shock_params.dust_image_particles - image_right_1d_p_num;

    std::vector<double> coord(real_1d_particles);
    std::vector<double> image_coord(all_1d_particles);

    fill_initial_sod_coord(coord, image_coord, real_left_1d_p_num, real_right_1d_p_num, image_left_1d_p_num,
                           image_right_1d_p_num);

    double yz_length = shock_params.yz_right - shock_params.yz_left;
    double step_yz = yz_length / (double)shock_params.dust_yz_particles;

    double image_left_length = image_coord.at(image_left_1d_p_num + real_left_1d_p_num) - image_coord.at(0);
    double mass = image_left_length * yz_length * yz_length * shock_params.dust_dens_left /
                  (image_left_1d_p_num + real_left_1d_p_num) / (double)shock_params.dust_yz_particles
                  / (double)shock_params.dust_yz_particles;

    for(int i = 0; i < (int) all_1d_particles; ++i)
    {
        for(int j = 0; j < (int) shock_params.dust_yz_particles; ++j)
        {
            for(int k = 0; k < (int) shock_params.dust_yz_particles; ++k)
            {
                Particle particle_main(Particle::Kind::Dust, image_coord.at(i), step_yz / 2. + j * step_yz,
                                       step_yz / 2. + k * step_yz);

                particle_main.mass = mass;
                particle_main.set_velocities(0, 0, 0);
                //делим объем плоскостью x = membrane
                if(particle_main.x <= shock_params.membrane)
                {
                    particle_main.density = shock_params.dust_dens_left;
                }
                else
                {
                    particle_main.density = shock_params.dust_dens_right;
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

                for(Particle & particle : cell.particles_of_kind(Particle::Kind::Dust))
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

                cell.dust_particles.erase(cell.dust_particles.begin(), cell.dust_particles.end());
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

                cell.dust_particles.erase(cell.dust_particles.begin(), cell.dust_particles.end());
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

                cell.dust_particles.erase(cell.dust_particles.begin(), cell.dust_particles.end());
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

                cell.dust_particles.erase(cell.dust_particles.begin(), cell.dust_particles.end());
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

                for(Particle & particle : cell.particles_of_kind(Particle::Kind::Dust))
                {
                    Particle ym_particle(particle);
                    ym_particle.set_coordinates(particle.x, particle.y - yz_length, particle.z);
                    init.with_copy_of(ym_particle);

                    Particle yp_particle(particle);
                    yp_particle.set_coordinates(particle.x, particle.y + yz_length, particle.z);
                    init.with_copy_of(yp_particle);

                    Particle zm_particle(particle);
                    zm_particle.set_coordinates(particle.x, particle.y, particle.z - yz_length);
                    init.with_copy_of(zm_particle);

                    Particle zp_particle(particle);
                    zp_particle.set_coordinates(particle.x, particle.y, particle.z + yz_length);
                    init.with_copy_of(zp_particle);

                    Particle particle_yp_zp(particle);
                    particle_yp_zp.set_coordinates(particle.x, particle.y + yz_length, particle.z + yz_length);
                    init.with_copy_of(particle_yp_zp);

                    Particle particle_ym_zp(particle);
                    particle_ym_zp.set_coordinates(particle.x, particle.y - yz_length, particle.z + yz_length);
                    init.with_copy_of(particle_ym_zp);

                    Particle particle_yp_zm(particle);
                    particle_yp_zm.set_coordinates(particle.x, particle.y + yz_length, particle.z - yz_length);
                    init.with_copy_of(particle_yp_zm);

                    Particle particle_ym_zm(particle);
                    particle_ym_zm.set_coordinates(particle.x, particle.y - yz_length, particle.z - yz_length);
                    init.with_copy_of(particle_ym_zm);
                }
            }
        }
    }

    std::ofstream dust_shock_file;
    std::ofstream gas_shock_file;
    if (PRINT_FILE)
    {
        dust_shock_file.open(OUTPUT_PATH + "dust_shock.dat");
        gas_shock_file.open(OUTPUT_PATH + "gas_shock.dat");
    }

    for(int i = 0; i < init.x_size; ++i)
    {
        for(int j = 0; j < init.y_size; ++j)
        {
            for(int k = 0; k < init.z_size; ++k)
            {
                Cell & cell = init.cells[i][j][k];

                if(PRINT_FILE)
                {
                    for (Particle &particle : cell.dust_particles)
                    {
                        dust_shock_file << particle.x << " " << particle.y << " " << particle.z << " "
                                        << particle.vx << " " << particle.vy << " " << particle.vz << " "
                                        << particle.density
                                        << std::endl;
                    }

                    for (Particle &particle : cell.gas_particles)
                    {
                        gas_shock_file << particle.x << " " << particle.y << " " << particle.z << " "
                                       << particle.vx << " " << particle.vy << " " << particle.vz << " "
                                       << particle.density << " " << particle.pressure << " " << particle.energy
                                       << std::endl;
                    }
                }
            }
        }
    }

    return init;
}

Grid Dusty_shock_3d::do_time_step(Grid & old_grid, int step_num, bool /*isIDIC*/)
{
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

                for(Particle & particle : cell.gas_particles)
                {
                    Particle new_particle(particle);
                    Point new_coords = find_new_coordinates(particle);

                    Point new_vel = {NAN, NAN, NAN};

                    if(IS_IDIC)
                    {
                        new_vel = idic::find_gas_velocity(&particle, &cell);
                    }
                    else
                    {
                        new_vel = monaghan::find_gas_velocity(&particle, &cell);
                    }

                    new_particle.density = NAN;
                    new_particle.dbg_state = 2;

                    new_particle.set_coordinates(new_coords.x, new_coords.y, new_coords.z);
                    new_particle.set_velocities(new_vel.x, new_vel.y, new_vel.z);

                    new_particle.energy = Sod_tube_3d::find_new_energy(&particle, &cell);

                    next_grid.with_copy_of(new_particle);

                    Particle ym_particle(particle);
                    ym_particle.set_coordinates(particle.x, particle.y - yz_length, particle.z);
                    next_grid.with_copy_of(ym_particle);

                    Particle yp_particle(particle);
                    yp_particle.set_coordinates(particle.x, particle.y + yz_length, particle.z);
                    next_grid.with_copy_of(yp_particle);

                    Particle zm_particle(particle);
                    zm_particle.set_coordinates(particle.x, particle.y, particle.z - yz_length);
                    next_grid.with_copy_of(zm_particle);

                    Particle zp_particle(particle);
                    zp_particle.set_coordinates(particle.x, particle.y, particle.z + yz_length);
                    next_grid.with_copy_of(zp_particle);

                    Particle particle_yp_zp(particle);
                    particle_yp_zp.set_coordinates(particle.x, particle.y + yz_length, particle.z + yz_length);
                    next_grid.with_copy_of(particle_yp_zp);

                    Particle particle_ym_zp(particle);
                    particle_ym_zp.set_coordinates(particle.x, particle.y - yz_length, particle.z + yz_length);
                    next_grid.with_copy_of(particle_ym_zp);

                    Particle particle_yp_zm(particle);
                    particle_yp_zm.set_coordinates(particle.x, particle.y + yz_length, particle.z - yz_length);
                    next_grid.with_copy_of(particle_yp_zm);

                    Particle particle_ym_zm(particle);
                    particle_ym_zm.set_coordinates(particle.x, particle.y - yz_length, particle.z - yz_length);
                    next_grid.with_copy_of(particle_ym_zm);
                }

                for(Particle & particle : cell.dust_particles)
                {
                    Particle new_particle(particle);
                    Point new_coords = find_new_coordinates(particle);

                    Point new_vel = {NAN, NAN, NAN};

                    if(IS_IDIC)
                    {
                        new_vel = idic::find_dust_velocity(&particle, &cell);
                    }
                    else
                    {
                        new_vel = monaghan::find_dust_velocity(&particle, &cell);
                    }

                    new_particle.density = NAN;
                    new_particle.dbg_state = 2;

                    new_particle.set_coordinates(new_coords.x, new_coords.y, new_coords.z);
                    new_particle.set_velocities(new_vel.x, new_vel.y, new_vel.z);

                    next_grid.with_copy_of(new_particle);

                    Particle ym_particle(particle);
                    ym_particle.set_coordinates(particle.x, particle.y - yz_length, particle.z);
                    next_grid.with_copy_of(ym_particle);

                    Particle yp_particle(particle);
                    yp_particle.set_coordinates(particle.x, particle.y + yz_length, particle.z);
                    next_grid.with_copy_of(yp_particle);

                    Particle zm_particle(particle);
                    zm_particle.set_coordinates(particle.x, particle.y, particle.z - yz_length);
                    next_grid.with_copy_of(zm_particle);

                    Particle zp_particle(particle);
                    zp_particle.set_coordinates(particle.x, particle.y, particle.z + yz_length);
                    next_grid.with_copy_of(zp_particle);

                    Particle particle_yp_zp(particle);
                    particle_yp_zp.set_coordinates(particle.x, particle.y + yz_length, particle.z + yz_length);
                    next_grid.with_copy_of(particle_yp_zp);

                    Particle particle_ym_zp(particle);
                    particle_ym_zp.set_coordinates(particle.x, particle.y - yz_length, particle.z + yz_length);
                    next_grid.with_copy_of(particle_ym_zp);

                    Particle particle_yp_zm(particle);
                    particle_yp_zm.set_coordinates(particle.x, particle.y + yz_length, particle.z - yz_length);
                    next_grid.with_copy_of(particle_yp_zm);

                    Particle particle_ym_zm(particle);
                    particle_ym_zm.set_coordinates(particle.x, particle.y - yz_length, particle.z - yz_length);
                    next_grid.with_copy_of(particle_ym_zm);
                }

                /*
                for(Particle * particle : cell.get_all_particles())
                {

                    Particle new_particle(*particle);
                    Point new_coords = find_new_coordinates(*particle);

                    Point new_vel = {NAN, NAN, NAN};

                    if(IS_IDIC)
                    {
                        if(particle->kind == Particle::Kind::Gas)
                        {
                            new_vel = idic::find_gas_velocity(particle, &cell);
                            new_particle.energy = Sod_tube_3d::find_new_energy(particle, &cell);
                        }
                        else if(particle->kind == Particle::Kind::Dust)
                        {
                            new_vel = idic::find_dust_velocity(particle, &cell);
                        }
                        else
                        {
                            assert(false);
                        }
                    }
                    else
                        {
                            if(particle->kind == Particle::Kind::Gas)
                            {
                                new_vel = monaghan::find_gas_velocity(particle, &cell);
                                new_particle.energy = Sod_tube_3d::find_new_energy(particle, &cell);
                            }
                            else if(particle->kind == Particle::Kind::Dust)
                            {
                                new_vel = monaghan::find_dust_velocity(particle, &cell);
                            }
                            else
                            {
                                assert(false);
                            }
                        }

                    new_particle.density = NAN;
                    new_particle.dbg_state = 2;

                    new_particle.set_coordinates(new_coords.x, new_coords.y, new_coords.z);
                    new_particle.set_velocities(new_vel.x, new_vel.y, new_vel.z);

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
                 */
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

                for(Particle & particle : cell.gas_particles)
                {
                    particle.density = find_density(&particle, &cell);
                    particle.pressure = Sod_tube_1d::find_new_pressure(&particle);
                    particle.dbg_state = 1;
                }
                for(Particle & particle : cell.dust_particles)
                {
                    particle.density = find_density(&particle, &cell);
                    particle.dbg_state = 1;
                }
                /*
                for(Particle * particle : cell.get_all_particles())
                {
                    particle->density = find_density(particle, &cell);
                    if(particle->kind == Particle::Kind::Gas)
                    {
                        particle->pressure = Sod_tube_1d::find_new_pressure(particle);
                    }
                    particle->dbg_state = 1;
                }
                 */
            }
        }
    }

    for(int i = x_left_border; i < x_right_border; ++i)
    {
        for(int j = 0; j < next_grid.y_size; ++j)
        {
            for(int k = yz_right_border; k < next_grid.z_size; ++k)
            {
                Cell & cell = next_grid.cells[i][j][k];

                cell.gas_particles.erase(cell.gas_particles.begin(), cell.gas_particles.end());
                cell.dust_particles.erase(cell.dust_particles.begin(), cell.dust_particles.end());
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
                cell.dust_particles.erase(cell.dust_particles.begin(), cell.dust_particles.end());
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
                cell.dust_particles.erase(cell.dust_particles.begin(), cell.dust_particles.end());
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
                cell.dust_particles.erase(cell.dust_particles.begin(), cell.dust_particles.end());
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

    std::ofstream f_dust;
    std::ofstream f_gas;

    if (PRINT_FILE)
    {
        f_dust.open(OUTPUT_PATH + "dust_part_" + std::to_string(step_num) + ".dat");
        f_gas.open(OUTPUT_PATH + "gas_part_" + std::to_string(step_num) + ".dat");
    }

    for(int i = x_left_border; i < x_right_border; ++i)
    {
        for(int j = yz_left_border; j < yz_right_border; ++j)
        {
            for(int k = yz_left_border; k < yz_right_border; ++k)
            {
                Cell & cell = next_grid.cells[i][j][k];

                if(PRINT_FILE) {
                    for (Particle &particle : cell.dust_particles)
                    {
                        f_dust << particle.x << " " << particle.y << " " << particle.z << " "
                               << particle.vx << " " << particle.vy << " " << particle.vz << " "
                               << particle.density
                               << std::endl;
                    }

                    for (Particle &particle : cell.gas_particles)
                    {
                        f_gas << particle.x << " " << particle.y << " " << particle.z << " "
                              << particle.vx << " " << particle.vy << " " << particle.vz << " "
                              << particle.density << " " << particle.pressure << " " << particle.energy
                              << std::endl;
                    }
                }
            }
        }
    }

    return next_grid;
}