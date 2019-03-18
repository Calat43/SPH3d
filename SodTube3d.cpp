#include "SodTube3d.h"

/*
Grid Sod_tube_3d::init()
{
    Params & params = Params::get_instance();

    Grid init(params.grid_step_x, params.grid_step_y, params.grid_step_z, params.border1, params.border2);

    uint real_1d_particles = params.real_particles;
    uint all_1d_particles = params.real_particles + params.image_particles;

    double gas_l2r = params.dens_left / params.dens_right;

    // считаем количества реальных частиц
    uint real_right_1d_p_num = (uint) round((double) params.real_particles / (gas_l2r + 1));
    uint real_left_1d_p_num = params.real_particles - real_right_1d_p_num;

    // аналогично для виртуальных
    uint image_right_1d_p_num = (uint) round((double) params.image_particles / (gas_l2r + 1));
    uint image_left_1d_p_num = params.image_particles - image_right_1d_p_num;

    std::vector<double> coord(real_1d_particles);
    std::vector<double> image_coord(all_1d_particles);

    fill_initial_sod_coord(coord, image_coord, real_left_1d_p_num, real_right_1d_p_num, image_left_1d_p_num,
                           image_right_1d_p_num);

    //TODO mass
    //По осям z и у ставим по yz_particles частиц на длину 0.1
    int yz_particles = 20;
    double step_yz = 0.1 / 20.;

    double image_left_lenght = image_coord.at(image_left_1d_p_num + real_left_1d_p_num) - image_coord.at(0);
    double mass = image_left_lenght * 0.1 * 0.1 * Params::get_instance().dens_left /
                  (image_left_1d_p_num + real_left_1d_p_num) / (double)yz_particles / (double)yz_particles;

    for(int i = 0; i < all_1d_particles; ++i)
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
                if(particle_main.x <= params.membrane)
                {
                    particle_main.density = params.dens_left;
                    particle_main.pressure = params.press_left;
                    particle_main.energy = params.ener_left;
                }
                else
                {
                    particle_main.density = params.dens_right;
                    particle_main.pressure = params.press_right;
                    particle_main.energy = params.ener_right;
                }
                init.with_copy_of(particle_main);

                Particle particle_zp(particle_main);
                particle_zp.set_coordinates(particle_main.x, particle_main.y, particle_main.z + 0.1);
                init.with_copy_of(particle_zp);

                Particle particle_zm(particle_main);
                particle_zm.set_coordinates(particle_main.x, particle_main.y, particle_main.z - 0.1);
                init.with_copy_of(particle_zm);

                Particle particle_yp(particle_main);
                particle_yp.set_coordinates(particle_main.x, particle_main.y + 0.1, particle_main.z);
                init.with_copy_of(particle_yp);

                Particle particle_ym(particle_main);
                particle_ym.set_coordinates(particle_main.x, particle_main.y - 0.1, particle_main.z);
                init.with_copy_of(particle_ym);
            }
        }
    }

    //частицы в расчетной области
    for(int i = 10; i < 40; ++i)
    {
        for(int j = 5; j < 10; ++j)
        {
            for(int k = 5; k < 10; ++k)
            {
                Cell & cell = init.cells[i][j][k];

                for(Particle * particle : cell.get_all_particles())
                {
                    particle->density = find_density(particle, &cell);
                    particle->pressure = Sod_tube_1d::find_new_pressure(particle);
                }
            }
        }
    }

    for(int i = 10; i < 40; ++i)
    {
        for(int j = 0; j < 5; ++j)
        {
            for(int k = 0; k < 5; ++k)
            {
                Cell & cell = init.cells[i][j][k];

                cell.gas_particles.erase(cell.gas_particles.begin(), cell.gas_particles.end());
            }
        }
    }
    for(int i = 10; i < 40; ++i)
    {
        for(int j = 0; j < 5; ++j)
        {
            for(int k = 10; k < 15; ++k)
            {
                Cell & cell = init.cells[i][j][k];

                cell.gas_particles.erase(cell.gas_particles.begin(), cell.gas_particles.end());
            }
        }
    }
    for(int i = 10; i < 40; ++i)
    {
        for(int j = 10; j < 15; ++j)
        {
            for(int k = 0; k < 5; ++k)
            {
                Cell & cell = init.cells[i][j][k];

                cell.gas_particles.erase(cell.gas_particles.begin(), cell.gas_particles.end());
            }
        }
    }
    for(int i = 10; i < 40; ++i)
    {
        for(int j = 10; j < 15; ++j)
        {
            for(int k = 10; k < 15; ++k)
            {
                Cell & cell = init.cells[i][j][k];

                cell.gas_particles.erase(cell.gas_particles.begin(), cell.gas_particles.end());
            }
        }
    }

    //частицы в расчетной области
    for(int i = 10; i < 40; ++i)
    {
        for(int j = 5; j < 10; ++j)
        {
            for(int k = 5; k < 10; ++k)
            {
                Cell & cell = init.cells[i][j][k];

                for(Particle * particle : cell.get_all_particles())
                {
                    Particle particle_zp(*particle);
                    particle_zp.set_coordinates(particle->x, particle->y, particle->z + 0.1);
                    init.with_copy_of(particle_zp);

                    Particle particle_zm(*particle);
                    particle_zm.set_coordinates(particle->x, particle->y, particle->z - 0.1);
                    init.with_copy_of(particle_zm);

                    Particle particle_yp(*particle);
                    particle_yp.set_coordinates(particle->x, particle->y + 0.1, particle->z);
                    init.with_copy_of(particle_yp);

                    Particle particle_ym(*particle);
                    particle_ym.set_coordinates(particle->x, particle->y - 0.1, particle->z);
                    init.with_copy_of(particle_ym);
                }
            }
        }
    }

    //TODO PRINT DENSITIES
    char init_dens[] = "init_dens";
    print_grid_dens(init, init_dens);

    return init;
}
*/

Grid Sod_tube_3d::init()
{
    Params & params = Params::get_instance();

    Grid init(params.grid_step_x, params.grid_step_y, params.grid_step_z, params.border1, params.border2);

    uint real_1d_particles = params.real_particles;
    uint all_1d_particles = params.real_particles + params.image_particles;

    double gas_l2r = params.dens_left / params.dens_right;

    // считаем количества реальных частиц
    uint real_right_1d_p_num = (uint) round((double) params.real_particles / (gas_l2r + 1));
    uint real_left_1d_p_num = params.real_particles - real_right_1d_p_num;

    // аналогично для виртуальных
    uint image_right_1d_p_num = (uint) round((double) params.image_particles / (gas_l2r + 1));
    uint image_left_1d_p_num = params.image_particles - image_right_1d_p_num;

    std::vector<double> coord(real_1d_particles);
    std::vector<double> image_coord(all_1d_particles);

    fill_initial_sod_coord(coord, image_coord, real_left_1d_p_num, real_right_1d_p_num, image_left_1d_p_num,
                           image_right_1d_p_num);

    //TODO mass
    //По осям z и у ставим по yz_particles частиц на длину 0.1
    int yz_particles = 20;
    double step_yz = 0.1 / 20.;

    double image_left_lenght = image_coord.at(image_left_1d_p_num + real_left_1d_p_num) - image_coord.at(0);
    double mass = image_left_lenght * 0.1 * 0.1 * Params::get_instance().dens_left /
                  (image_left_1d_p_num + real_left_1d_p_num) / (double)yz_particles / (double)yz_particles;

    for(int i = 0; i < all_1d_particles; ++i)
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
                if(particle_main.x <= params.membrane)
                {
                    particle_main.density = params.dens_left;
                    particle_main.pressure = params.press_left;
                    particle_main.energy = params.ener_left;
                }
                else
                {
                    particle_main.density = params.dens_right;
                    particle_main.pressure = params.press_right;
                    particle_main.energy = params.ener_right;
                }
                init.with_copy_of(particle_main);
            }
        }
    }
    recalc_density(init);
    Sod_tube_1d::recalc_pressure(init);

    if(PRINT_STUFF)
    {
        char init_dens[] = "init_dens";
        print_grid_dens(init, init_dens);
    }

    return init;
}

namespace viscosity_3d
{
    double find_sound_speed(Particle * particle)
    {
        assert(!__isnan(particle->pressure));
        assert(!__isnan(particle->density));
        assert(!__isnan(sqrt(Params::get_instance().gamma * particle->pressure / particle->density)));
        return sqrt(Params::get_instance().gamma * particle->pressure / particle->density);
    }

    double find_mu(Point vel_ab, Point coord_ab)
    {
        return (Params::get_instance().h * dot_product(vel_ab, coord_ab)) /
               (pow(magnitude(coord_ab), 2) + pow(Params::get_instance().nu, 2));
    }

    double find_viscosity(Particle * p1, Particle * p2)
    {
        Params & params = Params::get_instance();

        if(params.have_viscosity)
        {
            double vel_x_ab = p1->vx - p2->vx;
            double vel_y_ab = p1->vy - p2->vy;
            double vel_z_ab = p1->vz - p2->vz;

            double coord_x_ab = p1->x - p2->x;
            double coord_y_ab = p1->y - p2->y;
            double coord_z_ab = p1->z - p2->z;

            Point vel_ab(vel_x_ab, vel_y_ab, vel_z_ab);
            Point coord_ab(coord_x_ab, coord_y_ab, coord_z_ab);

            if(dot_product(vel_ab, coord_ab) >= 0)
            {
                return 0;
            }
            if(dot_product(vel_ab, coord_ab) < 0)
            {
                double c_a = find_sound_speed(p1);
                double c_b = find_sound_speed(p2);

                double dens_ab = 1. / 2. * (p1->density + p2->density);
                double c_ab = 1. / 2. * (c_a + c_b);

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

Point Sod_tube_3d::find_new_velocity(Particle * particle, Cell * cell)
{
    double sum_x = 0;
    double sum_y = 0;
    double sum_z = 0;

    double vx = 0;
    double vy = 0;
    double vz = 0;

    Params & params = Params::get_instance();

    assert(particle->density != 0);
    assert(!__isnan(particle->density));

    for (Cell * neighbour : cell->get_neighbours())
    {
        for (Particle & p : particle->kind == Particle::Kind::Gas ? neighbour->gas_particles : neighbour->dust_particles)
        { // TODO remove direct access
            assert(!__isnan(p.density));
            assert(!__isnan(kernel_gradient_x(*(particle), p, params.dimensions)));
            double viscosity = viscosity_3d::find_viscosity(particle, &p);

            sum_x += (particle->pressure / pow(particle->density, 2) + p.pressure / pow(p.density, 2) + viscosity)
                     * kernel_gradient_x(*(particle), p, params.dimensions);

            sum_y += 0;//(particle->pressure / pow(particle->density, 2) + p.pressure / pow(p.density, 2) + viscosity)
                     //* kernel_gradient_y(*(particle), p, params.dimensions);

            sum_z += 0;//(particle->pressure / pow(particle->density, 2) + p.pressure / pow(p.density, 2) + viscosity)
                     //* kernel_gradient_z(*(particle), p, params.dimensions);
        }
    }

    vx = particle->vx - params.tau * particle->mass * sum_x;
    vy = particle->vy - params.tau * particle->mass * sum_y;
    vz = particle->vz - params.tau * particle->mass * sum_z;

    return {vx, vy, vz};
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
            pres_member += vel_dot_grad_kernel(*particle, p);
            visc_member += vel_dot_grad_kernel(*particle, p) * viscosity_3d::find_viscosity(particle, &p);
        }
    }

    result = particle->pressure /  pow(particle->density, 2) * particle->mass * Params::get_instance().tau * pres_member +
             Params::get_instance().tau * particle->mass / 2. * visc_member + particle->energy;
    assert(visc_member >= 0);
    assert(result >= 0);

    return result;
}

/*
TODO This is hard dependence of the parameters
Grid Sod_tube_3d::do_time_step(Grid & old_grid, int step_num)
{
    char filename[512];
    FILE * f;
    if (PRINT_STUFF)
    {
        sprintf(filename, "/home/calat/tmp/part_%0d.dat", step_num);
        f = fopen(filename, "w");
    }

    Params & params = Params::get_instance();

    Grid next_grid(params.grid_step_x, params.grid_step_y, params.grid_step_z, params.border1, params.border2);


    clock_t coord_vel_start = clock();
    //частицы в расчетной области
    for(int i = 10; i < 40; ++i)
    {
        for(int j = 5; j < 10; ++j)
        {
            for(int k = 5; k < 10; ++k)
            {
                Cell & cell = old_grid.cells[i][j][k];

                for(Particle * particle : cell.get_all_particles())
                {
                    if (PRINT_STUFF)
                    {
                        fprintf(f, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", particle->x, particle->y, particle->z,
                                particle->vx, particle->vy, particle->vz, particle->density, particle->pressure,
                                particle->energy);
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

                    Particle particle_zp(*particle);
                    particle_zp.set_coordinates(particle->x, particle->y, particle->z + 0.1);
                    next_grid.with_copy_of(particle_zp);

                    Particle particle_zm(*particle);
                    particle_zm.set_coordinates(particle->x, particle->y, particle->z - 0.1);
                    next_grid.with_copy_of(particle_zm);

                    Particle particle_yp(*particle);
                    particle_yp.set_coordinates(particle->x, particle->y + 0.1, particle->z);
                    next_grid.with_copy_of(particle_yp);

                    Particle particle_ym(*particle);
                    particle_ym.set_coordinates(particle->x, particle->y - 0.1, particle->z);
                    next_grid.with_copy_of(particle_ym);
                }
            }
        }
    }

    clock_t coord_vel_fin = clock();
    std::cout << "coord and vel: " << (double)(coord_vel_fin - coord_vel_start) / CLOCKS_PER_SEC << std::endl;

    clock_t copy_for_dens_start = clock();
    //частицы в расчетной области
    //TODO remove this part?
    for(int i = 10; i < 40; ++i)
    {
        for(int j = 5; j < 10; ++j)
        {
            for(int k = 5; k < 10; ++k)
            {
                Cell & cell = old_grid.cells[i][j][k];

                for(Particle * particle : cell.get_all_particles())
                {
                    Particle particle_zp(*particle);
                    particle_zp.set_coordinates(particle->x, particle->y, particle->z + 0.1);
                    next_grid.with_copy_of(particle_zp);

                    Particle particle_zm(*particle);
                    particle_zm.set_coordinates(particle->x, particle->y, particle->z - 0.1);
                    next_grid.with_copy_of(particle_zm);

                    Particle particle_yp(*particle);
                    particle_yp.set_coordinates(particle->x, particle->y + 0.1, particle->z);
                    next_grid.with_copy_of(particle_yp);

                    Particle particle_ym(*particle);
                    particle_ym.set_coordinates(particle->x, particle->y - 0.1, particle->z);
                    next_grid.with_copy_of(particle_ym);
                }
            }
        }
    }

    //Частицы слева и справа (по х) от области не трогаем
    for(int i = 0; i < 10; ++i)
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
    for(int i = 40; i < old_grid.x_size; ++i)
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
    clock_t copy_for_dens_fin = clock();
    std::cout << "copy for dens: " << (double)(copy_for_dens_fin - copy_for_dens_start) / CLOCKS_PER_SEC << std::endl;

    clock_t dens_start = clock();

    //частицы в расчетной области
    for(int i = 10; i < 40; ++i)
    {
        for(int j = 5; j < 10; ++j)
        {
            for(int k = 5; k < 10; ++k)
            {
                Cell & cell = next_grid.cells[i][j][k];

                for(Particle * particle : cell.get_all_particles())
                {
                    particle->density = find_density(particle, &cell);
                    if (PRINT_STUFF)
                    {
                        std::cout << particle->get_id() << ", my dens: " << particle->density << std::endl;
                    }
                    assert(!__isnan(particle->density));
                    particle->dbg_state = 1;
                    particle->pressure = Sod_tube_1d::find_new_pressure(particle);
                }
            }
        }
    }
    //TODO PRINT DENSITIES
    char before_copy[] = "before_copy";
    print_grid_dens(next_grid, before_copy);

    clock_t dens_fin = clock();
    if (PRINT_STUFF)
    {
        std::cout << "dens: " << (double) (dens_fin - dens_start) / CLOCKS_PER_SEC << std::endl;
    }

    clock_t copy_dens_start = clock();

    for(int i = 10; i < 40; ++i)
    {
        for(int j = 0; j < 5; ++j)
        {
            for(int k = 0; k < 5; ++k)
            {
                Cell & cell = next_grid.cells[i][j][k];

                cell.gas_particles.erase(cell.gas_particles.begin(), cell.gas_particles.end());
            }
        }
    }
    for(int i = 10; i < 40; ++i)
    {
        for(int j = 0; j < 5; ++j)
        {
            for(int k = 10; k < 15; ++k)
            {
                Cell & cell = next_grid.cells[i][j][k];

                cell.gas_particles.erase(cell.gas_particles.begin(), cell.gas_particles.end());
            }
        }
    }
    for(int i = 10; i < 40; ++i)
    {
        for(int j = 10; j < 15; ++j)
        {
            for(int k = 0; k < 5; ++k)
            {
                Cell & cell = next_grid.cells[i][j][k];

                cell.gas_particles.erase(cell.gas_particles.begin(), cell.gas_particles.end());
            }
        }
    }
    for(int i = 10; i < 40; ++i)
    {
        for(int j = 10; j < 15; ++j)
        {
            for(int k = 10; k < 15; ++k)
            {
                Cell & cell = next_grid.cells[i][j][k];

                cell.gas_particles.erase(cell.gas_particles.begin(), cell.gas_particles.end());
            }
        }
    }

    частицы в расчетной области
    for(int i = 10; i < 40; ++i)
    {
        for(int j = 5; j < 10; ++j)
        {
            for(int k = 5; k < 10; ++k)
            {
                Cell & cell = next_grid.cells[i][j][k];

                for(Particle * particle : cell.get_all_particles())
                {
                    Particle particle_zp(*particle);
                    particle_zp.set_coordinates(particle->x, particle->y, particle->z + 0.1);
                    next_grid.with_copy_of(particle_zp);

                    Particle particle_zm(*particle);
                    particle_zm.set_coordinates(particle->x, particle->y, particle->z - 0.1);
                    next_grid.with_copy_of(particle_zm);

                    Particle particle_yp(*particle);
                    particle_yp.set_coordinates(particle->x, particle->y + 0.1, particle->z);
                    next_grid.with_copy_of(particle_yp);

                    Particle particle_ym(*particle);
                    particle_ym.set_coordinates(particle->x, particle->y - 0.1, particle->z);
                    next_grid.with_copy_of(particle_ym);
                }
            }
        }
    }
    clock_t copy_dens_fin = clock();

    std::cout << "copy dens: " << (double)(copy_dens_fin - copy_dens_start) / CLOCKS_PER_SEC << std::endl;

    if (PRINT_STUFF)
    {
        fclose(f);
    }

    //TODO PRINT DENSITIES
    char after_copy[] = "after_copy";
    print_grid_dens(next_grid, after_copy);

    return next_grid;
}
*/

Grid Sod_tube_3d::do_time_step(Grid & old_grid, int step_num)
{
    char filename[512];
    FILE * f;
    if(PRINT_FILE)
    {
        sprintf(filename, "/home/calat/tmp/part_%0d.dat", step_num);
        f = fopen(filename, "w");
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
                        fprintf(f, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", particle->x, particle->y, particle->z,
                                particle->vx, particle->vy, particle->vz, particle->density, particle->pressure,
                                particle->energy);
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

    if(PRINT_FILE)
    {
        fclose(f);
    }

    recalc_density(next_grid);
    Sod_tube_1d::recalc_pressure(next_grid);

    return next_grid;
}