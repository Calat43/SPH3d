
#include <fstream>
#include <iostream>
#include <string>

#include "Params.h"
#include "Cell.h"
#include "Grid.h"

#include "SodTube1d.h"
#include "BallInVacuum.h"
#include "SodTube3d.h"

double ball_analytic::time_through_R(double R)
{
    Params & params = Params::get_instance();

    //TODO something
    if(params.gamma != 4. / 3.)
    {
        return -1000;
    }

    double r0 = params.ball_init_radius;

    double c0 = 4. * params.A * pow(3. * params.ball_mass / 4. / PI, 1./3.);

    double b = c0 / r0 + params.ball_init_velocity * params.ball_init_velocity;

    double c1 = c0 * log(2. * sqrt(b) * r0 * sqrt(b - c0 / r0) + 2. * b * r0 - c0) / 2. / pow(b, 3. / 2.) +
                r0 * sqrt(b - c0 / r0) / b;

    double result = c0 * log(2. * sqrt(b) * R * sqrt(b - c0 / R) + 2. * b * R - c0) / 2. / pow(b, 3. / 2.) +
                    R * sqrt(b - c0 / R) / b - c1;

    return result;
}

double ball_analytic::R_bisection(double exact_time, double step, double defect)
{
    double r0 = Params::get_instance().ball_init_radius;
    double r = r0;

    double time = ball_analytic::time_through_R(r0);

    while(fabs(time - exact_time) > defect)
    {
        r += step;
        time = ball_analytic::time_through_R(r);
    }

    return r;
}

double ball_analytic::density(double exact_time, double step, double defect)
{
    Params & params = Params::get_instance();

    double R = R_bisection(exact_time, step, defect);

    return 3. * params.ball_mass / 4. / PI / pow(R, 3);
}

double ball_analytic::velocity(double r, double f)
{
    return r * f;
}

double ball_analytic::pressure(double r, double R, double ddot_R)
{
    return 3. * Params::get_instance().ball_mass / 4. / PI * ddot_R / 2. / R / R * (1. - r * r / R / R);
}

void ball_analytic::print_solution(double exact_time, double step_r, double step_R, double defect)
{
    Params & params = Params::get_instance();

    double density = ball_analytic::density(exact_time, step_R, defect);

    double R = ball_analytic::R_bisection(exact_time, step_R, defect);

    double dot_R = sqrt(4. * params.A * pow(3. * params.ball_mass / 4. / PI, 1. / 3.) *
                   (1. / params.ball_init_radius - 1. / R) + pow(params.ball_init_velocity, 2));

    double ddot_R = 2. * params.A * pow(3. * params.ball_mass / 4. / PI, 1. / 3.) / R / R;

    double f = dot_R / R;

    std::ofstream f_solv;

    f_solv.open(OUTPUT_PATH + "ball_solution_" + std::to_string(exact_time) + ".dat");

    for(double radius = 0; radius <= R; radius += step_r)
    {
        f_solv << radius << " " << ball_analytic::velocity(radius, f) << " "
               << ball_analytic::pressure(radius, R, ddot_R) << " " << density
               << std::endl;
    }

}

double ball_in_vacuum::find_pressure(Particle * particle, Cell * cell)
{
    Params & params = Params::get_instance();

    double sum = 0;

    for (Cell * neighbour : cell->get_neighbours())
    {
        for (Particle & p : particle->kind == Particle::Kind::Gas ? neighbour->gas_particles : neighbour->dust_particles)
        { // TODO remove direct access//
            assert(!__isnan(p.density));

            Point kernel_grad = kernel_gradient(*(particle), p, params.dimensions);
           sum +=  (p.vx - particle->vx) * kernel_grad.x +
                   (p.vy - particle->vy) * kernel_grad.y +
                   (p.vz - particle->vz) * kernel_grad.z;
        }
    }

    return - params.tau * params.gamma * particle->pressure / particle->density * sum + particle->pressure;
}

Point ball_in_vacuum::init_velocity(Point center, Point particle_coordinates)
{
    Point vel = particle_coordinates - center;

    double magn = magnitude(vel);

    return vel / magn * Params::get_instance().ball_init_velocity;
}

Grid ball_in_vacuum::init(double step)
{
    Params & params = Params::get_instance();

    double radius = params.ball_init_radius;

    Grid init(params.grid_step_x, params.grid_step_y, params.grid_step_z, params.border1, params.border2);

    double center_x = (params.border2.x - params.border1.x) / 2;
    double center_y = (params.border2.y - params.border1.y) / 2;
    double center_z = (params.border2.z - params.border1.z) / 2;

    Point center(center_x, center_y, center_z);

    int count = 0;

    for(double x = params.border1.x; x < params.border2.x; x += step)
    {
        for(double y = params.border1.y; y < params.border2.y; y += step)
        {
            for(double z = params.border1.z; z < params.border2.z; z += step)
            {
                if (distance(Point(x, y, z), center) <= radius)
                {
                    Particle particle(Particle::Kind::Gas, x, y, z);

                    particle.vx = ball_in_vacuum::init_velocity(center, Point(particle.x, particle.y, particle.z)).x;
                    particle.vy = ball_in_vacuum::init_velocity(center, Point(particle.x, particle.y, particle.z)).y;
                    particle.vz = ball_in_vacuum::init_velocity(center, Point(particle.x, particle.y, particle.z)).z;

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

    for (int i = 0; i < init.x_size; ++i)
    { // TODO foreach
        for(int j = 0; j < init.y_size; ++j)
        {
            for(int k = 0; k < init.z_size; ++k)
            {
                Cell & cell = init.cells[i][j][k];

                for(Particle * particle : cell.get_all_particles())
                {
                    particle->pressure = params.A * pow(particle->density, params.gamma);
                }
            }
        }
    }

    return init;
}

Grid ball_in_vacuum::do_time_step(Grid & old_grid, int step_num)
{
    std::ofstream part_file;
    if(PRINT_FILE)
    {
        part_file.open(OUTPUT_PATH + "part_" + std::to_string(step_num) + ".dat");
    }

    Params & params = Params::get_instance();

    Grid next_grid(params.grid_step_x, params.grid_step_y, params.grid_step_z, params.border1, params.border2);

    for(int i = 0; i < old_grid.x_size; ++i)
    { // TODO foreach
        for(int j = 0; j < old_grid.y_size; ++j)
        {
            for(int k = 0; k < old_grid.z_size; ++k)
            {
                Cell & cell = old_grid.cells[i][j][k];

                for(Particle * particle : cell.get_all_particles())
                {
                    Particle new_particle(*particle);
                    Point new_coords = find_new_coordinates(*particle);

                    Point new_vel = Sod_tube_3d::find_new_velocity(particle, &cell);
                    new_particle.density = NAN;
                    new_particle.dbg_state = 2;

                    new_particle.set_coordinates(new_coords.x, new_coords.y, new_coords.z);
                    new_particle.set_velocities(new_vel.x, new_vel.y, new_vel.z);

                    new_particle.pressure = ball_in_vacuum::find_pressure(particle, &cell);

                    next_grid.with_copy_of(new_particle);
                }
            }
        }
    }

    recalc_density(next_grid, Particle::Kind::Gas);

    for(int i = 0; i < old_grid.x_size; ++i)
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
                                  << particle->density << " " << particle->pressure << std::endl;
                    }
                }
            }
        }
    }

    return next_grid;
}
