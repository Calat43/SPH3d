#include <fstream>
#include <iostream>
#include <string>

#include "Params.h"
#include "Cell.h"
#include "Grid.h"

#include "SPHSolver.h"
#include "BallInVacuum.h"
#include "SodTube1d.h"
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

    std::cout << "c0: " << c0 << std::endl;
    std::cout << "c1: " << c1 << std::endl;
    std::cout << "b: " << b << std::endl;

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

    f_solv.open(OUTPUT_PATH + "ball_solution_" + std::to_string(exact_time) + "_"
                + std::to_string(params.ball_init_radius) + ".dat");

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

    double result = - params.tau * params.gamma * particle->mass * particle->pressure / particle->density * sum + particle->pressure;

    //assert(result >= 0);
    return result;
}

Point ball_in_vacuum::init_velocity(Point center, Point particle_coordinates)
{
    double radius = distance(particle_coordinates, center);

    double velocity = Params::get_instance().ball_init_velocity / Params::get_instance().ball_init_radius * radius;

    Point tmp_vel = particle_coordinates - center;

    return tmp_vel / magnitude(tmp_vel) * velocity;
}

Grid ball_in_vacuum::init(double step)
{
    Params & params = Params::get_instance();

    double radius = params.ball_init_radius;

    Grid init(params.grid_step_x, params.grid_step_y, params.grid_step_z, params.border1, params.border2);

    std::ofstream part_file;

    if(PRINT_FILE)
    {
        part_file.open(OUTPUT_PATH + "ball_init" + ".dat");
    }

    double center_x = (params.border2.x + params.border1.x) / 2;
    double center_y = (params.border2.y + params.border1.y) / 2;
    double center_z = (params.border2.z + params.border1.z) / 2;

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
                    particle->mass = params.ball_mass / count;
                }
            }
        }
    }
    std::cout << "------particles:" << count << std::endl;

    recalc_density(init, Particle::Kind::Gas);


    //костыль из аналитики
    double ddot_R = 2. * params.A * pow(3. * params.ball_mass / 4. / PI, 1. / 3.) / params.ball_init_radius
                    / params.ball_init_radius;

    for (int i = 0; i < init.x_size; ++i)
    { // TODO foreach
        for(int j = 0; j < init.y_size; ++j)
        {
            for(int k = 0; k < init.z_size; ++k)
            {
                Cell & cell = init.cells[i][j][k];

                for(Particle * particle : cell.get_all_particles())
                {
                    double r = distance(Point(particle->x, particle->y, particle->z), center);
                    particle->pressure = ball_analytic::pressure(r, params.ball_init_radius, ddot_R);
                }
            }
        }
    }

    for(int i = 0; i < init.x_size; ++i)
    { // TODO foreach
        for(int j = 0; j < init.y_size; ++j)
        {
            for(int k = 0; k < init.z_size; ++k)
            {
                Cell & cell = init.cells[i][j][k];

                for(Particle * particle : cell.get_all_particles())
                {
                    if(PRINT_FILE)
                    {
                        Point particle_point(particle->x, particle->y, particle->z);

                        double r = distance(particle_point, center);

                        part_file << r << " " << magnitude(particle->vx, particle->vy, particle->vz) <<
                                   " " << particle->pressure << " " << particle->density << std::endl;
                    }
                }
            }
        }
    }

    return init;
}

//analytic velocity with p squared, equation only for velocity
double vel_p_squared(double t, double r)
{
    Params & params = Params::get_instance();

    double init_density = 3. * params.ball_mass / 4. / PI / pow(params.ball_init_radius, 3);
    double beta = sqrt(200. / init_density);
    double alpha = 1.;

    double coeff = (beta - alpha) / (beta + alpha);

    double power = - 2. * beta * t;
    double expon = exp(- 2. * beta * t);

    return //beta * r * (1 - coeff * expon) / sqrt(beta * expon);
    //Euler velocity
    beta * r * (1. - coeff * expon) / (1. + coeff * expon);
}

void print_vel_p_squared(double t, double r_max, double r_step) {
    std::ofstream file;
    file.open(OUTPUT_PATH + "vel_p_squared" + ".dat");
    for (double r = 0; r < r_max; r += r_step) {
        file << r << " " << vel_p_squared(t, r) << std::endl;
    }
    file.close();
}

Grid ball_in_vacuum::do_time_step(Grid & old_grid, int step_num)
{
    std::ofstream part_file;
    if(PRINT_FILE)
    {
        part_file.open(OUTPUT_PATH + "part_" + std::to_string(step_num) + ".dat");
    }

    //std::ofstream diff_part;

    //diff_part.open(OUTPUT_PATH + "diff_" + std::to_string(step_num) + ".dat");

    Params & params = Params::get_instance();

    Grid next_grid(params.grid_step_x, params.grid_step_y, params.grid_step_z, params.border1, params.border2);

    double center_x = (params.border2.x + params.border1.x) / 2;
    double center_y = (params.border2.y + params.border1.y) / 2;
    double center_z = (params.border2.z + params.border1.z) / 2;

    Point center(center_x, center_y, center_z);


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
                    Point new_coords = find_new_coordinates_(*particle);

                    Point new_vel = Sod_tube_3d::find_new_velocity(particle, &cell, step_num);
                    new_particle.density = NAN;
                    new_particle.dbg_state = 2;

                    new_particle.set_coordinates(new_coords.x, new_coords.y, new_coords.z);
                    //new_particle.set_velocities(particle->vx, particle->vy, particle->vz);
                    new_particle.set_velocities(new_vel.x, new_vel.y, new_vel.z);

                    //double r = distance(new_coords, center);

                    //new_particle.pressure = 1. - 100. * r * r;
                    //particle->pressure;//ball_in_vacuum::find_pressure(particle, &cell);
                    //new_particle.density = particle->density;

                    next_grid.with_copy_of(new_particle);
                }
            }
        }
    }


    recalc_density(next_grid, Particle::Kind::Gas);

    for(int i = 0; i < next_grid.x_size; ++i)
    { // TODO foreach
        for (int j = 0; j < next_grid.y_size; ++j)
        {
            for (int k = 0; k < next_grid.z_size; ++k)
            {
                Cell &cell = next_grid.cells[i][j][k];

                for (Particle *particle : cell.get_all_particles())
                {
                    particle->pressure = pow(params.gas_sound_speed, 2) * particle->density;
                }
            }
        }
    }

    for(int i = 0; i < next_grid.x_size; ++i)
    { // TODO foreach
        for(int j = 0; j < next_grid.y_size; ++j)
        {
            for(int k = 0; k < next_grid.z_size; ++k)
            {
                Cell & cell = next_grid.cells[i][j][k];

                for(Particle * particle : cell.get_all_particles())
                {
                    if(PRINT_FILE)
                    {
                        Point particle_point(particle->x, particle->y, particle->z);

                        double radius = distance(particle_point, center);
                        double v_magn = magnitude(particle->vx, particle->vy, particle->vz);

                        part_file << radius << " " << particle->x << " " << particle->y << " " << particle->z << " "
                                  << v_magn << " " << particle->vx<< " " << particle->vy << " " << particle->vz
                                  << " " << particle->pressure << " " << particle->density << std::endl;
                    }
                }
            }
        }
    }

    return next_grid;
}
