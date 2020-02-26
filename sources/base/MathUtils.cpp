#include <cassert>
#include <iostream>

#include "Particle.h"
#include "Point.h"
#include "Cell.h"
#include "Grid.h"
#include "Params.h"
#include "MathUtils.h"


int max(int a, int b)
{
    return a > b ? a : b;
}

int min(int a, int b)
{
    return a < b ? a : b;
}

double magnitude(double x, double y, double z)
{
    return sqrt(x * x + y * y + z * z);
}

double magnitude(Point p)
{
    return sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
}

//dimensions = 1, 2, 3
double find_sigma(int dimensions)
{
    double sigma = 0;
    if(dimensions == 1)
    {
        sigma = 2. / 3.;
    }
    else if(dimensions == 2)
    {
        sigma = 10. / 7. / PI;
    }
    else if(dimensions == 3)
    {
        sigma = 1. / PI;
    }
    else
    {
        assert(false);
    }
    return sigma;
}

//cubic
double kernel(Particle const part1, Particle const part2, int dimensions)
{
    double sigma = find_sigma(dimensions);
    double h = Params::get_instance().h;
    double mag_r = magnitude(part1.x - part2.x, part1.y - part2.y, part1.z - part2.z);
    double q = mag_r / h;
    double result = 0;

    if (q >= 0 && q <= 1)
    {
        result = 1 - 3. / 2. * pow(q, 2) + 3. / 4. * pow(q, 3);
        return sigma / pow(h, dimensions) * result;
    }
    if (q > 1 && q <= 2)
    {
        result = 1. / 4. * pow((2. - q), 3);
        return sigma / pow(h, dimensions) * result;
    }
    if (q > 2)
    {
        return 0;
    }
    assert(false);
    return 0;
}

Point kernel_gradient(Particle const & part1, Particle const & part2, int dimensions)
{
    double sigma = find_sigma(dimensions);
    double h = Params::get_instance().h;
    double mag_r = magnitude(part1.x - part2.x, part1.y - part2.y, part1.z - part2.z);
    double q = mag_r / h;
    double result = 0;

    if(q >= 0 && q <= 1)
    {
        result = - 3. * q + 9. / 4. * q *  q;
    }
    if( q > 1 && q <= 2)
    {
        result = - 3. / 4. * pow((2 - q), 2);
    }

    if(part1.x == part2.x && part1.y == part2.y && part1.z == part2.z)
    {
        return {0, 0, 0};
    }
    else
    {
        double r_x = part1.x - part2.x;
        double r_y = part1.y - part2.y;
        double r_z = part1.z - part2.z;

        double grad_x = sigma / pow(h, dimensions) * r_x / h / mag_r * result;
        double grad_y = sigma / pow(h, dimensions) * r_y / h / mag_r * result;
        double grad_z = sigma / pow(h, dimensions) * r_z / h / mag_r * result;

        return {grad_x, grad_y, grad_z};
    }

}

/*
double kernel_gradient_x(Particle const & part1, Particle const & part2, int dimensions)
{
    double sigma = find_sigma(dimensions);
    double h = Params::get_instance().h;
    double r = part1.x - part2.x;
    double mag_r = magnitude(part1.x - part2.x, part1.y - part2.y, part1.z - part2.z);
    double q = mag_r / h;
    double result = 0;

    if(q >= 0 && q <= 1)
    {
        result = - 3 * q + 9. / 4. * q *  q;
    }
    if( q > 1 && q <= 2)
    {
        result = - 3. / 4. * pow((2 - q), 2);
    }

    if(part1.x == part2.x && part1.y == part2.y && part1.z == part2.z)
    {
        return 0;
    }
    else
    {
        return sigma / pow(h, dimensions) * r / h / mag_r * result;
    }
}

double kernel_gradient_y(Particle const & part1, Particle const & part2, int dimensions)
{
    double sigma = find_sigma(dimensions);

    double h = Params::get_instance().h;
    double r = part1.y - part2.y;
    double mag_r = magnitude(part1.x - part2.x, part1.y - part2.y, part1.z - part2.z);
    double q = mag_r / h;
    double result = 0;

    if(q >= 0 && q <= 1)
    {
        result = - 3 * q + 9. / 4. * q *  q;
    }
    if( q > 1 && q <= 2)
    {
        result = - 3. / 4. * pow((2 - q), 2);
    }

    if(part1.x == part2.x && part1.y == part2.y && part1.z == part2.z)
    {
        return 0;
    }
    else
    {
        return sigma / pow(h, dimensions) * r / h / mag_r * result;
    }
}

double kernel_gradient_z(Particle const & part1, Particle const & part2, int dimensions)
{
    double sigma = find_sigma(dimensions);

    double h = Params::get_instance().h;
    double r = part1.z - part2.z;
    double mag_r = magnitude(part1.x - part2.x, part1.y - part2.y, part1.z - part2.z);
    double q = mag_r / h;
    double result = 0;

    if(q >= 0 && q <= 1)
    {
        result = - 3 * q + 9. / 4. * q *  q;
    }
    if( q > 1 && q <= 2)
    {
        result = - 3. / 4. * pow((2 - q), 2);
    }

    if(part1.x == part2.x && part1.y == part2.y && part1.z == part2.z)
    {
        return 0;
    }
    else
    {
        return sigma / pow(h, dimensions) * r / h / mag_r * result;
    }
}
 */

double random_double(double from, double to)
{
    return from + (rand() / (double) RAND_MAX * (to - from));
}

double distance(Point const & p1, Point const & p2) {
    return sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2) + pow(p1.z - p2.z, 2));
}

double dot_product(Point const point1, Point const point2)
{
    return point1.x * point2.x + point1.y * point2.y + point1.z * point2.z;
}

double vel_dot_grad_kernel(Particle const part1, Particle const part2)
{
    int dimensions = Params::get_instance().dimensions;

    Point kernel_grad = kernel_gradient(part1, part2, dimensions);

    return (part1.vx - part2.vx) * kernel_grad.x +
           (part1.vy - part2.vy) * kernel_grad.y +
           (part1.vz - part2.vz) * kernel_grad.z;
}

Point find_new_coordinates(Particle const & particle)
{
    Params & params = Params::get_instance();
    double x = particle.x + params.tau * particle.vx;
    double y = particle.y + params.tau * particle.vy;
    double z = particle.z + params.tau * particle.vz;
    return {x, y, z};
}


Point MathUtils::find_new_velocity(Particle * particle, Cell * cell)
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
            double viscosity = MathUtils::find_viscosity(particle, &p); // was got from viscosity_3d

            double term1 = (particle->pressure / pow(particle->density, 2) + p.pressure / pow(p.density, 2) + viscosity);

            sum = sum + kernel_grad * term1;
        }
    }

    double term2 = params.tau * particle->mass;

    vel = prev_vel - sum * term2;

    return vel;
}

double MathUtils::find_new_energy(Particle * particle, Cell * cell)
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
            visc_member += vel_dot_gr_k * MathUtils::find_viscosity(particle, &p);
        }
    }

    result = particle->pressure /  pow(particle->density, 2) * particle->mass * Params::get_instance().tau * pres_member +
             Params::get_instance().tau * particle->mass / 2. * visc_member + particle->energy;
    assert(visc_member >= 0);
    assert(result >= 0);

    return result;
}

void MathUtils::recalc_density(Grid & grid, Particle::Kind kind)
{
    for(int i = 0; i < grid.x_size; ++i)
    {
        for(int j = 0; j < grid.y_size; ++j)
        {
            for(int k = 0; k < grid.z_size; ++k)
            {
                Cell & cell = grid.cells[i][j][k];
                for(Particle & particle : cell.particles_of_kind(kind))
                {
                    particle.density = MathUtils::find_density(&particle, &(cell));
                    //particle->density = find_density_no_sort(*particle, state.grid);
                    assert(!__isnan(particle.density));
                }
            }
        }
    }
}

void MathUtils::recalc_pressure(Grid & grid, Particle::Kind kind)
{
    for(int i = 0; i < grid.x_size; ++i)
    {
        for(int j = 0; j < grid.y_size; ++j)
        {
            for(int k = 0; k < grid.z_size; ++k)
            {
                Cell & cell = grid.cells[i][j][k];
                for(Particle & particle : cell.particles_of_kind(kind))
                {
                    particle.pressure = MathUtils::find_new_pressure(&particle); // was got from Sod_tube_1d
                    //particle->density = find_density_no_sort(*particle, state.grid);
                    assert(!__isnan(particle.pressure));
                }
            }
        }
    }
}


double MathUtils::find_viscosity(Particle * p1, Particle * p2)
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
            double c_a = MathUtils::find_sound_speed(p1); // was got from viscosity
            double c_b = MathUtils::find_sound_speed(p2); // was got from viscosity

            double dens_ab = 1. / 2. * (p1->density + p2->density);
            double c_ab = 1. / 2. * (c_a + c_b);

            double mu_ab = MathUtils::find_mu(vel_ab, coord_ab); // was got from viscosity

            assert(!__isnan(c_a));
            assert(!__isnan(c_b));
            assert(!__isnan(mu_ab));

            return (- params.alpha * c_ab * mu_ab + params.beta * pow(mu_ab, 2)) / dens_ab;
        }
    }
    return 0;
}

double MathUtils::find_density(Particle * particle, Cell * cell)
{
    double density = 0;
    int nears = 0;

    //Params & params = Params::get_instance();

    for (Cell * neighbour : cell->get_neighbours())
    {
        // std::cout << *neighbour << std::endl;
        for (Particle & p : particle->kind == Particle::Kind::Gas ? neighbour->gas_particles : neighbour->dust_particles)
        { // TODO remove direct access
            density += p.mass * kernel(*particle, p, Params::get_instance().dimensions);

            if(PRINT_STUFF)
            {
                if(kernel(*particle, p, Params::get_instance().dimensions) != 0)
                {
                    nears += 1;
                }
            }
        }
    }

    if(PRINT_STUFF)
    {
        std::cout << particle->get_id() << " " << nears << std::endl;
    }
    assert(!__isnan(density));

    return density;
}

double MathUtils::find_new_pressure(Particle * particle)
{
    return particle->density * particle->energy * (Params::get_instance().gamma - 1);
}

double MathUtils::find_sound_speed(Particle * particle)
{
    assert(!__isnan(particle->pressure));
    assert(!__isnan(particle->density));
    assert(!__isnan(sqrt(Params::get_instance().gamma * particle->pressure / particle->density)));
    return sqrt(Params::get_instance().gamma * particle->pressure / particle->density);
}

double MathUtils::find_mu(Point vel_ab, Point coord_ab)
{
    return (Params::get_instance().h * dot_product(vel_ab, coord_ab)) /
           (pow(magnitude(coord_ab), 2) + pow(Params::get_instance().nu, 2));
}

double MathUtils::find_mu(double vel_ab, double coord_ab)
{
    return (Params::get_instance().h * vel_ab * coord_ab) / (pow(coord_ab, 2) + pow(Params::get_instance().nu, 2));
}

