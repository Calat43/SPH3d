#include <cassert>

#include "Particle.h"
#include "Point.h"
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
