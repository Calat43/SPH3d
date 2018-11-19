#include "Common.h"
#include "Params.h"
#include "Numbered.h"
#include <cassert>


double magnitude(double x, double y, double z)
{
    return sqrt(x * x + y * y + z * z);
}

//dimensions = 1, 2, 3
double find_sigma(int dimensions)
{
    double sigma;
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
double kernel(Particle part1, Particle part2, int dimensions)
{
    double sigma = find_sigma(dimensions);
    double h = Params::get_instance().h;
    double r = magnitude(part1.x - part2.x, part1.y - part2.y, part1.z - part2.z);
    double q = r / h;
    double result = 0;

    if (q >= 0 && q <= 1)
    {
        result = 1 - 3. / 2 * pow(q, 2) + 3. / 4 * pow(q, 3);
        return sigma / pow(h, dimensions) * result;
    }
    if (q > 1 && q <= 2)
    {
        result = 1. / 4 * pow((2. - q), 3);
        return sigma / pow(h, dimensions) * result;
    }
    if (q > 2)
    {
        return 0;
    }
    assert(false);
}

double kernel_gradient_x(Particle part1, Particle part2, int dimensions)
{
    double sigma = find_sigma(dimensions);

    double h = Params::get_instance().h;
    double r = fabs(part1.x - part2.x);
    double q = r / h;
    double result = 0;

    if(q >= 0 && q <= 1)
    {
        result = - 3 * q + 9. / 4. * q *  q;
    }
    if( q > 1 && q <= 2)
    {
        result = - 3. / 4. * pow((2 - q), 2);
    }
    if (part1.x > part2.x)
    {
        return sigma / h / pow(h, dimensions) * result;
    }
    if (part1.x == part2.x)
    {
        return 0;
    }
    if (part1.x < part2.x)
    {
        return - sigma / h / pow(h, dimensions) * result;
    }

    assert(false);
}

double kernel_gradient_y(Particle part1, Particle part2, int dimensions)
{
    double sigma = find_sigma(dimensions);

    double h = Params::get_instance().h;
    double r = fabs(part1.y - part2.y);
    double q = r / h;
    double result = 0;

    if(q >= 0 && q <= 1)
    {
        result = - 3 * q + 9. / 4. * q *  q;
    }
    if( q > 1 && q <= 2)
    {
        result = - 3. / 4. * pow((2 - q), 2);
    }
    if (part1.y > part2.y)
    {
        return sigma / h / pow(h, dimensions) * result;
    }
    if (part1.y == part2.y)
    {
        return 0;
    }
    if (part1.y < part2.y)
    {
        return - sigma / h / pow(h, dimensions) * result;
    }

    assert(false);
}

double kernel_gradient_z(Particle part1, Particle part2, int dimensions)
{
    double sigma = find_sigma(dimensions);

    double h = Params::get_instance().h;
    double r = fabs(part1.z - part2.z);
    double q = r / h;
    double result = 0;

    if(q >= 0 && q <= 1)
    {
        result = - 3 * q + 9. / 4. * q *  q;
    }
    if( q > 1 && q <= 2)
    {
        result = - 3. / 4. * pow((2 - q), 2);
    }
    if (part1.z > part2.z)
    {
        return sigma / h / pow(h, dimensions) * result;
    }
    if (part1.z == part2.z)
    {
        return 0;
    }
    if (part1.z < part2.z)
    {
        return - sigma / h / pow(h, dimensions) * result;
    }

    assert(false);
}

double random_double(double from, double to)
{
    return from + (rand() / (double) RAND_MAX * (to - from));
}

