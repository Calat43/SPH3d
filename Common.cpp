#include "Common.h"
#include "Params.h"
#include "Numbered.h"
#include <cassert>
#include <iostream>


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

double max0(double value)
{
    return (value > 0) ? value : 0;
}

double cubic_3dkernel(Particle const part1, Particle const part2)
{
    double result = 0;
    double h = Params::get_instance().h;
    double r = magnitude(part1.x - part2.x, part1.y - part2.y, part1.z - part2.z);
    double q = r / h;

    result = 1. / h / h / h * 16. / PI * (pow(max0(1 - q), 3) - 4 * pow(max0(0.5 - q), 3));

    return result;
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
}

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
    /*
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
     */

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

    /*
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
    */

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

    /*
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
    */
    if(part1.x == part2.x && part1.y == part2.y && part1.z == part2.z)
    {
        return 0;
    }
    else
    {
        return sigma / pow(h, dimensions) * r / h / mag_r * result;
    }
}

double random_double(double from, double to)
{
    return from + (rand() / (double) RAND_MAX * (to - from));
}

double distance(Point const & p1, Point const & p2) {
    return sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2) + pow(p1.z - p2.z, 2));
}

