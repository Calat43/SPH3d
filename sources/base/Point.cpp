
#include <cassert>

#include "Point.h"

Point Point::operator +(const Point & p)
{
    return {this->x + p.x, this->y + p.y, this->z + p.z};
}

Point Point::operator -(const Point & p)
{
    return {this->x - p.x, this->y - p.y, this->z - p.z};
}

Point Point::operator *(double term)
{
    return {this->x * term, this->y * term, this->z * term};
}

Point Point::operator /(double term)
{
    assert(term != 0);
    return {this->x / term, this->y / term, this->z / term};
}

Point & Point::operator =(const Point & p)
{
    this->x = p.x;
    this->y = p.y;
    this->z = p.z;

    return *this;
}

Point::Point(const Point & p) : x(p.x), y(p.y), z(p.z)
{}

Point::Point(double x, double y, double z) : x(x), y(y), z(z) {}
