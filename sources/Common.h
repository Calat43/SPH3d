#pragma once

#include "Particle.h"
#include "Point.h"
#include <cmath>

#ifdef _WINDOWS
  #define __isnan std::isnan
  #define uint unsigned int
#endif

const double PI = acos(-1);

class Particle;

double magnitude(double x, double y, double z);
double magnitude(Point p);

double random_double(double from, double to);

double kernel(Particle const part1, Particle const part2, int dimensions);

Point kernel_gradient(Particle const & part1, Particle const & part2, int dimensions);

/*
double kernel_gradient_x(Particle const & part1, Particle const & part2, int dimensions);

double kernel_gradient_y(Particle const & part1, Particle const & part2, int dimensions);

double kernel_gradient_z(Particle const & part1, Particle const & part2, int dimensions);
*/

double distance(Point const & p1, Point const & p2);

double dot_product(Point const point1, Point const point2);

double vel_dot_grad_kernel(Particle const part1, Particle const part2);
