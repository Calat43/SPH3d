#pragma once

#include "Particle.h"
#include "Point.h"
#include <cmath>

const double PI = acos(-1);

const bool PRINT_DENSITY = false;

const bool PRINT_STUFF = false;

const bool PRINT_FILE = false;

class Particle;

double magnitude(double x, double y, double z);
double magnitude(Point p);

double random_double(double from, double to);

double kernel(Particle const part1, Particle const part2, int dimensions);

double kernel_gradient_x(Particle const & part1, Particle const & part2, int dimensions);

double kernel_gradient_y(Particle const & part1, Particle const & part2, int dimensions);

double kernel_gradient_z(Particle const & part1, Particle const & part2, int dimensions);

double distance(Point const & p1, Point const & p2);

double dot_product(Point const point1, Point const point2);

double vel_dot_grad_kernel(Particle const part1, Particle const part2);
