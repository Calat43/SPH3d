#pragma once

#include "Particle.h"
#include "Point.h"
#include <cmath>

const double PI = acos(-1);

const bool PRINT_DENSITY = false;

class Particle;

double random_double(double from, double to);

double kernel(Particle const part1, Particle const part2, int dimensions);

double kernel_gradient_x(Particle const & part1, Particle const & part2, int dimensions);

double kernel_gradient_y(Particle const & part1, Particle const & part2, int dimensions);

double kernel_gradient_z(Particle const & part1, Particle const & part2, int dimensions);

double distance(Point const & p1, Point const & p2);
