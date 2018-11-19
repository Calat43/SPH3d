#pragma once

#include "Particle.h"
#include <cmath>

const double PI = acos(-1);

class Particle;

double random_double(double from, double to);

double kernel(Particle part1, Particle part2, int dimensions);

double kernel_gradient_x(Particle part1, Particle part2, int dimensions);

double kernel_gradient_y(Particle part1, Particle part2, int dimensions);

double kernel_gradient_z(Particle part1, Particle part2, int dimensions);

