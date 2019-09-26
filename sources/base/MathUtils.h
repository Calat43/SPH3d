#pragma once

#ifdef _WINDOWS
  #define __isnan std::isnan
  #define uint unsigned int
#endif

#include <cmath>

#include "Particle.h"

class Point;
class Grid;
class Cell;

const double PI = acos(-1);

int max(int a, int b);

int min(int a, int b);

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

Point find_new_coordinates(Particle const & particle);

namespace MathUtils {

    Point find_new_velocity(Particle * particle, Cell * cell);

    double find_new_energy(Particle * particle, Cell * cell);
    
    void recalc_density(Grid & grid, Particle::Kind kind);

    void recalc_pressure(Grid & grid, Particle::Kind kind);

    double find_viscosity(Particle * p1, Particle * p2);

    double find_density(Particle * particle, Cell * cell);

    double find_new_pressure(Particle * particle);

    double find_sound_speed(Particle * particle);
    double find_mu(Point vel_ab, Point coord_ab);
    double find_mu(double vel_ab, double coord_ab);
}



