#pragma once

#include <iostream>

#include "Params.h"
#include "Common.h"
#include "Cell.h"
#include "ParticlesState.h"

Point find_new_coordinates(Particle const & particle);

double find_density_no_sort(Particle const & particle, Grid const & grid);

void recalc_density(ParticlesState & state, double radius);

Point find_new_velocity(Particle * particle, Cell * cell);

Point find_new_velocity_no_sort(Particle * particle, Grid * grid);