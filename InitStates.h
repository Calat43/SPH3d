#pragma once

#include "Particle.h"
#include "Grid.h"
#include "Common.h"
#include "Solver.h"

Grid squared_ball_init_state(double radius, double step);

Grid sod_tube();
