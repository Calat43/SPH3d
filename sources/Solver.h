#pragma once

class Cell;
class Grid;

Point find_new_coordinates(Particle const & particle);

double find_density(Particle * particle, Cell * cell);

void recalc_density(Grid & grid, Particle::Kind kind);
