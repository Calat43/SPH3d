#pragma once
#include <vector>
#include <cassert>
#include <cmath>
#include "Particle.h"
#include "Params.h"
#include "Grid.h"
#include "utils.h"

class Grid;
class Particle;

class Cell : public Numbered
{
public:
    static const Cell NO_CELL;

    std::vector<Particle> gas_particles;
    std::vector<Particle> dust_particles;

private:
    Grid * grid;
    int i;
    int j;
    int k;

public:
    Cell(Grid * grid, int i, int j, int k);

    Cell(Cell const & that);

    Cell & operator=(Cell const & that);

    void add_copy_of_particle(Particle const & particle);

    friend std::ostream & operator<<(std::ostream & stream, Cell const & c);

private:
    std::vector<Particle> & particles_of_kind(Particle::Kind kind) {
        switch(kind) {
            case Particle::Kind::Gas: return gas_particles;
            case Particle::Kind::Dust: return dust_particles;
        }
	    assert(false);
    }

public:
    std::vector<Particle *> get_all_particles();

    std::vector<Cell *> get_neighbours();

    /*void for_each_particle(Consumer<Particle *> f) {
        for (Particle * particle : get_all_particles()) {
            f.accept(particle);
        }
    }*/
};
