#pragma once
#include <vector>
#include <cassert>

#include "Numbered.h"
#include "Particle.h"

class Grid;


class Cell : public Numbered
{
    friend Grid;
public:
    static const Cell NO_CELL;

    std::vector<Particle> gas_particles;
    std::vector<Particle> dust_particles;

    Cell(Grid * grid, int i, int j, int k);

    Cell(Cell const & that);

    Cell & operator=(Cell const & that);

    void add_copy_of_particle(Particle const & particle);

    void remove_particle(Particle & particle);

    friend std::ostream & operator<<(std::ostream & stream, Cell const & c);

    std::vector<Particle *> get_all_particles();

    std::vector<Cell *> get_neighbours();

    /*void for_each_particle(Consumer<Particle *> f) {
        for (Particle * particle : get_all_particles()) {
            f.accept(particle);
        }
    }*/

    int get_i() const;

    int get_j() const;

    int get_k() const;

    int get_drag_mode() const;

    void set_drag_mode(int num);

    std::vector<Particle> & particles_of_kind(Particle::Kind kind) {
        switch(kind) {
            case Particle::Kind::Gas: return gas_particles;
            case Particle::Kind::Dust: return dust_particles;
        }
        assert(false);
        return gas_particles; // just to avoid warning
    }

    void clear_arrays();


private:
    Grid * grid;
    int i;
    int j;
    int k;
    int drag_mode;
};
