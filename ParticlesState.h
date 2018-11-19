#pragma once

#include "Grid.h"
#include "Params.h"

class ParticlesState {
public:
    Grid grid;

    ParticlesState();

    ParticlesState(ParticlesState const & that);

    ParticlesState & operator=(ParticlesState const & that);

    void with_copy_of(Particle const & particle);

    /*void print(PrintWriter writer) {
        grid.for_each_cell(
                (c) -> c->for_each_particle(
                        (p) -> writer.println(p->x + " " + p->y + " " + p->z)
                )
        );
    }*/

};
