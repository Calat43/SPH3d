#include "ParticlesState.h"

ParticlesState::ParticlesState() : grid(
        Params::get_instance().grid_step_x,
        Params::get_instance().grid_step_y,
        Params::get_instance().grid_step_z,
        Params::get_instance().border1,
        Params::get_instance().border2
) {}

ParticlesState::ParticlesState(ParticlesState const & that) : grid(that.grid) {}

ParticlesState & ParticlesState::operator=(ParticlesState const & that) {
    grid = that.grid;
    return *this;
}

void ParticlesState::with_copy_of(Particle const & particle) {
    Cell * cell = grid.find_cell(particle);
    if (cell == nullptr)
    {
        return;
    }
    cell->add_copy_of_particle(particle);
}
