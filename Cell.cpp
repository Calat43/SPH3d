#include "Cell.h"
#include <iostream>

const Cell Cell::NO_CELL = Cell(nullptr, -1, -1, -1);

Cell::Cell(Grid * grid, int i, int j, int k) : Numbered(), grid(grid), i(i), j(j), k(k) {}

Cell::Cell(Cell const & that) : Numbered(that)
{
    grid = that.grid;
    i = that.i;
    j = that.j;
    k = that.k;
    gas_particles = that.gas_particles;
    dust_particles = that.dust_particles;
}

Cell & Cell::operator=(Cell const & that)
{
    Numbered::operator=(that);
    grid = that.grid;
    i = that.i;
    j = that.j;
    k = that.k;
    gas_particles = that.gas_particles;
    dust_particles = that.dust_particles;
    return *this;
}

void Cell::add_copy_of_particle(Particle const &particle) {
    particles_of_kind(particle.kind).push_back(particle);
}

std::vector<Particle *> Cell::get_all_particles() {
    std::vector<Particle *> all_particles;
    for (Particle & p: gas_particles) {
        all_particles.push_back(&p);
    }
    for (Particle & p: dust_particles) {
        all_particles.push_back(&p);
    }
    return all_particles;
}

std::vector<Cell *> Cell::get_neighbours() {
    int offset_x = (int) ceil(Params::get_instance().smooth_radius / grid->step_x);
    int offset_y = (int) ceil(Params::get_instance().smooth_radius / grid->step_y);
    int offset_z = (int) ceil(Params::get_instance().smooth_radius / grid->step_z);

    std::vector<Cell *> neighbours;

    int form_x = max(0, i - offset_x);
    int to_x = min(grid->x_size, i + offset_x);

    int form_y = max(0, j - offset_x);
    int to_y = min(grid->y_size, j + offset_y);

    int form_z = max(0, k - offset_x);
    int to_z = min(grid->z_size, k + offset_z);

    for (int id_x = form_x; id_x < to_x; ++id_x) {
        for (int id_y = form_y; id_y < to_y; ++id_y) {
            for (int id_z = form_z; id_z < to_z; ++id_z) {
                neighbours.push_back(&(grid->cells[id_x][id_y][id_z]));
            }
        }
    }

    return neighbours;
}

std::ostream & operator<<(std::ostream & stream, Cell const & c)
{
    stream << "Cell#" << c.get_id() << "(" << c.i << ", " << c.j << ", " << c.k << ")";
    return stream;
}
