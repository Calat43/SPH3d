#include <iostream>

#include "Params.h"
#include "MathUtils.h"
#include "Grid.h"
#include "Cell.h"


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

    int from_x = max(0, i - offset_x);
    int to_x = min(grid->x_size - 1, i + offset_x);

    int from_y = max(0, j - offset_x);
    int to_y = min(grid->y_size - 1, j + offset_y);

    int from_z = max(0, k - offset_x);
    int to_z = min(grid->z_size - 1, k + offset_z);

    for (int id_x = from_x; id_x <= to_x; ++id_x) {
        for (int id_y = from_y; id_y <= to_y; ++id_y) {
            for (int id_z = from_z; id_z <= to_z; ++id_z) {
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

void Cell::remove_particle(Particle & /*particle*/)
{
    assert(false);
}

int Cell::get_i() const
{
    return i;
}

int Cell::get_j() const
{
    return j;
}

int Cell::get_k() const
{
    return k;
}

int Cell::get_drag_mode() const
{
    return drag_mode;
}

void Cell::set_drag_mode(int num)
{
    this->drag_mode = num;
}


void Cell::clear_arrays()
{
    gas_particles.clear();
    dust_particles.clear();
}
