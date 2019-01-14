#include <iostream>
#include "Grid.h"

Grid::Grid(double step_x, double step_y, double step_z, Point const & border1, Point const & border2)
    : cells(),
    border1(border1), border2(border2), step_x(step_x), step_y(step_y), step_z(step_z),
    x_size((int) ceil((border2.x - border1.x) / step_x)), // TODO maybe support last cell of different size?
    y_size((int) ceil((border2.y - border1.y) / step_y)),
    z_size((int) ceil((border2.z - border1.z) / step_z))
{
    assert(this->border1.x <= this->border2.x);
    assert(this->border1.y <= this->border2.y);
    assert(this->border1.z <= this->border2.z);

    for (int i = 0; i < x_size; ++i) {
        cells.emplace_back();
        for (int j = 0; j < y_size; ++j) {
            cells[i].emplace_back();
            for (int k = 0; k < z_size; ++k) {
                cells[i][j].push_back(Cell(this, i, j, k));
            }
        }
    }
}

Cell * Grid::find_cell(const Particle & particle) {
    if (!encloses_point(Point(particle.x, particle.y, particle.z))) {
        return nullptr; // TODO &(Cell::NO_CELL);
    }

    double x_relative = particle.x - border1.x;
    double y_relative = particle.y - border1.y;
    double z_relative = particle.z - border1.z;

    int i = (int) floor(x_relative / step_x);
    int j = (int) floor(y_relative / step_y);
    int k = (int) floor(z_relative / step_z);

    return & (cells[i][j][k]);
}

bool Grid::encloses_point(Point point) {
    if (point.x < border1.x || border2.x <= point.x) return false;
    if (point.y < border1.y || border2.y <= point.y) return false;
    if (point.z < border1.z || border2.z <= point.z) return false;
    return true;
}

Grid::Grid(Grid const & that) : Grid(that.step_x, that.step_y, that.step_z, that.border1, that.border2) {
    cells = that.cells;
    for_each_cell([this] (Cell & cell)
    {
        cell.grid = this;
    });
    x_size = that.x_size;
    y_size = that.y_size;
    z_size = that.z_size;
}

Grid & Grid::operator=(Grid const & that) {
    step_x = that.step_x;
    step_y = that.step_y;
    step_z = that.step_z;
    x_size = that.x_size;
    y_size = that.y_size;
    z_size = that.z_size;
    border1 = that.border1;
    border2 = that.border2;
    cells = that.cells;
    for_each_cell([this] (Cell & cell)
    {
        cell.grid = this;
    });
    return *this;
}

void Grid::for_each_cell(std::function<void(Cell &)> const & f)
{
    for (int i = 0; i < x_size; ++i) {
        for (int j = 0; j < y_size; ++j) {
            for (int k = 0; k < z_size; ++k) {
                Cell & cell = cells[i][j][k];
                f(cell);
            }
        }
    }
}

void Grid::for_each_cell(std::function<void(Cell *)> const & f)
{
    for (int i = 0; i < x_size; ++i) {
        for (int j = 0; j < y_size; ++j) {
            for (int k = 0; k < z_size; ++k) {
                Cell & cell = cells[i][j][k];
                f(&cell);
            }
        }
    }
}
