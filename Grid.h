#pragma once

#include <vector>
#include <cmath>
#include <cassert>
#include <functional>
#include "Cell.h"
#include "Point.h"


class Cell;

class Grid {

public:
    std::vector<std::vector<std::vector<Cell>>> cells;


    Point border1;
    Point border2;

    double step_x;
    double step_y;
    double step_z;

    int x_size;
    int y_size;
    int z_size;

    Grid(double step_x, double step_y, double step_z, Point const & border1, Point const & border2);

    Grid(Grid const & that);

    Grid & operator=(Grid const & that);

    Cell * find_cell(const Particle &particle);

    void for_each_cell(std::function<void(Cell &)> const & f);

    void for_each_cell(std::function<void(Cell *)> const & f);

private:
    bool encloses_point(Point point);
};
