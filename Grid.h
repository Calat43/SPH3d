#pragma once

#include <vector>
#include <cmath>
#include <cassert>
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

private:
    bool encloses_point(Point point);

public:
    /*void for_each_cell(Consumer<Cell *> f) {
        for (int i = 0; i < x_size; ++i) {
            for (int j = 0; j < y_size; ++j) {
                for (int k = 0; k < z_size; ++k) {
                    Cell cell = cells[i][j][k];
                    f.accept(&(cell));
                }
            }
        }
    }*/

};
