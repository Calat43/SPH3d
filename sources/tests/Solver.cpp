#include <fstream>
#include <iostream>

#include "Params.h"

#include "Grid.h"
#include "Cell.h"
#include "Viscosity.h"

#include "Solver.h"

Point find_new_coordinates(Particle const & particle)
{
    Params & params = Params::get_instance();
    double x = particle.x + params.tau * particle.vx;
    double y = particle.y + params.tau * particle.vy;
    double z = particle.z + params.tau * particle.vz;
    return {x, y, z};
}

double find_density(Particle * particle, Cell * cell)
{
    double density = 0;
    int nears = 0;

    //Params & params = Params::get_instance();

    for (Cell * neighbour : cell->get_neighbours())
    {
        // std::cout << *neighbour << std::endl;
        for (Particle & p : particle->kind == Particle::Kind::Gas ? neighbour->gas_particles : neighbour->dust_particles)
        { // TODO remove direct access
            density += p.mass * kernel(*particle, p, Params::get_instance().dimensions);

            if(PRINT_STUFF)
            {
                if(kernel(*particle, p, Params::get_instance().dimensions) != 0)
                {
                    nears += 1;
                }
            }
        }
    }

    if(PRINT_STUFF)
    {
        std::cout << particle->get_id() << " " << nears << std::endl;
    }
    assert(!__isnan(density));

    return density;
}

void recalc_density(Grid & grid, Particle::Kind kind)
{
    for(int i = 0; i < grid.x_size; ++i)
    {
        for(int j = 0; j < grid.y_size; ++j)
        {
            for(int k = 0; k < grid.z_size; ++k)
            {
                Cell & cell = grid.cells[i][j][k];
                for(Particle & particle : cell.particles_of_kind(kind))
                {
                    particle.density = find_density(&particle, &(cell));
                    assert(!__isnan(particle.density));
                }
            }
        }
    }
}