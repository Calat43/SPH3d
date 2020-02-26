#include "Params.h"
#include "Grid.h"

#include "SPHSolver.h"

SPHSolver::SPHSolver()
{
    // initialize grid
    Params & params = Params::get_instance();
    current_grid = new Grid(params.grid_step_x, params.grid_step_y, params.grid_step_z, params.border1, params.border2);
    next_grid = new Grid(params.grid_step_x, params.grid_step_y, params.grid_step_z, params.border1, params.border2);
}

SPHSolver::~SPHSolver()
{
    delete current_grid;
    delete next_grid;
}

int SPHSolver::generate_initial_distribution()
{
    Params & params = Params::get_instance();

    switch (params.dt) {
        case Params::dtUniform: 
        { 
            return generate_uniform_distribution();
        }
        case Params::dtBall: 
        { 
            return generate_ball_distribution();
        }
        default:
            return 1;
    }

}

int SPHSolver::generate_uniform_distribution()
{
    //Params & params = Params::get_instance();
    // params.sz_gas_particles_total;
    // params.sz_dust_particles_total;

    // calculate coordinates of dust / gas particles

    // put them into cells

    // calculate density and pressure for each particle
    return 0;
}


int SPHSolver::generate_ball_distribution()
{
    return 0;
}


int SPHSolver::do_time_step()
{
    Params & params = Params::get_instance();

    // iterate over all cells in grid and calculate new coordinates of particles

    for (int i = 0; i < current_grid->x_size; ++i) 
    {
        for (int j = 0; j < current_grid->y_size; ++j)
        {
            for (int k = 0; k < current_grid->z_size; ++k)
            {
                Cell & cell = current_grid->cells[i][j][k];
                for (Particle * particle : cell.get_all_particles())
                {
                    Particle new_particle(*particle);
                    Point new_coords = find_new_coordinates(*particle);

                    Point new_vel = MathUtils::find_new_velocity(particle, &cell); // was got from Sod_tube_3d
                    new_particle.density = NAN;
                    new_particle.dbg_state = 2;

                    new_particle.set_coordinates(new_coords.x, new_coords.y, new_coords.z);
                    new_particle.set_velocities(new_vel.x, new_vel.y, new_vel.z);

                    new_particle.energy = MathUtils::find_new_energy(particle, &cell); // was got from Sod_tube_3d

                    next_grid->with_copy_of(new_particle);
                }
            }
        }
    }

    MathUtils::recalc_density(*next_grid, Particle::Kind::Gas); // was got from file Sod_tube_1d (not its namespace)
    MathUtils::recalc_pressure(*next_grid, Particle::Kind::Gas); // was got from Sod_tube_1d

    // swap grids and clear next grid (prepare for the next timestep)
    swap_grids_and_clear_next();
    return 0;
}


void SPHSolver::swap_grids_and_clear_next()
{
    Grid* grid = current_grid;
    current_grid = next_grid;
    next_grid = grid;
    next_grid->clear_all_particles();
}

