#include <iostream>
//#include <string>
//#include <fstream>
//#include <time.h>

#include "Params.h"
#include "SPHSolver.h"


int main()
{
    Params & params = Params::get_instance();

    // 0. Solver initialization (grid, arrays, etc.)
    SPHSolver solver;

    // 1. Generate initial coniditions (values from global parameters are used)
    solver.generate_initial_conditions();

    // 2. Main Loop with time step increment
    int nb_time_steps = (int) floor(params.t / params.tau);
    //int nb_time_steps_diag = floor(params.t / params.t_diagnostics);
    for (int time_step_id = 0; time_step_id < nb_time_steps; ++time_step_id)
    {
        solver.do_time_step();

       // 3. Run intermediate diagnostics and store output
       // if (time_step_id % nb_time_steps_diag)
       // {
       //     run diagnostics
       // }
    }

    return 0;
}
