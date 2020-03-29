#include <iostream>
#include <string>
#include <fstream>
#include <time.h>

#include "Particle.h"
#include "MathUtils.h"
#include "Grid.h"
#include "Cell.h"
#include "Params.h"
#include "BallInVacuum.h"
#include "DustyWave.h"


void check_output_dir()
{
    std::ofstream test_file(OUTPUT_PATH + "test");
    if (!test_file.is_open())
    {
        std::cout << "Cant write to the output directory";
        std::exit(1);
    }
    else
    {
        test_file.close();
    }
}

int main()
{
    check_output_dir();

    Params & params = Params::get_instance();

/*
    std::cout << ball_analytic::time_through_R(0.9) << std::endl;
    std::cout << ball_analytic::time_through_R(0.87714) << std::endl;

    std::cout << ball_analytic::R_bisection(0.2, 0.00001, 0.000001) << std::endl;
    std::cout << ball_analytic::time_through_R(ball_analytic::R_bisection(0.2, 0.00001, 0.000001)) << std::endl;

 */
    //ball_analytic::print_solution(0., 0.01, 0.00001, 0.00001);

    /*
    Grid grid = ball_in_vacuum::init(0.05);

    clock_t startTime = clock();
    for (int frameId = 0; frameId < floor(params.t / params.tau); ++frameId)
    {
        clock_t step_start = clock();

        grid = ball_in_vacuum::do_time_step(grid, frameId);
        clock_t step_fin = clock();

        double step_time = (double)(step_fin - step_start) / CLOCKS_PER_SEC;
        std::cout << frameId << " " << step_time << " s" << std::endl;
        //break;
    }

    std::cout << "Done!" << std::endl;
    clock_t finishTime = clock();

    double executionTime = (double)(finishTime - startTime) / CLOCKS_PER_SEC;
    printf("Finished in %lf seconds.\n", executionTime);

    //print_vel_p_squared(0.4, 0.14, 0.001);
*/

    dusty_wave_init();

    return 0;
}
