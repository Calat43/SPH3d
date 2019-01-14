#pragma once

#include <iostream>
#include "Params.h"

void compare_third_column(std::string read_1, std::string read_2)
{
    Params & params = Params::get_instance();

    FILE * f_1 = fopen("/home/calat/tmp/no_sort.dat", "r");
    FILE * f_2 = fopen("/home/calat/tmp/sort.dat", "r");

    double density_1[params.n_gas];
    double density_2[params.n_gas];

    double none_1[params.n_gas];

    for(int i = 0; i < params.n_gas; ++i)
    {
        density_1[i] = 0;
        density_2[i] = 0;
    }

    int elem = 0;
    while(fscanf (f_1, "%lf %lf %lf", &(none_1[elem]), &(none_1[elem]), &(density_1[elem])) != EOF)
    {
        elem++;
    }

    elem = 0;
    while (fscanf (f_2, "%lf %lf %lf", &(none_1[elem]), &(none_1[elem]), &(density_2[elem])) != EOF)
    {
        elem++;
    }

    char filename[512];
    sprintf(filename, "/home/calat/tmp/compare_dens.dat");
    FILE * f = fopen(filename, "w");

    for(int i = 0; i < params.n_gas; ++i)
    {
        fprintf(f, "%lf \n", (density_1[i] - density_2[i]));
    }
}