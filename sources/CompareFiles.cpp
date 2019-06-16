#include "CompareFiles.h"

void compare_third_column(std::string /* read_1 */, std::string /* read_2 */)
{
    Params & params = Params::get_instance();

    FILE * f_1 = fopen((OUTPUT_PATH + "no_sort.dat").c_str(), "r");
    FILE * f_2 = fopen((OUTPUT_PATH + "sort.dat").c_str(), "r");

    double * density_1 = new double[params.n_gas];
    double * density_2 = new double[params.n_gas];

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

    std::string filename = OUTPUT_PATH + "compare_dens.dat";
    FILE * f = fopen(filename.c_str(), "w");


    for(int i = 0; i < params.n_gas; ++i)
    {
        fprintf(f, "%lf \n", (density_1[i] - density_2[i]));
    }
    fclose(f);
}

void print_grid_dens(Grid grid, char * name)
{
    char filename[512];
    sprintf(filename, (OUTPUT_PATH + "%s.dat").c_str(), name);
    FILE * f = fopen(filename, "w");

    Cell & cell = grid.cells[20][10][10];

    fprintf(f, "[%i][%i][%i]\n", 20, 10, 10);

    for(Particle * particle : cell.get_all_particles())
    {
        fprintf(f, "%lf %d\n", particle->density, particle->dbg_state);
    }

    fclose(f);
}

void centering()
{
    FILE * f_1 = fopen((OUTPUT_PATH + "part_99.dat").c_str(), "r");

    double coord_x[40000];
    double coord_y[40000];
    double coord_z[40000];
    double vel_x[40000];
    double vel_y[40000];
    double vel_z[40000];
    double dens[40000];
    double pres[40000];
    double ener[40000];

    int elem = 0;
    while(fscanf (f_1, "%lf %lf %lf %lf %lf %lf %lf %lf %lf", &(coord_x[elem]), &(coord_y[elem]), &(coord_z[elem]),
            &(vel_x[elem]), &(vel_y[elem]), &(vel_z[elem]), &(dens[elem]), &(pres[elem]), &(ener[elem])) != EOF)
    {
        elem++;
    }

    std::string filename = OUTPUT_PATH + "center_tube.dat";
    FILE * f = fopen(filename.c_str(), "w");

    for(int i = 0; i < elem; ++i)
    {
        if(coord_y[i] > 0.03 && coord_y[i] < 0.07 && coord_z[i] > 0.03 && coord_z[i] < 0.07)
        {
            fprintf(f, "%lf %lf %lf %lf %lf\n", coord_x[i], vel_x[i], dens[i], pres[i], ener[i]);
        }
    }

    fclose(f_1);
    fclose(f);
}