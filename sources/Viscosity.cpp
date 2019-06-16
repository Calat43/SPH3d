#include "Viscosity.h"

double viscosity::find_sound_speed(Particle * particle)
{
    assert(!__isnan(particle->pressure));
    assert(!__isnan(particle->density));
    assert(!__isnan(sqrt(Params::get_instance().gamma * particle->pressure / particle->density)));
    return sqrt(Params::get_instance().gamma * particle->pressure / particle->density);
}

double viscosity::find_mu(Point vel_ab, Point coord_ab)
{
    return (Params::get_instance().h * dot_product(vel_ab, coord_ab)) /
           (pow(magnitude(coord_ab), 2) + pow(Params::get_instance().nu, 2));
}

double viscosity::find_mu(double vel_ab, double coord_ab)
{
    return (Params::get_instance().h * vel_ab * coord_ab) / (pow(coord_ab, 2) + pow(Params::get_instance().nu, 2));
}

double viscosity::find_viscosity(Particle * p1, Particle * p2)
{
    Params & params = Params::get_instance();
    int dim = params.dimensions;
}

namespace viscosity_3d
{
    double find_viscosity(Particle * p1, Particle * p2)
    {
        Params & params = Params::get_instance();

        if(params.have_viscosity)
        {
            double vel_x_ab = p1->vx - p2->vx;
            double vel_y_ab = p1->vy - p2->vy;
            double vel_z_ab = p1->vz - p2->vz;

            double coord_x_ab = p1->x - p2->x;
            double coord_y_ab = p1->y - p2->y;
            double coord_z_ab = p1->z - p2->z;

            Point vel_ab(vel_x_ab, vel_y_ab, vel_z_ab);
            Point coord_ab(coord_x_ab, coord_y_ab, coord_z_ab);

            if(dot_product(vel_ab, coord_ab) >= 0)
            {
                return 0;
            }
            if(dot_product(vel_ab, coord_ab) < 0)
            {
                double c_a = viscosity::find_sound_speed(p1);
                double c_b = viscosity::find_sound_speed(p2);

                double dens_ab = 1. / 2. * (p1->density + p2->density);
                double c_ab = 1. / 2. * (c_a + c_b);

                double mu_ab = viscosity::find_mu(vel_ab, coord_ab);

                assert(!__isnan(c_a));
                assert(!__isnan(c_b));
                assert(!__isnan(mu_ab));

                return (- params.alpha * c_ab * mu_ab + params.beta * pow(mu_ab, 2)) / dens_ab;
            }
        }
        return 0;
    }
}

double viscosity_1d::find_viscosity(Particle * p1, Particle * p2)
{
    Params & params = Params::get_instance();

    if(params.have_viscosity)
    {
        double vel_ab = p1->vx - p2->vx;
        double coord_ab = p1->x - p2->x;

        if(vel_ab * coord_ab >= 0)
        {
            return 0;
        }
        if(vel_ab * coord_ab < 0)
        {
            double c_a = viscosity::find_sound_speed(p1);
            double c_b = viscosity::find_sound_speed(p2);
            double dens_ab = 1. / 2 * (p1->density + p2->density);
            double c_ab = 1. / 2 * (c_a + c_b);

            double mu_ab = viscosity::find_mu(vel_ab, coord_ab);

            assert(!__isnan(c_a));
            assert(!__isnan(c_b));
            assert(!__isnan(mu_ab));

            return (- params.alpha * c_ab * mu_ab + params.beta * pow(mu_ab, 2)) / dens_ab;
        }
    }
    return 0;
}