#include <fstream>
#include <iostream>

#include "Params.h"

#include "Grid.h"
#include "Cell.h"
#include "Viscosity.h"
#include "DustyShock1d.h"
#include "NonLinear.h"
#include "MathUtils.h"

#include "Solver.h"

double sph_dens_approx(double coord, Particle::Kind approx_kind, Grid & grid)
{
    Particle tmp_particle(approx_kind, coord, 0, 0);

    Cell * cell = grid.find_cell(tmp_particle);

    double result = 0;

    for (Cell * neighbour : cell->get_neighbours())
    {
        for (Particle & p : approx_kind == Particle::Kind::Gas ? neighbour->gas_particles : neighbour->dust_particles)
        {
            result += p.mass * kernel(tmp_particle, p, Params::get_instance().dimensions);
        }
    }

    return result;
}

double theta(Particle * particle, Grid & grid)
{
    double intrinstic_dens = Params::get_instance().particle_density;
    double result = 0;

    if(particle->kind == Particle::Kind::Dust)
    {
        result = particle->density / intrinstic_dens;
    } else if(particle->kind == Particle::Kind::Gas)
    {
        result = sph_dens_approx(particle->x, Particle::Kind::Dust, grid) / intrinstic_dens;
    }

    return result;
}

double gas_pressure_term(Particle * gas_p, Cell * cell) {

    Params & params = Params::get_instance();
    int dim = params.dimensions;
    double c_s = params.gas_sound_speed;

    double result = 0;
    double gas_term = 0;
    double dust_term = 0;

    for(Cell *neighbour : cell->get_neighbours())
    {
        for(Particle &pg : neighbour->gas_particles)
        {
            gas_term += (1. / gas_p->density + 1. / pg.density) * kernel_gradient(*gas_p, pg, dim).x;
        }
        for(Particle &pd : neighbour->dust_particles)
        {
            dust_term += pd.mass * kernel_gradient(*gas_p, pd, dim).x;
        }
    }

    if(params.particle_density == 0)
    {
        result = -gas_term * c_s * c_s * gas_p->mass;
    }
    else
    {
        result = -gas_term * c_s * c_s * gas_p->mass - dust_term * c_s * c_s / params.particle_density / (1. - gas_p->theta);
    }

    assert(!__isnan(result));

    return result;
}

double dust_pressure_term(Particle * dust_p, Cell * cell)
{
    Params & params = Params::get_instance();
    int dim = params.dimensions;

    double gas_term = 0;

    if(params.particle_density == 0)
        return 0;

    for (Cell *neighbour : cell->get_neighbours())
    {
        for (Particle &p : neighbour->gas_particles)
        {
            gas_term += p.mass / (1. - p.theta) * kernel_gradient(*dust_p, p, dim).x;
        }
    }

    return pow(params.gas_sound_speed, 2) / params.particle_density * gas_term;
}

double gas_pres_term_astr(Cell * cell)
{
    size_t gas_size = cell->gas_particles.size();

    double a_p_x_astr = 0;

    if(gas_size == 0)
    {
        return 0;
    }
    else
    {
        for (Particle & p : cell->gas_particles)
        {
            a_p_x_astr += gas_pressure_term(&p, cell);
        }
    }

    return a_p_x_astr / gas_size;
}

double dust_pres_term_astr(Cell * cell)
{
    size_t dust_size = cell->dust_particles.size();

    double a_p_x_astr = 0;

    if(dust_size == 0)
    {
        return 0;
    }
    else
    {
        for (Particle & p : cell->dust_particles)
        {
            a_p_x_astr += dust_pressure_term(&p, cell);
        }
    }

    return a_p_x_astr / dust_size;
}


double idic_1d::vel_asterisk(Cell * cell, Particle::Kind kind)
{
    int size = (int)cell->particles_of_kind(kind).size();

    double vel_x = 0;

    if(size == 0)
    {
        return vel_x;
    }

    for(Particle & particle : cell->particles_of_kind(kind))
    {
        vel_x += particle.vx;
    }

    return vel_x / (double)size;
}

double idic_1d::eps_asterisk(Cell * cell)
{
    //TODO int or not?
    size_t dust_size = cell->dust_particles.size();
    size_t gas_size = cell->gas_particles.size();

    double dust_mass = 0;

    if(gas_size == 0)
    {
        return 0;
    }
    else
    {
        if(dust_size != 0)
        {
            dust_mass = cell->dust_particles.at(0).mass;
        }

        double gas_mass = cell->gas_particles.at(0).mass;

        double tmp = dust_mass * dust_size / (gas_mass * gas_size);

        return tmp;
    }
}

double idic_1d::x_through_vel(Cell * cell)
{
    return idic_1d::vel_asterisk(cell, Particle::Kind::Gas) - idic_1d::vel_asterisk(cell, Particle::Kind::Dust);
}

double idic_1d::y_through_vel(Cell * cell)
{
    return idic_1d::vel_asterisk(cell, Particle::Kind::Gas) + idic_1d::vel_asterisk(cell, Particle::Kind::Dust)
                                                              * idic_1d::eps_asterisk(cell);
}

double idic_1d::t_stop_asterisk(Cell *cell)
{
    size_t dust_size = cell->dust_particles.size();

    double sum = 0;

    if(dust_size == 0)
    {
        return 0;
    }

    for(Particle & particle : cell->dust_particles)
    {
        sum += particle.density;
    }

    return Params::get_instance().t_stop;
    //return sum / dust_size / Params::get_instance().K;
}

double idic_1d::find_x(Cell * cell)
{
    double tau = Params::get_instance().tau;

    double x_prev = idic_1d::x_through_vel(cell);
    //double a_p_astr = idic_1d::pressure_term_asterisk(cell);
    double a_p_astr = gas_pres_term_astr(cell) - dust_pres_term_astr(cell);
    double t_stop_astr = idic_1d::t_stop_asterisk(cell); //nl_idic::t_stop_asterisk(cell);

    double front_next_x = 1. / tau + (idic_1d::eps_asterisk(cell) + 1.) / t_stop_astr;


    return (x_prev / tau + a_p_astr) / front_next_x;
}

double idic_1d::find_y(Cell * cell)
{
    double tau = Params::get_instance().tau;

    double y_prev = idic_1d::y_through_vel(cell);
    double a_p_asrt = eps_asterisk(cell) * dust_pres_term_astr(cell) + gas_pres_term_astr(cell);
            //idic_1d::pressure_term_asterisk(cell);

    return y_prev + a_p_asrt * tau;
}

double idic_1d::find_gas_vel_asterisk(Cell * cell)
{
    return (idic_1d::find_y(cell) + idic_1d::find_x(cell) * idic_1d::eps_asterisk(cell)) / (1. + idic_1d::eps_asterisk(cell));
}

double idic_1d::find_dust_vel_asterisk(Cell * cell)
{
    return (idic_1d::find_y(cell) - idic_1d::find_x(cell)) / (1. + idic_1d::eps_asterisk(cell));
}

double idic_1d::find_gas_velocity(Particle * particle, Cell * cell)
{
    double tau = Params::get_instance().tau;

    double prev_vel = particle->vx;

    double result = 0;

    double front_next_vel = 0;

    double pres_term = gas_pressure_term(particle, cell); //idic_1d::pressure_term(particle, cell);

    double t_stop_astr = idic_1d::t_stop_asterisk(cell); //nl_idic::t_stop_asterisk(cell);

    if(t_stop_astr == 0)
    {
        front_next_vel = 1. / tau;

        result = (prev_vel / tau + pres_term) / front_next_vel;
    }

    if(t_stop_astr != 0)
    {
        front_next_vel = 1. / tau + idic_1d::eps_asterisk(cell) / t_stop_astr;

        result = (prev_vel / tau + idic_1d::find_dust_vel_asterisk(cell) * (idic_1d::eps_asterisk(cell)
                                                                            / t_stop_astr)
                 + pres_term) / front_next_vel;
    }

    assert(!__isnan(result));

    return result;
}

//TODO assert Kind
double idic_1d::find_dust_velocity(Particle * particle, Cell * cell)
{
    double tau = Params::get_instance().tau;

    double t_stop_astr = idic_1d::t_stop_asterisk(cell); //nl_idic::t_stop_asterisk(cell);

    double front_next_vel = 1. / tau + 1. / t_stop_astr;

    double pres_term = dust_pressure_term(particle, cell);

    double prev_vel = particle->vx;
    double result = 0;

    if(t_stop_astr != 0)
    {
        result = (prev_vel / tau + idic_1d::find_gas_vel_asterisk(cell) / t_stop_astr + pres_term) / front_next_vel;
    }

    assert(!__isnan(result));

    return result;
}

Point find_new_coordinates(Particle const &particle)
{
    Params & params = Params::get_instance();
    double x = particle.x + params.tau * particle.vx;
    double y = particle.y + params.tau * particle.vy;
    double z = particle.z + params.tau * particle.vz;
    return {x, y, z};
}

//functions for DustyWave with non-zero theta

double gas_velocity_MK(Particle * particle, Cell * cell)
{
    Params & params = Params::get_instance();

    double gas_pres_sum = 0;
    double dust_pres_sum = 0;
    double drag_sum = 0;
    double result = 0;
    double rja = 0;

    for(Cell * neighbour : cell->get_neighbours())
    {
        for(Particle &gas_p : neighbour->gas_particles)
        {
            gas_pres_sum += (1. / particle->density + 1. / gas_p.density) * kernel_gradient(*particle, gas_p, params.dimensions).x;
        }

        for(Particle &dust_p : neighbour->dust_particles)
        {
            dust_pres_sum += dust_p.mass  / (1 - particle->theta) * kernel_gradient(*particle, dust_p, params.dimensions).x;

            rja = dust_p.x - particle->x;

            drag_sum += dust_p.mass * 1. / (particle->density * dust_p.density)
                     * ((particle->vx - dust_p.vx) * rja / (rja * rja + params.eta_squared)) * rja
                     * kernel(dust_p, *particle, params.dimensions);
        }
    }

    gas_pres_sum *= pow(params.c_s, 2) * particle->mass;

    if(params.particle_density != 0) {
        dust_pres_sum *= 1. / params.particle_density * pow(params.c_s, 2);
    } else dust_pres_sum = 0;

    drag_sum *= params.sigma * params.K;

    result = params.tau * (- gas_pres_sum - dust_pres_sum - drag_sum) + particle->vx;

    return result;
}

double dust_velocity_MK(Particle * particle, Cell * cell)
{
    Params & params = Params::get_instance();

    double pres_sum = 0;
    double drag_sum = 0;
    double result = 0;
    double rja = 0;

    for(Cell * neighbour : cell->get_neighbours())
    {
        for(Particle &gas_p : neighbour->gas_particles)
        {
            pres_sum += gas_p.mass / (1 - gas_p.theta) * kernel_gradient(*particle, gas_p, params.dimensions).x;

            rja = particle->x - gas_p.x;

            drag_sum += gas_p.mass / (particle->density * gas_p.density)
                        * ((gas_p.vx - particle->vx) * rja / (rja * rja + params.eta_squared)) * rja
                        * kernel(*particle, gas_p, params.dimensions);
        }
    }

    if(params.particle_density != 0) {
        pres_sum *= pow(params.c_s, 2) / params.particle_density;
    }else pres_sum = 0;

    drag_sum *= params.sigma * params.K;

    result = params.tau * (-pres_sum + drag_sum) + particle->vx;

    return result;
}

double gas_velocity_shock_MK(Particle * particle, Cell * cell)
{
    Params & params = Params::get_instance();

    double gas_pres_sum = 0;
    double dust_pres_sum = 0;
    double drag_sum = 0;
    double result = 0;
    double rja = 0;
    double viscosity = 0;
    double theta = params.theta_const;
    double rho_s = params.particle_density;

    for(Cell * neighbour : cell->get_neighbours())
    {
        for(Particle &gas_p : neighbour->gas_particles)
        {
            viscosity = viscosity_1d::find_viscosity(particle, &gas_p);
            gas_pres_sum += (gas_p.pressure / pow(gas_p.density, 2) + particle->pressure / pow(particle->density, 2)
                         + viscosity) * kernel_gradient(*particle, gas_p, params.dimensions).x;
        }

        for(Particle &dust_p : neighbour->dust_particles)
        {
            rja = dust_p.x - particle->x;

            drag_sum += dust_p.mass * 1. / (particle->density * dust_p.density)
                        * ((particle->vx - dust_p.vx) * rja / (rja * rja + params.eta_squared)) * rja
                        * kernel(dust_p, *particle, params.dimensions);

            dust_pres_sum += dust_p.mass * kernel_gradient(*particle, dust_p, params.dimensions).x;
        }
    }

    gas_pres_sum = - (1 - theta) * gas_pres_sum * particle->mass;

    if(rho_s != 0) {
        dust_pres_sum = -particle->pressure / particle->density / rho_s * dust_pres_sum;
    } else dust_pres_sum = 0;

    drag_sum = - drag_sum * params.sigma * params.K;

    result = params.tau * (gas_pres_sum + dust_pres_sum + drag_sum) + particle->vx;

    return result;
}

double dust_velocity_shock_MK(Particle * particle, Cell * cell)
{
    Params & params = Params::get_instance();

    double pres_sum = 0;
    double drag_sum = 0;
    double result = 0;
    double rja = 0;

    for(Cell * neighbour : cell->get_neighbours())
    {
        for(Particle &gas_p : neighbour->gas_particles)
        {
            rja = particle->x - gas_p.x;

            drag_sum += gas_p.mass / (particle->density * gas_p.density)
                        * ((gas_p.vx - particle->vx) * rja / (rja * rja + params.eta_squared)) * rja
                        * kernel(*particle, gas_p, params.dimensions);

            pres_sum += gas_p.mass * gas_p.pressure / gas_p.density * kernel_gradient(gas_p, *particle, params.dimensions).x;
        }
    }

    drag_sum *= params.sigma * params.K;

    if(params.particle_density != 0) {
        pres_sum = pres_sum / params.particle_density;
    } else pres_sum = 0;

    result = params.tau * (pres_sum + drag_sum) + particle->vx;

    return result;
}
