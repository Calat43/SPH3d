
#include "Params.h"
#include "Cell.h"
#include "NonLinear.h"

#include "DustyShock.h"


/*
double find_lambda(double gas_density)
{
    return 5. * pow(10, -6) / gas_density;
}

double sound_speed(Particle * gas_particle)
{
    return sqrt(Params::get_instance().gamma * gas_particle->pressure / gas_particle->density);
}

double mach_number(Particle * gas_particle, Particle * dust_particle)
{
    return fabs(gas_particle->vx - dust_particle->vx) / sound_speed(gas_particle);
}

double knudsen_number(double gas_density)
{
    return find_lambda(gas_density) / Params::get_instance().a;
}

double reynolds_number(Particle * gas_particle, Particle * dust_particle)
{
    return 4. * mach_number(gas_particle, dust_particle) / knudsen_number(gas_particle->density);
}

double find_t_stop(Particle * gas_particle, Particle * dust_particle)
{
    double lambda = find_lambda(gas_particle->density);
    double Ma = mach_number(gas_particle, dust_particle);
    double Re = reynolds_number(gas_particle, dust_particle);
    double sound_speed = sound_speed(gas_particle);

    double vel_ratio = fabs(gas_particle->vx - dust_particle->vx);

    double a = Params::get_instance().a;
    double particle_density = Params::get_instance().particle_density;

    double result = 0;

    if(a < 9. / 4. * lambda || Ma == 0)
    {
        result = a * particle_density / sound_speed / gas_particle->density;
    }
    else if(Re <= 1)
    {
        result = 8. / 3. * 24. * 4. * a / lambda * a * particle_density / sound_speed / gas_particle->density;
                // particle_density * Re / gas_particle->density / vel_ratio;
    }
    else if(Re <= 800)
    {
        result = 8. / 3. * 24 * a * particle_density * pow(Re, 0.6) / gas_particle->density / vel_ratio;
    }
    else if(Re > 800)
    {
        result = 8. / 3. * a * particle_density / gas_particle->density / vel_ratio * 0.44;
    }

    return result;
}
 */

double nl_idic::lambda_asterisk(Cell * cell)
{
    if(cell->gas_particles.size() == 0)
    {
        return 0;
    }
    double gas_density_astr = idic::density_asterisk(cell, Particle::Kind::Gas);

    return 5. * pow(10, -6) / gas_density_astr;
}

double nl_idic::sound_speed_asterisk(Cell * cell)
{
    double result = 0;

    int size = (int)cell->gas_particles.size();

    if(size == 0)
    {
        return result;
    }

    for(Particle & particle : cell->gas_particles)
    {
        result += sqrt(particle.pressure / particle.density);
        //result += sqrt(particle.energy);
    }

    return sqrt(Params::get_instance().gamma) * result / size;

}

double nl_idic::mach_number_asterisk(Cell * cell)
{
    double gas_vel_astr = idic::vel_asterisk(cell, Particle::Kind::Gas).x;
    double dust_vel_astr = idic::vel_asterisk(cell, Particle::Kind::Dust).x;

    if(sound_speed_asterisk(cell) == 0)
    {
        return 0;
    }
    return fabs(gas_vel_astr - dust_vel_astr) / sound_speed_asterisk(cell);
}

double nl_idic::knudsen_number_asterisk(Cell * cell)
{
    return nl_idic::lambda_asterisk(cell) / Params::get_instance().a;
}

double nl_idic::reynolds_number_asterisk(Cell * cell)
{
    if(nl_idic::knudsen_number_asterisk(cell) == 0)
    {
        return 0;
    }

    return  4. * nl_idic::mach_number_asterisk(cell) / nl_idic::knudsen_number_asterisk(cell);
}

double nl_idic::t_stop_asterisk(Cell * cell)
{
    double lambda_astr = nl_idic::lambda_asterisk(cell);
    double Ma_astr = nl_idic::mach_number_asterisk(cell);
    double Re_astr = nl_idic::reynolds_number_asterisk(cell);
    double c_s_astr = nl_idic::sound_speed_asterisk(cell);

    double gas_vel_astr = idic::vel_asterisk(cell, Particle::Kind::Gas).x;
    double dust_vel_astr = idic::vel_asterisk(cell, Particle::Kind::Dust).x;

    double vel_ratio = fabs(gas_vel_astr - dust_vel_astr);

    double a = Params::get_instance().a;
    double particle_density = Params::get_instance().particle_density;
    double gas_density_astr = idic::density_asterisk(cell, Particle::Kind::Gas);

    double result = 0;

    if(cell->dust_particles.size() == 0 || cell->gas_particles.size() == 0)
    {
        cell->set_drag_mode(-2);
        return 0;
    }
    if(a < 9. / 4. * lambda_astr || Ma_astr == 0)
    {
        result = a * particle_density / c_s_astr / gas_density_astr;
        cell->set_drag_mode(1);
    }
    else if(Re_astr <= 1)
    {
        result = 8. / 3. / 24. * 4. * a / lambda_astr * a * particle_density / c_s_astr / gas_density_astr;
        cell->set_drag_mode(2);
    }
    else if(Re_astr <= 800)
    {
        result = 8. / 3. / 24. * a * particle_density * pow(Re_astr, 0.6) / gas_density_astr / vel_ratio;
                cell->set_drag_mode(3);
    }
    else if(Re_astr > 800)
    {
        result = 8. / 3. * a * particle_density / gas_density_astr / vel_ratio / 0.44;
        cell->set_drag_mode(4);
    }

    return result;
}


