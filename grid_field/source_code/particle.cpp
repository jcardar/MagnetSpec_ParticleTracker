#include "particle.h"

#include "threevector.h"
#include "threematrix.h"
#include "my_functions.h"
#include <cmath>
#include <iostream>


Particle::Particle(ThreeVec pos, ThreeVec momentum, int charge, double& time, std::ofstream& OUT_time, 
            std::ofstream& OUT_posx, std::ofstream& OUT_posy, std::ofstream& OUT_posz,
            std::ofstream& OUT_velx, std::ofstream& OUT_vely, std::ofstream& OUT_velz, std::ofstream& OUT_energy)
                : m_time(&time), m_out_time(&OUT_time), m_out_posx(&OUT_posx), m_out_posy(&OUT_posy), m_out_posz(&OUT_posz),
                    m_out_px(&OUT_velx), m_out_py(&OUT_vely), m_out_pz(&OUT_velz), m_out_energy(&OUT_energy)
    {
        set_pos(pos);
        set_p(momentum);
        m_charge = charge;
        m_vel.setX(momentum.getX()/m_energy);
        m_vel.setY(momentum.getY()/m_energy);
        m_vel.setZ(momentum.getZ()/m_energy);
    }

void Particle::set_pos(ThreeVec pos)
{
    m_pos = pos;
}

void Particle::set_pos(int index, double value)
{
    m_pos.set(index, value);
}

void Particle::set_pos(double x, double y, double z)
{
    m_pos.setX(x);
    m_pos.setY(y);
    m_pos.setZ(z);
}

void Particle::set_p(ThreeVec momentum)
{
    m_p = momentum;
    set_energy(sqrt(momentum.mag()*momentum.mag() + 1));
    m_vel = momentum/m_energy;
} 

void Particle::set_p(int index, double value)
{
    m_p.set(index, value);
    set_energy(sqrt(m_p.mag()*m_p.mag() + 1));
    m_vel.set(index, value/m_energy);
}

void Particle::set_p(double x, double y, double z)
{
    m_p.setX(x);
    m_p.setY(y);
    m_p.setZ(z);
    set_energy(sqrt(m_p.mag()*m_p.mag() + 1));
    m_vel = m_p/m_energy;
}

void Particle::set_energy(double t_energy)
{
    m_energy = t_energy;
}

void Particle::set_time(double& time)
{
    m_time = &time;
}

void Particle::set_vel(ThreeVec vel)
{
    m_vel = vel;
}

void Particle::set_vel(int index, double value)
{
    m_vel.set(index, value);
}



/*Particle::Particle(ThreeVec pos, int charge, ThreeVec momentum, double& time, std::ofstream& OUT_time, 
            std::ofstream& OUT_posx, std::ofstream& OUT_posy, std::ofstream& OUT_posz,
            std::ofstream& OUT_velx, std::ofstream& OUT_vely, std::ofstream& OUT_velz, std::ofstream& OUT_energy)
                : m_time(&time), m_out_time(&OUT_time), m_out_posx(&OUT_posx), m_out_posy(&OUT_posy), m_out_posz(&OUT_posz),
                    m_out_px(&OUT_velx), m_out_py(&OUT_vely), m_out_pz(&OUT_velz), m_out_energy(&OUT_energy)
    {
        set_pos(pos);
        if(energy.getX() >= 0.0)
        {
            set_energy(0, (( energy.getX() )) );
        }
            else
            {
                m_energy.setX( (((-1)*energy.getX()) ) );
                double rel_energy_x = ((-1)*energy.getX() + m_rest_en)/m_rest_en;
                m_vel.setX( 0.0-sqrt(1 - ( 1/(rel_energy_x*rel_energy_x) )) );
            }
        
        if(energy.getY() >= 0.0)
        {
            set_energy(1, (( energy.getY() )) );
        }
            else
            {
                m_energy.setY( (((-1)*energy.getY()) ) );
                double rel_energy_y = ((-1)*energy.getY() + m_rest_en)/m_rest_en;
                m_vel.setY( 0.0-sqrt(1 - ( 1/(rel_energy_y*rel_energy_y) )) );
            }
        if(energy.getZ() >= 0.0)
        {
            set_energy(2, (( energy.getZ() )) );
        }
            else
            {
                m_energy.setZ( (((-1)*energy.getZ() )));
                double rel_energy_z = ((-1)*energy.getZ() + m_rest_en)/m_rest_en;
                m_vel.setZ( 0.0-sqrt(1 - ( 1/(rel_energy_z*rel_energy_z) )) );
            }
        m_charge = charge;
    }

    */

/*
void Particle::set_energy(ThreeVec energy_t)
{
    
    if(energy_t.getX() >= 0.0)
    {
        m_energy.setX(energy_t.getX());
        double rel_energy_x = (energy_t.getX() + m_rest_en)/m_rest_en;
        set_vel(0, sqrt( 1 - (1/((rel_energy_x)*(rel_energy_x))) ) ); 
    }
        else
        { 
            m_energy.setX((energy_t.getX()*(-1)) );
            double rel_energy_x = (-1*energy_t.getX() + m_rest_en)/m_rest_en;
            set_vel(0, 0.0-sqrt( 1 - (1/((rel_energy_x)*(rel_energy_x))) ) ); 
        }
    if(energy_t.getY() >= 0.0)
    {
        m_energy.setY(energy_t.getY());
        double rel_energy_y = (energy_t.getY() + m_rest_en)/m_rest_en;
        set_vel(1, sqrt( 1 - (1/((rel_energy_y)*(rel_energy_y))) ) );
        //std::cerr << m_energy.getY() << "  ";
    }
        else
        {
            m_energy.setY(-1.0*energy_t.getY() );
            double rel_energy_y = (-1*energy_t.getY() + m_rest_en)/m_rest_en;
            set_vel(1, 0.0-sqrt( 1 - (1/((rel_energy_y)*(rel_energy_y))) ) );
        }
    if(energy_t.getZ() >= 0.0)
    {
        m_energy.setZ(energy_t.getZ());
        double rel_energy_z = (energy_t.getZ() + m_rest_en)/m_rest_en;
        set_vel(2, sqrt( 1 - (1/((rel_energy_z)*(rel_energy_z))) ) );
        //std::cerr << m_energy.getZ() << '\n';
    }
        else
        {
            m_energy.setZ(-1*energy_t.getZ());
            double rel_energy_z = (-1*energy_t.getZ() + m_rest_en)/m_rest_en;
            set_vel(2, 0.0-sqrt( 1 - (1/((rel_energy_z)*(rel_energy_z))) ) );
        }
    //std::cerr << m_vel.getX() << '  ' << m_vel.getY() << '  ' << m_vel.getZ() << '\n';
    //std::cerr << "energy in x is " << get_energy(0) << " and velocity in x is " << get_vel(0) << "\n"; 
    //std::cerr << "energy in y is " << get_energy(1) << " and velocity in y is " << get_vel(1) << "\n"; 
    //std::cerr << "energy in z is " << get_energy(2) << " and velocity in z is " << get_vel(2) << "\n"; 
}

void Particle::set_momentum(int index, double value)
{
    if(value >= 0)
    {
        m_energy.set(index, value);
        double rel_gamma = (value+m_rest_en)/m_rest_en;
        set_vel(index, sqrt( 1 - (1/(rel_gamma*rel_gamma)) ));
    }
}
*/