#include "beam.h"

#include "threevector.h"
#include "threematrix.h"
#include "my_functions.h"
#include <cmath>
#include <iostream>

Beam::Beam(int num_particle, double particle_charge, double particle_mass, double energy_central, double energy_spread, ThreeVec position_initial,
                ThreeVec position_spread, ThreeVec angle_initial, ThreeVec angle_spread,
                PositionInitializationTypes pos_init, EnergyInitializationTypes energy_init, DivergenceInitializationTypes diverge_init)
                : m_position_initialization_type(pos_init), m_energy_initialization_type(energy_init), m_divergence_initialization_type(diverge_init)
{
        m_num_particles                  = num_particle;
        m_particle_charge                = particle_charge;
        m_energy_central                 = energy_central;
        m_energy_spread                  = energy_spread;
        m_position_central               = position_initial;
        m_position_spread                = position_spread;
        m_angle_central                  = angle_initial;
        m_angle_spread                   = angle_spread;
//        m_position_initialization_type   = pos_init;
//        m_energy_initialization_type     = energy_init;
//        m_divergence_initialization_type = diverge_init;
        m_particle_counter               = 0;

        m_particle.set_charge(m_particle_charge);
        m_particle.set_mass(particle_mass);

        
        
        switch(pos_init)
        {
            case PositionInitializationTypes::INITIALIZE_GAUSSIAN_POS:
                double x_init = gaussian_init(position_initial.getX(), position_spread.getX());
                double y_init = gaussian_init(position_initial.getY(), position_spread.getY());
                double z_init = gaussian_init(position_initial.getZ(), position_spread.getZ());
                m_particle.set_pos(x_init, y_init, z_init);
                break;

        }

        switch(energy_init)
        {
            case EnergyInitializationTypes::INITIALIZE_GAUSSIAN_EN:
                //m_particle.set_energy(gaussian_init(m_energy_central, m_energy_spread));
                /*
                 * First particle will have central energy)
                 */
                m_particle.set_energy(m_energy_central);
                break;
        }

        switch(diverge_init)
        {
            case DivergenceInitializationTypes::INITIALIZE_GAUSSIAN_DIV:
                double angle_x = gaussian_init(m_angle_central.getX(), m_angle_spread.getX());
                double angle_y = gaussian_init(m_angle_central.getY(), m_angle_spread.getY());
                double angle_z = gaussian_init(m_angle_central.getZ(), m_angle_spread.getZ());
                
                double p_mag   = sqrt(m_particle.get_energy()*m_particle.get_energy() - 1.0);
                double pz      = p_mag*cos(angle_z);
                double py      = p_mag*cos(angle_y);
                double px      = sqrt((p_mag*p_mag) - (py*py) - (pz*pz) );
                //double pz      = p_mag*cos(angle_z);
                m_particle.set_p(px, py, pz);
                break;
        }
}









double Beam::gaussian_init( double initializiation_central_value, double initialization_radius_of_values )
{
    double value = gaussian()*initialization_radius_of_values + initializiation_central_value;
    return value;
}









void Beam::next_particle(int& particle_counter,
                PositionInitializationTypes pos_init, 
                EnergyInitializationTypes energy_init, 
                DivergenceInitializationTypes diverge_init)//Default init values defined in beam.h
{
    particle_counter = particle_counter + 1;


    switch(pos_init)
        {
            case PositionInitializationTypes::INITIALIZE_GAUSSIAN_POS:
                double x_init = gaussian_init(m_position_central.getX(), m_position_spread.getX());
                double y_init = gaussian_init(m_position_central.getY(), m_position_spread.getY());
                double z_init = gaussian_init(m_position_central.getZ(), m_position_spread.getZ());
                m_particle.set_pos(x_init, y_init, z_init);
                break;

        }

        switch(energy_init)
        {
            case EnergyInitializationTypes::INITIALIZE_GAUSSIAN_EN:
                m_particle.set_energy(gaussian_init(m_energy_central, m_energy_spread));
                break;
        }

        switch(diverge_init)
        {
            case DivergenceInitializationTypes::INITIALIZE_GAUSSIAN_DIV:
                double angle_x = gaussian_init(m_angle_central.getX(), m_angle_spread.getX());
                double angle_y = gaussian_init(m_angle_central.getY(), m_angle_spread.getY());
                double angle_z = gaussian_init(m_angle_central.getZ(), m_angle_spread.getZ());
                
                double p_mag   = sqrt(m_particle.get_energy()*m_particle.get_energy() - 1.0);
                //double px      = p_mag*cos(angle_x);
                double py      = p_mag*cos(angle_y);
                double pz      = p_mag*cos(angle_z);
                double px      = sqrt((p_mag*p_mag) - (py*py) - (pz*pz) );
                m_particle.set_p(px, py, pz);
                break;
        }
}