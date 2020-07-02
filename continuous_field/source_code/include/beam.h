#ifndef BEAM_H
#define BEAM_H

#include "threematrix.h"
#include "threevector.h"
#include <fstream>
#include "particle.h"

class Beam
{
private:
    int      m_num_particles;
    double   m_particle_charge;
    double   m_energy_central; 
    double   m_energy_spread;
    ThreeVec m_position_central;
    ThreeVec m_position_spread;
    ThreeVec m_angle_central;      //radians from central axis
    ThreeVec m_angle_spread;
 //   double   m_divergence; //in milliradians
    Particle m_particle;

public:
    Beam (){};

    int m_particle_counter = 0;
    
    enum EnergyInitializationTypes
    {
        INITIALIZE_GAUSSIAN_EN = 0,
        INITIALIZE_UNIFORM_EN  = 1,
        INITIALIZE_LOG_EN      = 2,
    };

    enum PositionInitializationTypes
    {
        INITIALIZE_GAUSSIAN_POS = 0,
        INITIALIZE_UNIFORM_POS  = 1,
        INITIALIZE_SCAN_POS  = 2,
    };

    enum DivergenceInitializationTypes
    {
        INITIALIZE_GAUSSIAN_DIV = 0,
        INITIALIZE_UNIFORM_DIV  = 1,
        INITIALIZE_SCAN_DIV     = 2,
    };

    EnergyInitializationTypes     m_energy_initialization_type     = INITIALIZE_GAUSSIAN_EN;
    PositionInitializationTypes   m_position_initialization_type   = INITIALIZE_GAUSSIAN_POS;
    DivergenceInitializationTypes m_divergence_initialization_type = INITIALIZE_GAUSSIAN_DIV;

    Beam (int num_particle, double particle_charge, double particle_mass, double central_energy, double energy_spread,
                ThreeVec position_initial, ThreeVec position_spread, ThreeVec angle_initial, ThreeVec angle_spread, 
                PositionInitializationTypes pos_init       = PositionInitializationTypes::INITIALIZE_GAUSSIAN_POS, 
                EnergyInitializationTypes energy_init      = EnergyInitializationTypes::INITIALIZE_GAUSSIAN_EN, 
                DivergenceInitializationTypes diverge_init = DivergenceInitializationTypes::INITIALIZE_GAUSSIAN_DIV);

    
    

    double gaussian_init( double initializiation_central_value, double initialization_radius_of_values );

    void next_particle()
    {
        next_particle(m_particle_counter, m_position_initialization_type, m_energy_initialization_type, m_divergence_initialization_type);
    }
    void next_particle( int& particle_counter,
                PositionInitializationTypes pos_init, 
                EnergyInitializationTypes energy_init, 
                DivergenceInitializationTypes diverge_init);

    Particle& get_particle()
    {
        return m_particle;
    }

};



#endif