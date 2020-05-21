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

    int m_particle_counter;
    
    enum EnergyInitializationTypes
    {
        INITIALIZE_GAUSSIAN_EN = 0,

    };

    enum PositionInitializationTypes
    {
        INITIALIZE_GAUSSIAN_POS = 0,

    };

    enum DivergenceInitializationTypes
    {
        INITIALIZE_GAUSSIAN_DIV = 0,

    };

    EnergyInitializationTypes     m_energy_initialization_type;
    PositionInitializationTypes   m_position_initialization_type;
    DivergenceInitializationTypes m_divergence_initialization_type;

    Beam (int num_particle, double particle_charge, double particle_mass, double central_energy, double energy_spread,
                ThreeVec position_initial, ThreeVec position_spread, ThreeVec angle_initial, ThreeVec angle_spread, 
                PositionInitializationTypes pos_init       = PositionInitializationTypes::INITIALIZE_GAUSSIAN_POS, 
                EnergyInitializationTypes energy_init      = EnergyInitializationTypes::INITIALIZE_GAUSSIAN_EN, 
                DivergenceInitializationTypes diverge_init = DivergenceInitializationTypes::INITIALIZE_GAUSSIAN_DIV);

    
    

    double gaussian_init( double initializiation_central_value, double initialization_radius_of_values );

    void next_particle( int& particle_counter,
                PositionInitializationTypes pos_init       = PositionInitializationTypes::INITIALIZE_GAUSSIAN_POS, 
                EnergyInitializationTypes energy_init      = EnergyInitializationTypes::INITIALIZE_GAUSSIAN_EN, 
                DivergenceInitializationTypes diverge_init = DivergenceInitializationTypes::INITIALIZE_GAUSSIAN_DIV);

    Particle& get_particle()
    {
        return m_particle;
    }

};



#endif