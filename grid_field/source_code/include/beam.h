#ifndef BEAM_H
#define BEAM_H

#include "threematrix.h"
#include "threevector.h"
#include <fstream>
#include "particle.h"

class Beam
{
    int m_num_particles;
    double m_energy_spread;
    

public:
    Beam() {};

    Beam(int num_particles, double energy_spread);

    enum InitializationTypes
    {
        INITIALIZE_GAUSSIAN,
        INITIALIZE_UNIFORM_POS_DIST,
        INITIALIZE_UNIFORM_EN_DIST,
        INITIALIZE_POINT_SOURCE_GAUS,
        INITIALIZE_POINT_SOURCE_UNI,
    };

};

#endif