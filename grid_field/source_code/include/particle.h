#ifndef PARTICLE_H
#define PARTICLE_H

#include "threematrix.h"
#include "threevector.h"

#include <fstream>

class Particle
{
    ThreeVec m_pos;
    ThreeVec m_p;                  //in p/mc
    ThreeVec m_vel;
    double m_energy;
    double m_charge;
    double m_rest_mass;
    double *m_time;
    const double m_rest_en = 1;     //E0/E0

public:
    Particle() {};

    Particle(ThreeVec pos, ThreeVec momentum, int charge, double mass, double& time, std::ofstream& out_time, 
            std::ofstream& out_posx, std::ofstream& out_posy, std::ofstream& out_posz,
            std::ofstream& out_velx, std::ofstream& out_vely, std::ofstream& out_velz, std::ofstream& out_energy);

    /*Particle(ThreeVec pos, int charge, ThreeVec momentum, double& time, std::ofstream& out_time, 
            std::ofstream& out_posx, std::ofstream& out_posy, std::ofstream& out_posz,
            std::ofstream& out_velx, std::ofstream& out_vely, std::ofstream& out_velz, std::ofstream& out_energy);
    */
    ThreeVec get_pos()
    {
        return m_pos;
    }

    double get_charge()
    {
        return m_charge;
    }

    double get_mass()
    {
        return m_rest_mass;
    }

    double get_pos(int i)
    {
        return m_pos.get(i);
    }

    ThreeVec get_p()
    {
        return m_p;
    }

    double get_p(int i)
    {
        return m_p.get(i);
    }

    double get_energy()
    {
        return m_energy;
    }

    double get_time()
    {
        return *m_time;
    }

    ThreeVec get_vel()
    {
        return m_vel;
    }

    double get_vel(int ii)
    {
        return m_vel.get(ii);
    }

    void set_pos(ThreeVec pos);
    void set_pos(int index, double value);
    void set_pos(double x, double y, double z);
    void set_energy(double energy_t);
    //void set_energy(int index_t, double value_t);
    void set_p(ThreeVec momentum);
    void set_p(int index, double value);
    void set_p(double x, double y, double z);
    void set_time(double& time);

    void set_vel(ThreeVec vel);
    void set_vel(int index, double value);

    void set_charge(double charge = -1.0);
    void set_mass(double mass = 1.0);

    

    std::ofstream *m_out_time;
    std::ofstream *m_out_posx;
    std::ofstream *m_out_posy;
    std::ofstream *m_out_posz;
    std::ofstream *m_out_px;
    std::ofstream *m_out_py;
    std::ofstream *m_out_pz; 
    std::ofstream *m_out_energy;

    void set_outfiles( std::ofstream& out_time, 
            std::ofstream& out_posx, std::ofstream& out_posy, std::ofstream& out_posz,
            std::ofstream& out_velx, std::ofstream& out_vely, std::ofstream& out_velz, std::ofstream& out_energy)
        {
            m_out_time   = &out_time;
            m_out_posx   = &out_posx;
            m_out_posy   = &out_posy;
            m_out_posz   = &out_posz;
            m_out_px     = &out_velx;
            m_out_py     = &out_vely;
            m_out_pz     = &out_velz;
            m_out_energy = &out_energy;
        }

};

#endif // !PARTICLE_H
