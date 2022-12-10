#include "beam.h"

#include "threevector.h"
#include "threematrix.h"
#include "my_functions.h"
#include <cmath>
#include <iostream>
//#include <fstream>
#include <limits>

double num_par_gaussian_multiplier = 100;
int file_downsample_factor = 1000;

std::ifstream& GotoLine(std::ifstream& file, unsigned int num){
    file.seekg(std::ios::beg);
    for(int i=0; i < num - 1; ++i){
        file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    }
    return file;
}

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
            {
                double x_init = gaussian_init(position_initial.getX(), position_spread.getX());
                double y_init = gaussian_init(position_initial.getY(), position_spread.getY());
                double z_init = gaussian_init(position_initial.getZ(), position_spread.getZ());
                m_particle.set_pos(x_init, y_init, z_init);
                break;
            }
            case PositionInitializationTypes::INITIALIZE_UNIFORM_POS:
            {
                double x_init = position_initial.getX()-position_spread.getX();
                double y_init = position_initial.getY()-position_spread.getY();
                double z_init = position_initial.getZ()-position_spread.getZ();
                m_particle.set_pos(x_init, y_init, z_init);
                break;
            }
            case PositionInitializationTypes::INITIALIZE_SCAN_POS:
            {
                //central, left, right (Y), bottom, top (Z), back, front (X)
                double x_init = position_initial.getX();
                double y_init = position_initial.getY();
                double z_init = position_initial.getZ();
                m_particle.set_pos(x_init, y_init, z_init);
                break;
            }
            case PositionInitializationTypes::INITIALIZE_INPUT_FILE_POS:
            {
                std::ifstream xpos ("../data/analysis/beaminput_position_x.txt");
                double x_init;
                xpos >> x_init;
                std::ifstream ypos ("../data/analysis/beaminput_position_y.txt");
                double y_init;
                ypos >> y_init;
                std::ifstream zpos ("../data/analysis/beaminput_position_z.txt");
                double z_init;
                zpos >> z_init;
                m_particle.set_pos(x_init, y_init, z_init);
                xpos.close();
                ypos.close();
                zpos.close();
                break;
            }
        }

        switch(energy_init)
        {
            case EnergyInitializationTypes::INITIALIZE_GAUSSIAN_EN:
            {
                /*
                 * First particle will have central energy
                 */
                // do{
                    m_particle.set_energy(m_energy_central);
                // }
                // while(m_particle.get_energy() > 1);
                break;
            }

            case EnergyInitializationTypes::INITIALIZE_UNIFORM_EN:
            {
                /*
                 * First particle will have central - FWHM energy
                 */
                m_particle.set_energy(abs(m_energy_central - m_energy_spread));
                if(m_particle.get_energy() < 1.0)
                {
                    m_particle.set_energy(1.0);
                }
                break;
            }
            case EnergyInitializationTypes::INITIALIZE_LOG_EN:
            {
                /*
                 * First particle will have central - FWHM energy
                 */
                m_particle.set_energy(abs(m_energy_central - m_energy_spread));
                break;
            }
            case EnergyInitializationTypes::INITIALIZE_INPUT_FILE_EN:
            {
                std::ifstream gamma_input ("../data/analysis/beaminput_gamma.txt");
                double file_gamma;
                gamma_input >> file_gamma;
                m_particle.set_energy(file_gamma);
                gamma_input.close();
                break;
            }
        }

        if(pos_init == PositionInitializationTypes::INITIALIZE_SCAN_POS)
        {
            diverge_init = DivergenceInitializationTypes::INITIALIZE_SCAN_DIV;
        }
        switch(diverge_init)
        {
            case DivergenceInitializationTypes::INITIALIZE_GAUSSIAN_DIV:
            {
                m_num_particles = m_num_particles*num_par_gaussian_multiplier;
                std::cerr << "num_par is now" << m_num_particles << '\n';
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

            case DivergenceInitializationTypes::INITIALIZE_UNIFORM_DIV:
                {
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

            case DivergenceInitializationTypes::INITIALIZE_SCAN_DIV:
            {
                m_num_particles = m_num_particles*7.0;
                double angle_x = m_angle_central.getX();
                double angle_y = m_angle_central.getY();
                double angle_z = m_angle_central.getZ();
                
                double p_mag   = sqrt(m_particle.get_energy()*m_particle.get_energy() - 1.0);
                double px, py, pz;
                px      = p_mag*cos(angle_x);
                py      = p_mag*cos(angle_y);
                pz      = p_mag*cos(angle_z);
                if(px == p_mag)
                {
                    py = 0.0;
                    pz = 0.0;
                }
                else if(py == p_mag)
                {
                    px = 0.0;
                    pz = 0.0;
                }
                else if(pz == p_mag)
                {
                    px = 0.0;
                    pz = 0.0;
                }
                else
                {
                    px      = sqrt((p_mag*p_mag) - (py*py) - (pz*pz) );
                }

                if((angle_x>1.57079632679 && angle_x<4.71238898038) || (angle_x>-4.71238898038 && angle_x<-1.57079632679))
                {
                    px = -px;
                }
                m_particle.set_p(px, py, pz);
                break;
            }

            case DivergenceInitializationTypes::INITIALIZE_GAMMA_SCAN_DIV:
            {
                m_num_particles = m_num_particles*7.0;
                double angle_x  = m_angle_central.getX();
                double angle_y  = m_angle_central.getY();
                double angle_z  = m_angle_central.getZ();
                
                double p_mag    = sqrt(m_particle.get_energy()*m_particle.get_energy() - 1.0);
                double py       = p_mag*cos(angle_y);
                double pz       = p_mag*cos(angle_z);
                double px       = sqrt((p_mag*p_mag) - (py*py) - (pz*pz) );
                if((angle_x>1.57079632679 && angle_x<4.71238898038) || (angle_x>-4.71238898038 && angle_x<-1.57079632679))
                {
                    px = -px;
                }
                m_particle.set_p(px, py, pz);
                break;
            }

            case DivergenceInitializationTypes::INITIALIZE_GAMMA_GAUSSIAN_DIV:
            {
                m_num_particles = m_num_particles*num_par_gaussian_multiplier;
                //FIRST ASSUMING TRAVEL ALONG X AXIS, WILL CHANGE
                m_angle_spread.setX(1/m_particle.get_energy());
                m_angle_spread.setY(1/m_particle.get_energy());
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

            case DivergenceInitializationTypes::INITIALIZE_INPUT_FILE_DIV:
            {
                std::ifstream px_input ("../data/analysis/beaminput_p_x.txt");
                double px;
                px_input >> px;
                std::ifstream py_input ("../data/analysis/beaminput_p_y.txt");
                double py;
                py_input >> py;
                std::ifstream pz_input ("../data/analysis/beaminput_p_z.txt");
                double pz;
                pz_input >> pz;
                m_particle.set_p(px, py, pz);
                px_input.close();
                py_input.close();
                pz_input.close();
                break;
            }
}
}









double Beam::gaussian_init( double initializiation_central_value, double initialization_fwhm )
{
    double value = gaussian()*initialization_fwhm + initializiation_central_value;
    return value;
}









void Beam::next_particle(int& particle_counter,
                PositionInitializationTypes pos_init, 
                EnergyInitializationTypes energy_init, 
                DivergenceInitializationTypes diverge_init)//Default values defined in beam.h
{
    particle_counter = particle_counter + 1;


    switch(pos_init)
        {
            case PositionInitializationTypes::INITIALIZE_GAUSSIAN_POS:
            {
                double x_init = gaussian_init(m_position_central.getX(), m_position_spread.getX());
                double y_init = gaussian_init(m_position_central.getY(), m_position_spread.getY());
                double z_init = gaussian_init(m_position_central.getZ(), m_position_spread.getZ());
                m_particle.set_pos(x_init, y_init, z_init);
                break;
            }

            case PositionInitializationTypes::INITIALIZE_UNIFORM_POS:
            {
                break;
            }

            case PositionInitializationTypes::INITIALIZE_SCAN_POS:
            {
                if(particle_counter%7 == 0)
                {
                    m_particle.set_pos(m_position_central.getX(), m_position_central.getY(), m_position_central.getZ());
                }
                else if(particle_counter%7 == 1)
                {
                    m_particle.set_pos(m_position_central.getX(), m_position_central.getY()-m_position_spread.getY(), m_position_central.getZ());
                }
                else if(particle_counter%7 == 2)
                {
                    m_particle.set_pos(m_position_central.getX(), m_position_central.getY()+m_position_spread.getY(), m_position_central.getZ());
                }
                else if(particle_counter%7 == 3)
                {
                    m_particle.set_pos(m_position_central.getX(), m_position_central.getY(), m_position_central.getZ()-m_position_spread.getZ());
                }
                else if(particle_counter%7 == 4)
                {
                    m_particle.set_pos(m_position_central.getX(), m_position_central.getY(), m_position_central.getZ()+m_position_spread.getZ());
                }
                else if(particle_counter%7 == 5)
                {
                    m_particle.set_pos(m_position_central.getX()-m_position_spread.getX(), m_position_central.getY(), m_position_central.getZ());
                }
                else if(particle_counter%7 == 6)
                {
                    m_particle.set_pos(m_position_central.getX()+m_position_spread.getX(), m_position_central.getY(), m_position_central.getZ());
                }
                break;
            }
            case PositionInitializationTypes::INITIALIZE_INPUT_FILE_POS:
            {
                std::ifstream xpos ("../data/analysis/beaminput_position_x.txt");
                double x_init;
                GotoLine(xpos, (particle_counter+1)*file_downsample_factor);
                xpos >> x_init;
                std::ifstream ypos ("../data/analysis/beaminput_position_y.txt");
                double y_init;
                GotoLine(ypos, (particle_counter+1)*file_downsample_factor);
                ypos >> y_init;
                std::ifstream zpos ("../data/analysis/beaminput_position_z.txt");
                double z_init;
                GotoLine(zpos, (particle_counter+1)*file_downsample_factor);
                zpos >> z_init;
                m_particle.set_pos(x_init, y_init, z_init);
                xpos.close();
                ypos.close();
                zpos.close();
                break;
            }
        }

        switch(energy_init)
        {
            case EnergyInitializationTypes::INITIALIZE_GAUSSIAN_EN:
            {  
                do{
                    m_particle.set_energy(gaussian_init(m_energy_central, m_energy_spread));
                }
                while(m_particle.get_energy() > 1);
                break;
            }

            case EnergyInitializationTypes::INITIALIZE_UNIFORM_EN:
            {
                if(pos_init != PositionInitializationTypes::INITIALIZE_SCAN_POS && diverge_init != DivergenceInitializationTypes::INITIALIZE_SCAN_DIV && diverge_init != DivergenceInitializationTypes::INITIALIZE_GAUSSIAN_DIV)
                    {
                        m_particle.set_energy(m_energy_central - m_energy_spread + particle_counter*2*m_energy_spread/m_num_particles);
                    }
                else if(diverge_init == DivergenceInitializationTypes::INITIALIZE_SCAN_DIV)
                    {
                        if((particle_counter) % 7 == 0)
                        {
                            m_particle.set_energy(m_energy_central - m_energy_spread + (particle_counter/7)*2*m_energy_spread/((m_num_particles-7)/7));
                            //std::cout << "Energy of this particle is " << m_particle.get_energy() << '\n';
                        }
                    }
                else if(diverge_init == DivergenceInitializationTypes::INITIALIZE_GAUSSIAN_DIV)
                    {
                        if((particle_counter) % int(num_par_gaussian_multiplier) == 0)
                        {
                            std::cerr << "NEXT ENERGY!\n";
                            m_particle.set_energy(m_energy_central - m_energy_spread + (particle_counter/10)*2*m_energy_spread/((m_num_particles-10)/10));
                        }
                    }
                break;
            }

            case EnergyInitializationTypes::INITIALIZE_LOG_EN:
            {
                if(diverge_init != DivergenceInitializationTypes::INITIALIZE_SCAN_DIV)
                    {
                        m_particle.set_energy(m_energy_central - m_energy_spread + particle_counter+m_energy_spread/m_num_particles);
                    }
                else if(diverge_init == DivergenceInitializationTypes::INITIALIZE_SCAN_DIV)
                    {
                        if((particle_counter) % 7 == 0)
                        {
                            m_particle.set_energy((m_energy_central - m_energy_spread)*(particle_counter/3));
                        }
                    }
                break;
            }

            case EnergyInitializationTypes::INITIALIZE_INPUT_FILE_EN:
            {
                std::ifstream gamma_input ("../data/analysis/beaminput_gamma.txt");
                double file_gamma;
                GotoLine(gamma_input, (particle_counter+1)*file_downsample_factor);
                gamma_input >> file_gamma;
                m_particle.set_energy(file_gamma);
                gamma_input.close();
                break;
            }
        }

        switch(diverge_init)
        {
            case DivergenceInitializationTypes::INITIALIZE_GAUSSIAN_DIV:
            {
                double angle_x = gaussian_init(m_angle_central.getX(), m_angle_spread.getX());
                double angle_y = gaussian_init(m_angle_central.getY(), m_angle_spread.getY());
                double angle_z = gaussian_init(m_angle_central.getZ(), m_angle_spread.getZ());
                
                double p_mag   = sqrt(m_particle.get_energy()*m_particle.get_energy() - 1.0);
                //double px      = p_mag*cos(angle_x);
                double py      = p_mag*cos(angle_y);
                double pz      = p_mag*cos(angle_z);
                double px      = sqrt((p_mag*p_mag) - (py*py) - (pz*pz) );
                if((angle_x>1.57079632679 && angle_x<4.71238898038) || (angle_x>-4.71238898038 && angle_x<-1.57079632679))
                {
                    px = -px;
                }
                m_particle.set_p(px, py, pz);
                break;
            }

            case DivergenceInitializationTypes::INITIALIZE_UNIFORM_DIV:
                {break;}

            case DivergenceInitializationTypes::INITIALIZE_SCAN_DIV:
            {
                //std::cerr << "In scan divergence initialization\n";
                //scans central, left, right, top, bottom of divergence extremes
                //std::cerr << "(particle_counter)%7 = " << (particle_counter)%7 << "\n";
                //std::cerr << "Particle energy is " << m_particle.get_energy() << "\n";
                if((particle_counter)%7 == 0)
                {//CENTRAL TRAJECTORY
                        double angle_x = m_angle_central.getX();
                        double angle_y = m_angle_central.getY();
                        double angle_z = m_angle_central.getZ();
                        
                        double p_mag   = sqrt(m_particle.get_energy()*m_particle.get_energy() - 1.0);
                        double px, py, pz;
                        px      = p_mag*cos(angle_x);
                        py      = p_mag*cos(angle_y);
                        pz      = p_mag*cos(angle_z);
                        if(px == p_mag)
                        {
                            py = 0.0;
                            pz = 0.0;
                        }
                        else if(py == p_mag)
                        {
                            px = 0.0;
                            pz = 0.0;
                        }
                        else if(pz == p_mag)
                        {
                            px = 0.0;
                            pz = 0.0;
                        }
                        else
                        {
                            px      = sqrt((p_mag*p_mag) - (py*py) - (pz*pz) );
                        }


                        if((angle_x>1.57079632679 && angle_x<4.71238898038) || (angle_x>-4.71238898038 && angle_x<-1.57079632679))
                        {
                            px = -px;
                        }
                        if(angle_y < 0)
                        {
                            py = -py;
                        }
                        if(angle_z < 0)
                        {
                            pz = -pz;
                        }
                        m_particle.set_p(px, py, pz);
                        //std::cerr << "Energy after setting p is " << m_particle.get_energy() << std::endl;
                }
                else if((particle_counter) % 7 == 1)// && particle_counter+1 % 4 != 0)
                    {//TO THE SIDE DIVERGENCE TRAJECTORY
                    //subtracting angle from 2 dimensions b/c it is a 3D grid that we're splitting 1 angle into
                        double angle_x = m_angle_central.getX() - m_angle_spread.getY()/2.0;
                        double angle_y = m_angle_central.getY() + m_angle_spread.getY()/2.0;
                        double angle_z = m_angle_central.getZ();
                        
                        double p_mag   = sqrt(m_particle.get_energy()*m_particle.get_energy() - 1.0);
                        double px, py, pz;
                        px      = p_mag*cos(angle_x);
                        py      = p_mag*cos(angle_y);
                        pz      = p_mag*cos(angle_z);
                        if(px == p_mag)
                        {
                            py = 0.0;
                            pz = 0.0;
                        }
                        else if(py == p_mag)
                        {
                            px = 0.0;
                            pz = 0.0;
                        }
                        else if(pz == p_mag)
                        {
                            px = 0.0;
                            pz = 0.0;
                        }
                        else
                        {
                            px      = sqrt((p_mag*p_mag) - (py*py) - (pz*pz) );
                        }


                        if((angle_x>1.57079632679 && angle_x<4.71238898038) || (angle_x>-4.71238898038 && angle_x<-1.57079632679))
                        {
                            px = -px;
                        }
                        if(angle_y < 0)
                        {
                            py = -py;
                        }
                        if(angle_z < 0)
                        {
                            pz = -pz;
                        }
                        m_particle.set_p(px, py, pz);
                    }
                else if(particle_counter % 7 == 2)// && particle_counter+1 % 6 != 0)
                    {
                        double angle_x = m_angle_central.getX() + m_angle_spread.getY()/2.0;
                        double angle_y = m_angle_central.getY() - m_angle_spread.getY()/2.0;
                        double angle_z = m_angle_central.getZ();
                        
                        double p_mag   = sqrt(m_particle.get_energy()*m_particle.get_energy() - 1.0);
                        double px, py, pz;
                        px      = p_mag*cos(angle_x);
                        py      = p_mag*cos(angle_y);
                        pz      = p_mag*cos(angle_z);
                        if(px == p_mag)
                        {
                            py = 0.0;
                            pz = 0.0;
                        }
                        else if(py == p_mag)
                        {
                            px = 0.0;
                            pz = 0.0;
                        }
                        else if(pz == p_mag)
                        {
                            px = 0.0;
                            pz = 0.0;
                        }
                        else{
                            px      = sqrt((p_mag*p_mag) - (py*py) - (pz*pz) );
                        }


                        if((angle_x>1.57079632679 && angle_x<4.71238898038) || (angle_x>-4.71238898038 && angle_x<-1.57079632679))
                        {
                            px = -px;
                        }
                        if(angle_y < 0)
                        {
                            py = -py;
                        }
                        if(angle_z < 0)
                        {
                            pz = -pz;
                        }
                        m_particle.set_p(px, py, pz);
                    }
                else if((particle_counter) % 7 == 3)
                    {
                        double angle_x = m_angle_central.getX() - m_angle_spread.getZ()/2.0;
                        double angle_y = m_angle_central.getY();
                        double angle_z = m_angle_central.getZ() + m_angle_spread.getZ()/2.0;
                     
                        double p_mag   = sqrt(m_particle.get_energy()*m_particle.get_energy() - 1.0);
                        double px, py, pz;
                        px      = p_mag*cos(angle_x);
                        py      = p_mag*cos(angle_y);
                        pz      = p_mag*cos(angle_z);
                        if(px == p_mag)
                        {
                            py = 0.0;
                            pz = 0.0;
                        }
                        else if(py == p_mag)
                        {
                            px = 0.0;
                            pz = 0.0;
                        }
                        else if(pz == p_mag)
                        {
                            px = 0.0;
                            pz = 0.0;
                        }
                        else
                        {
                            px      = sqrt((p_mag*p_mag) - (py*py) - (pz*pz) );
                        }
                        
                        
                        if((angle_x>1.57079632679 && angle_x<4.71238898038) || (angle_x>-4.71238898038 && angle_x<-1.57079632679))
                        {
                            px = -px;
                        }
                        if(angle_y < 0)
                        {
                            py = -py;
                        }
                        if(angle_z < 0)
                        {
                            pz = -pz;
                        }
                        m_particle.set_p(px, py, pz);
                    }
                else if((particle_counter) % 7 == 4)
                    {
                        double angle_x = m_angle_central.getX() + m_angle_spread.getZ()/2.0;
                        double angle_y = m_angle_central.getY();
                        double angle_z = m_angle_central.getZ() - m_angle_spread.getZ()/2.0;
                     
                        double p_mag   = sqrt(m_particle.get_energy()*m_particle.get_energy() - 1.0);
                        double px, py, pz;
                        px      = p_mag*cos(angle_x);
                        py      = p_mag*cos(angle_y);
                        pz      = p_mag*cos(angle_z);
                        if(px == p_mag)
                        {
                            py = 0.0;
                            pz = 0.0;
                        }
                        else if(py == p_mag)
                        {
                            px = 0.0;
                            pz = 0.0;
                        }
                        else if(pz == p_mag)
                        {
                            px = 0.0;
                            pz = 0.0;
                        }
                        else
                        {
                            px      = sqrt((p_mag*p_mag) - (py*py) - (pz*pz) );
                        }


                        if((angle_x>1.57079632679 && angle_x<4.71238898038) || (angle_x>-4.71238898038 && angle_x<-1.57079632679))
                        {
                            px = -px;
                        }
                        if(angle_y < 0)
                        {
                            py = -py;
                        }
                        if(angle_z < 0)
                        {
                            pz = -pz;
                        }
                        //std::cerr << "(px, py, pz) = (" << ThreeVec(px,py,pz) << ")\n";
                        m_particle.set_p(px, py, pz);
                    }
                else if((particle_counter) % 7 == 5)
                    {
                        double angle_x = m_angle_central.getX() + m_angle_spread.getX();
                        double angle_y = m_angle_central.getY() - m_angle_spread.getX();
                        double angle_z = m_angle_central.getZ();
                     
                        double p_mag   = sqrt(m_particle.get_energy()*m_particle.get_energy() - 1.0);
                        //std::cerr << "p_mag is " << p_mag << "\n";
                        double px, py, pz;
                        px      = p_mag*cos(angle_x);
                        py      = p_mag*cos(angle_y);
                        pz      = p_mag*cos(angle_z);
                        if(px == p_mag)
                        {
                            py = 0.0f;
                            pz = 0.0f;
                        }
                        else if(py == p_mag)
                        {
                            px = 0.0f;
                            pz = 0.0f;
                        }
                        else if(pz == p_mag)
                        {
                            px = 0.0f;
                            pz = 0.0f;
                        }
                        else
                        {
                            pz      = sqrt((p_mag*p_mag) - (py*py) - (px*px) );
                        }
                        if((px*px + py*py + pz*pz) > p_mag*p_mag)
                        {
                            std::cout << "momenta greater than p_mag: " << (px*px + py*py + pz*pz) << " > " << p_mag*p_mag << '\n';
                            double diff = (px*px + py*py + pz*pz) - p_mag*p_mag;
                            std::cout << "    Greater by difference of " << diff << '\n';
                            std::cout << "        0.005*p_mag*p_mag = " << 0.005*p_mag*p_mag;
                            if(diff < 0.005*p_mag*p_mag){
                                std::cerr << "Correction to pz";
                                pz      = sqrt((p_mag*p_mag) - (py*py) - (px*px) + diff);
                            }
                        }
                        //double pz      = p_mag*cos(angle_z);
                        //std::cerr << "px = " << px << ", py = " << py << "\n"; 
                        //std::cerr << "p_mag^2 = " << p_mag*p_mag << "\n";
                        
                        if(angle_x < 0)
                        {
                            px = -px;
                        }
                        if(angle_y < 0)
                        {
                            py = -py;
                        }
                        if((angle_z>1.57079632679 && angle_z<4.71238898038) || (angle_z>-4.71238898038 && angle_z<-1.57079632679))
                        {
                            pz = -pz;
                        }
                        //std::cerr << "(px, py, pz) = " << ThreeVec(px,py,pz) << "\n";
                        m_particle.set_p(px, py, pz);
                    }
                else if((particle_counter) % 7 == 6)
                    {
                        double angle_x = m_angle_central.getX() - m_angle_spread.getX()/2.0;
                        double angle_y = m_angle_central.getY() + m_angle_spread.getX();
                        double angle_z = m_angle_central.getZ();
                     
                        double p_mag   = sqrt(m_particle.get_energy()*m_particle.get_energy() - 1.0);
                        double px, py, pz;
                        px      = p_mag*cos(angle_x);
                        py      = p_mag*cos(angle_y);
                        pz      = p_mag*cos(angle_z);
                        if(px == p_mag)
                        {
                            py = 0.0;
                            pz = 0.0;
                        }
                        else if(py == p_mag)
                        {
                            px = 0.0;
                            pz = 0.0;
                        }
                        else if(pz == p_mag)
                        {
                            px = 0.0;
                            pz = 0.0;
                        }
                        else
                        {
                            pz = sqrt((p_mag*p_mag) - (py*py) - (px*px) );
                        }
                        if((px*px + py*py + pz*pz) > p_mag*p_mag)
                        {
                            //std::cout << "momenta greater than p_mag: " << (px*px + py*py + pz*pz) << " > " << p_mag*p_mag << '\n';
                            double diff = (px*px + py*py + pz*pz) - p_mag*p_mag;
                            //std::cout << "    Greater by difference of " << diff << '\n';
                            if(diff < 0.005*p_mag*p_mag){
                                pz      = sqrt((p_mag*p_mag) - (py*py) - (px*px) + diff);
                            }
                        }


                        if(angle_x < 0)
                        {
                            px = -px;
                        }
                        if(angle_y < 0)
                        {
                            py = -py;
                        }
                        if((angle_z>1.57079632679 && angle_z<4.71238898038) || (angle_z>-4.71238898038 && angle_z<-1.57079632679))
                        {
                            pz = -pz;
                        }
                        m_particle.set_p(px, py, pz);
                    }
                    break;
            }

            case DivergenceInitializationTypes::INITIALIZE_GAMMA_SCAN_DIV:
            {
                if((particle_counter)%7 == 0)
                {
                        double angle_x = m_angle_central.getX();
                        double angle_y = m_angle_central.getY();
                        double angle_z = m_angle_central.getZ();
                        
                        double p_mag   = sqrt(m_particle.get_energy()*m_particle.get_energy() - 1.0);
                        //std::cerr << "Energy to set p is " << m_particle.get_energy() << std::endl;
                        double py      = p_mag*cos(angle_y);
                        double pz      = p_mag*cos(angle_z);
                        double px      = sqrt((p_mag*p_mag) - (py*py) - (pz*pz) );
                        if(angle_x < 0)
                        {
                            px = -px;
                        }
                        if(angle_y < 0)
                        {
                            py = -py;
                        }
                        if(angle_z < 0)
                        {
                            pz = -pz;
                        }
                        m_particle.set_p(px, py, pz);
                        //std::cerr << "Energy after setting p is " << m_particle.get_energy() << std::endl;
                }
                else if((particle_counter) % 7 == 1)// && particle_counter+1 % 4 != 0)
                    {
                        double angle_x = m_angle_central.getX();
                        double angle_y = m_angle_central.getY() - m_particle.get_energy()/2.0;
                        double angle_z = m_angle_central.getZ();
                        
                        double p_mag   = sqrt(m_particle.get_energy()*m_particle.get_energy() - 1.0);
                        double py      = p_mag*cos(angle_y);
                        double pz      = p_mag*cos(angle_z);
                        double px      = sqrt((p_mag*p_mag) - (py*py) - (pz*pz) );
                        if(angle_x < 0)
                        {
                            px = -px;
                        }
                        if(angle_y < 0)
                        {
                            py = -py;
                        }
                        if(angle_z < 0)
                        {
                            pz = -pz;
                        }
                        m_particle.set_p(px, py, pz);
                    }
                else if(particle_counter % 7 == 2)// && particle_counter+1 % 6 != 0)
                    {
                        double angle_x = m_angle_central.getX();
                        double angle_y = m_angle_central.getY() + m_particle.get_energy()/2.0;
                        double angle_z = m_angle_central.getZ();
                        
                        double p_mag   = sqrt(m_particle.get_energy()*m_particle.get_energy() - 1.0);
                        double py      = p_mag*cos(angle_y);
                        double pz      = p_mag*cos(angle_z);
                        double px      = sqrt((p_mag*p_mag) - (py*py) - (pz*pz) );
                        if(angle_x < 0)
                        {
                            px = -px;
                        }
                        if(angle_y < 0)
                        {
                            py = -py;
                        }
                        if(angle_z < 0)
                        {
                            pz = -pz;
                        }
                        m_particle.set_p(px, py, pz);
                    }
                else if((particle_counter) % 7 == 3)
                    {
                        double angle_x = m_angle_central.getX();
                        double angle_y = m_angle_central.getY();
                        double angle_z = m_angle_central.getZ() - m_particle.get_energy()/2.0;
                     
                        double p_mag   = sqrt(m_particle.get_energy()*m_particle.get_energy() - 1.0);
                        //double px      = p_mag*cos(angle_x);
                        double py      = p_mag*cos(angle_y);
                        double pz      = p_mag*cos(angle_z);
                        double px      = sqrt((p_mag*p_mag) - (py*py) - (pz*pz) );
                        if(angle_x < 0)
                        {
                            px = -px;
                        }
                        if(angle_y < 0)
                        {
                            py = -py;
                        }
                        if(angle_z < 0)
                        {
                            pz = -pz;
                        }
                        m_particle.set_p(px, py, pz);
                    }
                else if((particle_counter) % 7 == 4)
                    {
                        double angle_x = m_angle_central.getX();
                        double angle_y = m_angle_central.getY();
                        double angle_z = m_angle_central.getZ() + m_particle.get_energy()/2.0;
                     
                        double p_mag   = sqrt(m_particle.get_energy()*m_particle.get_energy() - 1.0);
                        //double px      = p_mag*cos(angle_x);
                        double py      = p_mag*cos(angle_y);
                        double pz      = p_mag*cos(angle_z);
                        double px      = sqrt((p_mag*p_mag) - (py*py) - (pz*pz) );
                        if(angle_x < 0)
                        {
                            px = -px;
                        }
                        if(angle_y < 0)
                        {
                            py = -py;
                        }
                        if(angle_z < 0)
                        {
                            pz = -pz;
                        }
                        m_particle.set_p(px, py, pz);
                    }
                else if((particle_counter) % 7 == 5)
                    {
                        double angle_x = m_angle_central.getX() - m_particle.get_energy()/2.0;
                        double angle_y = m_angle_central.getY();
                        double angle_z = m_angle_central.getZ();
                     
                        double p_mag   = sqrt(m_particle.get_energy()*m_particle.get_energy() - 1.0);
                        double px      = p_mag*cos(angle_x);
                        double py      = p_mag*cos(angle_y);
                        //double pz      = p_mag*cos(angle_z);
                        double pz      = sqrt((p_mag*p_mag) - (py*py) - (px*px) );
                        if(angle_x < 0)
                        {
                            px = -px;
                        }
                        if(angle_y < 0)
                        {
                            py = -py;
                        }
                        if(angle_z < 0)
                        {
                            pz = -pz;
                        }
                        m_particle.set_p(px, py, pz);
                    }
                else if((particle_counter) % 7 == 6)
                    {
                        double angle_x = m_angle_central.getX() + m_particle.get_energy()/2.0;
                        double angle_y = m_angle_central.getY();
                        double angle_z = m_angle_central.getZ();
                     
                        double p_mag   = sqrt(m_particle.get_energy()*m_particle.get_energy() - 1.0);
                        double px      = p_mag*cos(angle_x);
                        double py      = p_mag*cos(angle_y);
                        //double pz      = p_mag*cos(angle_z);
                        double pz      = sqrt((p_mag*p_mag) - (py*py) - (px*px) );
                        if(angle_x < 0)
                        {
                            px = -px;
                        }
                        if(angle_y < 0)
                        {
                            py = -py;
                        }
                        if(angle_z < 0)
                        {
                            pz = -pz;
                        }
                        m_particle.set_p(px, py, pz);
                    }
                break;
            }

            case DivergenceInitializationTypes::INITIALIZE_GAMMA_GAUSSIAN_DIV:
            {
                double angle_x = gaussian_init(m_angle_central.getX(), m_particle.get_energy());
                double angle_y = gaussian_init(m_angle_central.getY(), m_particle.get_energy());
                double angle_z = gaussian_init(m_angle_central.getZ(), m_particle.get_energy());
                
                double p_mag   = sqrt(m_particle.get_energy()*m_particle.get_energy() - 1.0);
                //double px      = p_mag*cos(angle_x);
                double py      = p_mag*cos(angle_y);
                double pz      = p_mag*cos(angle_z);
                double px      = sqrt((p_mag*p_mag) - (py*py) - (pz*pz) );
                m_particle.set_p(px, py, pz);
                break;
            }

            case DivergenceInitializationTypes::INITIALIZE_INPUT_FILE_DIV:
            {
                std::ifstream px_input ("../data/analysis/beaminput_p_x.txt");
                double px;
                GotoLine(px_input, (particle_counter+1)*file_downsample_factor);
                px_input >> px;
                std::ifstream py_input ("../data/analysis/beaminput_p_y.txt");
                double py;
                GotoLine(py_input, (particle_counter+1)*file_downsample_factor);
                py_input >> py;
                std::ifstream pz_input ("../data/analysis/beaminput_p_z.txt");
                double pz;
                GotoLine(pz_input, (particle_counter+1)*file_downsample_factor);
                pz_input >> pz;
                m_particle.set_p(px, py, pz);
                px_input.close();
                py_input.close();
                pz_input.close();
                break;
            }
                
        }
}