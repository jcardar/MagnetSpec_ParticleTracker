#include "my_functions.h"

#include <iostream>
#include <fstream>
#include "threevector.h"
#include "threematrix.h"
#include "magnet.h"
#include "screen.h"
#include "particle.h"
#include <cmath>
#include <cfloat>
#include <string>
#include <vector>
#include "math.h"
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void outfile_tab(double& time, std::ofstream& out_time, ThreeVec pos, std::ofstream& out_xpos, std::ofstream& out_ypos, std::ofstream& out_zpos, ThreeVec vel, std::ofstream& out_vx, std::ofstream& out_vy, std::ofstream& out_vz)
{
    out_time << time << "\t";
    out_xpos << pos.get(0) << "\t";
    out_ypos << pos.get(1) << "\t";
    out_zpos << pos.get(2) << "\t";
    out_vx << vel.get(0) << "\t";
    out_vy << vel.get(1) << "\t";
    out_vz << vel.get(2) << "\t";
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void outfile_part_writeAndTab(Particle& particle)
{
    *(particle.m_out_time) << (particle.get_time()) << "\t";
    *(particle.m_out_posx) << (particle.get_pos(0)) << "\t";
    *(particle.m_out_posy) << (particle.get_pos(1)) << "\t";
    *(particle.m_out_posz) << (particle.get_pos(2)) << "\t";
    *(particle.m_out_px) << (particle.get_p(0)) << "\t";
    *(particle.m_out_py) << (particle.get_p(1)) << "\t";
    *(particle.m_out_pz) << (particle.get_p(2)) << "\t";
    *(particle.m_out_energy) << (particle.get_energy()) << "\t";
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void outfile_part_writeAndComma(Particle& particle)
{
//        std::cerr << "outfile prompt entered.\n";
    *(particle.m_out_time) << (particle.get_time()) << ",";
//        std::cerr << "particle time written\n";
    *(particle.m_out_posx) << (particle.get_pos(0)) << ",";
//       std::cerr << "particle pos0 written\n";
    *(particle.m_out_posy) << (particle.get_pos(1)) << ",";
//       std::cerr << "particle pos1 written\n";
    *(particle.m_out_posz) << (particle.get_pos(2)) << ",";
//       std::cerr << "particle pos2 written\n";
    *(particle.m_out_px) << (particle.get_p(0)) << ",";
//       std::cerr << "particle mom0 written\n";
    *(particle.m_out_py) << (particle.get_p(1)) << ",";
//       std::cerr << "particle mom1 written\n";
    *(particle.m_out_pz) << (particle.get_p(2)) << ",";
//       std::cerr << "particle mom2 written\n";
    *(particle.m_out_energy) << particle.get_energy() << ",";
//       std::cerr << "particle energy written\n";
}

void outfile_part_commaAndWrite(Particle& particle)
{
//        std::cerr << "outfile prompt entered.\n";
    *(particle.m_out_time) << "," << (particle.get_time());
//        std::cerr << "particle time written\n";
    *(particle.m_out_posx)  << "," << (particle.get_pos(0));
//       std::cerr << "particle pos0 written\n";
    *(particle.m_out_posy) << "," << (particle.get_pos(1));
//       std::cerr << "particle pos1 written\n";
    *(particle.m_out_posz) << "," << (particle.get_pos(2));
//       std::cerr << "particle pos2 written\n";
    *(particle.m_out_px) << "," << (particle.get_p(0));
//       std::cerr << "particle mom0 written\n";
    *(particle.m_out_py) << "," << (particle.get_p(1));
//       std::cerr << "particle mom1 written\n";
    *(particle.m_out_pz) << "," << (particle.get_p(2)) ;
//       std::cerr << "particle mom2 written\n";
    *(particle.m_out_energy) << "," << particle.get_energy();
//       std::cerr << "particle energy written\n";
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void outfile_part_write(Particle& particle)
{
    *(particle.m_out_time) << (particle.get_time());
    *(particle.m_out_posx) << (particle.get_pos(0));
    *(particle.m_out_posy) << (particle.get_pos(1));
    *(particle.m_out_posz) << (particle.get_pos(2));
    *(particle.m_out_px) << (particle.get_p(0));
    *(particle.m_out_py) << (particle.get_p(1));
    *(particle.m_out_pz) << (particle.get_p(2));
    *(particle.m_out_energy) << particle.get_energy();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void outfile_newline(std::ofstream& out_time, std::ofstream& out_xpos, std::ofstream& out_ypos, std::ofstream& out_zpos, std::ofstream& out_vx, std::ofstream& out_vy, std::ofstream& out_vz)
{
    out_time << "\n";
    out_xpos << "\n";
    out_ypos << "\n";
    out_zpos << "\n";
    out_vx << "\n";
    out_vy << "\n";
    out_vz << "\n";
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void outfile_part_newline(Particle& particle)
{
    *(particle.m_out_time) << "\n";
    *(particle.m_out_posx) << "\n";
    *(particle.m_out_posy) << "\n";
    *(particle.m_out_posz) << "\n";
    *(particle.m_out_px) << "\n";
    *(particle.m_out_py) << "\n";
    *(particle.m_out_pz) << "\n";
    *(particle.m_out_energy) << "\n";
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void outfile_part_comma(Particle& particle)
{
    *(particle.m_out_time) << ","; 
    *(particle.m_out_posx) << ",";
    //std::cerr << "wrote comma on x\n";
    *(particle.m_out_posy) << ",";
    //std::cerr << "wrote comma on y\n";
    *(particle.m_out_posz) << ",";
    //std::cerr << "wrote comma on z\n";
    *(particle.m_out_px) << ",";
    *(particle.m_out_py) << ",";
    *(particle.m_out_pz) << ",";
    *(particle.m_out_energy) << ",";
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void outfile_uniform_magnet(Magnet& magnet, int counter)
{
    if(counter==0)
    {
        *(magnet.m_out_magnet) << "Num ," << "Bz(norm)," << "mag_posx," << "mag_posy," << "mag_posz," << "length," << "width," << "height" << "\n";
    }
    *(magnet.m_out_magnet) << (++counter) << "," << magnet.get_B0(2) << "," 
                            << magnet.get_pos(0) << "," << magnet.get_pos(1) << "," << magnet.get_pos(2) << ","
                            << magnet.get_length() << "," << magnet.get_width() << "," << magnet.get_height() << '\n';
} 



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void outfile_screen_single(Screen& screen, int counter)
{
    if(counter==0)
    {
        *(screen.m_out_screen) << "Num," << "screen_low_energy_edgex," << "screen_low_energy_edgey," << "screen_low_energy_edgez," << "degrees about x-axis," << "degrees about y-axis," << "degrees about z-axis," << "length," << "height" << "\n";
    }
    *(screen.m_out_screen) << (++counter) << ","
                            << screen.get_pos(0) << "," << screen.get_pos(1) << "," << screen.get_pos(2) << ","
                            << screen.get_angle_about_x('d') << "," << screen.get_angle_about_y('d') << "," << screen.get_angle_about_z('d') << "," << screen.get_length() << "," << screen.get_height() << '\n';
} 


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Gaussian random distribution function
double gaussian()
{
	double gaussian_width=3.0; // +/- 3 sigma range
	double x, y;
	do
	{
		x = (2.0*rand()/RAND_MAX-1.0)*gaussian_width;
		y = exp(-x*x*0.5);
	}
	while (1.0*rand()/RAND_MAX > y);
	return x;
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double uniform_dist_single(int num_par_t, double vel0_t, double radius_v0_t, int &counter_t)
{
    //From v_min to v_max in evenly divided intervals, iterates for each particle
    double del_vel = (2*radius_v0_t)/(static_cast<double>(num_par_t));

    double counter_d = static_cast<double>(counter_t);
    counter_t++;

    return (vel0_t-radius_v0_t + (del_vel * counter_d));
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void uniform_en_dist(double& initial_x, double& initial_y, double& initial_z, double& initial_enx, double& initial_eny, double& initial_enz, double length_before, int* posy_counter, int* posz_counter, ThreeVec r0, ThreeVec radius_r0, ThreeVec energy0, int num_par)
{
    initial_x = 0.0;
    initial_y = uniform_dist_single(num_par, r0.getY(), radius_r0.getY(), *posy_counter);
    initial_z = uniform_dist_single(num_par, r0.getZ(), radius_r0.getZ(), *posz_counter);
    double tot_energy = sqrt(energy0.getX()*energy0.getX() + energy0.getY()*energy0.getY() + energy0.getZ()*energy0.getZ());
    double length_on_screen = sqrt(initial_y*initial_y + initial_z*initial_z);
    double half_divergence;
    if(length_on_screen != 0)
    { half_divergence = atan(length_before/length_on_screen);}
    else { half_divergence = 0; }
    double angle_on_yz_plane;
    if(initial_y != 0)
    { angle_on_yz_plane = atan(initial_z/initial_y); }
    else { angle_on_yz_plane = M_PI/2.0; }
    initial_enx = tot_energy*cos(half_divergence);
    double remaining_energy = sqrt(tot_energy*tot_energy - initial_enx*initial_enx);
    if(initial_y >= 0.0)
    { initial_eny = remaining_energy*cos(angle_on_yz_plane);}
    else
    { initial_eny = -remaining_energy*cos(angle_on_yz_plane);}
    if(initial_z >= 0.0)
    { initial_enz = remaining_energy*sin(angle_on_yz_plane);}
    else
    { initial_enz = -remaining_energy*sin(angle_on_yz_plane); }            
}


void uniform_pos_dist(double& initial_x, double& initial_y, double& initial_z, double& initial_enx, double& initial_eny, double& initial_enz, double length_before,int* posx_counter, int* posy_counter, int* posz_counter, ThreeVec p0, ThreeVec radius_p0, ThreeVec r0, int num_par, bool point_source = false, bool dist_x = true, bool dist_y = false, bool dist_z = false)
{
    initial_x = r0.getX();
    initial_y = r0.getY();
    initial_z = r0.getZ();
    initial_enx = 0.0;
    initial_eny = 0.0;
    initial_enz = 0.0;
    if(dist_x) { initial_enx = uniform_dist_single(num_par, p0.getY(), radius_p0.getY(), *posx_counter); }
    if(dist_y) { initial_eny = uniform_dist_single(num_par, p0.getZ(), radius_p0.getZ(), *posy_counter); }
    if(dist_z) { initial_enz = uniform_dist_single(num_par, p0.getZ(), radius_p0.getZ(), *posz_counter); }
    if(point_source) { num_par = num_par/5; }
    double tot_energy = sqrt(p0.getX()*p0.getX() + p0.getY()*p0.getY() + p0.getZ()*p0.getZ());
    double length_on_screen = sqrt(initial_y*initial_y + initial_z*initial_z);
    double half_divergence;
    if(length_on_screen != 0)
    { half_divergence = atan(length_before/length_on_screen);}
    else { half_divergence = 0; }
    double angle_on_yz_plane;
    if(initial_y != 0)
    { angle_on_yz_plane = atan(initial_z/initial_y); }
    else { angle_on_yz_plane = M_PI/2.0; }
    initial_enx = tot_energy*cos(half_divergence);
    double remaining_energy = sqrt(tot_energy*tot_energy - initial_enx*initial_enx);
    if(initial_y >= 0.0)
    { initial_eny = remaining_energy*cos(angle_on_yz_plane);}
    else
    { initial_eny = -remaining_energy*cos(angle_on_yz_plane);}
    if(initial_z >= 0.0)
    { initial_enz = remaining_energy*sin(angle_on_yz_plane);}
    else
    { initial_enz = -remaining_energy*sin(angle_on_yz_plane); }            
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void stepThroughMagnet_LeapUniformB(Particle &electron, ThreeVec vn_plus, ThreeVec vn_minus, ThreeVec B0, ThreeVec rn_plus, ThreeVec rn_minus, double& time, const double& del_time, double& width_bmap, double& length_bmap)
{
    bool check;
    vn_minus = (electron.get_vel() ) - (( (electron.get_vel())^B0 )*del_time);
    //vn_plus;

    rn_minus = (electron.get_pos() ) - ((electron.get_vel())*del_time);
    //rn_plus;
    do
        {
            vn_plus = vn_minus + ((electron.get_vel())^B0)*2*del_time;    //leapfrog algarithm vn for uniform B0
            rn_plus = rn_minus + (electron.get_vel())*2*del_time;         //leapfrog algarithm rn for uniform B0

            rn_minus = electron.get_pos();                          //iterate xn_minus up a timestep
            electron.set_pos(rn_plus);                              //do the same to xn

            vn_minus = (electron.get_vel());                        //iterate vn_minus up a timestep
            electron.set_vel(vn_plus);                              //do the same to vn

            time += del_time;
            
            check = (electron.get_pos(0) >= 0.0) && (electron.get_pos(0) <= length_bmap) && (electron.get_pos(1) >= (-width_bmap/2.0)) && (electron.get_pos(1) <= (width_bmap/2.0));
            if(check)
                { outfile_part_writeAndComma(electron); }
            else if(!(check))
                { outfile_part_write(electron); }            
            
        } while(check);
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void stepThroughMagnet_Leap(Particle &electron, Magnet &magnet, double& time, const double &del_time)
{
    bool check;
    ThreeVec vn_minus = (electron.get_vel() ) - (( (electron.get_vel())^(magnet.get_B0()) )*del_time);
    ThreeVec vn_plus;

    ThreeVec rn_minus = (electron.get_pos() ) - ((electron.get_vel())*del_time);
    ThreeVec rn_plus;
    do
        {
            vn_plus = vn_minus + ((electron.get_vel())^(magnet.get_B0()))*2*del_time;    //leapfrog algarithm vn for uniform B0
            rn_plus = rn_minus + (electron.get_vel())*2*del_time;         //leapfrog algarithm rn for uniform B0

            rn_minus = electron.get_pos();                          //iterate xn_minus up a timestep
            electron.set_pos(rn_plus);                              //do the same to xn

            vn_minus = (electron.get_vel());                        //iterate vn_minus up a timestep
            electron.set_vel(vn_plus);                              //do the same to vn

            time += del_time;

            electron.set_time(time);

            check = ((electron.get_pos(0) >= (magnet.get_pos(0))) && (electron.get_pos(0) <= (magnet.get_length()+(magnet.get_pos(0)))) && (electron.get_pos(1) >= ((magnet.get_pos(1))-((magnet.get_width())/2.0))) && (electron.get_pos(1) <= ((magnet.get_pos(1))+(magnet.get_width())/2.0)));
            
            if(check)
                { outfile_part_writeAndComma(electron); }
            else if(!(check))
                { outfile_part_write(electron); }            
            
        } while(check);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void stepThroughMagnet_Leap(Particle *electron, Magnet &magnet, double& time, const double &del_time)
{
    bool check;
    ThreeVec vn_minus = (electron->get_vel() ) - (( (electron->get_vel())^(magnet.get_B0()) )*del_time);
    ThreeVec vn_plus;

    ThreeVec rn_minus = (electron->get_pos() ) - ((electron->get_vel())*del_time);
    ThreeVec rn_plus;
    do
        {
            vn_plus = vn_minus + ((electron->get_vel())^(magnet.get_B0()))*2*del_time;    //leapfrog algarithm vn for uniform B0
            rn_plus = rn_minus + (electron->get_vel())*2*del_time;         //leapfrog algarithm rn for uniform B0

            rn_minus = electron->get_pos();                          //iterate xn_minus up a timestep
            electron->set_pos(rn_plus);                              //do the same to xn

            vn_minus = (electron->get_vel());                        //iterate vn_minus up a timestep
            electron->set_vel(vn_plus);                              //do the same to vn

            time += del_time;

            check = ((electron->get_pos(0) >= (magnet.get_pos(0))) && (electron->get_pos(0) <= (magnet.get_length()+(magnet.get_pos(0)))) && (electron->get_pos(1) >= ((magnet.get_pos(1))-((magnet.get_width())/2.0))) && (electron->get_pos(1) <= ((magnet.get_pos(1))+(magnet.get_width())/2.0)));
            
            if(check)
                { outfile_part_writeAndComma(*electron); }
            else if(!(check))
                { outfile_part_write(*electron); }            
            
        } while(check);
        //std::cerr << electron->get_pos() << '\n' << '~' << '\n';
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void boris(Particle &electron_t, Magnet &magnet_t, double del_t, int counter)
{
    //Rotation of momentum vector in uniform field
    // p+ = p- + (p- + (p- x t)) x s
    // t = B q del_t / (2 m gamma)
    // s = 2 t / (1 + |t|^2)
  unsigned long parnum,x1;
  double igamma, psquared, Bsquared;

  // t vector in Boris method
  double tt[3],ttsquared;
  // s vector in Boris method
  double ss[3];
  // vstar vector in Boris method
  double vstar[3];
  // perpendicular component of v in Boris method
  double vperp[3];

  //initial half-position update
    psquared  = ((electron_t.get_p(0) * electron_t.get_p(0)) + (electron_t.get_p(1) * electron_t.get_p(1)) 
                + (electron_t.get_p(2) * electron_t.get_p(2)));
    igamma    = 1.0/(sqrt(1.0+psquared)); 

        electron_t.set_pos( 0, electron_t.get_pos(0) + (electron_t.get_p(0) * del_t * igamma * 0.5) );
        electron_t.set_pos( 1, electron_t.get_pos(1) + (electron_t.get_p(1) * del_t * igamma * 0.5) );
        electron_t.set_pos( 2, electron_t.get_pos(2) + (electron_t.get_p(2) * del_t * igamma * 0.5) ); 

    Bsquared  = ( ( (magnet_t.get_B0(0))*(magnet_t.get_B0(0)) ) + ( (magnet_t.get_B0(1))*(magnet_t.get_B0(1)) )
                + ( (magnet_t.get_B0(2))*(magnet_t.get_B0(2)) ) );
    
    double old_px, old_py, old_pz;

        ttsquared = 0.0;
    //if(counter!=0)
    //{
        for (x1=0; x1<3; ++x1)
            {
            tt[x1]     = (magnet_t.get_B0(x1) )*del_t*0.5*igamma;
            ttsquared += (tt[x1]*tt[x1]);
            }
        
        for (x1=0; x1<3; ++x1)
            {
            ss[x1]     = 2.0*tt[x1]/(1.0+ttsquared);
            }
        double u_cross_tt_plus_u[3];
        u_cross_tt_plus_u[0] = electron_t.get_p(0) + ((electron_t.get_p(1) * tt[2]) - (electron_t.get_p(2) * tt[1]));
        u_cross_tt_plus_u[1] = electron_t.get_p(1) + ((electron_t.get_p(2) * tt[0]) - (electron_t.get_p(0) * tt[2]));
        u_cross_tt_plus_u[2] = electron_t.get_p(2) + ((electron_t.get_p(0) * tt[1]) - (electron_t.get_p(1) * tt[0]));

        double all_cross_ss[3];
        all_cross_ss[0] = (u_cross_tt_plus_u[1]*ss[2]) - (u_cross_tt_plus_u[2]*ss[1]);
        all_cross_ss[1] = (u_cross_tt_plus_u[2]*ss[0]) - (u_cross_tt_plus_u[0]*ss[2]);
        all_cross_ss[2] = (u_cross_tt_plus_u[0]*ss[1]) - (u_cross_tt_plus_u[1]*ss[0]);

        old_px = electron_t.get_p(0);
        old_py = electron_t.get_p(1);
        old_pz = electron_t.get_p(2);


        electron_t.set_p( 0, ( (electron_t.get_p(0)) + (all_cross_ss[0])));
        electron_t.set_p( 1, ( (electron_t.get_p(1)) + (all_cross_ss[1])));
        electron_t.set_p( 2, ( (electron_t.get_p(2)) + (all_cross_ss[2])));

    //Do average velocity to update position

    ThreeVec average_vel((old_px+electron_t.get_p(0))*igamma*0.5, (old_py+electron_t.get_p(1))*igamma*0.5, (old_pz+electron_t.get_p(2))*igamma*0.5);
    
    //ThreeVec half_pos((electron_t.get_pos(0) + (old_px*igamma*del_t*0.5)),(electron_t.get_pos(1) + (old_py*igamma*del_t*0.5)),(electron_t.get_pos(2) + (old_pz*igamma*del_t*0.5)));

    electron_t.set_pos( 0, electron_t.get_pos(0) + (electron_t.get_p(0) * del_t * igamma * 0.5) );
    electron_t.set_pos( 1, electron_t.get_pos(1) + (electron_t.get_p(1) * del_t * igamma * 0.5) );
    electron_t.set_pos( 2, electron_t.get_pos(2) + (electron_t.get_p(2) * del_t * igamma * 0.5) );
    
    //electron_t.set_pos( 0, half_pos.getX() + (average_vel.getX() * del_t) );
    //electron_t.set_pos( 1, half_pos.getY() + (average_vel.getY() * del_t) );
    //electron_t.set_pos( 2, half_pos.getZ() + (average_vel.getZ() * del_t) );
    //electron_t.set_pos( 0, electron_t.get_pos(0) + (average_vel.getX() * del_t) );
    //electron_t.set_pos( 1, electron_t.get_pos(1) + (average_vel.getY() * del_t) );
    //electron_t.set_pos( 2, electron_t.get_pos(2) + (average_vel.getZ() * del_t) );
}




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void step_through_magnet_mag_boris(Particle &electron, Magnet &magnet, double& time, const double &del_time, double time_out)
{
    //TESTING CHANING TIME STEP:
    //const double del_time = del_time_1 * 0.5;

    bool check_x, check_y, check_z;
    double psquared;
    int counter = 0;
    do
        {
            boris(electron, magnet, del_time, counter); //updates particle velocity & position in magnetic field 
            counter++;

            time += del_time;
//            std::cerr << del_time << std::endl;

            electron.set_time(time);

            
            check_x = (electron.get_pos(0) >= (magnet.get_pos(0))) && (electron.get_pos(0) <= (magnet.get_length()+(magnet.get_pos(0))));
            check_y = (electron.get_pos(1) >= ((magnet.get_pos(1))-((magnet.get_width())/2.0)))  && (electron.get_pos(1) <= ((magnet.get_pos(1))+(magnet.get_width())/2.0));
            check_z = (electron.get_pos(2) >= ((magnet.get_pos(2))-((magnet.get_height())/2.0))) && (electron.get_pos(2) <= ((magnet.get_pos(2))+(magnet.get_height())/2.0));

            if((check_x && check_y && check_z) && !(time>=time_out))
                { outfile_part_commaAndWrite(electron); }
            else if( (!( check_x && check_y && check_z )) || (time >= time_out) )
                {
                    outfile_part_commaAndWrite(electron); 
                }
            if(time >= time_out)
                { std::cout << "Particle Timed-Out.\n"; }
            if( (check_x && check_y && check_z) == false )
                {
                    std::cout << "Particle Out of Magnet.\n";
                    if(check_x == false)
                        { std::cout << "Out of bounds in x.\n"; } 
                    if(check_y == false)
                        { std::cout << "Out of bounds in y.\n"; } 
                    if(check_z == false)
                        { std::cout << "Out of bounds in z.\n"; } 
                }
            
        } while((check_x && check_y && check_z) && (time < time_out));
}
 




void first_half_position_step(Particle &electron_t, const double del_t)
{
    double psquared  = ((electron_t.get_p(0) * electron_t.get_p(0)) + (electron_t.get_p(1) * electron_t.get_p(1)) 
                + (electron_t.get_p(2) * electron_t.get_p(2)));
    double igamma = 1.0/(sqrt(1.0+psquared));

    electron_t.set_pos(0, electron_t.get_pos(0) + (electron_t.get_p(0)*igamma*del_t*0.5));
    electron_t.set_pos(1, electron_t.get_pos(1) + (electron_t.get_p(1)*igamma*del_t*0.5));
    electron_t.set_pos(2, electron_t.get_pos(2) + (electron_t.get_p(2)*igamma*del_t*0.5));
}






bool inside_of_mag(Magnet magnet_t, Particle particle_t)
{
    bool inside_of_mag;

    bool inside_x_limits = (particle_t.get_pos(0) >= magnet_t.get_pos(0)) && (particle_t.get_pos(0) <= (magnet_t.get_pos(0)+magnet_t.get_length()));
    bool inside_y_limits = (particle_t.get_pos(1) >= (magnet_t.get_pos(1) - (magnet_t.get_width()/2.0)))  && (particle_t.get_pos(1) <= (magnet_t.get_pos(1) + (magnet_t.get_width()/2.0)));
    bool inside_z_limits = (particle_t.get_pos(2) >= (magnet_t.get_pos(2) - (magnet_t.get_height()/2.0))) && (particle_t.get_pos(2) <= (magnet_t.get_pos(2) + (magnet_t.get_height()/2.0)));
    if(inside_x_limits && inside_y_limits && inside_z_limits)
    {
        inside_of_mag = true;
    }
    else
    {
        inside_of_mag = false;
    }
    return inside_of_mag;
}



double time_to_magnet_boundary(Magnet magnet_t, Particle particle_t)
{
    //Assumes that particle has already been checked if it is at-or-within magnet's boundary

    double time_btwn_mags_x_front;
    double time_btwn_mags_x_back;
    if(particle_t.get_vel(0) !=0)
        { 
            time_btwn_mags_x_front = (magnet_t.get_pos(0) - particle_t.get_pos(0))/particle_t.get_vel(0); 
            time_btwn_mags_x_back  = (magnet_t.get_pos(0)+magnet_t.get_length() - particle_t.get_pos(0))/particle_t.get_vel(0);
        }
    else
        {
            time_btwn_mags_x_back  = DBL_MAX;
            time_btwn_mags_x_front = DBL_MAX;
        }

    double time_btwn_mags_y_top;
    double time_btwn_mags_y_bottom;
    if(particle_t.get_vel(1) !=0)
        { 
            time_btwn_mags_y_top     = (magnet_t.get_pos(1)+(magnet_t.get_width()/2.0) - particle_t.get_pos(1))/particle_t.get_vel(1); 
            time_btwn_mags_y_bottom  = (magnet_t.get_pos(1)-(magnet_t.get_width()/2.0) - particle_t.get_pos(1))/particle_t.get_vel(1);
        }
    else
        {
            time_btwn_mags_y_top    = DBL_MAX;
            time_btwn_mags_y_bottom = DBL_MAX;
        }

    double time_btwn_mags_z_top;
    double time_btwn_mags_z_bottom;
    if(particle_t.get_vel(2) !=0)
        { 
            time_btwn_mags_z_top     = (magnet_t.get_pos(2)+(magnet_t.get_height()/2.0) - particle_t.get_pos(2))/particle_t.get_vel(2); 
            time_btwn_mags_z_bottom  = (magnet_t.get_pos(2)-(magnet_t.get_height()/2.0) - particle_t.get_pos(2))/particle_t.get_vel(2);
        }
    else
        {
            time_btwn_mags_z_top    = DBL_MAX;
            time_btwn_mags_z_bottom = DBL_MAX;
        }

    double shortest_time = DBL_MAX;
    if((time_btwn_mags_x_front < shortest_time) && (time_btwn_mags_x_front>0.0))
        { shortest_time = time_btwn_mags_x_front;  }
    if((time_btwn_mags_x_back < shortest_time) && (time_btwn_mags_x_back>0.0))
        { shortest_time = time_btwn_mags_x_back;   }
    if((time_btwn_mags_y_top < shortest_time) && (time_btwn_mags_y_top>0.0))
        { shortest_time = time_btwn_mags_y_top;    }
    if((time_btwn_mags_y_bottom < shortest_time) && (time_btwn_mags_y_bottom>0.0))
        { shortest_time = time_btwn_mags_y_bottom; }
    if((time_btwn_mags_z_top < shortest_time) && (time_btwn_mags_z_top>0.0))
        { shortest_time = time_btwn_mags_z_top;    }
    if((time_btwn_mags_z_bottom < shortest_time) && (time_btwn_mags_z_bottom>0.0))
        { shortest_time = time_btwn_mags_z_bottom; }

    if(shortest_time == DBL_MAX)
    {
        return -1.0;
    }
    return shortest_time;
}



bool intersect_mag(Magnet magnet_t, Particle particle_t)
{
    /* Will check if particle intersects with a magnet
     * Will NOT update particle's positions! We'll want to check if there's another
     * magnet that is closer than the one currently being checked before updating position.
     */
    bool intersect;

    double time_to_magnet = time_to_magnet_boundary(magnet_t, particle_t);
    if(time_to_magnet < 0.0)
    {
        return false;
    }

    ThreeVec pos_at_shortest_time;
    pos_at_shortest_time.setX(particle_t.get_pos(0) + (time_to_magnet * particle_t.get_vel(0)));
    pos_at_shortest_time.setY(particle_t.get_pos(1) + (time_to_magnet * particle_t.get_vel(1)));
    pos_at_shortest_time.setZ(particle_t.get_pos(2) + (time_to_magnet * particle_t.get_vel(2)));

    bool within_x_bounds = (pos_at_shortest_time.getX() >= magnet_t.get_pos(0)) && (pos_at_shortest_time.getX() <= (magnet_t.get_pos(0)+magnet_t.get_length()));
    bool within_y_bounds = (pos_at_shortest_time.getY() >= (magnet_t.get_pos(1)-(magnet_t.get_width()/2.0)))  && (pos_at_shortest_time.getY() <= (magnet_t.get_pos(1)+(magnet_t.get_width()/2.0)));
    bool within_z_bounds = (pos_at_shortest_time.getZ() >= (magnet_t.get_pos(2)-(magnet_t.get_height()/2.0))) && (pos_at_shortest_time.getZ() <= (magnet_t.get_pos(2)+(magnet_t.get_height()/2.0)));

    if(within_x_bounds && within_y_bounds && within_z_bounds)
    { intersect = true; }
    else
    { intersect = false; }
 
    return intersect;
}




double dist_to_mag(Magnet magnet_t, Particle particle_t)
{
    double dist_to_mag;
    double time_to_mag = time_to_magnet_boundary(magnet_t, particle_t);

    if(time_to_mag <= 0.0)
    {
        return 0.0;
    }

    ThreeVec pos_at_shortest_time;
    pos_at_shortest_time.setX(particle_t.get_pos(0) + (time_to_mag * particle_t.get_vel(0)));
    pos_at_shortest_time.setY(particle_t.get_pos(1) + (time_to_mag * particle_t.get_vel(1)));
    pos_at_shortest_time.setZ(particle_t.get_pos(2) + (time_to_mag * particle_t.get_vel(2)));

    dist_to_mag = sqrt((pos_at_shortest_time.getX() - particle_t.get_pos(0))*(pos_at_shortest_time.getX() - particle_t.get_pos(0)) + (pos_at_shortest_time.getY() - particle_t.get_pos(1))*(pos_at_shortest_time.getY() - particle_t.get_pos(1)) + (pos_at_shortest_time.getZ() - particle_t.get_pos(2))*(pos_at_shortest_time.getZ() - particle_t.get_pos(2)));
    return dist_to_mag;
}





void move_particle_to_magnet(Magnet magnet_t, Particle &particle_t)
{
    double time_to_mag = time_to_magnet_boundary(magnet_t, particle_t);

    ThreeVec pos_at_shortest_time;
    pos_at_shortest_time.setX(particle_t.get_pos(0) + (time_to_mag * particle_t.get_vel(0)));
    pos_at_shortest_time.setY(particle_t.get_pos(1) + (time_to_mag * particle_t.get_vel(1)));
    pos_at_shortest_time.setZ(particle_t.get_pos(2) + (time_to_mag * particle_t.get_vel(2)));

    particle_t.set_pos(0, pos_at_shortest_time.getX());
    particle_t.set_pos(1, pos_at_shortest_time.getY());
    particle_t.set_pos(2, pos_at_shortest_time.getZ());

    double time_at_mag = particle_t.get_time() + time_to_mag;
    
    particle_t.set_time(time_at_mag);
}





void move_through_magnets(Magnet magnet_t[], int num_mags, Particle &particle_t, double &time, double del_time, double time_limit)
{
    bool check_inside_magnet = false;
    bool check_intersect_magnet = false;
    double distance_to_mag_ii[num_mags];
    int boris_counter = 0;
    for(int ii=0; ii<num_mags; ii++)
    {
        check_inside_magnet = inside_of_mag(magnet_t[ii], particle_t);
        if(check_inside_magnet==true)
        {
           //outfile_part_comma(particle_t);
            step_through_magnet_mag_boris(particle_t, magnet_t[ii], time, del_time, time_limit);
            boris_counter++;
            ii = -1; //restart loop (which will iterate ii by 1, to zero)
            continue;
        }

        check_intersect_magnet = intersect_mag(magnet_t[ii], particle_t);
        if(check_intersect_magnet==true)
        {
            distance_to_mag_ii[ii] = dist_to_mag(magnet_t[ii], particle_t);
        }
        else
        {
            distance_to_mag_ii[ii] = 0.0;
        }

        if(ii == (num_mags - 1))
        {
            double shortest_distance = DBL_MAX;
            int index_of_shortest = -1;
            for(int jj = 0; jj<num_mags; jj++)
            {
                if(distance_to_mag_ii[jj] != 0.0 && distance_to_mag_ii[jj] < shortest_distance)
                {
                    ii = -1;
                    shortest_distance = distance_to_mag_ii[jj];
                    index_of_shortest = jj;
                }
            }

            if(index_of_shortest != -1)
            {                
                move_particle_to_magnet(magnet_t[index_of_shortest], particle_t);
                //std::cerr << "Sending to next magnet" << std::endl;
                outfile_part_commaAndWrite(particle_t);

                double particle_time = particle_t.get_time();
                time = time + particle_time;
                continue;
            }
        }
    }
}

void half_time_step(double &time_step)
{
    time_step = time_step * 0.5;
    //time_step = time_step - 0.1;
    //return time_step;
}




double screen_plane_half_equation(Screen screen_t, double x_pos, double y_pos, double z_pos)
{
    double alpha = screen_t.get_angle_about_z();
    double beta  = screen_t.get_angle_about_y();
    double gamma = screen_t.get_angle_about_x();
    double i_hat = x_pos*(sin(alpha)*cos(beta)*cos(beta)*cos(gamma) + (sin(beta)*(sin(alpha)*sin(beta)*cos(gamma) - cos(alpha)*sin(gamma))));
    double j_hat = y_pos*(cos(alpha)*cos(beta)*cos(beta)*cos(gamma) + (sin(beta)*(cos(alpha)*sin(beta)*cos(gamma) + sin(alpha)*sin(gamma))));
    double k_hat = z_pos*(cos(alpha)*cos(beta)*(sin(alpha)*sin(beta)*cos(gamma) - cos(alpha)*sin(gamma)) - sin(alpha)*cos(beta)*(cos(alpha)*sin(beta)*cos(gamma) + sin(alpha)*sin(gamma)));
    return (i_hat - j_hat + k_hat);
}




double t_line_particle(Screen screen_t, Particle particle_t)
{
    double t_numerator = screen_plane_half_equation(screen_t, (screen_t.get_pos(0) - particle_t.get_pos(0)), (screen_t.get_pos(1) - particle_t.get_pos(1)), (screen_t.get_pos(2) - particle_t.get_pos(2)));
    double t_denom     = screen_plane_half_equation(screen_t, particle_t.get_vel(0), particle_t.get_vel(1), particle_t.get_vel(2));
    // double t_numerator = (sin(screen_t.get_angle_x())*cos(screen_t.get_angle_z())*(screen_t.get_pos(0) - particle_t.get_pos(0)) - (((cos(screen_t.get_angle_x())*cos(screen_t.get_angle_z())*(1 - sin(screen_t.get_angle_z())) - cos(screen_t.get_angle_x())*sin(screen_t.get_angle_z())*sin(screen_t.get_angle_z()))*(screen_t.get_pos(1) - particle_t.get_pos(1)))) - (sin(screen_t.get_angle_x())*sin(screen_t.get_angle_z())*(screen_t.get_pos(2) - particle_t.get_pos(2))) );
    // double t_denom     = particle_t.get_vel(0)*sin(screen_t.get_angle_x())*cos(screen_t.get_angle_z()) - particle_t.get_vel(1)*(cos(screen_t.get_angle_x())*cos(screen_t.get_angle_z())*(1 - sin(screen_t.get_angle_z())) - cos(screen_t.get_angle_x())*sin(screen_t.get_angle_z())*sin(screen_t.get_angle_z())) - particle_t.get_vel(2)*sin(screen_t.get_angle_x())*sin(screen_t.get_angle_z());
    return (t_numerator / t_denom);
}



void min_max_x(double corner1_x, double corner2_x, double corner3_x, double corner4_x, double min_max_array[])
{
    double max = corner1_x; /* assume x is the largest */
	if (corner2_x > max) { /* if y is larger than max, assign y to max */
		max = corner2_x;
	} /* end if */
	if (corner3_x > max) { /* if z is larger than max, assign z to max */
		max = corner3_x;
	} /* end if */
    if (corner4_x > max) { /* if z is larger than max, assign z to max */
		max = corner4_x;
	} /* end if */

    double min = corner1_x; /* assume x is the smallest */
	if (corner2_x < min) { /* if y is larger than max, assign y to max */
		min = corner2_x;
	} /* end if */
	if (corner3_x < min) { /* if z is larger than max, assign z to max */
		min = corner3_x;
	} /* end if */
    if (corner4_x < min) { /* if z is larger than max, assign z to max */
		min = corner4_x;
	} /* end if */

    min_max_array[0] = min;
    min_max_array[1] = max;
    //return min_max;
}



bool check_if_intersect_screen(Screen screen_t, Particle particle_t)
{
    bool intersect = false;
    double t_line  = t_line_particle(screen_t, particle_t);

    if(t_line < 0.0)
    {
        return false;
    }
    double alpha = screen_t.get_angle_about_z();
    double beta  = screen_t.get_angle_about_y();
    double gamma = screen_t.get_angle_about_x();
    ThreeVec corner1(screen_t.get_pos(0) + (screen_t.get_height()*0.5)*(cos(alpha)*sin(beta)*cos(gamma) + sin(alpha)*sin(gamma)), screen_t.get_pos(1) + (screen_t.get_height()*0.5)*(sin(alpha)*sin(beta)*cos(gamma) - cos(alpha)*sin(gamma)), screen_t.get_pos(2) + (screen_t.get_height()*0.5)*(cos(beta)*cos(gamma)));
    ThreeVec corner2(screen_t.get_pos(0) - (screen_t.get_height()*0.5)*(cos(alpha)*sin(beta)*cos(gamma) + sin(alpha)*sin(gamma)), screen_t.get_pos(1) - (screen_t.get_height()*0.5)*(sin(alpha)*sin(beta)*cos(gamma) - cos(alpha)*sin(gamma)), screen_t.get_pos(2) - (screen_t.get_height()*0.5)*(cos(beta)*cos(gamma)));
    ThreeVec corner3(screen_t.get_pos(0) + screen_t.get_length()*cos(alpha)*cos(beta) + (screen_t.get_height()*0.5)*(cos(alpha)*sin(beta)*cos(gamma) + sin(alpha)*sin(gamma)), screen_t.get_pos(1) + screen_t.get_length()*sin(alpha)*cos(beta) + (screen_t.get_height()*0.5)*(sin(alpha)*sin(beta)*cos(gamma) - cos(alpha)*sin(gamma)), screen_t.get_pos(2) - screen_t.get_length()*sin(beta) + (screen_t.get_height()*0.5)*(cos(beta)*cos(gamma)));
    ThreeVec corner4(screen_t.get_pos(0) + screen_t.get_length()*cos(alpha)*cos(beta) - (screen_t.get_height()*0.5)*(cos(alpha)*sin(beta)*cos(gamma) + sin(alpha)*sin(gamma)), screen_t.get_pos(1) + screen_t.get_length()*sin(alpha)*cos(beta) - (screen_t.get_height()*0.5)*(sin(alpha)*sin(beta)*cos(gamma) - cos(alpha)*sin(gamma)), screen_t.get_pos(2) - screen_t.get_length()*sin(beta) - (screen_t.get_height()*0.5)*(cos(beta)*cos(gamma)));
    ThreeVec intersection(particle_t.get_pos(0) + particle_t.get_vel(0)*t_line, particle_t.get_pos(1) + particle_t.get_vel(1)*t_line, particle_t.get_pos(2) + particle_t.get_vel(2)*t_line);
    
    //check areas of 4 triangles made by conecting corners to intersect point; if area equals area of rectangle, then its bounded and intersects
    //method from https://math.stackexchange.com/questions/190111/how-to-check-if-a-point-is-inside-a-rectangle and https://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
    double area0 = screen_t.get_length()*screen_t.get_height();
    double tmin, dist_to_intersect;
    ThreeVec closest_point(0.0,0.0,0.0);
    
    tmin = -(((corner1 - intersection) * (corner2 - corner1))/abs((corner2 - corner1)*(corner2 - corner1)));
    closest_point.setX(corner1.getX() + (corner2.getX() - corner1.getX())*tmin);
    closest_point.setY(corner1.getY() + (corner2.getY() - corner1.getY())*tmin);
    closest_point.setZ(corner1.getZ() + (corner2.getZ() - corner1.getZ())*tmin);
    dist_to_intersect = sqrt(((closest_point.getX()-intersection.getX())*(closest_point.getX()-intersection.getX())) + ((closest_point.getY()-intersection.getY())*(closest_point.getY()-intersection.getY())) + ((closest_point.getZ()-intersection.getZ())*(closest_point.getZ()-intersection.getZ())));
    double area12 = 0.5*screen_t.get_height()*dist_to_intersect;

    tmin = -(((corner3 - intersection) * (corner1 - corner3))/abs((corner1 - corner3)*(corner1 - corner3)));
    closest_point.setX(corner3.getX() + (corner1.getX() - corner3.getX())*tmin);
    closest_point.setY(corner3.getY() + (corner1.getY() - corner3.getY())*tmin);
    closest_point.setZ(corner3.getZ() + (corner1.getZ() - corner3.getZ())*tmin);
    dist_to_intersect = sqrt(((closest_point.getX()-intersection.getX())*(closest_point.getX()-intersection.getX())) + ((closest_point.getY()-intersection.getY())*(closest_point.getY()-intersection.getY())) + ((closest_point.getZ()-intersection.getZ())*(closest_point.getZ()-intersection.getZ())));
    double area13 = 0.5*screen_t.get_length()*dist_to_intersect;

    tmin = -(((corner4 - intersection) * (corner2 - corner4))/abs((corner2 - corner4)*(corner2 - corner4)));
    closest_point.setX(corner4.getX() + (corner2.getX() - corner4.getX())*tmin);
    closest_point.setY(corner4.getY() + (corner2.getY() - corner4.getY())*tmin);
    closest_point.setZ(corner4.getZ() + (corner2.getZ() - corner4.getZ())*tmin);
    dist_to_intersect = sqrt(((closest_point.getX()-intersection.getX())*(closest_point.getX()-intersection.getX())) + ((closest_point.getY()-intersection.getY())*(closest_point.getY()-intersection.getY())) + ((closest_point.getZ()-intersection.getZ())*(closest_point.getZ()-intersection.getZ())));
    double area24 = 0.5*screen_t.get_length()*dist_to_intersect;

    tmin = -(((corner3 - intersection) * (corner4 - corner3))/abs((corner4 - corner3)*(corner4 - corner3)));
    closest_point.setX(corner3.getX() + (corner4.getX() - corner3.getX())*tmin);
    closest_point.setY(corner3.getY() + (corner4.getY() - corner3.getY())*tmin);
    closest_point.setZ(corner3.getZ() + (corner4.getZ() - corner3.getZ())*tmin);
    dist_to_intersect = sqrt(((closest_point.getX()-intersection.getX())*(closest_point.getX()-intersection.getX())) + ((closest_point.getY()-intersection.getY())*(closest_point.getY()-intersection.getY())) + ((closest_point.getZ()-intersection.getZ())*(closest_point.getZ()-intersection.getZ())));
    double area34 = 0.5*screen_t.get_height()*dist_to_intersect;

    //std::cout << "area 0 = " << area0 << " and areas combined are " << (area12 + area13 + area24 + area34) << std::endl;
    if((area0 + area0*0.01) >= (area12 + area13 + area24 + area34))
    {
        intersect = true;
    }

    return intersect;
}


void move_to_screens(Screen screen_t[], int num_screen, Particle particle_t)
{
    double dist_to_screen_ii[num_screen];
    int jj = -1;
    bool intersected = false;
    double t_numerator, t_denom, t_line;

    for(int ii = 0; ii < num_screen; ii++)
    {
        bool intersect_check = check_if_intersect_screen(screen_t[ii], particle_t);
        if(intersect_check && (jj != ii))
        {
            intersected = true;
            t_line      = t_line_particle(screen_t[ii], particle_t);

            double x_intersect = particle_t.get_pos(0) + particle_t.get_vel(0)*t_line;
            double y_intersect = particle_t.get_pos(1) + particle_t.get_vel(1)*t_line;
            double z_intersect = particle_t.get_pos(2) + particle_t.get_vel(2)*t_line;

            dist_to_screen_ii[ii] = sqrt((x_intersect - particle_t.get_pos(0))*(x_intersect - particle_t.get_pos(0)) + (y_intersect - particle_t.get_pos(1))*(y_intersect - particle_t.get_pos(1)) + (z_intersect - particle_t.get_pos(2))*(z_intersect - particle_t.get_pos(2)));
        }
        else
        {
            dist_to_screen_ii[ii] = -1;
        }


        if( ii == (num_screen - 1) && intersected)
        {
            double shortest_distance = DBL_MAX;
            int index_of_shortest = -1;
            for(int kk = 0; kk < num_screen; kk++)
            {
                if(dist_to_screen_ii[kk] != -1 && dist_to_screen_ii[kk] < shortest_distance)
                {
                    ii = -1;
                    intersected = false;
                    jj = kk;
                    
                    shortest_distance = dist_to_screen_ii[kk];
                    index_of_shortest = kk;
                }
            }

            if(index_of_shortest != -1)
            {                
                t_line             = t_line_particle(screen_t[index_of_shortest], particle_t);
                double x_intersect = particle_t.get_pos(0) + particle_t.get_vel(0)*t_line;
                double y_intersect = particle_t.get_pos(1) + particle_t.get_vel(1)*t_line;
                double z_intersect = particle_t.get_pos(2) + particle_t.get_vel(2)*t_line;
                particle_t.set_pos(x_intersect, y_intersect, z_intersect);

                double particle_time = particle_t.get_time();
                particle_time += t_line;
                particle_t.set_time(particle_time);
                outfile_part_commaAndWrite(particle_t);
                continue;
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void readUnits(std::ifstream &input_stream, std::vector<std::string> &desired_units) {
    std::string length_unit;
    std::string energy_unit;
    std::string angle_unit;
    std::string magnetic_field_unit;

    input_stream >> length_unit >> energy_unit >> angle_unit >> magnetic_field_unit;

    desired_units = {length_unit, energy_unit, angle_unit, magnetic_field_unit};
}

int readNumOf(std::ifstream &input_stream) {
    std::string numStr;
    input_stream >> numStr;
    int numInt = std::stoi(numStr);
    
    return numInt;
}

void readMagnet(std::ifstream &input_stream, int &magNum, std::vector<std::vector<std::vector<double>>> &magInfo) {
    magNum = readNumOf(input_stream);
    std::string tempStr;

    //outside index goes dimensions, position, magnetic field value
    //dimensions go width, length, height
    for(int i=0; i<3; ++i) {
        std::vector<std::vector<double>> tempMagSet;
        
        if(i!=2) {
            for(int j=0; j<magNum; ++j) {
                std::vector<double> tempInfoBits;

                for(int k=0; k<3; ++k) {
                    input_stream >> tempStr;
                    double tempDbl = std::stod(tempStr);
                    tempInfoBits.push_back(tempDbl);
                }
                tempMagSet.push_back(tempInfoBits);
            }
            magInfo.push_back(tempMagSet);
        }
        else {
            for(int j=0; j<magNum; ++j) {
                std::vector<double> tempInfoBits;
                
                input_stream >> tempStr;
                double tempDbl = std::stod(tempStr);
                tempInfoBits.push_back(tempDbl);
                tempMagSet.push_back(tempInfoBits);
            }
            magInfo.push_back(tempMagSet);
        }
    }
}

void readPermanentMagDim(std::ifstream &input_stream, int magNum, std::vector<double> &PmagDim) {
    std::string tempStr;
    
    for(int i=0; i<magNum; ++i) {
        input_stream >> tempStr;
        double tempDbl = std::stod(tempStr);
        PmagDim.push_back(tempDbl);
    }
}

void readMagAxis(std::ifstream &input_stream, int magNum, std::vector<char> &axisInfo) {
    char tempChar;
    
    for(int i=0; i<magNum; ++i) {
        input_stream >> tempChar;
        axisInfo.push_back(tempChar);
    }
}

void readBeam(std::ifstream &input_stream, std::vector<std::vector<double>> &beamInfo) {
    std::string tempStr;

    //outside index goes position, energy, direction
    for(int i=0; i<3; ++i) {
        std::vector<double> tempInfoBits;

        if(i==1) {
            input_stream >> tempStr;
            double tempDbl = std::stod(tempStr);
            tempInfoBits.push_back(tempDbl);
        }
        else {
            for(int j=0; j<3; ++j) {
                input_stream >> tempStr;
                double tempDbl = std::stod(tempStr);
                tempInfoBits.push_back(tempDbl);
            }
        }
        beamInfo.push_back(tempInfoBits);
    }
}

void readSpread(std::ifstream &input_stream, std::vector<std::vector<double>> &spreadInfo) {
    std::string tempStr;

    //outside index goes position, energy, divergence
    for(int i=0; i<3; ++i) {
        std::vector<double> tempInfoBits;

        if(i==1) {
            input_stream >> tempStr;
            double tempDbl = std::stod(tempStr);
            tempInfoBits.push_back(tempDbl);
        }
        else {
            for(int k=0; k<3; ++k) {
                input_stream >> tempStr;
                double tempDbl = std::stod(tempStr);
                tempInfoBits.push_back(tempDbl);
            }
        }
        spreadInfo.push_back(tempInfoBits);
    }
}

void readScreen(std::ifstream &input_stream, int &screenNum, std::vector<std::vector<std::vector<double>>> &screenInfo) {
    screenNum = readNumOf(input_stream);
    std::string tempStr;

    //outside index goes dimensions, position, angles
    //dimensions go length, height
    //angles go z-axis, y-axis, x-axis
    for(int i=0; i<3; ++i) {
        std::vector<std::vector<double>> tempScreenSet;

        for(int j=0; j<screenNum; ++j) {
            std::vector<double> tempInfoBits;

            if(i==0) {
                for(int k=0; k<2; ++k) {
                    input_stream >> tempStr;
                    double tempDbl = std::stod(tempStr);
                    tempInfoBits.push_back(tempDbl);
                }
            }
            else {
                for(int k=0; k<3; ++k) {
                    input_stream >> tempStr;
                    double tempDbl = std::stod(tempStr);
                    tempInfoBits.push_back(tempDbl);
                }
            }
            tempScreenSet.push_back(tempInfoBits);
        }
        screenInfo.push_back(tempScreenSet);
    }
}

void ReadInitTypes(std::ifstream &input_stream, std::vector<int> &init_types) {
    std::string pos_type_str;
    std::string energy_type_str;
    std::string div_type_str;

    input_stream >> pos_type_str >> energy_type_str >> div_type_str;
    int pos_type, energy_type, div_type;
    pos_type = std::stoi(pos_type_str);
    energy_type = std::stoi(energy_type_str);
    div_type = std::stoi(div_type_str);

    init_types = {pos_type, energy_type, div_type};
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double find_magnetization(Magnet &magnet, double mag_height) {
    // dimensions of permanent magnet
    double b = 0.5*magnet.get_width();   // half width
    double a = 0.5*magnet.get_length();  // half length
    double c = 0.5*mag_height;           // half height
    
    double z = 0.5*magnet.get_height() + c;  // distance from magnet center to center between magnets
    double Bz = magnet.get_B0(2);
    
    double magnetization = (Bz*4*pow(10.,7))/(2*(atan((a*b)/((z-c)*sqrt(pow(b,2)+pow(a,2)+pow(z-c,2))))-
                                                            atan((a*b)/((z+c)*sqrt(pow(a,2)+pow(b,2)+pow(z+c,2))))));
    // divided by 2 because of the contributions of both magnets
    
    return magnetization;
}

void calc_grid_B_comps(double factor, double a, double b, double c, double x, double y, double z) {
    
   double temp_B1 = factor*log( (( (sqrt(pow(-x+a,2)+pow(y-b,2)+pow(-z+c,2))+b-y)/(sqrt(pow(-x+a,2)+pow(y+b,2)+pow(-z+c,2))-b-y) )*( (sqrt(pow(x+a,2)+pow(y-b,2)+pow(z+c,2))+b-y)/(sqrt(pow(x+a,2)+pow(y+b,2)+pow(z+c,2))-b-y) ))/(( (sqrt(pow(x+a,2)+pow(y-b,2)+pow(-z+c,2))+b-y)/(sqrt(pow(x+a,2)+pow(y+b,2)+pow(-z+c,2))-b-y) )*( (sqrt(pow(-x+a,2)+pow(y-b,2)+pow(z+c,2))+b-y)/(sqrt(pow(-x+a,2)+pow(y+b,2)+pow(z+c,2))-b-y) )) );
    
   double temp_B2 = factor*log( (( (sqrt(pow(x-a,2)+pow(-y+b,2)+pow(-z+c,2))+a-x)/(sqrt(pow(x+a,2)+pow(-y+b,2)+pow(-z+c,2))-a-x) )*( (sqrt(pow(x-a,2)+pow(y+b,2)+pow(z+c,2))+a-x)/(sqrt(pow(x+a,2)+pow(y+b,2)+pow(z+c,2))-a-x) ))/(( (sqrt(pow(x-a,2)+pow(y+b,2)+pow(-z+c,2))+a-x)/(sqrt(pow(x+a,2)+pow(y+b,2)+pow(-z+c,2))-a-x) )*( (sqrt(pow(x-a,2)+pow(-y+b,2)+pow(z+c,2))+a-x)/(sqrt(pow(x+a,2)+pow(-y+b,2)+pow(z+c,2))-a-x) )) );
    
   double temp_B3 = -factor*( atan(((-x+a)*(y+b))/((z+c)*sqrt(pow(-x+a,2)+pow(y+b,2)+pow(z+c,2)))) + atan(((-x+a)*(y+b))/((-z+c)*sqrt(pow(-x+a,2)+pow(y+b,2)+pow(-z+c,2)))) + atan(((-x+a)*(-y+b))/((z+c)*sqrt(pow(-x+a,2)+pow(-y+b,2)+pow(z+c,2)))) + atan(((-x+a)*(-y+b))/((-z+c)*sqrt(pow(-x+a,2)+pow(-y+b,2)+pow(-z+c,2)))) + atan(((x+a)*(y+b))/((z+c)*sqrt(pow(x+a,2)+pow(y+b,2)+pow(z+c,2)))) + atan(((x+a)*(y+b))/((-z+c)*sqrt(pow(x+a,2)+pow(y+b,2)+pow(-z+c,2)))) + atan(((x+a)*(-y+b))/((z+c)*sqrt(pow(x+a,2)+pow(-y+b,2)+pow(z+c,2)))) + atan(((x+a)*(-y+b))/((-z+c)*sqrt(pow(x+a,2)+pow(-y+b,2)+pow(-z+c,2)))) );
}

bool B_within_margin(double magnetization, double B1, double B2, double B3) {
    bool in_margin = true;
    
    // magnitude of B field from 1 magnet
    double magnitude = sqrt( pow(B1,2) + pow(B2,2) + pow(B3,2) );
    double percent = 0.01;
    double cutoff_value = fabs(percent * magnetization);
    
    if(magnitude < cutoff_value) {
        in_margin = false;
    }
    
    return in_margin;
}

ThreeVec calc_grid_point_B(ThreeVec &grid_point, Magnet &magnet, double mag_dim, double magnetization, char axis) {
    ThreeVec grid_point_B;

    if(axis=='z') {
        // dimensions of permanent magnet
        double b = 0.5*magnet.get_width();   // half width for permanent magnet
        double a = 0.5*magnet.get_length();  // half length for permanent magnet
        double c = 0.5*mag_dim;           // half height for permanent magnet
        
        double offset = 0.5*magnet.get_height() + c;  // moves user defined point for magnet space to the magnet center
        
        double grid_B1 = 0.0;  //B1 and B2 are components perpendicular to magnetization axis
        double grid_B2 = 0.0;  // for z-axis, B1 = x comp and B2 = y comp
        double grid_B3 = 0.0;  // B3 is along axis
        
        int i = 1;
        while (i > -2) {
            double mag_center_x = magnet.get_pos(0) + a;
            double mag_center_y = magnet.get_pos(1);
            double mag_center_z = magnet.get_pos(2) - (i*offset);
            
            double x = grid_point.getX() - mag_center_x;
            double y = grid_point.getY() - mag_center_y;
            double z = grid_point.getZ() - mag_center_z;  // distance from magnet to grid point
            
            double factor = magnetization * pow(10.,-7) * i;
            double temp_B1;
            double temp_B2;
            double temp_B3;
            calc_grid_B_comps(factor, a, b, c, x, y, z);
            
            bool in_margin = B_within_margin(magnetization, temp_B1, temp_B2, temp_B3);
            if(in_margin) {
                grid_B1 += temp_B1;
                grid_B2 += temp_B2;
                grid_B3 += temp_B3;
            }
            i -= 2;
        }
        grid_point_B.setX(grid_B1);
        grid_point_B.setY(grid_B2);
        grid_point_B.setZ(grid_B3);
    }
    else if(axis=='y') {
        double b = 0.5*magnet.get_height();
        double a = 0.5*magnet.get_length();
        double c = 0.5*mag_dim;
        
        double offset = 0.5*magnet.get_width() + c;
        
        double grid_B1 = 0.0;
        double grid_B2 = 0.0;
        double grid_B3 = 0.0;
        
        int i = 1;
        while (i > -2) {
            double mag_center_x = magnet.get_pos(0) + a;
            double mag_center_y = magnet.get_pos(1) - (i*offset);
            double mag_center_z = magnet.get_pos(2);
            
            double x = grid_point.getX() - mag_center_x;
            double y = grid_point.getY() - mag_center_y;
            double z = grid_point.getZ() - mag_center_z;
            
            double factor = magnetization * pow(10.,-7) * i;
            double temp_B1;
            double temp_B2;
            double temp_B3;
            calc_grid_B_comps(factor, a, b, c, x, y, z);
            
            bool in_margin = B_within_margin(magnetization, temp_B1, temp_B2, temp_B3);
            if(in_margin) {
                grid_B1 += temp_B1;
                grid_B2 += temp_B2;
                grid_B3 += temp_B3;
            }
            i -= 2;
        }
        grid_point_B.setX(grid_B1);
        grid_point_B.setY(grid_B3);
        grid_point_B.setZ(grid_B2);
    }
    else {
        // axis = x
        double b = 0.5*magnet.get_width();
        double a = 0.5*magnet.get_height();
        double c = 0.5*mag_dim;
        
        double grid_B1 = 0.0;
        double grid_B2 = 0.0;
        double grid_B3 = 0.0;
        
        int i = 1;
        while (i > -2) {
            // user defined point is not at halfway of the magnet space length
            if(i > 0) {
                double offset = -1*c;
            }
            else {
                double offset = magnet.get_length() + c;
            }
            
            double mag_center_x = magnet.get_pos(0) + offset;
            double mag_center_y = magnet.get_pos(1);
            double mag_center_z = magnet.get_pos(2);
            
            double x = grid_point.getX() - mag_center_x;
            double y = grid_point.getY() - mag_center_y;
            double z = grid_point.getZ() - mag_center_z;
            
            double factor = magnetization * pow(10.,-7) * i;
            double temp_B1;
            double temp_B2;
            double temp_B3;
            calc_grid_B_comps(factor, a, b, c, x, y, z);
            
            bool in_margin = B_within_margin(magnetization, temp_B1, temp_B2, temp_B3);
            if(in_margin) {
                grid_B1 += temp_B1;
                grid_B2 += temp_B2;
                grid_B3 += temp_B3;
            }
            i -= 2;
        }
        grid_point_B.setX(grid_B3);
        grid_point_B.setY(grid_B2);
        grid_point_B.setZ(grid_B1);
    }
    
    return grid_point_B;
}