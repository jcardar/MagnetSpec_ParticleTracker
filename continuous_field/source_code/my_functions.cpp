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
        *(magnet.m_out_magnet) << "Num ," << "Bz(norm) ," << "mag_posx ," << "mag_posy ," << "mag_posz ," << "length ," << "width ," << "height" << "\n";
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
        *(screen.m_out_screen) << "Num ," << "screen_low_energy_edgex ," << "screen_low_energy_edgey ," << "screen_low_energy_edgez ," << "angle ," << "length ," << "height ," << "\n";
    }
    *(screen.m_out_screen) << (++counter) << ","
                            << screen.get_pos(0) << "," << screen.get_pos(1) << "," << screen.get_pos(2) << ","
                            << screen.get_angle('d') << "," << screen.get_length() << "," << screen.get_height() << '\n';
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
                //if(boris_counter != 0)
                //{
                    //outfile_part_comma(particle_t);
                    //std::cerr << "Writing Comma" << std::endl;
                //}
                
                move_particle_to_magnet(magnet_t[index_of_shortest], particle_t);
                std::cerr << "Sending to next magnet" << std::endl;
                outfile_part_commaAndWrite(particle_t);
                //boris_counter++;

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