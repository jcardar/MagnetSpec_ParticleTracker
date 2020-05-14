#include "my_functions.h"

#include <iostream>
#include <fstream>
#include "threevector.h"
#include "threematrix.h"
#include "magnet.h"
#include "screen.h"
#include "particle.h"
#include <cmath>

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
    *(particle.m_out_time) << (particle.get_time()) << ",";
    *(particle.m_out_posx) << (particle.get_pos(0)) << ",";
    *(particle.m_out_posy) << (particle.get_pos(1)) << ",";
    *(particle.m_out_posz) << (particle.get_pos(2)) << ",";
    *(particle.m_out_px) << (particle.get_p(0)) << ",";
    *(particle.m_out_py) << (particle.get_p(1)) << ",";
    *(particle.m_out_pz) << (particle.get_p(2)) << ",";
    *(particle.m_out_energy) << particle.get_energy() << ",";
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
    *(particle.m_out_posy) << ",";
    *(particle.m_out_posz) << ",";
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
	double gaussian_width=1.0; // +/- 1 sigma range
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
    //From v_min to v_max, iterates for each particle
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
void boris(Particle &electron_t, Magnet &magnet_t, const double del_t)
{
  long unsigned parnum,x1;
  double igamma, psquared, Bsquared, alpha;

  // t vector in Boris method
  double tt[3],ttsquared;
  // s vector in Boris method
  double ss[3];
  // vstar vector in Boris method
  double vstar[3];
  // perpendicular component of v in Boris method
  double vperp[3];
      
    Bsquared  = ( ( (magnet_t.get_B0(0))*(magnet_t.get_B0(0)) ) + ( (magnet_t.get_B0(1))*(magnet_t.get_B0(1)) )
                + ( (magnet_t.get_B0(2))*(magnet_t.get_B0(2)) ) );
    
    psquared  = ((electron_t.get_p(0) * electron_t.get_p(0)) + (electron_t.get_p(1) * electron_t.get_p(1)) 
                + (electron_t.get_p(2) * electron_t.get_p(2)));
    igamma    = 1.0/(sqrt(1.0+psquared)); 
    alpha     = sqrt( 1.0/(1 + ((del_t/2)*(del_t/2))) ); //SHOULD DEL_T/2 BE MULTIPLIED BY GAMMA???
    ttsquared = 0.0;
    for (x1=0;x1<3;x1++)
        {
        tt[x1]     = (magnet_t.get_B0(x1) )*del_t*0.5;
        ttsquared += (tt[x1]*tt[x1]);
        }
    
    for (x1=0;x1<3;x1++)
        {
        ss[x1]     = 2.0*tt[x1]/(1.0+ttsquared);
        }
    
    // calculated vperp
    for (x1=0;x1<3;++x1)
        {
        vperp[x1]  = (electron_t.get_p(x1)) * (1.0 - (magnet_t.get_B0(x1))/sqrt(Bsquared));
        }
    
    //calculate vstar
    // for efficiency, component by component
    vstar[0] = vperp[0]+(vperp[1]*tt[2]-vperp[2]*tt[1])*igamma;
    vstar[1] = vperp[1]+(vperp[2]*tt[0]-vperp[0]*tt[2])*igamma;
    vstar[2] = vperp[2]+(vperp[0]*tt[1]-vperp[1]*tt[0])*igamma;

    //Finally update momentum
    /*
    electron_t.set_p(0, ( (electron_t.get_p(0)) + (( (vstar[1]*ss[2]) - (vstar[2]*ss[1]) ) * igamma) ));
    electron_t.set_p(1, ( (electron_t.get_p(1)) + (( (vstar[2]*ss[0]) - (vstar[0]*ss[2]) ) * igamma) ));
    electron_t.set_p(2, ( (electron_t.get_p(2)) + (( (vstar[0]*ss[1]) - (vstar[1]*ss[0]) ) * igamma) ));
    electron_t.set_pos( 0, electron_t.get_pos(0) + (electron_t.get_p(0) * igamma * del_t) );
    electron_t.set_pos( 1, electron_t.get_pos(1) + (electron_t.get_p(1) * igamma * del_t) );
    electron_t.set_pos( 2, electron_t.get_pos(2) + (electron_t.get_p(2) * igamma * del_t) );
    */

    // /*
    double v_fin[3];
    v_fin[0] = ((( (vstar[1]*ss[2]) - (vstar[2]*ss[1]) ) ) + ((1- alpha)*(0) )); //UPDATE 0 IF GRAD B IS PRESENT
    v_fin[1] = ((( (vstar[2]*ss[0]) - (vstar[0]*ss[2]) ) ) + ((1- alpha)*(0) ));
    v_fin[2] = ((( (vstar[0]*ss[1]) - (vstar[1]*ss[0]) ) ) + ((1- alpha)*(0) ));

    electron_t.set_p(0, ( (electron_t.get_p(0)) + (v_fin[0] * igamma))); 
    electron_t.set_p(1, ( (electron_t.get_p(1)) + (v_fin[1] * igamma)));
    electron_t.set_p(2, ( (electron_t.get_p(2)) + (v_fin[2] * igamma)));

    double vpar[3];
    double v_drift[3];
    double v_eff[3];
    for(x1=0;x1<3;x1++)
    {
        v_drift[x1] = 0; //correct if drift is added (non-uniform B)
        //vperp_eff[x1]   = (alpha * ) + ((1-alpha) * v_drift[x1]);
        vperp[x1]   = (electron_t.get_p(x1)) * (1.0 - (magnet_t.get_B0(x1))/sqrt(Bsquared))*igamma;
        vpar[x1]    = (electron_t.get_p(x1)) * ((magnet_t.get_B0(x1))/sqrt(Bsquared))*igamma;
        v_eff[x1]   = vpar[x1] + (alpha*vperp[x1]) + ((1-alpha)*v_drift[x1]);
    }
    electron_t.set_pos( 0, electron_t.get_pos(0) + (v_eff[0] * del_t) );
    electron_t.set_pos( 1, electron_t.get_pos(1) + (v_eff[1] * del_t) );
    electron_t.set_pos( 2, electron_t.get_pos(2) + (v_eff[2] * del_t) );
    // */
    
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void step_through_magnet_mag_boris(Particle &electron, Magnet &magnet, double& time, const double &del_time)
{
    bool check;
    double psquared;
    do
        {
            boris(electron, magnet, del_time); //updates particle velocity in magnetic field 

            time += del_time;

            electron.set_time(time);

            
            check = ((electron.get_pos(0) >= (magnet.get_pos(0))) && (electron.get_pos(0) <= (magnet.get_length()+(magnet.get_pos(0)))) && (electron.get_pos(1) >= ((magnet.get_pos(1))-((magnet.get_width())/2.0))) && (electron.get_pos(1) <= ((magnet.get_pos(1))+(magnet.get_width())/2.0)));
            //std::cerr << "Check Boris: " << check << '\n';
            if(check)
                { outfile_part_writeAndComma(electron); }
            else if(!(check))
                { outfile_part_write(electron); }    
            
            
        } while(check);
}
 