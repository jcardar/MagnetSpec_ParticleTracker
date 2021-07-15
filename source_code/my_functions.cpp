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



void outfile_part_on_screen_first_line(Screen& screen_t)
{
    *(screen_t.m_out_particle_on_screen) << "screen_number," << "particle_number," << "(part_x-screen_x)," << "(part_y-screen_y)," << "(part_z-screen_z)" << "\n";
}

void outfile_part_on_screen(Screen& screen_t, int particle_number_t, Particle& particle_t)
{
    *(screen_t.m_out_particle_on_screen) << screen_t.get_index() << ","
                                        << particle_number_t << ","
                                        << (particle_t.get_pos(0) - screen_t.get_pos(0)) << ","
                                        << (particle_t.get_pos(1) - screen_t.get_pos(1)) << ","
                                        << (particle_t.get_pos(2) - screen_t.get_pos(2)) << "\n";
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Gaussian random distribution function
double gaussian()
{
    //Produces sample of value between -gaussian_width and + gaussian width with FWHM of 1.0
    double c = 1.0/2.35485;
	double gaussian_width=0.12; // +/- 3 sigma range
	double x, y;
	do
	{
		x = ((2.0*rand()/RAND_MAX)-1.0)*gaussian_width;
		y = exp(-x*x*0.5/(c*c));
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
    else 
        { half_divergence = 0; }
    double angle_on_yz_plane;
    if(initial_y != 0)
        { angle_on_yz_plane = atan(initial_z/initial_y); }
    else
        { angle_on_yz_plane = M_PI/2.0; }
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
            rn_plus = rn_minus + (electron->get_vel())*2*del_time;                        //leapfrog algarithm rn for uniform B0

            rn_minus = electron->get_pos();                                               //iterate xn_minus up a timestep
            electron->set_pos(rn_plus);                                                   //do the same to xn

            vn_minus = (electron->get_vel());                                             //iterate vn_minus up a timestep
            electron->set_vel(vn_plus);                                                   //do the same to vn
      
            time += del_time;

            check = ((electron->get_pos(0) >= (magnet.get_pos(0))) && (electron->get_pos(0) <= (magnet.get_length()+(magnet.get_pos(0)))) && (electron->get_pos(1) >= ((magnet.get_pos(1))-((magnet.get_width())/2.0))) && (electron->get_pos(1) <= ((magnet.get_pos(1))+(magnet.get_width())/2.0)));
            
            if(check)
                { outfile_part_writeAndComma(*electron); }
            else if(!(check))
                { outfile_part_write(*electron); }            
            
        } while(check);
        //std::cerr << electron->get_pos() << '\n' << '~' << '\n';
}

double ReadMu0(std::ifstream &input_stream);
double find_magnetization(Magnet &magnet, double mag_dim, double mu_0, char axis);
void calc_grid_B_comps(double factor, double a, double b, double c, double x, double y, double z, double &temp_B1, double &temp_B2, double &temp_B3, bool iterative_add);
ThreeVec calc_dipole_B(ThreeVec &grid_point, Magnet &magnet, double mag_dim, double mu_0, char axis, double magnetization);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void boris_analytic(Particle &electron_t, Magnet magnet_t[], double del_t, int counter, double mu_0, bool inside_array[], int num_mags)
{
    //std::cerr << "Magnet x-Position " << magnet_t.get_pos(0) << '\n';
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
    
    //double B1_atPosition_temp = 0.0;
    //double B2_atPosition_temp = 0.0;
    //double B3_atPosition_temp = 0.0;
    double B1_atPosition = 0.0;
    double B2_atPosition = 0.0;
    double B3_atPosition = 0.0;
    ThreeVec B_field_temp;
    ThreeVec particle_position = electron_t.get_pos();
    double magnetization, factor, x_pos, y_pos, z_pos;
    for(int jj = 0; jj < num_mags; jj++)
    {
        char mag_type = magnet_t[jj].get_type();
        if(inside_array[jj] == true)
        {
            if(mag_type == 'd')
            {
                magnetization = find_magnetization(magnet_t[jj], magnet_t[jj].get_height_of_dipole_block(), mu_0, magnet_t[jj].get_axis_of_magnetization());
                //std::cerr << "Magnitization is " << magnetization << '\n';
                factor = (mu_0 * magnetization)/(4 * M_PI);
                //std::cerr << "factor is " << factor << '\n';
                //Magnitization and factor are staying the same across magnets: they're not the problem
                //x_pos = electron_t.get_pos(0) - (magnet_t[jj].get_pos(0)+0.5*magnet_t[jj].get_length());
                //y_pos = electron_t.get_pos(1) - magnet_t[jj].get_pos(1);
                //z_pos = electron_t.get_pos(2) - magnet_t[jj].get_pos(2);
                //calc_grid_B_comps(factor, 0.5*magnet_t[jj].get_length(), 0.5*magnet_t[jj].get_width(), 0.5*magnet_t[jj].get_height(), x_pos, y_pos, z_pos, B1_atPosition, B2_atPosition, B3_atPosition);
                //ThreeVec calc_grid_point_B(ThreeVec &grid_point, Magnet &magnet, double mag_dim, double mu_0, char axis, double magnetization);
                B_field_temp = calc_dipole_B(particle_position, magnet_t[jj], magnet_t[jj].get_height_of_dipole_block(), mu_0, magnet_t[jj].get_axis_of_magnetization(), magnetization);
                B1_atPosition = B1_atPosition + B_field_temp.getX();
                B2_atPosition = B2_atPosition + B_field_temp.getY();
                B3_atPosition = B3_atPosition + B_field_temp.getZ();
                //std::cerr << "The B1 field at this position was updated as " << B1_atPosition << ", from magnet number " << jj << "\n";
                //std::cerr << "The magnitude B field at this position was updated as " << sqrt(B1_atPosition*B1_atPosition + B2_atPosition*B2_atPosition + B3_atPosition*B3_atPosition) << ", and jj is " << jj << "\n";
            }
            else if(mag_type == 'u')
            {
                if(magnet_t[jj].get_axis_of_magnetization() == 'z')
                {
                    B3_atPosition = B3_atPosition + magnet_t[jj].get_B0(2);
                }
                else if(magnet_t[jj].get_axis_of_magnetization() == 'y')
                {
                    B2_atPosition = B2_atPosition + magnet_t[jj].get_B0(1);
                }
                else if(magnet_t[jj].get_axis_of_magnetization() == 'x')
                {
                    B1_atPosition = B1_atPosition + magnet_t[jj].get_B0(0);
                }
            }
            else if(mag_type == 'q')
            {
                //magnetization = find_magnetization(magnet_t[jj], magnet_t[jj].get_height_of_dipole_block(), mu_0, magnet_t[jj].get_axis_of_magnetization());
                //std::cerr << "Magnitization is " << magnetization << '\n';
                //factor = (mu_0 * magnetization)/(4 * M_PI);
                //std::cerr << "factor is " << factor << '\n';
                //Magnitization and factor are staying the same across magnets: they're not the problem
                //x_pos = electron_t.get_pos(0) - (magnet_t[jj].get_pos(0)+0.5*magnet_t[jj].get_length());
                //y_pos = electron_t.get_pos(1) - magnet_t[jj].get_pos(1);
                //z_pos = electron_t.get_pos(2) - magnet_t[jj].get_pos(2);
                
                //calc_grid_B_comps(factor, 0.5*magnet_t[jj].get_length(), 0.5*magnet_t[jj].get_width(), 0.5*magnet_t[jj].get_height(), x_pos, y_pos, z_pos, B1_atPosition, B2_atPosition, B3_atPosition);
                //ThreeVec calc_grid_point_B(ThreeVec &grid_point, Magnet &magnet, double mag_dim, double mu_0, char axis, double magnetization);
                B_field_temp = calc_quadrupole_B(particle_position, magnet_t[jj], electron_t.get_charge());
                //B_field_temp = calc_dipole_B(particle_position, magnet_t[jj], magnet_t[jj].get_height_of_dipole_block(), mu_0, magnet_t[jj].get_axis_of_magnetization(), magnetization);
                B1_atPosition = B1_atPosition + B_field_temp.getX();
                B2_atPosition = B2_atPosition + B_field_temp.getY();
                B3_atPosition = B3_atPosition + B_field_temp.getZ();
            }
            else if(mag_type == 'h')
            {
                //magnetization = find_magnetization(magnet_t[jj], magnet_t[jj].get_height_of_dipole_block(), mu_0, magnet_t[jj].get_axis_of_magnetization());
                //std::cerr << "Magnitization is " << magnetization << '\n';
                //factor = (mu_0 * magnetization)/(4 * M_PI);
                //std::cerr << "factor is " << factor << '\n';
                //Magnitization and factor are staying the same across magnets: they're not the problem
                //x_pos = electron_t.get_pos(0) - (magnet_t[jj].get_pos(0)+0.5*magnet_t[jj].get_length());
                //y_pos = electron_t.get_pos(1) - magnet_t[jj].get_pos(1);
                //z_pos = electron_t.get_pos(2) - magnet_t[jj].get_pos(2);
                //calc_grid_B_comps(factor, 0.5*magnet_t[jj].get_length(), 0.5*magnet_t[jj].get_width(), 0.5*magnet_t[jj].get_height(), x_pos, y_pos, z_pos, B1_atPosition, B2_atPosition, B3_atPosition);
                //ThreeVec calc_grid_point_B(ThreeVec &grid_point, Magnet &magnet, double mag_dim, double mu_0, char axis, double magnetization);
                B_field_temp = calc_halbach_B(particle_position, magnet_t[jj], electron_t.get_charge());
                //B_field_temp = calc_dipole_B(particle_position, magnet_t[jj], magnet_t[jj].get_height_of_dipole_block(), mu_0, magnet_t[jj].get_axis_of_magnetization(), magnetization);
                B1_atPosition = B1_atPosition + B_field_temp.getX();
                B2_atPosition = B2_atPosition + B_field_temp.getY();
                B3_atPosition = B3_atPosition + B_field_temp.getZ();
                
            }
        }
    }
    Bsquared = ( B1_atPosition*B1_atPosition ) + (B2_atPosition * B2_atPosition) + (B3_atPosition * B3_atPosition);
    double old_px, old_py, old_pz;
    ThreeVec B_at_position(B1_atPosition, B2_atPosition, B3_atPosition);
    //std::cerr << electron_t.get_pos() << '\n';
    //std::cerr << B_at_position << '\n';

        ttsquared = 0.0;
    //if(counter!=0)
    //{
        for (x1=0; x1<3; ++x1)
            {
            tt[x1]     = electron_t.get_charge()*(B_at_position.get(x1) )*del_t*0.5*igamma;
            //tt[x1]     = 1.0*(B_at_position.get(x1) )*del_t*0.5*igamma;
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

    electron_t.set_pos( 0, electron_t.get_pos(0) + (electron_t.get_p(0) * del_t * igamma * 0.5) );
    electron_t.set_pos( 1, electron_t.get_pos(1) + (electron_t.get_p(1) * del_t * igamma * 0.5) );
    electron_t.set_pos( 2, electron_t.get_pos(2) + (electron_t.get_p(2) * del_t * igamma * 0.5) );
}

void boris_analytic_uniform_field(Particle &electron_t, Magnet magnet_t[], double del_t, int counter, double mu_0, bool inside_array[], int num_mags)
{
    //std::cerr << "Magnet x-Position " << magnet_t.get_pos(0) << '\n';
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
    
    //double B1_atPosition_temp = 0.0;
    //double B2_atPosition_temp = 0.0;
    //double B3_atPosition_temp = 0.0;
    double B1_atPosition = 0.0;
    double B2_atPosition = 0.0;
    double B3_atPosition = 0.0;
    ThreeVec B_field_temp;
    ThreeVec particle_position = electron_t.get_pos();
    for(int jj = 0; jj < num_mags; jj++)
    {
        if(inside_array[jj] == true)
        {
            if(magnet_t[jj].get_axis_of_magnetization() == 'z')
            {
                B3_atPosition = B3_atPosition + magnet_t[jj].get_B0(2);
            }
            else if(magnet_t[jj].get_axis_of_magnetization() == 'y')
            {
                B2_atPosition = B2_atPosition + magnet_t[jj].get_B0(1);
            }
            else if(magnet_t[jj].get_axis_of_magnetization() == 'x')
            {
                B1_atPosition = B1_atPosition + magnet_t[jj].get_B0(0);
            }
        }
    }
    //std::cerr << "Magnetic field uniform is (" << B1_atPosition << ", " << B2_atPosition << ", " << B3_atPosition << ")\n";

    Bsquared = ( B1_atPosition*B1_atPosition ) + (B2_atPosition * B2_atPosition) + (B3_atPosition * B3_atPosition);
    double old_px, old_py, old_pz;
    ThreeVec B_at_position(B1_atPosition, B2_atPosition, B3_atPosition);

        ttsquared = 0.0;

        for (x1=0; x1<3; ++x1)
            {
            tt[x1]     = electron_t.get_charge()*(B_at_position.get(x1) )*del_t*0.5*igamma;
            //tt[x1]     = 1.0*(B_at_position.get(x1) )*del_t*0.5*igamma;
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

    electron_t.set_pos( 0, electron_t.get_pos(0) + (electron_t.get_p(0) * del_t * igamma * 0.5) );
    electron_t.set_pos( 1, electron_t.get_pos(1) + (electron_t.get_p(1) * del_t * igamma * 0.5) );
    electron_t.set_pos( 2, electron_t.get_pos(2) + (electron_t.get_p(2) * del_t * igamma * 0.5) );
}


void boris_analytic_ONEMAGNET(Particle &electron_t, Magnet &magnet_t, double del_t, int counter, double mu_0)
{
    //std::cerr << "Magnet x-Position " << magnet_t.get_pos(0) << '\n';
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

    //Bsquared  = ( ( (magnet_t.get_B0(0))*(magnet_t.get_B0(0)) ) + ( (magnet_t.get_B0(1))*(magnet_t.get_B0(1)) )
    //           + ( (magnet_t.get_B0(2))*(magnet_t.get_B0(2)) ) );
    
    double B1_atPosition, B2_atPosition, B3_atPosition;
    double magnetization = find_magnetization(magnet_t, magnet_t.get_height_of_dipole_block(), mu_0, magnet_t.get_axis_of_magnetization());
    //std::cerr << "Magnitization is " << magnetization << '\n';
    double factor = (mu_0 * magnetization)/(4 * M_PI);
    //std::cerr << "factor is " << factor << '\n';
    //Magnitization and factor are staying the same across magnets: they're not the problem
    double x_pos = electron_t.get_pos(0) - (magnet_t.get_pos(0)+0.5*magnet_t.get_length());
    double y_pos = electron_t.get_pos(1) - magnet_t.get_pos(1);
    double z_pos = electron_t.get_pos(2) - magnet_t.get_pos(2);
    calc_grid_B_comps(factor, 0.5*magnet_t.get_length(), 0.5*magnet_t.get_width(), 0.5*magnet_t.get_height(), x_pos, y_pos, z_pos, B1_atPosition, B2_atPosition, B3_atPosition);

    Bsquared = ( B1_atPosition*B1_atPosition ) + (B2_atPosition * B2_atPosition) + (B3_atPosition * B3_atPosition);
    //std::cerr << "B-squared is "<< Bsquared << '\n';
    //Bsquared is constantly decreasing!!!
    double old_px, old_py, old_pz;
    ThreeVec B_at_position(B1_atPosition, B2_atPosition, B3_atPosition);

        ttsquared = 0.0;
    //if(counter!=0)
    //{
        for (x1=0; x1<3; ++x1)
            {
            tt[x1]     = electron_t.get_charge()*(B_at_position.get(x1) )*del_t*0.5*igamma;
            //tt[x1]     = 1.0*(B_at_position.get(x1) )*del_t*0.5*igamma;
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

    //ThreeVec average_vel((old_px+electron_t.get_p(0))*igamma*0.5, (old_py+electron_t.get_p(1))*igamma*0.5, (old_pz+electron_t.get_p(2))*igamma*0.5);
    
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
bool inside_of_mag_dipole(Magnet magnet_t, Particle particle_t);
bool is_this_array_only_zeros(bool array[], int num_elements);
bool check_if_intersect_screen_in_next_dt(Screen screen_t, Particle particle_t, double del_t);


bool step_through_magnet_mag_boris_analytic_general(Particle &electron, Magnet magnet[], double& time, const double &del_time, double mu_0, double time_out, bool inside_of_magnet[], int num_mags, Screen screen[], int num_screens, bool supress_output=false)
{
    //RETURNS TRUE IF PARTICLE EXITS MAGNET WITHOUT INTERSECTING A SCREEN
    //RETURNS FALSE IF PARTICLE INTERSECTS SCREEN WHILE INSIDE MAGNETIC REGION
    //TESTING CHANING TIME STEP:
    //const double del_time = del_time_1 * 0.5;

    //bool check_x, check_y, check_z;
    double psquared;
    int counter = 0;
    bool inside_of_magnet_check = true;
    bool intersect_screen_array[num_screens];
    bool intersect_screen_check = false;
    do
        {
            //check if intersects screen in next dt first, then push particle, then determine if particle is still inside the magnetic field
            for(int jj = 0; jj < num_screens; jj++)
            {
                intersect_screen_array[jj] = check_if_intersect_screen_in_next_dt(screen[jj], electron, del_time);
                if(intersect_screen_array[jj] == true)
                {
                    intersect_screen_check = true;
                    break;
                }
            }
            if(intersect_screen_check == true)
            {
                break;
            }

            boris_analytic(electron, magnet, del_time, counter, mu_0, inside_of_magnet, num_mags); //updates particle velocity & position in magnetic field 
            
            counter++;
            time += del_time;
            electron.set_time(time);

            for(int ii=0; ii<num_mags; ii++)
            {
                inside_of_magnet[ii] = inside_of_mag_general(magnet[ii], electron);
            }
            
            inside_of_magnet_check = !(is_this_array_only_zeros(inside_of_magnet, num_mags));

            //if((inside_of_magnet_check) && !(time>=time_out))
             //   { 
            outfile_part_commaAndWrite(electron); 
            //    }
            //else if( ( !inside_of_magnet_check ) || (time >= time_out) )
            //    {
            //        outfile_part_commaAndWrite(electron); 
            //    }
            if(supress_output)
            {
                supress_output = true;
            }
            else
            {
                if(time >= time_out)
                    { std::cout << "Particle Timed-Out.\n"; }
                if( (inside_of_magnet_check) == false )
                    {
                        std::cout << "Particle Out of Magnet.\n";
                        //if(check_x == false)
                        //    { std::cout << "Out of bounds in x.\n"; } 
                        //if(check_y == false)
                        //    { std::cout << "Out of bounds in y.\n"; } 
                        //if(check_z == false)
                        //   { std::cout << "Out of bounds in z.\n"; } 
                    }
            }
            
        } while((inside_of_magnet_check) && (time < time_out));
    if(inside_of_magnet_check == false)
    {
        for(int ii = 0; ii < num_mags; ii++)
        {
            char mag_type = magnet[ii].get_type();
            if(mag_type == 'd')
            {
                bool inside_x_limits = (electron.get_pos(0) >= (magnet[ii].get_pos(0)- magnet[ii].get_length() )) && (electron.get_pos(0) <= (magnet[ii].get_pos(0)+2.0*magnet[ii].get_length()));
                bool inside_y_limits = (electron.get_pos(1) >= (magnet[ii].get_pos(1) - (magnet[ii].get_width()*5.0)))  && (electron.get_pos(1) <= (magnet[ii].get_pos(1) + (magnet[ii].get_width()*5.0)));
                bool inside_z_limits = (electron.get_pos(2) >= (magnet[ii].get_pos(2) - (magnet[ii].get_height()/2.0))) && (electron.get_pos(2) <= (magnet[ii].get_pos(2) + (magnet[ii].get_height()/2.0)));
                if(inside_x_limits && inside_y_limits && (!inside_z_limits))
                {
                    return false;
                }
            }
            else if(mag_type == 'u')
            {
                bool inside_x_limits = (electron.get_pos(0) >= (magnet[ii].get_pos(0) )) && (electron.get_pos(0) <= (magnet[ii].get_pos(0)+magnet[ii].get_length()));
                bool inside_y_limits = (electron.get_pos(1) >= (magnet[ii].get_pos(1) - (magnet[ii].get_width()/2.0)))  && (electron.get_pos(1) <= (magnet[ii].get_pos(1) + (magnet[ii].get_width()*0.5)));
                bool inside_z_limits = (electron.get_pos(2) >= (magnet[ii].get_pos(2) - (magnet[ii].get_height()/2.0))) && (electron.get_pos(2) <= (magnet[ii].get_pos(2) + (magnet[ii].get_height()/2.0)));
                if(inside_x_limits && inside_y_limits && (!inside_z_limits))
                {
                    return false;
                }
            }
            else if(mag_type == 'q' || mag_type == 'h')
            {
                bool inside_x_limits = (electron.get_pos(0) >= (magnet[ii].get_pos(0) - magnet[ii].get_length() )) && (electron.get_pos(0) <= (magnet[ii].get_pos(0)+2.0*magnet[ii].get_length()));
                bool inside_y_limits = (electron.get_pos(1) >= (magnet[ii].get_pos(1) - (magnet[ii].get_width())))  && (electron.get_pos(1) <= (magnet[ii].get_pos(1) + (magnet[ii].get_width())));
                bool inside_z_limits = (electron.get_pos(2) >= (magnet[ii].get_pos(2) - (magnet[ii].get_width()))) && (electron.get_pos(2) <= (magnet[ii].get_pos(2) + (magnet[ii].get_width())));
                if(inside_x_limits && inside_y_limits && (!inside_z_limits))
                {
                    return false;
                }
            }
        }
    }
    return !(intersect_screen_check);
}



bool step_through_magnet_mag_boris_analytic_dipole(Particle &electron, Magnet magnet[], double& time, const double &del_time, double mu_0, double time_out, bool inside_of_magnet[], int num_mags, Screen screen[], int num_screens, bool supress_output=false)
{
    //RETURNS TRUE IF PARTICLE EXITS MAGNET WITHOUT INTERSECTING A SCREEN
    //RETURNS FALSE IF PARTICLE INTERSECTS SCREEN WHILE INSIDE MAGNETIC REGION
    //TESTING CHANING TIME STEP:
    //const double del_time = del_time_1 * 0.5;

    //bool check_x, check_y, check_z;
    double psquared;
    int counter = 0;
    bool inside_of_magnet_check = true;
    bool intersect_screen_array[num_screens];
    bool intersect_screen_check = false;
    do
        {
            //check if intersects screen in next dt first, then push particle, then determine if particle is still inside the magnetic field
            for(int jj = 0; jj < num_screens; jj++)
            {
                intersect_screen_array[jj] = check_if_intersect_screen_in_next_dt(screen[jj], electron, del_time);
                if(intersect_screen_array[jj] == true)
                {
                    intersect_screen_check = true;
                    break;
                }
            }
            if(intersect_screen_check == true)
            {
                break;
            }

            boris_analytic(electron, magnet, del_time, counter, mu_0, inside_of_magnet, num_mags); //updates particle velocity & position in magnetic field 
            
            counter++;
            time += del_time;
            electron.set_time(time);

            for(int ii=0; ii<num_mags; ii++)
            {
                inside_of_magnet[ii] = inside_of_mag_dipole(magnet[ii], electron);
            }
            
            inside_of_magnet_check = !(is_this_array_only_zeros(inside_of_magnet, num_mags));

            //if((inside_of_magnet_check) && !(time>=time_out))
             //   { 
            outfile_part_commaAndWrite(electron); 
            //    }
            //else if( ( !inside_of_magnet_check ) || (time >= time_out) )
            //    {
            //        outfile_part_commaAndWrite(electron); 
            //    }
            if(supress_output)
            {
                supress_output = true;
            }
            else
            {
                if(time >= time_out)
                    { std::cout << "Particle Timed-Out.\n"; }
                if( (inside_of_magnet_check) == false )
                    {
                        std::cout << "Particle Out of Magnet.\n";
                        //if(check_x == false)
                        //    { std::cout << "Out of bounds in x.\n"; } 
                        //if(check_y == false)
                        //    { std::cout << "Out of bounds in y.\n"; } 
                        //if(check_z == false)
                        //   { std::cout << "Out of bounds in z.\n"; } 
                    }
            }
            
        } while((inside_of_magnet_check) && (time < time_out));
    if(inside_of_magnet_check == false)
    {
        for(int ii = 0; ii < num_mags; ii++)
        {
            bool inside_x_limits = (electron.get_pos(0) >= (magnet[ii].get_pos(0)- magnet[ii].get_length() )) && (electron.get_pos(0) <= (magnet[ii].get_pos(0)+2.0*magnet[ii].get_length()));
            bool inside_y_limits = (electron.get_pos(1) >= (magnet[ii].get_pos(1) - (magnet[ii].get_width()*5.0)))  && (electron.get_pos(1) <= (magnet[ii].get_pos(1) + (magnet[ii].get_width()*5.0)));
            bool inside_z_limits = (electron.get_pos(2) >= (magnet[ii].get_pos(2) - (magnet[ii].get_height()/2.0))) && (electron.get_pos(2) <= (magnet[ii].get_pos(2) + (magnet[ii].get_height()/2.0)));
            if(inside_x_limits && inside_y_limits && (!inside_z_limits))
            {
                return false;
            }
        }
    }
    return !(intersect_screen_check);
}

bool inside_of_mag_uniform(Magnet magnet_t, Particle particle_t);
bool step_through_magnet_mag_boris_analytic_uniform_field(Particle &electron, Magnet magnet[], double& time, const double &del_time, double mu_0, double time_out, bool inside_of_magnet[], int num_mags, Screen screen[], int num_screens, bool supress_output)
{
    //RETURNS TRUE IF PARTICLE EXITS MAGNET WITHOUT INTERSECTING A SCREEN
    //RETURNS FALSE IF PARTICLE INTERSECTS SCREEN WHILE INSIDE MAGNETIC REGION
    //TESTING CHANING TIME STEP:
    //const double del_time = del_time_1 * 0.5;
    //std::cerr << "Inside a uniform field magnet\n";
    //bool check_x, check_y, check_z;
    double psquared;
    int counter = 0;
    bool inside_of_magnet_check = true;
    bool intersect_screen_array[num_screens];
    bool intersect_screen_check = false;
    do
        {
            //check if intersects screen in next dt first, then push particle, then determine if particle is still inside the magnetic field
            for(int jj = 0; jj < num_screens; jj++)
            {
                intersect_screen_array[jj] = check_if_intersect_screen_in_next_dt(screen[jj], electron, del_time);
                if(intersect_screen_array[jj] == true)
                {
                    intersect_screen_check = true;
                    break;
                }
            }
            if(intersect_screen_check == true)
            {
                break;
            }

            boris_analytic_uniform_field(electron, magnet, del_time, counter, mu_0, inside_of_magnet, num_mags); //updates particle velocity & position in magnetic field 
            
            counter++;
            time += del_time;
            electron.set_time(time);

            for(int ii=0; ii<num_mags; ii++)
            {
                inside_of_magnet[ii] = inside_of_mag_uniform(magnet[ii],electron);
            }
            
            inside_of_magnet_check = !(is_this_array_only_zeros(inside_of_magnet, num_mags));

            outfile_part_commaAndWrite(electron); 

            if(supress_output)
            {
                supress_output = true;
            }
            else
            {
                if(time >= time_out)
                    { std::cout << "Particle Timed-Out.\n"; }
                if( (inside_of_magnet_check) == false )
                    {
                        std::cout << "Particle Out of Magnet.\n";
                        //if(check_x == false)
                        //    { std::cout << "Out of bounds in x.\n"; } 
                        //if(check_y == false)
                        //    { std::cout << "Out of bounds in y.\n"; } 
                        //if(check_z == false)
                        //   { std::cout << "Out of bounds in z.\n"; } 
                    }
            }
            
        } while((inside_of_magnet_check) && (time < time_out));
    if(inside_of_magnet_check == false)
    {
        for(int ii = 0; ii < num_mags; ii++)
        {
            bool inside_x_limits = (electron.get_pos(0) >= (magnet[ii].get_pos(0) )) && (electron.get_pos(0) <= (magnet[ii].get_pos(0)+magnet[ii].get_length()));
            bool inside_y_limits = (electron.get_pos(1) >= (magnet[ii].get_pos(1) - (magnet[ii].get_width()*0.5)))  && (electron.get_pos(1) <= (magnet[ii].get_pos(1) + (magnet[ii].get_width()*0.5)));
            bool inside_z_limits = (electron.get_pos(2) >= (magnet[ii].get_pos(2) - (magnet[ii].get_height()/2.0))) && (electron.get_pos(2) <= (magnet[ii].get_pos(2) + (magnet[ii].get_height()/2.0)));
            if(inside_x_limits && inside_y_limits && (!inside_z_limits))
            {
                return false;
            }
        }
    }
    return !(intersect_screen_check);
}

void step_through_magnet_mag_boris_analytic_ONEMAGNET(Particle &electron, Magnet magnet[], double& time, const double &del_time, double mu_0, double time_out, int triggered_mag_number, int num_mags, bool supress_output = true)
{
    //TESTING CHANING TIME STEP:
    //const double del_time = del_time_1 * 0.5;

    bool check_x, check_y, check_z;
    double psquared;
    int counter = 0;
    do
        {
            boris_analytic_ONEMAGNET(electron, magnet[triggered_mag_number], del_time, counter, mu_0); //updates particle velocity & position in magnetic field 
            //for(int ii=0; ii<num_mags; ii++)
            //{
                //if(ii != triggered_mag_number)
            //    {
            //        bool inside_of_magnet = inside_of_mag(magnet[ii], electron);
            //        if(inside_of_mag)
            //        {
            //            boris_analytic(electron, magnet[ii], del_time, counter, mu_0);
            //        }
            //    }
            //}
            counter++;

            time += del_time;
//            std::cerr << del_time << std::endl;

            electron.set_time(time);

            
            //check_x = (electron.get_pos(0) >= (magnet.get_pos(0))) && (electron.get_pos(0) <= (magnet.get_length()+(magnet.get_pos(0))));
            //check_y = (electron.get_pos(1) >= ((magnet.get_pos(1))-((magnet.get_width())/2.0)))  && (electron.get_pos(1) <= ((magnet.get_pos(1))+(magnet.get_width())/2.0));
            //check_z = (electron.get_pos(2) >= ((magnet.get_pos(2))-((magnet.get_height())/2.0))) && (electron.get_pos(2) <= ((magnet.get_pos(2))+(magnet.get_height())/2.0));
            check_x = (electron.get_pos(0) >= (magnet[triggered_mag_number].get_pos(0)-1.0*magnet[triggered_mag_number].get_length())) && (electron.get_pos(0) <= (2.0*magnet[triggered_mag_number].get_length()+(magnet[triggered_mag_number].get_pos(0))));
            //^ LENGTH TAG
            check_y = (electron.get_pos(1) >= ((magnet[triggered_mag_number].get_pos(1))-((100.0*magnet[triggered_mag_number].get_width())/2.0)))  && (electron.get_pos(1) <= ((magnet[triggered_mag_number].get_pos(1))+(100.0*magnet[triggered_mag_number].get_width())/2.0));
            //^ WIDTH TAG
            check_z = (electron.get_pos(2) >= ((magnet[triggered_mag_number].get_pos(2))-((magnet[triggered_mag_number].get_height())/2.0))) && (electron.get_pos(2) <= ((magnet[triggered_mag_number].get_pos(2))+(magnet[triggered_mag_number].get_height())/2.0));
            //if(counter%10==0)
            //    {std::cerr << counter << '\n'; std::cerr << electron.get_p(0) << '\n';}
            if((check_x && check_y && check_z) && !(time>=time_out))
                { outfile_part_commaAndWrite(electron); }
            else if( (!( check_x && check_y && check_z )) || (time >= time_out) )
                {
                    outfile_part_commaAndWrite(electron); 
                }
            if(supress_output)
            {
                supress_output = true;
            }
            else
            {
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


bool inside_of_mag_general(Magnet magnet_t, Particle particle_t)
{
    char mag_type = magnet_t.get_type();
    if(mag_type == 'd')
    {
        bool inside_of_mag;

        bool inside_x_limits = (particle_t.get_pos(0) >= (magnet_t.get_pos(0) - 1.0*magnet_t.get_length())) && (particle_t.get_pos(0) <= (magnet_t.get_pos(0)+2.0*magnet_t.get_length()));
        //^ LENGTH TAG
        bool inside_y_limits = (particle_t.get_pos(1) >= (magnet_t.get_pos(1) - (magnet_t.get_width()*5.0)))  && (particle_t.get_pos(1) <= (magnet_t.get_pos(1) + (magnet_t.get_width()*5.0)));
        //^ WIDTH TAG
        bool inside_z_limits = (particle_t.get_pos(2) >= (magnet_t.get_pos(2) - (magnet_t.get_height()/2.0))) && (particle_t.get_pos(2) <= (magnet_t.get_pos(2) + (magnet_t.get_height()/2.0)));
        if(inside_x_limits && inside_y_limits && inside_z_limits)
        {
            inside_of_mag = true;
            //std::cerr << "Inside of magnet of type " << magnet_t.get_type() << std::endl; 
        }
        else
        {
            inside_of_mag = false;
        }
        
        return inside_of_mag;
    }
    else if(mag_type == 'u')
    {
        bool inside_of_mag;

        bool inside_x_limits = (particle_t.get_pos(0) >= (magnet_t.get_pos(0) )) && (particle_t.get_pos(0) <= (magnet_t.get_pos(0)+magnet_t.get_length()));
        //^ LENGTH TAG
        bool inside_y_limits = (particle_t.get_pos(1) >= (magnet_t.get_pos(1) - (magnet_t.get_width()*0.5)))  && (particle_t.get_pos(1) <= (magnet_t.get_pos(1) + (magnet_t.get_width()*0.5)));
        //^ WIDTH TAG
        bool inside_z_limits = (particle_t.get_pos(2) >= (magnet_t.get_pos(2) - (magnet_t.get_height()/2.0))) && (particle_t.get_pos(2) <= (magnet_t.get_pos(2) + (magnet_t.get_height()/2.0)));
        if(inside_x_limits && inside_y_limits && inside_z_limits)
        {
            inside_of_mag = true;
            //std::cerr << "Inside of magnet of type " << magnet_t.get_type() << std::endl; 
        }
        else
        {
            inside_of_mag = false;
        }
        
        return inside_of_mag;
    }
    else if(mag_type == 'q' || mag_type == 'h')
    {
        bool inside_of_mag;

        bool inside_x_limits = (particle_t.get_pos(0) >= (magnet_t.get_pos(0) - 1.0*magnet_t.get_length())) && (particle_t.get_pos(0) <= (magnet_t.get_pos(0)+2.0*magnet_t.get_length()));
        //^ LENGTH TAG
        bool inside_y_limits = (particle_t.get_pos(1) >= (magnet_t.get_pos(1) - (magnet_t.get_width())))  && (particle_t.get_pos(1) <= (magnet_t.get_pos(1) + (magnet_t.get_width())));
        //^ WIDTH TAG
        bool inside_z_limits = (particle_t.get_pos(2) >= (magnet_t.get_pos(2) - (magnet_t.get_width()))) && (particle_t.get_pos(2) <= (magnet_t.get_pos(2) + (magnet_t.get_width())));
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
    else
    {
        std::cout << "No valid magnet type detected!\n";
        return 0;
    }
}


bool inside_of_mag_dipole(Magnet magnet_t, Particle particle_t)
{
    bool inside_of_mag;
    
    bool inside_x_limits = (particle_t.get_pos(0) >= (magnet_t.get_pos(0) - 1.0*magnet_t.get_length())) && (particle_t.get_pos(0) <= (magnet_t.get_pos(0)+2.0*magnet_t.get_length()));
    //^ LENGTH TAG
    bool inside_y_limits = (particle_t.get_pos(1) >= (magnet_t.get_pos(1) - (magnet_t.get_width()*5.0)))  && (particle_t.get_pos(1) <= (magnet_t.get_pos(1) + (magnet_t.get_width()*5.0)));
    //^ WIDTH TAG
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

bool inside_of_mag_uniform(Magnet magnet_t, Particle particle_t)
{
    bool inside_of_mag;
    
    bool inside_x_limits = (particle_t.get_pos(0) >= (magnet_t.get_pos(0) )) && (particle_t.get_pos(0) <= (magnet_t.get_pos(0)+magnet_t.get_length()));
    //^ LENGTH TAG
    bool inside_y_limits = (particle_t.get_pos(1) >= (magnet_t.get_pos(1) - (magnet_t.get_width()*0.5)))  && (particle_t.get_pos(1) <= (magnet_t.get_pos(1) + (magnet_t.get_width()*0.5)));
    //^ WIDTH TAG
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


double time_to_magnet_boundary_general(Magnet magnet_t, Particle particle_t)
{
    //Assumes that particle has already been checked if it is at-or-within magnet's boundary
    char mag_type = magnet_t.get_type();
    if(mag_type == 'd')
    {
        double time_btwn_mags_x_front;
        double time_btwn_mags_x_back;
        if(particle_t.get_vel(0) !=0)
            { 
                //LENGTH TAG
                time_btwn_mags_x_front = (magnet_t.get_pos(0)-1.0*magnet_t.get_length() - particle_t.get_pos(0))/particle_t.get_vel(0); 
                time_btwn_mags_x_back  = (magnet_t.get_pos(0)+2.0*magnet_t.get_length() - particle_t.get_pos(0))/particle_t.get_vel(0);
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
                ///WIDTH TAG
                time_btwn_mags_y_top     = (magnet_t.get_pos(1)+((magnet_t.get_width()*5.0)) - particle_t.get_pos(1))/particle_t.get_vel(1); 
                time_btwn_mags_y_bottom  = (magnet_t.get_pos(1)-((magnet_t.get_width()*5.0)) - particle_t.get_pos(1))/particle_t.get_vel(1);
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
    else if(mag_type == 'u')
    {
        double time_btwn_mags_x_front;
        double time_btwn_mags_x_back;
        if(particle_t.get_vel(0) !=0)
            { 
                //LENGTH TAG
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
                ///WIDTH TAG
                time_btwn_mags_y_top     = (magnet_t.get_pos(1)+((magnet_t.get_width()*0.5)) - particle_t.get_pos(1))/particle_t.get_vel(1); 
                time_btwn_mags_y_bottom  = (magnet_t.get_pos(1)-((magnet_t.get_width()*0.5)) - particle_t.get_pos(1))/particle_t.get_vel(1);
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
    else if(mag_type == 'q')
    {
        double time_btwn_mags_x_front;
        double time_btwn_mags_x_back;
        if(particle_t.get_vel(0) !=0)
            { 
                //LENGTH TAG
                time_btwn_mags_x_front = (magnet_t.get_pos(0)-1.0*magnet_t.get_length() - particle_t.get_pos(0))/particle_t.get_vel(0); 
                time_btwn_mags_x_back  = (magnet_t.get_pos(0)+2.0*magnet_t.get_length() - particle_t.get_pos(0))/particle_t.get_vel(0);
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
                ///WIDTH TAG
                time_btwn_mags_y_top     = (magnet_t.get_pos(1)+((magnet_t.get_width())) - particle_t.get_pos(1))/particle_t.get_vel(1); 
                time_btwn_mags_y_bottom  = (magnet_t.get_pos(1)-((magnet_t.get_width())) - particle_t.get_pos(1))/particle_t.get_vel(1);
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
                time_btwn_mags_z_top     = (magnet_t.get_pos(2)+(magnet_t.get_width()) - particle_t.get_pos(2))/particle_t.get_vel(2); 
                time_btwn_mags_z_bottom  = (magnet_t.get_pos(2)-(magnet_t.get_width()) - particle_t.get_pos(2))/particle_t.get_vel(2);
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
    else if(mag_type == 'h')
    {
        double time_btwn_mags_x_front;
        double time_btwn_mags_x_back;
        if(particle_t.get_vel(0) !=0)
            { 
                //LENGTH TAG
                time_btwn_mags_x_front = (magnet_t.get_pos(0)-1.0*magnet_t.get_length() - particle_t.get_pos(0))/particle_t.get_vel(0); 
                time_btwn_mags_x_back  = (magnet_t.get_pos(0)+2.0*magnet_t.get_length() - particle_t.get_pos(0))/particle_t.get_vel(0);
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
                ///WIDTH TAG
                time_btwn_mags_y_top     = (magnet_t.get_pos(1)+((magnet_t.get_width())) - particle_t.get_pos(1))/particle_t.get_vel(1); 
                time_btwn_mags_y_bottom  = (magnet_t.get_pos(1)-((magnet_t.get_width())) - particle_t.get_pos(1))/particle_t.get_vel(1);
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
                time_btwn_mags_z_top     = (magnet_t.get_pos(2)+(magnet_t.get_width()) - particle_t.get_pos(2))/particle_t.get_vel(2); 
                time_btwn_mags_z_bottom  = (magnet_t.get_pos(2)-(magnet_t.get_width()) - particle_t.get_pos(2))/particle_t.get_vel(2);
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
    else 
    {std::cout << "No valid magnet type detected.\n"; return -1.0;}
}


double time_to_magnet_boundary_dipole(Magnet magnet_t, Particle particle_t)
{
    //Assumes that particle has already been checked if it is at-or-within magnet's boundary

    double time_btwn_mags_x_front;
    double time_btwn_mags_x_back;
    if(particle_t.get_vel(0) !=0)
        { 
            //LENGTH TAG
            time_btwn_mags_x_front = (magnet_t.get_pos(0)-1.0*magnet_t.get_length() - particle_t.get_pos(0))/particle_t.get_vel(0); 
            time_btwn_mags_x_back  = (magnet_t.get_pos(0)+2.0*magnet_t.get_length() - particle_t.get_pos(0))/particle_t.get_vel(0);
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
            ///WIDTH TAG
            time_btwn_mags_y_top     = (magnet_t.get_pos(1)+((magnet_t.get_width()*5.0)) - particle_t.get_pos(1))/particle_t.get_vel(1); 
            time_btwn_mags_y_bottom  = (magnet_t.get_pos(1)-((magnet_t.get_width()*5.0)) - particle_t.get_pos(1))/particle_t.get_vel(1);
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

double time_to_magnet_boundary_uniform(Magnet magnet_t, Particle particle_t)
{
    //Assumes that particle has already been checked if it is at-or-within magnet's boundary

    double time_btwn_mags_x_front;
    double time_btwn_mags_x_back;
    if(particle_t.get_vel(0) !=0)
        { 
            //LENGTH TAG
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
            ///WIDTH TAG
            time_btwn_mags_y_top     = (magnet_t.get_pos(1)+((magnet_t.get_width()*0.5)) - particle_t.get_pos(1))/particle_t.get_vel(1); 
            time_btwn_mags_y_bottom  = (magnet_t.get_pos(1)-((magnet_t.get_width()*0.5)) - particle_t.get_pos(1))/particle_t.get_vel(1);
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



bool intersect_mag_general(Magnet magnet_t, Particle particle_t)
{
    char mag_type = magnet_t.get_type();
    if(mag_type == 'd')
    {
        bool intersect;
        double time_to_magnet = time_to_magnet_boundary_dipole(magnet_t, particle_t);
        if(time_to_magnet < 0.0)
        {
            return false;
        }
        ThreeVec pos_at_shortest_time;
        pos_at_shortest_time.setX(particle_t.get_pos(0) + (time_to_magnet * particle_t.get_vel(0)));
        pos_at_shortest_time.setY(particle_t.get_pos(1) + (time_to_magnet * particle_t.get_vel(1)));
        pos_at_shortest_time.setZ(particle_t.get_pos(2) + (time_to_magnet * particle_t.get_vel(2)));

        bool within_x_bounds = (pos_at_shortest_time.getX() >= magnet_t.get_pos(0)-1.0*magnet_t.get_length()) && (pos_at_shortest_time.getX() <= (magnet_t.get_pos(0)+2.0*magnet_t.get_length()));
        //^ LENGTH TAG: edited to make width region 5*width away from center instead of 0.5*width.
        bool within_y_bounds = (pos_at_shortest_time.getY() >= (magnet_t.get_pos(1)-(magnet_t.get_width()*5.0)))  && (pos_at_shortest_time.getY() <= (magnet_t.get_pos(1)+(magnet_t.get_width()*5.0)));
        //^ WIDTH TAG
        bool within_z_bounds = (pos_at_shortest_time.getZ() >= (magnet_t.get_pos(2)-(magnet_t.get_height()/2.0))) && (pos_at_shortest_time.getZ() <= (magnet_t.get_pos(2)+(magnet_t.get_height()/2.0)));

        if(within_x_bounds && within_y_bounds && within_z_bounds)
        { intersect = true; }
        else
        { intersect = false; }
    
        return intersect;
    }
    else if(mag_type == 'u')
    {
        bool intersect;
        double time_to_magnet = time_to_magnet_boundary_uniform(magnet_t, particle_t);
        if(time_to_magnet < 0.0)
        {
            return false;
        }
        ThreeVec pos_at_shortest_time;
        pos_at_shortest_time.setX(particle_t.get_pos(0) + (time_to_magnet * particle_t.get_vel(0)));
        pos_at_shortest_time.setY(particle_t.get_pos(1) + (time_to_magnet * particle_t.get_vel(1)));
        pos_at_shortest_time.setZ(particle_t.get_pos(2) + (time_to_magnet * particle_t.get_vel(2)));

        bool within_x_bounds = (pos_at_shortest_time.getX() >= magnet_t.get_pos(0)) && (pos_at_shortest_time.getX() <= (magnet_t.get_pos(0)+magnet_t.get_length()));
        //^ LENGTH TAG: edited to make width region 5*width away from center instead of 0.5*width.
        bool within_y_bounds = (pos_at_shortest_time.getY() >= (magnet_t.get_pos(1)-(magnet_t.get_width()*0.5)))  && (pos_at_shortest_time.getY() <= (magnet_t.get_pos(1)+(magnet_t.get_width()*0.5)));
        //^ WIDTH TAG
        bool within_z_bounds = (pos_at_shortest_time.getZ() >= (magnet_t.get_pos(2)-(magnet_t.get_height()/2.0))) && (pos_at_shortest_time.getZ() <= (magnet_t.get_pos(2)+(magnet_t.get_height()/2.0)));

        if(within_x_bounds && within_y_bounds && within_z_bounds)
        { intersect = true; }
        else
        { intersect = false; }
    
        return intersect;
    }
    else if(mag_type == 'q' || mag_type =='h')
    {
        bool intersect;
        double time_to_magnet = time_to_magnet_boundary_dipole(magnet_t, particle_t);
        if(time_to_magnet < 0.0)
        {
            return false;
        }
        ThreeVec pos_at_shortest_time;
        pos_at_shortest_time.setX(particle_t.get_pos(0) + (time_to_magnet * particle_t.get_vel(0)));
        pos_at_shortest_time.setY(particle_t.get_pos(1) + (time_to_magnet * particle_t.get_vel(1)));
        pos_at_shortest_time.setZ(particle_t.get_pos(2) + (time_to_magnet * particle_t.get_vel(2)));

        bool within_x_bounds = (pos_at_shortest_time.getX() >= magnet_t.get_pos(0)-1.0*magnet_t.get_length()) && (pos_at_shortest_time.getX() <= (magnet_t.get_pos(0)+2.0*magnet_t.get_length()));
        //^ LENGTH TAG: edited to make width region 5*width away from center instead of 0.5*width.
        bool within_y_bounds = (pos_at_shortest_time.getY() >= (magnet_t.get_pos(1)-(magnet_t.get_width()))) && (pos_at_shortest_time.getY() <= (magnet_t.get_pos(1)+(magnet_t.get_width())));
        //^ WIDTH TAG
        bool within_z_bounds = (pos_at_shortest_time.getZ() >= (magnet_t.get_pos(2)-(magnet_t.get_width()))) && (pos_at_shortest_time.getZ() <= (magnet_t.get_pos(2)+(magnet_t.get_width())));

        if(within_x_bounds && within_y_bounds && within_z_bounds)
        { intersect = true; }
        else
        { intersect = false; }
    
        return intersect;
    }
    else
    {
        std::cout << "No valid magnet type detected in intersect mag function.\n";
        return false;
    }
}


bool intersect_mag_dipole(Magnet magnet_t, Particle particle_t)
{
    /* Will check if particle intersects with a magnet
     * Will NOT update particle's positions! We'll want to check if there's another
     * magnet that is closer than the one currently being checked before updating position.
     */
    bool intersect;

    double time_to_magnet = time_to_magnet_boundary_dipole(magnet_t, particle_t);
    if(time_to_magnet < 0.0)
    {
        return false;
    }

    ThreeVec pos_at_shortest_time;
    pos_at_shortest_time.setX(particle_t.get_pos(0) + (time_to_magnet * particle_t.get_vel(0)));
    pos_at_shortest_time.setY(particle_t.get_pos(1) + (time_to_magnet * particle_t.get_vel(1)));
    pos_at_shortest_time.setZ(particle_t.get_pos(2) + (time_to_magnet * particle_t.get_vel(2)));

    bool within_x_bounds = (pos_at_shortest_time.getX() >= magnet_t.get_pos(0)-1.0*magnet_t.get_length()) && (pos_at_shortest_time.getX() <= (magnet_t.get_pos(0)+2.0*magnet_t.get_length()));
    //^ LENGTH TAG: edited to make width region 5*width away from center instead of 0.5*width.
    bool within_y_bounds = (pos_at_shortest_time.getY() >= (magnet_t.get_pos(1)-(magnet_t.get_width()*5.0)))  && (pos_at_shortest_time.getY() <= (magnet_t.get_pos(1)+(magnet_t.get_width()*5.0)));
    //^ WIDTH TAG
    bool within_z_bounds = (pos_at_shortest_time.getZ() >= (magnet_t.get_pos(2)-(magnet_t.get_height()/2.0))) && (pos_at_shortest_time.getZ() <= (magnet_t.get_pos(2)+(magnet_t.get_height()/2.0)));

    if(within_x_bounds && within_y_bounds && within_z_bounds)
    { intersect = true; }
    else
    { intersect = false; }
 
    return intersect;
}

bool intersect_mag_uniform(Magnet magnet_t, Particle particle_t)
{
    /* Will check if particle intersects with a magnet
     * Will NOT update particle's positions! We'll want to check if there's another
     * magnet that is closer than the one currently being checked before updating position.
     */
    bool intersect;

    double time_to_magnet = time_to_magnet_boundary_uniform(magnet_t, particle_t);
    if(time_to_magnet < 0.0)
    {
        return false;
    }

    ThreeVec pos_at_shortest_time;
    pos_at_shortest_time.setX(particle_t.get_pos(0) + (time_to_magnet * particle_t.get_vel(0)));
    pos_at_shortest_time.setY(particle_t.get_pos(1) + (time_to_magnet * particle_t.get_vel(1)));
    pos_at_shortest_time.setZ(particle_t.get_pos(2) + (time_to_magnet * particle_t.get_vel(2)));

    bool within_x_bounds = (pos_at_shortest_time.getX() >= magnet_t.get_pos(0)) && (pos_at_shortest_time.getX() <= (magnet_t.get_pos(0)+magnet_t.get_length()));
    //^ LENGTH TAG: edited to make width region 5*width away from center instead of 0.5*width.
    bool within_y_bounds = (pos_at_shortest_time.getY() >= (magnet_t.get_pos(1)-(magnet_t.get_width()*0.5)))  && (pos_at_shortest_time.getY() <= (magnet_t.get_pos(1)+(magnet_t.get_width()*0.5)));
    //^ WIDTH TAG
    bool within_z_bounds = (pos_at_shortest_time.getZ() >= (magnet_t.get_pos(2)-(magnet_t.get_height()/2.0))) && (pos_at_shortest_time.getZ() <= (magnet_t.get_pos(2)+(magnet_t.get_height()/2.0)));

    if(within_x_bounds && within_y_bounds && within_z_bounds)
    { intersect = true; }
    else
    { intersect = false; }
 
    return intersect;
}


double dist_to_mag_general(Magnet magnet_t, Particle particle_t)
{
    char mag_type = magnet_t.get_type();
    if(mag_type == 'd')
    {
        double dist_to_mag;
        double time_to_mag = time_to_magnet_boundary_general(magnet_t, particle_t);
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
    else if(mag_type == 'u')
    {
        double dist_to_mag;
        double time_to_mag = time_to_magnet_boundary_general(magnet_t, particle_t);
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
    else if(mag_type == 'q' || mag_type =='h')
    {
        double dist_to_mag;
        double time_to_mag = time_to_magnet_boundary_general(magnet_t, particle_t);
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
    else
    {
        std::cout << "No valid magnet type detected in dist_to_mag function.\n";
        return DBL_MAX;
    }
}


double dist_to_mag_dipole(Magnet magnet_t, Particle particle_t)
{
    double dist_to_mag;
    double time_to_mag = time_to_magnet_boundary_dipole(magnet_t, particle_t);

    if(time_to_mag <= 0.0)
    {
        return 0.0;
    }

    ThreeVec pos_at_shortest_time;
    pos_at_shortest_time.setX(particle_t.get_pos(0) + (time_to_mag * particle_t.get_vel(0)));
    pos_at_shortest_time.setY(particle_t.get_pos(1) + (time_to_mag * particle_t.get_vel(1)));
    pos_at_shortest_time.setZ(particle_t.get_pos(2) + (time_to_mag * particle_t.get_vel(2)));

    dist_to_mag = sqrt((pos_at_shortest_time.getX() - particle_t.get_pos(0))*(pos_at_shortest_time.getX() - particle_t.get_pos(0)) + (pos_at_shortest_time.getY() - particle_t.get_pos(1))*(pos_at_shortest_time.getY() - particle_t.get_pos(1)) + (pos_at_shortest_time.getZ() - particle_t.get_pos(2))*(pos_at_shortest_time.getZ() - particle_t.get_pos(2)));
    //std::cerr<< dist_to_mag << '\n';
    return dist_to_mag;
}

double dist_to_mag_uniform(Magnet magnet_t, Particle particle_t)
{
    double dist_to_mag;
    double time_to_mag = time_to_magnet_boundary_uniform(magnet_t, particle_t);

    if(time_to_mag <= 0.0)
    {
        return 0.0;
    }

    ThreeVec pos_at_shortest_time;
    pos_at_shortest_time.setX(particle_t.get_pos(0) + (time_to_mag * particle_t.get_vel(0)));
    pos_at_shortest_time.setY(particle_t.get_pos(1) + (time_to_mag * particle_t.get_vel(1)));
    pos_at_shortest_time.setZ(particle_t.get_pos(2) + (time_to_mag * particle_t.get_vel(2)));

    dist_to_mag = sqrt((pos_at_shortest_time.getX() - particle_t.get_pos(0))*(pos_at_shortest_time.getX() - particle_t.get_pos(0)) + (pos_at_shortest_time.getY() - particle_t.get_pos(1))*(pos_at_shortest_time.getY() - particle_t.get_pos(1)) + (pos_at_shortest_time.getZ() - particle_t.get_pos(2))*(pos_at_shortest_time.getZ() - particle_t.get_pos(2)));
    //std::cerr<< dist_to_mag << '\n';
    return dist_to_mag;
}





void move_particle_to_magnet_general(Magnet magnet_t, Particle &particle_t)
{
    double time_to_mag = time_to_magnet_boundary_general(magnet_t, particle_t);

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

void move_particle_to_magnet_dipole(Magnet magnet_t, Particle &particle_t)
{
    double time_to_mag = time_to_magnet_boundary_dipole(magnet_t, particle_t);

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

void move_particle_to_magnet_uniform(Magnet magnet_t, Particle &particle_t)
{
    double time_to_mag = time_to_magnet_boundary_uniform(magnet_t, particle_t);

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


bool is_this_array_only_zeros(bool array[], int num_elements)
{
    //returns true if only zeros, returns false if any non-zero elements
    if(num_elements >= 1)
    {
        for(int ii = 0; ii< num_elements; ii++)
        {
            if(array[ii] != 0)
            {
                return false;
            }
        }
    }
    else
    {
        std::cout << "Error in 'is_this_array_only_zeros' function - number of elements is not at least 1.\n";
        return true;
    }
    return true;
}


double dist_to_screen_ii(Screen screen_t, Particle particle_t);
bool check_if_intersect_screen(Screen screen_t, Particle particle_t);





/////UNDER CONSTRUCTION
/////
/////
bool move_through_magnets_general(Magnet magnet_t[], int num_mags, Particle &particle_t, double &time, double del_time, double mu_0, double time_limit, Screen screen_t[], int num_screens)
{
    //if particle lands on screen in this process, return true. else, return false.
    
    bool check_inside_any_magnet = false;
    bool check_intersect_magnet = false;
    double distance_to_mag_ii[num_mags];
    double distance_to_screen_ii[num_screens];
    int boris_counter = 0;
    bool inside_magnet_array[num_mags];
    bool went_through_magnet_without_intersecting_screen = true;
    double shortest_magnet_distance, shortest_screen_distance;

    for(int ii=0; ii<num_mags; ii++)
    {
        inside_magnet_array[ii] = inside_of_mag_general(magnet_t[ii], particle_t);
    }
    check_inside_any_magnet = !(is_this_array_only_zeros(inside_magnet_array, num_mags));
    for(int ii=0; ii<num_mags; ii++)
    {
        check_intersect_magnet = intersect_mag_general(magnet_t[ii], particle_t);
        if(check_intersect_magnet == true)
        { break; }
    }
    //for(int ii=0; ii<num_mags; ii++)
    //{
    while(check_inside_any_magnet || check_intersect_magnet)
    {
        if(check_inside_any_magnet==true)
        {
            went_through_magnet_without_intersecting_screen = step_through_magnet_mag_boris_analytic_general(particle_t, magnet_t, time, del_time, mu_0, time_limit, inside_magnet_array, num_mags, screen_t, num_screens, true);
            boris_counter++;
            //ii = -1; //restart loop (which will iterate ii by 1, to zero)
        }
        if(went_through_magnet_without_intersecting_screen == false)
        { return true; }

        int index_of_shortest_magnet= -1;
        for(int ii = 0; ii < num_mags; ii++)
        {
            check_intersect_magnet = intersect_mag_general(magnet_t[ii], particle_t);
            if(check_intersect_magnet==true)
            {
                distance_to_mag_ii[ii] = dist_to_mag_general(magnet_t[ii], particle_t);
            }
            else
            {
                distance_to_mag_ii[ii] = 0.0;
            } 

            if(ii == (num_mags - 1))
            {
                shortest_magnet_distance = DBL_MAX;
                for(int jj = 0; jj<num_mags; jj++)
                {
                    if(distance_to_mag_ii[jj] != 0.0 && distance_to_mag_ii[jj] < shortest_magnet_distance)
                    {
                        //ii = -1;
                        shortest_magnet_distance = distance_to_mag_ii[jj];
                        index_of_shortest_magnet = jj;
                    }
                }
            }
        }

        bool check_intersect_screen = false;
        int index_of_shortest_screen = -1;
        for(int ii = 0; ii<num_screens; ii++)
        {
            check_intersect_screen = check_if_intersect_screen(screen_t[ii], particle_t);
            if(check_intersect_screen == true)
            {
                distance_to_screen_ii[ii] = dist_to_screen_ii(screen_t[ii], particle_t);
            }
            else
            {
                distance_to_screen_ii[ii] = 0.0;
            }

            if(ii == (num_screens - 1))
            {
                shortest_screen_distance = DBL_MAX;
                for(int jj = 0; jj<num_screens; jj++)
                {
                    if(distance_to_screen_ii[jj] != 0.0 && distance_to_screen_ii[jj] < shortest_screen_distance)
                    {
                        //ii = -1;
                        shortest_screen_distance = distance_to_screen_ii[jj];
                        index_of_shortest_screen = jj;
                    }
                }
            }
        }
            if(index_of_shortest_magnet != -1 && index_of_shortest_screen == -1)
            {
                move_particle_to_magnet_general(magnet_t[index_of_shortest_magnet], particle_t);
                outfile_part_commaAndWrite(particle_t);

                double particle_time = particle_t.get_time();
                time = time + particle_time;
                //continue;
            }
            else if(index_of_shortest_screen != -1 && index_of_shortest_magnet ==-1)
            {
                return true;
            }
            else if(index_of_shortest_magnet != -1 && index_of_shortest_screen != -1 && (shortest_screen_distance < shortest_magnet_distance))
            {
                return true;
            }
            else if(index_of_shortest_magnet != -1 && index_of_shortest_screen != -1 && (shortest_magnet_distance < shortest_screen_distance))
            {
                move_particle_to_magnet_general(magnet_t[index_of_shortest_magnet], particle_t);
                outfile_part_commaAndWrite(particle_t);

                double particle_time = particle_t.get_time();
                time = time + particle_time;
                //continue;
            }
        check_inside_any_magnet = false;
        check_intersect_magnet = false;
        
        for(int ii=0; ii<num_mags; ii++)
        {
            inside_magnet_array[ii] = inside_of_mag_general(magnet_t[ii], particle_t);
        }
        check_inside_any_magnet = !(is_this_array_only_zeros(inside_magnet_array, num_mags));
        for(int ii=0; ii<num_mags; ii++)
        {
            check_intersect_magnet = intersect_mag_general(magnet_t[ii], particle_t);
            if(check_intersect_magnet == true)
            { break; }
        }
    }
    //}
    return false;
}
/////
/////
/////UNDER CONSTRUCTION







bool move_through_magnets_dipole(Magnet magnet_t[], int num_mags, Particle &particle_t, double &time, double del_time, double mu_0, double time_limit, Screen screen_t[], int num_screens)
{
    //if particle lands on screen in this process, return true. else, return false.
    
    bool check_inside_any_magnet = false;
    bool check_intersect_magnet = false;
    double distance_to_mag_ii[num_mags];
    double distance_to_screen_ii[num_screens];
    int boris_counter = 0;
    bool inside_magnet_array[num_mags];
    bool went_through_magnet_without_intersecting_screen = true;
    double shortest_magnet_distance, shortest_screen_distance;

    for(int ii=0; ii<num_mags; ii++)
    {
        inside_magnet_array[ii] = inside_of_mag_dipole(magnet_t[ii], particle_t);
    }
    check_inside_any_magnet = !(is_this_array_only_zeros(inside_magnet_array, num_mags));
    for(int ii=0; ii<num_mags; ii++)
    {
        check_intersect_magnet = intersect_mag_dipole(magnet_t[ii], particle_t);
        if(check_intersect_magnet == true)
        { break; }
    }
    //for(int ii=0; ii<num_mags; ii++)
    //{
    while(check_inside_any_magnet || check_intersect_magnet)
    {
        if(check_inside_any_magnet==true)
        {
            went_through_magnet_without_intersecting_screen = step_through_magnet_mag_boris_analytic_dipole(particle_t, magnet_t, time, del_time, mu_0, time_limit, inside_magnet_array, num_mags, screen_t, num_screens, true);
            boris_counter++;
            //ii = -1; //restart loop (which will iterate ii by 1, to zero)
        }
        if(went_through_magnet_without_intersecting_screen == false)
        { return true; }

        int index_of_shortest_magnet= -1;
        for(int ii = 0; ii < num_mags; ii++)
        {
            check_intersect_magnet = intersect_mag_dipole(magnet_t[ii], particle_t);
            if(check_intersect_magnet==true)
            {
                distance_to_mag_ii[ii] = dist_to_mag_dipole(magnet_t[ii], particle_t);
            }
            else
            {
                distance_to_mag_ii[ii] = 0.0;
            } 

            if(ii == (num_mags - 1))
            {
                shortest_magnet_distance = DBL_MAX;
                for(int jj = 0; jj<num_mags; jj++)
                {
                    if(distance_to_mag_ii[jj] != 0.0 && distance_to_mag_ii[jj] < shortest_magnet_distance)
                    {
                        //ii = -1;
                        shortest_magnet_distance = distance_to_mag_ii[jj];
                        index_of_shortest_magnet = jj;
                    }
                }
            }
        }

        bool check_intersect_screen = false;
        int index_of_shortest_screen = -1;
        for(int ii = 0; ii<num_screens; ii++)
        {
            check_intersect_screen = check_if_intersect_screen(screen_t[ii], particle_t);
            if(check_intersect_screen == true)
            {
                distance_to_screen_ii[ii] = dist_to_screen_ii(screen_t[ii], particle_t);
            }
            else
            {
                distance_to_screen_ii[ii] = 0.0;
            }

            if(ii == (num_screens - 1))
            {
                shortest_screen_distance = DBL_MAX;
                for(int jj = 0; jj<num_screens; jj++)
                {
                    if(distance_to_screen_ii[jj] != 0.0 && distance_to_screen_ii[jj] < shortest_screen_distance)
                    {
                        //ii = -1;
                        shortest_screen_distance = distance_to_screen_ii[jj];
                        index_of_shortest_screen = jj;
                    }
                }
            }
        }
            if(index_of_shortest_magnet != -1 && index_of_shortest_screen == -1)
            {
                move_particle_to_magnet_dipole(magnet_t[index_of_shortest_magnet], particle_t);
                outfile_part_commaAndWrite(particle_t);

                double particle_time = particle_t.get_time();
                time = time + particle_time;
                //continue;
            }
            else if(index_of_shortest_screen != -1 && index_of_shortest_magnet ==-1)
            {
                return true;
            }
            else if(index_of_shortest_magnet != -1 && index_of_shortest_screen != -1 && (shortest_screen_distance < shortest_magnet_distance))
            {
                return true;
            }
            else if(index_of_shortest_magnet != -1 && index_of_shortest_screen != -1 && (shortest_magnet_distance < shortest_screen_distance))
            {
                move_particle_to_magnet_dipole(magnet_t[index_of_shortest_magnet], particle_t);
                outfile_part_commaAndWrite(particle_t);

                double particle_time = particle_t.get_time();
                time = time + particle_time;
                //continue;
            }
        check_inside_any_magnet = false;
        check_intersect_magnet = false;
        
        for(int ii=0; ii<num_mags; ii++)
        {
            inside_magnet_array[ii] = inside_of_mag_dipole(magnet_t[ii], particle_t);
        }
        check_inside_any_magnet = !(is_this_array_only_zeros(inside_magnet_array, num_mags));
        for(int ii=0; ii<num_mags; ii++)
        {
            check_intersect_magnet = intersect_mag_dipole(magnet_t[ii], particle_t);
            if(check_intersect_magnet == true)
            { break; }
        }
    }
    //}
    return false;
}

bool move_through_magnets_uniform(Magnet magnet_t[], int num_mags, Particle &particle_t, double &time, double del_time, double mu_0, double time_limit, Screen screen_t[], int num_screens)
{
    //std::cerr << "Inside of move through magnets uniform function\n";
    //if particle lands on screen in this process, return true. else, return false.
    
    bool check_inside_any_magnet = false;
    bool check_intersect_magnet = false;
    double distance_to_mag_ii[num_mags];
    double distance_to_screen_ii[num_screens];
    int boris_counter = 0;
    bool inside_magnet_array[num_mags];
    bool went_through_magnet_without_intersecting_screen = true;
    double shortest_magnet_distance, shortest_screen_distance;

    for(int ii=0; ii<num_mags; ii++)
    {
        inside_magnet_array[ii] = inside_of_mag_uniform(magnet_t[ii], particle_t);
    }
    check_inside_any_magnet = !(is_this_array_only_zeros(inside_magnet_array, num_mags));
    for(int ii=0; ii<num_mags; ii++)
    {
        check_intersect_magnet = intersect_mag_uniform(magnet_t[ii], particle_t);
        //std::cout << check_intersect_magnet << "\n";
        if(check_intersect_magnet == true)
        { break; }
    }
    //for(int ii=0; ii<num_mags; ii++)
    //{
        //std::cerr << "Instersect check is " << check_intersect_magnet << " and inside check is " << check_inside_any_magnet << "\n";
    while(check_inside_any_magnet || check_intersect_magnet)
    {
        if(check_inside_any_magnet==true)
        {
            //std::cout << "Inside of uniform magnet\n";
            went_through_magnet_without_intersecting_screen = step_through_magnet_mag_boris_analytic_uniform_field(particle_t, magnet_t, time, del_time, mu_0, time_limit, inside_magnet_array, num_mags, screen_t, num_screens, false);
            boris_counter++;
            //ii = -1; //restart loop (which will iterate ii by 1, to zero)
        }
        if(went_through_magnet_without_intersecting_screen == false)
        { return true; }

        int index_of_shortest_magnet= -1;
        for(int ii = 0; ii < num_mags; ii++)
        {
            check_intersect_magnet = intersect_mag_uniform(magnet_t[ii], particle_t);
            if(check_intersect_magnet==true)
            {
                distance_to_mag_ii[ii] = dist_to_mag_uniform(magnet_t[ii], particle_t);
            }
            else
            {
                distance_to_mag_ii[ii] = 0.0;
            } 

            if(ii == (num_mags - 1))
            {
                shortest_magnet_distance = DBL_MAX;
                for(int jj = 0; jj<num_mags; jj++)
                {
                    if(distance_to_mag_ii[jj] != 0.0 && distance_to_mag_ii[jj] < shortest_magnet_distance)
                    {
                        //ii = -1;
                        shortest_magnet_distance = distance_to_mag_ii[jj];
                        index_of_shortest_magnet = jj;
                    }
                }
            }
        }

        bool check_intersect_screen = false;
        int index_of_shortest_screen = -1;
        for(int ii = 0; ii<num_screens; ii++)
        {
            check_intersect_screen = check_if_intersect_screen(screen_t[ii], particle_t);
            if(check_intersect_screen == true)
            {
                distance_to_screen_ii[ii] = dist_to_screen_ii(screen_t[ii], particle_t);
            }
            else
            {
                distance_to_screen_ii[ii] = 0.0;
            }

            if(ii == (num_screens - 1))
            {
                shortest_screen_distance = DBL_MAX;
                for(int jj = 0; jj<num_screens; jj++)
                {
                    if(distance_to_screen_ii[jj] != 0.0 && distance_to_screen_ii[jj] < shortest_screen_distance)
                    {
                        //ii = -1;
                        shortest_screen_distance = distance_to_screen_ii[jj];
                        index_of_shortest_screen = jj;
                    }
                }
            }
        }
            if(index_of_shortest_magnet != -1 && index_of_shortest_screen == -1)
            {
                move_particle_to_magnet_uniform(magnet_t[index_of_shortest_magnet], particle_t);
                outfile_part_commaAndWrite(particle_t);

                double particle_time = particle_t.get_time();
                time = time + particle_time;
                //continue;
            }
            else if(index_of_shortest_screen != -1 && index_of_shortest_magnet ==-1)
            {
                return true;
            }
            else if(index_of_shortest_magnet != -1 && index_of_shortest_screen != -1 && (shortest_screen_distance < shortest_magnet_distance))
            {
                return true;
            }
            else if(index_of_shortest_magnet != -1 && index_of_shortest_screen != -1 && (shortest_magnet_distance < shortest_screen_distance))
            {
                //std::cout << "Moving to next uniform magnet\n";
                move_particle_to_magnet_uniform(magnet_t[index_of_shortest_magnet], particle_t);
                outfile_part_commaAndWrite(particle_t);

                double particle_time = particle_t.get_time();
                time = time + particle_time;
                //continue;
            }
        check_inside_any_magnet = false;
        check_intersect_magnet = false;
        
        for(int ii=0; ii<num_mags; ii++)
        {
            inside_magnet_array[ii] = inside_of_mag_uniform(magnet_t[ii], particle_t);
            //std::cout << "End-loop inside magnet check for magnet " << ii << " is " << inside_magnet_array[ii] << "\n";
        }
        check_inside_any_magnet = !(is_this_array_only_zeros(inside_magnet_array, num_mags));
        for(int ii=0; ii<num_mags; ii++)
        {
            check_intersect_magnet = intersect_mag_uniform(magnet_t[ii], particle_t);
            if(check_intersect_magnet == true)
            { break; }
        }
    }
    //}
    return false;
}

bool move_through_magnets_OLD(Magnet magnet_t[], int num_mags, Particle &particle_t, double &time, double del_time,double mu_0, double time_limit)
{
    //if particle lands on screen in this process, return true. else, return false.
    bool check_inside_magnet = false;
    bool check_intersect_magnet = false;
    double distance_to_mag_ii[num_mags];
    int boris_counter = 0;
    for(int ii=0; ii<num_mags; ii++)
    {
        check_inside_magnet = inside_of_mag_dipole(magnet_t[ii], particle_t);
        if(check_inside_magnet==true)
        {
            //std::cerr << "Inside of magnet " << ii+1 << '\n';
            //outfile_part_comma(particle_t);
            //step_through_magnet_mag_boris_analytic(particle_t, magnet_t[ii], time, del_time, time_limit);
            step_through_magnet_mag_boris_analytic_ONEMAGNET(particle_t, magnet_t, time, del_time, time_limit, ii, num_mags, false);
            boris_counter++;
            ii = -1; //restart loop (which will iterate ii by 1, to zero)
            continue;
        }

        check_intersect_magnet = intersect_mag_dipole(magnet_t[ii], particle_t);
        if(check_intersect_magnet==true)
        {
            distance_to_mag_ii[ii] = dist_to_mag_dipole(magnet_t[ii], particle_t);
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
                //std::cerr << "current particle X-location is " << particle_t.get_pos(0) << std::endl;
                //std::cerr << "distance to next magnet is " << shortest_distance << std::endl;
                move_particle_to_magnet_dipole(magnet_t[index_of_shortest], particle_t);
                //std::cerr << "Sending to magnet " << index_of_shortest+1 << std::endl;
                //std::cerr << "Sending to next magnet" << std::endl;
                //std::cerr << "New particle X-location is " << particle_t.get_pos(0) << std::endl;
                outfile_part_commaAndWrite(particle_t);

                double particle_time = particle_t.get_time();
                time = time + particle_time;
                continue;
            }
        }
    }
    return false;
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
    if((area0 + area0*0.00001) >= (area12 + area13 + area24 + area34))
    //if((area0) >= (area12 + area13 + area24 + area34))
    {
        intersect = true;
    }

    return intersect;
}

bool check_if_intersect_screen_in_next_dt(Screen screen_t, Particle particle_t, double del_t)
{
    bool intersect = false;
    double t_line  = t_line_particle(screen_t, particle_t);

    if(t_line < 0.0 || t_line > del_t)
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
    //if(area0 >= (area12 + area13 + area24 + area34))
    if((area0 + area0*0.00001) >= (area12 + area13 + area24 + area34))
    {
        intersect = true;
    }

    return intersect;
}

double dist_to_screen_ii(Screen screen_t, Particle particle_t)
{
    double dist_to_screen_ii;
    int jj = -1;
    double t_numerator, t_denom, t_line;
    t_line      = t_line_particle(screen_t, particle_t);

    if(t_line < 0.0)
    {return -1.0;}

    double x_intersect = particle_t.get_pos(0) + particle_t.get_vel(0)*t_line;
    double y_intersect = particle_t.get_pos(1) + particle_t.get_vel(1)*t_line;
    double z_intersect = particle_t.get_pos(2) + particle_t.get_vel(2)*t_line;

    dist_to_screen_ii = sqrt((x_intersect - particle_t.get_pos(0))*(x_intersect - particle_t.get_pos(0)) + (y_intersect - particle_t.get_pos(1))*(y_intersect - particle_t.get_pos(1)) + (z_intersect - particle_t.get_pos(2))*(z_intersect - particle_t.get_pos(2)));

    return dist_to_screen_ii;

}


void move_to_screens(Screen screen_t[], int num_screen, Particle particle_t, int particle_counter_t)
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
                //std::cout<< "Particle " << particle_counter_t <<" read as hitting screen " << screen_t[index_of_shortest].get_index() << "\n";
                outfile_part_on_screen(screen_t[index_of_shortest], particle_counter_t, particle_t);
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

void readMagnetType(std::ifstream &input_stream, int magNum, std::vector<char> &PmagType) {
    //std::string tempStr;
    char tempStr;
    for(int i=0; i<magNum; ++i) {
        input_stream >> tempStr;
        //double tempDbl = std::stod(tempStr);
        //std::cerr << "Magnet dim read\n";
        PmagType.push_back(tempStr);
    }
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
                    //std::cerr << "Magnet partially read\n";
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
                //std::cerr << "Magnet partially read\n";
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
        //std::cerr << "Magnet dim read\n";
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
            //std::cerr << "Beam partially read\n";
            tempInfoBits.push_back(tempDbl);
        }
        else {
            for(int j=0; j<3; ++j) {
                input_stream >> tempStr;
                double tempDbl = std::stod(tempStr);
                //std::cerr << "Beam partially read\n";
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
            //std::cerr << "beam spread partially read\n";
            tempInfoBits.push_back(tempDbl);
        }
        else {
            for(int k=0; k<3; ++k) {
                input_stream >> tempStr;
                double tempDbl = std::stod(tempStr);
                //std::cerr << "beam spread partially read\n";

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
                    //std::cerr << "screen partially read\n";

                    tempInfoBits.push_back(tempDbl);
                }
            }
            else {
                for(int k=0; k<3; ++k) {
                    input_stream >> tempStr;
                    double tempDbl = std::stod(tempStr);
                    //std::cerr << "screen partially read\n";
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

double ReadMu0(std::ifstream &input_stream) {
    std::string tempStr;
    input_stream >> tempStr;
    double mu_0 = std::stod(tempStr);
    //std::cerr << "mu0 read\n";

    return mu_0;
}

double ReadSpeciesCharge(std::ifstream &input_stream) {
    std::string tempStr;
    input_stream >> tempStr;
    double charge = std::stod(tempStr);

    return charge;
}

char ReadDipoleMagnetFieldType(std::ifstream &input_stream) {
    //std::string tempStr;
    char magtype;
    input_stream >> magtype;
    //char magtype = std::stod(tempStr);
    return magtype;
}

double F2(double a, double b, double c, double x, double y, double z)
{
    double numerator = sqrt( ((x+a)*(x+a)) + ((y-b)*(y-b)) + ((z+c)*(z+c)) ) + b - y;
    double denomenator = (sqrt( ((x+a)*(x+a)) + ((y+b)*(y+b)) + ((z+c)*(z+c)) ) - b - y);
    return (numerator/denomenator);
}

double F1(double a, double b, double c, double x, double y, double z)
{
    return atan( ( (x+a)*(y+b) )/( (z+c)*sqrt( ((x+a)*(x+a)) + ((y+b)*(y+b)) + ((z+c)*(z+c)) ) ) );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double find_magnetization(Magnet &magnet, double mag_dim, double mu_0, char axis) {
    //mag_dim is the 'height' of the dipole's magnetic blocks
    double B_center, a, b, c, d;
    if(axis=='z') {
        B_center = magnet.get_B0(2);
        //std::cerr << "Central Magnetic Field is " << B_center << "\n";

        // dimensions of permanent magnet
        a = 0.5*magnet.get_length();  // half length for permanent magnet
        b = 0.5*magnet.get_width();   // half width for permanent magnet
        c = 0.5*mag_dim;           // half height for permanent magnet
        
        d = 0.5*magnet.get_height() + c;  // distance from mag center to center between mags along the axis of magnetization
    }
    else if(axis=='y') {
        B_center = magnet.get_B0(1);

        a = 0.5*magnet.get_length();
        b = 0.5*magnet.get_height();
        c = 0.5*mag_dim;

        d = 0.5*magnet.get_width() + c;
    }
    else {
        // axis = x
        B_center = magnet.get_B0(0);

        a = 0.5*magnet.get_height();
        b = 0.5*magnet.get_width();
        c = 0.5*mag_dim;

        d = 0.5*magnet.get_length() + c;
    }
    //double magnetization = B_center/((mu_0/M_PI)*(atan((a*b)/((d-c)*sqrt(pow(a,2)+pow(b,2)+pow(d-c,2))))-atan((a*b)/((d+c)*sqrt(pow(a,2)+pow(b,2)+pow(d+c,2)))))*2);
    double magnetization = -0.5*B_center/((mu_0/M_PI)*(F1(a,b,c,0,0,d) + F1(a,b,c,0,0,-d)));

    // divided by 2 because of the contributions of both magnets
    
    return magnetization;
}

void calc_grid_B_comps(double factor, double a, double b, double c, double x, double y, double z, double &temp_B1, double &temp_B2, double &temp_B3, bool iterative_add) 
{
    if(iterative_add == true)
    {
        temp_B1 = temp_B1 + factor*log( ( F2(a,b,c,-x,y,-z)*F2(a,b,c,x,y,z) )/( F2(a,b,c,x,y,-z)*F2(a,b,c,-x,y,z) )  );

        temp_B2 = temp_B2 + factor*log( ( F2(b,a,c,-y,x,-z)*F2(b,a,c,y,x,z) )/( F2(b,a,c,y,x,-z)*F2(b,a,c,-y,x,z) )  );

        temp_B3 = temp_B3 - factor*( F1(a,b,c,-x,y,z) + F1(a,b,c,-x,y,-z) + F1(a,b,c,-x,-y,z) + F1(a,b,c,-x,-y,-z) + F1(a,b,c,x,y,z) + F1(a,b,c,x,y,-z) + F1(a,b,c,x,-y,z) + F1(a,b,c,x,-y,-z) );
    }
    else
    {
        temp_B1 = factor*log( ( F2(a,b,c,-x,y,-z)*F2(a,b,c,x,y,z) )/( F2(a,b,c,x,y,-z)*F2(a,b,c,-x,y,z) )  );

        temp_B2 = factor*log( ( F2(b,a,c,-y,x,-z)*F2(b,a,c,y,x,z) )/( F2(b,a,c,y,x,-z)*F2(b,a,c,-y,x,z) )  );

        temp_B3 = -factor*( F1(a,b,c,-x,y,z) + F1(a,b,c,-x,y,-z) + F1(a,b,c,-x,-y,z) + F1(a,b,c,-x,-y,-z) + F1(a,b,c,x,y,z) + F1(a,b,c,x,y,-z) + F1(a,b,c,x,-y,z) + F1(a,b,c,x,-y,-z) );
    }
}

bool B_within_margin(double B_center_val, double B1, double B2, double B3) 
{
    bool in_margin = true;
    
    // magnitude of B field from 1 magnet
    double magnitude = sqrt( (B1*B1) + (B2*B2) + (B3*B3) );
    double percent = 0.01;
    double cutoff_value = fabs(percent * B_center_val);

    if(magnitude < cutoff_value) {
        in_margin = false;
    }
    
    return in_margin;
}

ThreeVec calc_dipole_B(ThreeVec &grid_point, Magnet &magnet, double mag_dim, double mu_0, char axis, double magnetization) {
    ThreeVec grid_point_B;
    double grid_B1 = 0.0;  //B1 and B2 are components perpendicular to magnetization axis
    double grid_B2 = 0.0;  // for z-axis, B1 = x comp and B2 = y comp
    double grid_B3 = 0.0;  // B3 is along axis

    double factor = (mu_0 * magnetization)/(4 * M_PI);
    //std::cerr << "magnetization is " << magnetization << "\n";

    if(axis=='z') {

        // dimensions of permanent magnet
        double a = 0.5*magnet.get_length();  // half length for permanent magnet
        double b = 0.5*magnet.get_width();   // half width for permanent magnet
        double c = 0.5*mag_dim;           // half height for permanent magnet
        
        double offset = 0.5*magnet.get_height() + c;  // moves user defined point for magnet space to the magnet center
        
        double B_center_val = magnet.get_B0(2);
        
        double hold_B1 = 0.0;
        double hold_B2 = 0.0;
        double hold_B3 = 0.0;
        
        int i = 1;
        while(i > -2) {
            double mag_center_x = magnet.get_pos(0) + a;
            double mag_center_y = magnet.get_pos(1);
            double mag_center_z = magnet.get_pos(2) - (i*offset);
            
            double x = grid_point.getX() - mag_center_x;
            double y = grid_point.getY() - mag_center_y;
            double z = grid_point.getZ() - mag_center_z;  // distance from magnet to grid point
            
            double temp_B1 = 0.0;
            double temp_B2 = 0.0;
            double temp_B3 = 0.0;
            calc_grid_B_comps(factor, a, b, c, x, y, z, temp_B1, temp_B2, temp_B3);
            
            hold_B1 += temp_B1;
            hold_B2 += temp_B2;
            hold_B3 += temp_B3;
            //std::cerr << "B1 calculation is " << temp_B1 << "\n";
            //std::cerr << "B3 calculation is " << temp_B3 << "\n";
            //std::cerr << "hold B3 is " << hold_B3 << "\n";
            
            i -= 2;
            //if(i == -3)
            //{
            //    std::cerr << "Field at magnetic position (" << x << ", " << y <<", " << z << ") is (" << hold_B1 << ", " << hold_B2 << ", " << hold_B3 << ")\n";
            //}
        }
        //bool in_margin = B_within_margin(B_center_val, hold_B1, hold_B2, hold_B3);
        bool in_margin = true;
        if(in_margin) {
            grid_B1 = hold_B1;
            grid_B2 = hold_B2;
            grid_B3 = hold_B3;
        }

        grid_point_B.setX(grid_B1);
        grid_point_B.setY(grid_B2);
        grid_point_B.setZ(grid_B3);
    }
    else if(axis=='y') {

        double a = 0.5*magnet.get_length();
        double b = 0.5*magnet.get_height();
        double c = 0.5*mag_dim;
        
        double offset = 0.5*magnet.get_width() + c;
        
        double B_center_val = magnet.get_B0(1);
        
        double hold_B1 = 0.0;
        double hold_B2 = 0.0;
        double hold_B3 = 0.0;
        
        int i = 1;
        while (i > -2) {
            double mag_center_x = magnet.get_pos(0) + a;
            double mag_center_y = magnet.get_pos(1) - (i*offset);
            double mag_center_z = magnet.get_pos(2);
            
            double x = grid_point.getX() - mag_center_x;
            double y = grid_point.getY() - mag_center_y;
            double z = grid_point.getZ() - mag_center_z;
            
            double temp_B1 = 0.0;
            double temp_B2 = 0.0;
            double temp_B3 = 0.0;
            calc_grid_B_comps(factor, a, b, c, x, z, y, temp_B1, temp_B2, temp_B3);

            hold_B1 += temp_B1;
            hold_B2 += temp_B2;
            hold_B3 += temp_B3;
            
            i -= 2;
        }
        //bool in_margin = B_within_margin(B_center_val, hold_B1, hold_B2, hold_B3);
        bool in_margin = true;
        if(in_margin) {
            grid_B1 = hold_B1;
            grid_B2 = hold_B2;
            grid_B3 = hold_B3;
        }

        grid_point_B.setX(grid_B1);
        grid_point_B.setY(grid_B3);
        grid_point_B.setZ(grid_B2);
    }
    else 
    {
        // axis = x

        double a = 0.5*magnet.get_height();
        double b = 0.5*magnet.get_width();
        double c = 0.5*mag_dim;
        
        double B_center_val = magnet.get_B0(0);
        
        double hold_B1 = 0.0;
        double hold_B2 = 0.0;
        double hold_B3 = 0.0;
        
        int i = 1;
        while (i > -2) {
            // user defined point is not at halfway of the magnet space length
            double offset;
            if(i > 0) {
                offset = -1*c;
            }
            else {
                offset = magnet.get_length() + c;
            }
            
            double mag_center_x = magnet.get_pos(0) + offset;
            double mag_center_y = magnet.get_pos(1);
            double mag_center_z = magnet.get_pos(2);
            
            double x = grid_point.getX() - mag_center_x;
            double y = grid_point.getY() - mag_center_y;
            double z = grid_point.getZ() - mag_center_z;
            
            double temp_B1 = 0.0;
            double temp_B2 = 0.0;
            double temp_B3 = 0.0;
            calc_grid_B_comps(factor, a, b, c, z, y, x, temp_B1, temp_B2, temp_B3);

            hold_B1 += temp_B1;
            hold_B2 += temp_B2;
            hold_B3 += temp_B3;
            
            i -= 2;
        }
        //bool in_margin = B_within_margin(B_center_val, hold_B1, hold_B2, hold_B3);
        bool in_margin = true;
        if(in_margin) {
            grid_B1 = hold_B1;
            grid_B2 = hold_B2;
            grid_B3 = hold_B3;
        }

        grid_point_B.setX(grid_B3);
        grid_point_B.setY(grid_B2);
        grid_point_B.setZ(grid_B1);
    }
    
    return grid_point_B;
}


ThreeVec calc_quadrupole_Bq(ThreeVec &grid_point, Magnet &magnet, double charge) {
    ThreeVec grid_point_B;

    // rotate coordinates 45 degrees
    double rot_x = grid_point.getX();
    double rot_y = grid_point.getY()/sqrt(2.) + grid_point.getZ()/sqrt(2.);
    double rot_z = -grid_point.getY()/sqrt(2.) + grid_point.getZ()/sqrt(2.);
    //double rot_x = grid_point.getX();
    // double rot_y = grid_point.getY();
    // double rot_z = grid_point.getZ();

    // aperture = 2.0*magnet.get_width()
    double block_side = magnet.get_height() - magnet.get_width();

    double block_Bx = 0.0;
    double block_By = 0.0;
    double block_Bz = 0.0;

    double offset = magnet.get_width() + 0.5*block_side;

    double a = 0.5*magnet.get_length();  // half length for permanent magnet
    double b = 0.5*magnet.get_width();   // half width for permanent magnet
    double c = 0.5*block_side;  // half height for permanent magnet

    int focus_direction;
    if(magnet.get_axis_of_magnetization() == 'z') {
        focus_direction = 1;
    }
    else {
        focus_direction = -1;
    }

    double factor = charge * focus_direction * 1.333*magnet.get_Br()/(4 * M_PI);
    //std::cerr<<magnet.get_Br() << std::endl;

    // vertical mags [referenced from -z mag] (real z is formula z, real y is formula y, real x is formula x)
    int i = 1;
    while(i > -2) {
        double x1 = rot_x - a - magnet.get_pos(0);
        double y1 = i*(rot_y - magnet.get_pos(1));
        double z1 = i*(rot_z + i*offset - magnet.get_pos(2));
        
        double temp1_B1;
        double temp1_B2;
        double temp1_B3;
        //std::cerr << temp1_B1 << std::endl;
        calc_grid_B_comps(factor, a, b, c, x1, y1, z1, temp1_B1, temp1_B2, temp1_B3, false);
        
        // block_Bx += i*temp1_B1;
        block_Bx += i*temp1_B1;
        block_By += i*temp1_B2;
        block_Bz += i*temp1_B3;
        // temp1_B1 = 0.0;
        // temp1_B2 = 0.0;
        // temp1_B3 = 0.0;
        
        i -= 2;
    }

    // horizontal mags [referenced from -y mag] (real z is formula y, real y is formula -z, real x is formula x)
    int j = 1;
    while (j > -2) {
        double x2 = rot_x - a - magnet.get_pos(0);
        double z2 = j*(-rot_y - j*offset - magnet.get_pos(1));
        double y2 = j*(rot_z - magnet.get_pos(2));
        
        double temp2_B1;
        double temp2_B2;
        double temp2_B3;
        calc_grid_B_comps(factor, a, b, c, x2, y2, z2, temp2_B1, temp2_B2, temp2_B3, false);

        block_Bx += -j*temp2_B1;
        // block_Bx += -temp2_B1;
        block_By += -j*temp2_B3;
        block_Bz += j*temp2_B2;
        // temp2_B1 = 0.0;
        // temp2_B2 = 0.0;
        // temp2_B3 = 0.0;
        
        j -= 2;
    }

    // nearest corner of off-axis mags is a distance of radius, i.e. b

    // offset for off-axis mags
    //double offset_c30_s60 = 1.1*magnet.get_width()*(sqrt(3.)/2);
    //double offset_s30_c60 = 1.1*magnet.get_width()/2;
    double offset_c30_s60 = 1.2*magnet.get_width() * (1 + sqrt(5.))/4.;
    double offset_s30_c60 = 1.2*magnet.get_width() * sqrt(10 - 2*sqrt(5.))/4.;

    // z-axis rotated 30 mags [referenced to -z mag] (real z is formula y, real y is formula -z, real x is formula x)
    int k = 1;
    while(k > -2) {
        double x3 = rot_x - a - magnet.get_pos(0);
        double z3 = k*(-rot_y + k*(offset_s30_c60 + b) - magnet.get_pos(1));
        double y3 = k*(rot_z + k*(offset_c30_s60 + c) - magnet.get_pos(2));
        
        double temp3_B1;
        double temp3_B2;
        double temp3_B3;
        calc_grid_B_comps(factor, a, c, b, x3, y3, z3, temp3_B1, temp3_B2, temp3_B3, false);
        
        block_Bx += k*temp3_B1;
        // block_Bx += temp3_B1;
        block_By += -k*temp3_B3;
        block_Bz += k*temp3_B2;
        // temp3_B1 = 0.0;
        // temp3_B2 = 0.0;
        // temp3_B3 = 0.0;
        
        k -= 2;
    }

    // y-axis rotated 30 mags [referenced to -y mag] (real z is formula -z, real y is formula -y, real x is formula x)
    int l = 1;
    while(l > -2) {
        double x4 = rot_x - a - magnet.get_pos(0);
        double y4 = l*(-rot_y - l*(offset_c30_s60 + c) - magnet.get_pos(1));
        double z4 = l*(-rot_z - l*(offset_s30_c60 + b) - magnet.get_pos(2));
        
        double temp4_B1;
        double temp4_B2;
        double temp4_B3;
        calc_grid_B_comps(factor, a, c, b, x4, y4, z4, temp4_B1, temp4_B2, temp4_B3, false);
        
        block_Bx += -l*temp4_B1;
        // block_Bx += -temp4_B1;
        block_By += -l*temp4_B2;
        block_Bz += -l*temp4_B3;
        // temp4_B1 = 0.0;
        // temp4_B2 = 0.0;
        // temp4_B3 = 0.0;
        
        l -= 2;
    }

    // z-axis rotated 60 mags [referenced to -z mag] (real z is formula -z, real y is formula -y, real x is formula x)
    int m = 1;
    while(m > -2) {
        double x5 = rot_x - a - magnet.get_pos(0);
        double y5 = m*(-rot_y + m*(offset_c30_s60 + c) - magnet.get_pos(1));
        double z5 = m*(-rot_z - m*(offset_s30_c60 + b) - magnet.get_pos(2));
        
        double temp5_B1;
        double temp5_B2;
        double temp5_B3;
        calc_grid_B_comps(factor, a, c, b, x5, y5, z5, temp5_B1, temp5_B2, temp5_B3, false);
        
        block_Bx += -m*temp5_B1;
        // block_Bx += -temp5_B1;
        block_By += -m*temp5_B2;
        block_Bz += -m*temp5_B3;
        // temp5_B1 = 0.0;
        // temp5_B2 = 0.0;
        // temp5_B3 = 0.0;
        
        m -= 2;
    }

    // y-axis rotated 60 mags [referenced to -y mag] (real z is formula -y, real y is formula z, real x is formula x)
    int n = 1;
    while(n > -2) {
        double x6 = rot_x - a - magnet.get_pos(0);
        double z6 = n*(rot_y + n*(offset_s30_c60 + b) - magnet.get_pos(1));
        double y6 = n*(-rot_z - n*(offset_c30_s60 + c) - magnet.get_pos(2));
        
        double temp6_B1;
        double temp6_B2;
        double temp6_B3;
        calc_grid_B_comps(factor, a, c, b, x6, y6, z6, temp6_B1, temp6_B2, temp6_B3, false);
        
        block_Bx += n*temp6_B1;
        // block_Bx += temp6_B1;
        block_By += n*temp6_B3;
        block_Bz += -n*temp6_B2;
        // // temp6_B1 = 0.0;
        // temp6_B2 = 0.0;
        // temp6_B3 = 0.0;
        
        n -= 2;
    }

    // rotate B coordinates back 45 degrees
    double rot_By = -block_By/sqrt(2.) + block_Bz/sqrt(2.);
    double rot_Bz = block_By/sqrt(2.) + block_Bz/sqrt(2.);
    // rot_By = block_By;
    // rot_Bz = block_Bz;


    grid_point_B.setX(block_Bx);
    grid_point_B.setY(rot_By);
    grid_point_B.setZ(rot_Bz);

    return grid_point_B;
}

ThreeVec calc_halbach_B(ThreeVec &grid_point, Magnet &magnet, double charge) {
    ThreeVec grid_point_B;
    grid_point_B.setX(0.0);

    // rotated
    //double y = grid_point.getY()/sqrt(2.) + grid_point.getZ()/sqrt(2.) - magnet.get_pos(1);
    //double z = -grid_point.getY()/sqrt(2.) + grid_point.getZ()/sqrt(2.) - magnet.get_pos(2);
    // not rotated
    double y = grid_point.getY() - magnet.get_pos(1);
    double z = grid_point.getZ() - magnet.get_pos(2);

    int focus_direction;
    if(magnet.get_axis_of_magnetization() == 'z') {
        focus_direction = 1;
    }
    else {
        focus_direction = -1;
    }

    double k2 = charge * focus_direction * pow( cos(M_PI/12), 2.0 ) * sin(M_PI/6) * (6/M_PI);
    double By = (z/magnet.get_width()) * magnet.get_Br() * 2 * (1 - (magnet.get_width()/magnet.get_height())) * k2;
    double Bz = (y/magnet.get_width()) * magnet.get_Br() * 2 * (1 - (magnet.get_width()/magnet.get_height())) * k2;

    //double rot_By = By/sqrt(2.) - Bz/sqrt(2.);
    //double rot_Bz = By/sqrt(2.) + Bz/sqrt(2.);

    grid_point_B.setY(By);
    grid_point_B.setZ(Bz);

    return grid_point_B;
}