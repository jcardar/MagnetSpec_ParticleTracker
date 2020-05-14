/*Leap Frog Numerical Method Code
 *for tracking a particle in static B and possible E field
 *Author: Jason Cardarelli
 *NERS 574: Computational Plasma Physics
 *Prof. Alexander Thomas
 */

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include "my_functions.h"
#include "threevector.h"
#include "threematrix.h"
#include "particle.h"
#include "screen.h"
#include "magnet.h"


//double M_c = 299792458;           //SPEED OF LIGHT
//double globe_central_energy_x = 1000/0.511;     //MeV/E0
//double globe_central_gamma_x  = globe_central_energy_x/0.511;

//double globe_gamma = 90.0/0.511;

int main(int argc, char *argv[])
{
    int num_par     {5};
    double time     {0.0};                            //Define a variable time in which will be stepped over
    double del_time {0.25};                           //Define a time step

    // std::ifstream bmap( "bmap_in.txt" );   
    std::ofstream outfile_readme  ("../data/README.txt");
	std::ofstream outfile_time    ("../data/TIME.csv");
	std::ofstream outfile_xpos    ("../data/XPOS.csv");
    std::ofstream outfile_ypos    ("../data/YPOS.csv");
    std::ofstream outfile_zpos    ("../data/ZPOS.csv");
	std::ofstream outfile_px      ("../data/MOMENTUM_x.csv");
	std::ofstream outfile_py      ("../data/MOMENTUM_y.csv");
    std::ofstream outfile_pz      ("../data/MOMENTUM_z.csv");
    std::ofstream outfile_energy  ("../data/ENERGY.csv");
    std::ofstream outfile_magnets ("../data/MAGNETS.csv");
	std::ofstream outfile_screens ("../data/SCREENS.csv");
    std::ofstream outfile_array_size ("../data/ARRAY_SIZE.csv");
    std::ofstream outfile_del_t   ("../data/DEL_T.csv");
	std::ofstream outfile_dump;
    //std::ifstream infile;


/////////////////////
    //Define Magnets:
    const int num_magnets{1};
    Magnet magnet[num_magnets];
    //Have code read values from user, and then populate it here
    //Replace switch with just magnet[ii] and loop through user input values
        for(int ii{0}; ii<num_magnets; ii++)
        {
            switch(ii)
            {
            case 0:
                magnet[ii].set_B0(0, 0.0); magnet[ii].set_B0(1, 0.0); magnet[ii].set_B0(2, 1.0);
                magnet[ii].set_pos(0, -1958.0-2.0); magnet[ii].set_pos(1, 0.0); magnet[ii].set_pos(2, 0.0);
                magnet[ii].set_length(1958.0*2.1);
                magnet[ii].set_width(1957.95*8.0);
                magnet[ii].set_height(17.597*2.0);
                magnet[ii].set_outfile(outfile_magnets);
                outfile_uniform_magnet(magnet[ii], ii);
                break;
            
            case 1:
                {
                magnet[ii].set_B0(0, 0.0); magnet[ii].set_B0(1, 0.0); magnet[ii].set_B0(2, 1.0);
                double mag_x_pos = 17.597 + magnet[0].get_pos(0) + magnet[0].get_length();
                magnet[ii].set_pos(0, mag_x_pos); magnet[ii].set_pos(1, 0.0); magnet[ii].set_pos(2, 0.0);
                magnet[ii].set_length(175.972);
                magnet[ii].set_width(58.657);
                magnet[ii].set_height(17.597);
                magnet[ii].set_outfile(outfile_magnets);
                outfile_uniform_magnet(magnet[ii], ii);
                break;
                }   //<---END OF CASE 1//
            }   //<---END OF SWITCH STATEMENT//
        }   //<---END OF MAGNET 'FOR' LOOP//


////////////////////
    //Define Screens:
    const int num_screens = 1;
    Screen screen[num_screens];
    for(int ii{0}; ii < num_screens; ii++)
    {
        switch(ii)
        {
            case 0:
            {
                double srn_x_pos = 58.65 + magnet[num_magnets-1].get_pos(0)+magnet[num_magnets-1].get_length();
                double srn_y_pos = -17.597 + magnet[num_magnets-1].get_pos(1);
                double srn_z_pos = magnet[num_magnets-1].get_pos(2);
                screen[ii].set_pos(0, srn_x_pos); screen[ii].set_pos(1, srn_y_pos); screen[ii].set_pos(2, srn_z_pos);
                screen[ii].set_length(55);
                screen[ii].set_height(17.597);
                screen[ii].set_angle(90.0);
                screen[ii].set_outfile(outfile_screens);
                outfile_screen_single (screen[ii], ii);
                break;
            }   //<---END OF CASE 0//
            case 1:
                screen[ii].set_outfile(outfile_screens);
                outfile_screen_single (screen[ii], ii);
                break;
        }   //<---END OF SWITCH STATEMENT//
    }   //<---END OF SCREEN 'FOR' LOOP//


///////////////////////
    //Define Particles:

    double energy0 = (1000.0+0.511)/0.511;            //Normalized Central Energy (gamma)
    double p0_mag = sqrt(energy0*energy0 - 1);
    //std::cerr << p0_mag << '\n';
    ThreeVec p0(p0_mag, 0.0, 0.0);            //Kinetic energy (in MeV)
    ThreeVec radius_p0;
    {
    double percent_en_spread_x = 0.50;
    double percent_en_spread_y = 0.00;
    double percent_en_spread_z = 0.00;
    double p0_spread_x = energy0*(percent_en_spread_x*energy0)/p0_mag;
    double p0_spread_y = energy0*(percent_en_spread_y*energy0)/p0_mag;
    double p0_spread_z = energy0*(percent_en_spread_z*energy0)/p0_mag;
    radius_p0.setX(p0_spread_x);
    radius_p0.setY(p0_spread_y);
    radius_p0.setZ(p0_spread_z);
    }
    
    ThreeVec r0(magnet[0].get_pos(0),magnet[0].get_pos(1),magnet[0].get_pos(2));                   //INITIAL CENTRAL POSTIION OF // //
    ThreeVec radius_r0(0.0,0.0,0.0);            //RADIUS OF INT POSITIIONS IN PHASE SPACE
    int qe = -1;                                //CHARGE OF PARTICLE SPECIES (normalizd to charge of proton)

    //ThreeVec energy;

    //Particle electron(r0, v0, qe, time, outfile_time, outfile_xpos, outfile_ypos, outfile_zpos, outfile_px, outfile_py, outfile_pz, outfile_energy);
    Particle electron(r0, p0, qe, time, outfile_time, outfile_xpos, outfile_ypos, outfile_zpos, outfile_px, outfile_py, outfile_pz, outfile_energy);
    
    double initial_x;
    double initial_y;
    double initial_z;
    double initial_px;
    double initial_py;
    double initial_pz;

    Particle::InitializationTypes initialize = electron.INITIALIZE_UNIFORM_EN_DIST;
        if(initialize == electron.INITIALIZE_POINT_SOURCE_GAUS || initialize == electron.INITIALIZE_POINT_SOURCE_UNI)
                {num_par = num_par * 5;}

    //VALUES FOR UNIFORM DISTRIBUTION
    int posx_counter = 0;
    int posy_counter = 0;
    int posz_counter = 0;
    int velx_counter = 0;
    int vely_counter = 0;
    int velz_counter = 0;

    //VALUES FOR POINT SOURCE DISTRIBUTION
    const double length_before  = 17.5; //1cm
    double divergence           = 0.01; //10 mrad
    int energy_cycle            = 0;
    double tot_p;
    double p_x_val;
    double p_y_val;
    double p_z_val;





    int array_counter_max_cols=0;
    for(int ii{0}; ii < num_par; ii++)
    {
        int magnet_counter = 0;
        int screen_counter = 0;
        time = 0.0;
        std::cout << "particle number " << (ii+1) << '\n';
        
        switch (initialize)
        {
        case electron.INITIALIZE_GAUSSIAN:
            initial_x  = 0.0;
            initial_y  = gaussian()*(radius_r0.getY() ) + r0.getY();
            initial_z  = gaussian()*(radius_r0.getZ() ) + r0.getZ();
            initial_px = gaussian()*(radius_p0.getX() ) + p0.getX();
            initial_py = gaussian()*(radius_p0.getY() ) + p0.getY();
            initial_pz = gaussian()*(radius_p0.getZ() ) + p0.getZ();
            break;
        
        case electron.INITIALIZE_UNIFORM_EN_DIST:   //sweeps initial position from (0,-y_max,-z_max) to (0,y_max,z_max), keeping total energy the central energy
            uniform_en_dist(initial_x, initial_y, initial_z, initial_px, initial_py, initial_pz, length_before, &posy_counter, &posz_counter, r0, radius_r0, p0, num_par);
            break;    
        
        case electron.INITIALIZE_UNIFORM_POS_DIST:
            //add later
            break;

        case electron.INITIALIZE_POINT_SOURCE_GAUS:
            initial_x = 0.0;
            if(ii-(5*energy_cycle) == 0)
            {
                p_x_val = gaussian()*(radius_p0.getX() ) + p0.getX();
                p_y_val = gaussian()*(radius_p0.getY() ) + p0.getY();
                p_z_val = gaussian()*(radius_p0.getZ() ) + p0.getZ();
                tot_p       = sqrt(p_x_val*p_x_val + p_y_val*p_y_val + p_z_val*p_z_val);
                initial_y   = 0.0;
                initial_z   = 0.0;
                initial_px  = tot_p;
                initial_py  = 0.0;
                initial_pz  = 0.0;
            }
            else if(ii-(5*energy_cycle) == 1)
            {
                initial_y   = length_before*tan(divergence/2);
                initial_z   = 0.0;
                initial_px = tot_p*cos(divergence/2.0);
                initial_py = tot_p*sin(divergence/2.0);
                initial_pz = 0.0;
            }
            else if(ii-(5*energy_cycle) == 2)
            {
                initial_y = -length_before*tan(divergence/2);
                initial_z = 0.0;
                initial_px = tot_p*cos(divergence/2.0);
                initial_py = -tot_p*sin(divergence/2.0);
                initial_pz = 0.0;
            }
            else if(ii-(5*energy_cycle) == 3)
            {
                initial_y = 0.0;
                initial_z = length_before*tan(divergence/2);
                initial_px = tot_p*cos(divergence/2.0);
                initial_py = 0.0;
                initial_pz = tot_p*sin(divergence/2.0);
            }
            else if(ii-(5*energy_cycle) == 4)
            {
                initial_y = 0.0;
                initial_z = -length_before*tan(divergence/2);
                initial_px = tot_p*cos(divergence/2.0);
                initial_py = 0.0;
                initial_pz = -tot_p*sin(divergence/2.0);
                energy_cycle += 1;
            }
            electron.set_pos(-length_before,0.0,0.0);
            electron.set_p(initial_px, initial_py, initial_pz);
            electron.set_time(time);
            outfile_part_writeAndComma(electron);
            break;

        case Particle::INITIALIZE_POINT_SOURCE_UNI:
            initial_x   = 0.0;
            bool x_dist = true;
            bool y_dist = false;
            bool z_dist = false;
            if(ii-(5*energy_cycle) == 0)
            {
                p_x_val     = 0.0;
                p_y_val     = 0.0;
                p_z_val     = 0.0;
                if(x_dist) { p_x_val = uniform_dist_single(num_par/5, p0.getX(), radius_p0.getX(), velx_counter); }
                if(y_dist) { p_y_val = uniform_dist_single(num_par/5, p0.getY(), radius_p0.getY(), vely_counter); }
                if(z_dist) { p_z_val = uniform_dist_single(num_par/5, p0.getZ(), radius_p0.getZ(), velz_counter); }
                tot_p       = sqrt(p_x_val*p_x_val + p_y_val*p_y_val + p_z_val*p_z_val);
                initial_y   = 0.0;
                initial_z   = 0.0;
                initial_px  = tot_p;
                initial_py  = 0.0;
                initial_pz  = 0.0;
            }
            else if(ii-(5*energy_cycle) == 1)
            {
                initial_y   = length_before*tan(divergence/2);
                initial_z   = 0.0;
                initial_px = tot_p*cos(divergence/2.0);
                initial_py = tot_p*sin(divergence/2.0);
                initial_pz = 0.0;

                //del_time = del_time*0.5;
            }
            else if(ii-(5*energy_cycle) == 2)
            {
                initial_y = -length_before*tan(divergence/2);
                initial_z = 0.0;
                initial_px = tot_p*cos(divergence/2.0);
                initial_py = -tot_p*sin(divergence/2.0);
                initial_pz = 0.0;

                //del_time = del_time*0.5;
            }
            else if(ii-(5*energy_cycle) == 3)
            {
                initial_y = 0.0;
                initial_z = length_before*tan(divergence/2);
                initial_px = tot_p*cos(divergence/2.0);
                initial_py = 0.0;
                initial_pz = tot_p*sin(divergence/2.0);

                //del_time = 0.5*del_time;
            }
            else if(ii-(5*energy_cycle) == 4)
            {
                initial_y = 0.0;
                initial_z = -length_before*tan(divergence/2);
                initial_px = tot_p*cos(divergence/2.0);
                initial_py = 0.0;
                initial_pz = -tot_p*sin(divergence/2.0);

                //del_time = del_time*0.5;
                energy_cycle += 1;
            }
            electron.set_pos(-length_before,0.0,0.0);
            electron.set_p(initial_px, initial_py, initial_pz);
            electron.set_time(time);
            outfile_part_writeAndComma(electron);
            break;
        }

        outfile_del_t << del_time << '\n';
        if(initial_px == 0)
            { initial_px = 0.000000000001; }
        
        ThreeVec r_int(initial_x, initial_y, initial_z);
        ThreeVec p_init(initial_px, initial_py, initial_pz);
        electron.set_pos(r_int);
        electron.set_p(p_init);
        electron.set_time(time);
        outfile_part_writeAndComma(electron);        

        //step_through_magnet(electron, vn_plus, vn_minus, B0, rn_plus, rn_minus, time, del_time, width_bmap[magnet_counter], length_bmap[magnet_counter]);
        //step_through_magnet_mag_leap(electron, magnet[magnet_counter], vn_plus, vn_minus, rn_plus, rn_minus, time, del_time);
        step_through_magnet_mag_boris(electron,magnet[magnet_counter],time,del_time);
        
        while(++magnet_counter < num_magnets)
        {               //Loop through magnets after first magnet
            
            double dist_x_to_mag = magnet[(magnet_counter)].get_pos(0) - magnet[magnet_counter-1].get_pos(0);

            if( ( (dist_x_to_mag > 0) && ((electron.get_vel(0)) > 0) ) || ( (dist_x_to_mag < 0) && ((electron.get_vel(0)) < 0) ) )
            {           //First check that electron is moving toward next magnet
                double time_btwn_mags = 0.0;
                time_btwn_mags        = (((magnet[magnet_counter].get_pos(0))) - (electron.get_pos(0)))/(electron.get_vel(0)); 
                double y_at_time      = 0.0;
                y_at_time             = (electron.get_pos(1)) + ((electron.get_vel(1))*(time_btwn_mags));
                double z_at_time      = 0.0;
                z_at_time             = (electron.get_pos(2)) + ((electron.get_vel(2))*(time_btwn_mags));
                time                  = time + time_btwn_mags;

                if( ((y_at_time >= (magnet[magnet_counter].get_pos(1) - (magnet[magnet_counter].get_width()/2.0) ) ) && (y_at_time <= (magnet[magnet_counter].get_pos(1) + (magnet[magnet_counter].get_width()/2.0) ))))
                {       //Then check that particle ends up in magnet region
                //uncomment if using leapfrog method
                    //vn_minus = (electron.get_vel() ) - (((electron.get_vel())^B0)*del_time);
                    //vn_minus = electron.get_vel();
                    electron.set_pos(0 ,(magnet[magnet_counter].get_pos(0))); 
                    electron.set_pos(1, y_at_time);
                    electron.set_pos(2, z_at_time);
                    //rn_minus = (electron.get_pos() ) - ((electron.get_vel())*del_time);
                    //rn_minus = electron.get_pos();
                    
                    outfile_part_comma(electron);
                    //step_through_magnet_mag_leap(electron, magnet[magnet_counter], time, del_time);
                    step_through_magnet_mag_boris(electron, magnet[magnet_counter], time, del_time);
                    
                }
            }
        }   //<- end of magnetic while loop

        do
        {
            double x_intersect = (electron.get_pos(1) - (electron.get_pos(0)*electron.get_vel(1)/electron.get_vel(0)) - (screen[screen_counter].get_pos(1)) + (screen[screen_counter].get_pos(0)*tan(screen[screen_counter].get_angle('r'))))/(tan(screen[screen_counter].get_angle('r')) - (electron.get_vel(1)/electron.get_vel(0)));
            double y_intersect  = x_intersect*electron.get_vel(1)/electron.get_vel(0) + (electron.get_pos(1) - (electron.get_pos(0)*electron.get_vel(1)/electron.get_vel(0)));
            bool check1, check2;
            check1 = (screen[screen_counter].get_angle()==90 && y_intersect > screen[screen_counter].get_pos(1) && y_intersect < screen[screen_counter].get_pos(1)+screen[screen_counter].get_length() );
            check2 = ((x_intersect > screen[screen_counter].get_pos(0)) && (x_intersect < (screen[screen_counter].get_pos(0) + (screen[screen_counter].get_length()*cos(screen[screen_counter].get_angle('r'))))));
            if(check1 || check2)
            {
                
                double time_to_scrn = (x_intersect - electron.get_pos(0)) / (electron.get_vel(0));
                double z_intersect  = (electron.get_pos(2)) + ((electron.get_vel(2))*(time_to_scrn));
                time                = time + time_to_scrn;

                electron.set_pos(0, x_intersect);
                electron.set_pos(1, y_intersect);
                electron.set_pos(2, z_intersect);
                outfile_part_comma(electron);
                outfile_part_write(electron);
            }
        } while (++screen_counter < num_screens);
        

        outfile_part_newline(electron);
    }   //<-END OF PARTICLE 'FOR' LOOP

//outfile_array_size << array_counter_max_cols <<','<< num_par;

    outfile_readme. close();
    //outfile_grid. close();
    outfile_time.   close();
    outfile_xpos.   close();
    outfile_ypos.   close();
    outfile_zpos.   close();
    outfile_px.     close();
    outfile_py.     close();
    outfile_pz.     close();
    outfile_magnets.close();
    outfile_energy. close();
    outfile_array_size. close();
    outfile_del_t.  close();
    //infile.close();

    return 0;
}
