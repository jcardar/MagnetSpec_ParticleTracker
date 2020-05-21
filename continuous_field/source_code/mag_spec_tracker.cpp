/*Leap Frog Numerical Method Code
 *for tracking a particle in static B field
 *Author: Jason Cardarelli
 *Prof. Alexander Thomas
 */
/////Development notes: currently the code assumes that particle beam starts in the region of the magnet. Got to fix that!

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include "my_functions.h"
#include "threevector.h"
#include "threematrix.h"
#include "particle.h"
#include "beam.h"
#include "screen.h"
#include "magnet.h"



int main(int argc, char *argv[])
{
    double time     {0.0};                            //Define a variable time in which will be stepped over
    double del_time {2};                              //Define a time step
    bool time_step_test = false;

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
    std::ofstream outfile_del_t   ("../data/DEL_T.csv");
	std::ofstream outfile_dump;
    //std::ifstream infile;

///////////////////
    //Define Particle Beam:
    int num_par          {1};
    double charge        {-1.0};
    double mass          {1.0};
    double energy0       {10.0};              //Normalized Energy = gamma
    double energy_spread {0.0};
    ThreeVec initial_position(0.0, 0.0, 0.0);
    ThreeVec initial_position_spread(0.0, 0.0, 0.0);
    ThreeVec initial_angular_direction(0.0, M_PI/2.000000, M_PI/2.000000);
    ThreeVec initial_angular_spread(0.0, 0.0, 0.0);
    Beam electron_beam(num_par, charge, mass, energy0, energy_spread, initial_position, initial_position_spread, initial_angular_direction, initial_angular_spread);


/////////////////////
    //Define Magnets:
    const int num_magnets{1};
    Magnet magnet[num_magnets];
        for(int ii{0}; ii<num_magnets; ii++)
        {
            switch(ii)
            {
            case 0:
                magnet[ii].set_B0(0, 0.0); magnet[ii].set_B0(1, 0.0); magnet[ii].set_B0(2, 1.0);
                magnet[ii].set_pos(0, -1958.0*3.0); magnet[ii].set_pos(1, 0.0); magnet[ii].set_pos(2, 0.0);
                magnet[ii].set_length(1958.0*50.0);
                magnet[ii].set_width(1957.95*20.0);
                magnet[ii].set_height(1957.95*8.0);
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










    //MAIN LOOP FOR STEPPING PARTICLES THROUGH SYSTEM:
    for(int ii{0}; ii < num_par; ii++)
    {
        int magnet_counter = 0;
        int screen_counter = 0;
        time               = 0.0;
        std::cout << "particle number " << (ii+1) << '\n';

        if(ii>0)
            {electron_beam.next_particle(electron_beam.m_particle_counter);}
        
        Particle electron = electron_beam.get_particle();
        electron.set_outfiles(outfile_time, outfile_xpos, outfile_ypos, outfile_zpos, outfile_px, outfile_py, outfile_pz, outfile_energy);

        electron.set_time(time);

        outfile_del_t << del_time << '\n';

        outfile_part_writeAndComma(electron); 

        double particle_time_limit = 2*M_PI*10;
        step_through_magnet_mag_boris(electron,magnet[magnet_counter],time,del_time,particle_time_limit);
        
        while(++magnet_counter < num_magnets)
        {               
                //Loop through magnets after first magnet
                //First check that the magnet is heading toward the next magnet in the x-direction.
                //Determine if next magnet is in front of or behind the particle's position in the x-dimension:
            double dist_x_to_mag = magnet[(magnet_counter)].get_pos(0) - electron.get_pos(0);
                //Then check if particle is heading in that direection
            if( ( (dist_x_to_mag > 0) && ((electron.get_vel(0)) > 0) ) || ( (dist_x_to_mag < 0) && ((electron.get_vel(0)) < 0) ) || ( dist_x_to_mag == 0) )
            {           //First check that electron is moving toward next magnet
                double time_btwn_mags = 0.0;
                if(dist_x_to_mag !=0)
                    { time_btwn_mags        = (((magnet[magnet_counter].get_pos(0))) - (electron.get_pos(0)))/(electron.get_vel(0)); }
                double y_at_time      = 0.0;
                y_at_time             = (electron.get_pos(1)) + ((electron.get_vel(1))*(time_btwn_mags));
                double z_at_time      = 0.0;
                z_at_time             = (electron.get_pos(2)) + ((electron.get_vel(2))*(time_btwn_mags));
                time                  = time + time_btwn_mags;

                bool check_y_bound, check_z_bound;
                check_y_bound = ((y_at_time >= (magnet[magnet_counter].get_pos(1) - (magnet[magnet_counter].get_width()/2.0) ) ) && (y_at_time <= (magnet[magnet_counter].get_pos(1) + (magnet[magnet_counter].get_width()/2.0) )));
                check_z_bound = ((z_at_time >= (magnet[magnet_counter].get_pos(2) - (magnet[magnet_counter].get_height()/2.0) ) ) && (z_at_time <= (magnet[magnet_counter].get_pos(2) + (magnet[magnet_counter].get_height()/2.0) )));
                if( check_y_bound && check_z_bound )
                {       //Then check that particle ends up in magnet region

                    electron.set_pos(0 ,(magnet[magnet_counter].get_pos(0))); 
                    electron.set_pos(1, y_at_time);
                    electron.set_pos(2, z_at_time);
                    
                    outfile_part_comma(electron);
                    step_through_magnet_mag_boris(electron, magnet[magnet_counter], time, del_time, particle_time_limit);
                    
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
        if(time_step_test)
            {half_time_step(del_time);}
    }   //<-END OF PARTICLE STEPPING 'FOR' LOOP

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
    //outfile_array_size. close();
    outfile_del_t.  close();
    //infile.close();

    return 0;
}
