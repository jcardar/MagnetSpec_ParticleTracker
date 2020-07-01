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
#include <iomanip>



int main(int argc, char *argv[])
{
    double time     {0.0};                            //Define a variable time in which will be stepped over
    double del_time {0.05};                              //Define a time step
    //Hard-coded:
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
    std::ifstream infile          ("../data/analysis/input_deck.txt");

///////////////////
    //Read input deck values:
    //object numbers are defined here
    std::vector<std::string> desired_units;
    readUnits(infile, desired_units);

    int num_magnets;
    std::vector<std::vector<std::vector<double>>> magnet_info;
    readMagnet(infile, num_magnets, magnet_info);

    int num_par = readNumOf(infile);
    std::vector<std::vector<double>> beam_info;
    readBeam(infile, beam_info);

    std::vector<std::vector<double>> beam_spread_info;
    readSpread(infile, beam_spread_info);

    int num_screens;
    std::vector<std::vector<std::vector<double>>> screen_info;
    readScreen(infile, num_screens, screen_info);
    
///////////////////
    //Define Particle Beam:
    double charge     {-1.0};    //hard-coded, no input from user
    double mass       {1.0};    //hard-coded, no input from user
    double energy0 =         beam_info[1][0];   //Normalized Energy = gamma
    double energy_spread =   beam_spread_info[1][0];
    //double int_y = sqrt(-(1.0 - energy0*energy0)/(energy0*energy0))*energy0;
    //std::cout << std::setprecision(15);
    //del_time = 2*M_PI*energy0*100;
    //std::cout << int_y << std::endl;
    ThreeVec initial_position(beam_info[0][0], beam_info[0][1], beam_info[0][2]);
    ThreeVec initial_position_spread(beam_spread_info[0][0], beam_spread_info[0][1], beam_spread_info[0][2]);
    ThreeVec initial_angular_direction(beam_info[2][0], beam_info[2][1], beam_info[2][2]);
    ThreeVec initial_angular_spread(beam_spread_info[2][0], beam_spread_info[2][1], beam_spread_info[2][2]);
    
    Beam electron_beam(num_par, charge, mass, energy0, energy_spread, initial_position, initial_position_spread, initial_angular_direction, initial_angular_spread);


/////////////////////
    //Define Magnets:
    Magnet magnet[num_magnets];
    for(int ii=0; ii<num_magnets; ++ii) {
        
        magnet[ii].set_B0(0, magnet_info[2][ii][0]);
        magnet[ii].set_B0(1, magnet_info[2][ii][1]);
        magnet[ii].set_B0(2, magnet_info[2][ii][2]);
        //magnet[ii].set_pos(0, -1958.0*3000.0); magnet[ii].set_pos(1, 0.0); magnet[ii].set_pos(2, 0.0);
        magnet[ii].set_pos(0, magnet_info[1][ii][0]); 
        magnet[ii].set_pos(1, magnet_info[1][ii][1]); 
        magnet[ii].set_pos(2, magnet_info[1][ii][2]);
        //magnet[ii].set_length(1958.0*50000.0);
        magnet[ii].set_length(magnet_info[0][ii][1]);
        //magnet[ii].set_width(1957.95*20000.0);
        magnet[ii].set_width(magnet_info[0][ii][0]);
        magnet[ii].set_height(magnet_info[0][ii][2]);
        magnet[ii].set_outfile(outfile_magnets);
        outfile_uniform_magnet(magnet[ii], ii);
    }   //<---END OF MAGNET 'FOR' LOOP//


////////////////////
    //Define Screens:
    Screen screen[num_screens];
    
    if(num_screens <= 0) {
        
        Screen null_screen;
        null_screen.set_pos(0.0, 0.0, 0.0);
        null_screen.set_length(0.0);
        null_screen.set_height(0.0);
        null_screen.set_angle_about_x(0.0);
        null_screen.set_angle_about_y(0.0);
        null_screen.set_angle_about_z(0.0);
        null_screen.set_outfile(outfile_screens);
        outfile_screen_single(null_screen, 0);
    }
    else {
        for(int ii=0; ii<num_screens; ++ii) {

            screen[ii].set_pos(0, screen_info[1][ii][0]);
            screen[ii].set_pos(1, screen_info[1][ii][1]);
            screen[ii].set_pos(2, screen_info[1][ii][2]);
            screen[ii].set_length(screen_info[0][ii][0]);
            screen[ii].set_height(screen_info[0][ii][1]);
            screen[ii].set_angle_about_z(screen_info[2][ii][0], 'r');
            screen[ii].set_angle_about_y(screen_info[2][ii][1], 'r');
            screen[ii].set_angle_about_x(screen_info[2][ii][2], 'r');
            screen[ii].set_outfile(outfile_screens);
            outfile_screen_single (screen[ii], ii);
        }   //<---END OF SCREEN 'FOR' LOOP//
    }
    


    //MAIN LOOP FOR STEPPING PARTICLES THROUGH SYSTEM:
    for(int ii{0}; ii < num_par; ii++)
    {
        int magnet_counter = 0;
        int screen_counter = 0;
        time               = 0.0;
        std::cout << "particle number " << (ii+1) << '\n';

        if(ii>0)
            {
                electron_beam.next_particle();
            }
        
        Particle electron = electron_beam.get_particle();
        electron.set_outfiles(outfile_time, outfile_xpos, outfile_ypos, outfile_zpos, outfile_px, outfile_py, outfile_pz, outfile_energy);

        electron.set_time(time);

        outfile_del_t << del_time << '\n';

        outfile_part_write(electron);

        double particle_time_limit = (2*M_PI*energy0)*10.0;

        move_through_magnets(magnet, num_magnets, electron, time, del_time, particle_time_limit);
        move_to_screens(screen, num_screens, electron);

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
    infile.close();

    return 0;
}
