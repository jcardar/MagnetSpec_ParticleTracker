/*
 *Main code for tracking a particle in analytic dipole B field
 *Author: Jason Cardarelli, Evan Mahler
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
#include "beam.h"
#include "screen.h"
#include "magnet.h"
#include <iomanip>

//#define THREAD_NUM 5
#include <omp.h>

int main(int argc, char *argv[])
{
    
    double time     {0.0};                           //Define a variable time in which will be stepped over
    double del_time {0.07};                           //Define a time step
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
    std::ofstream outfile_particle_on_screen ("../data/PARTICLE_ON_SCREENS.csv");
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
    
    std::vector<double> PmagDim;
    readPermanentMagDim(infile, num_magnets, PmagDim);
    
    std::vector<char> magnet_axis_info;
    readMagAxis(infile, num_magnets, magnet_axis_info);

    std::vector<char> magnet_type;
    readMagnetType(infile, num_magnets, magnet_type);

    int num_par = readNumOf(infile);
    std::vector<std::vector<double>> beam_info;
    readBeam(infile, beam_info);

    std::vector<std::vector<double>> beam_spread_info;
    readSpread(infile, beam_spread_info);

    int num_screens;
    std::vector<std::vector<std::vector<double>>> screen_info;
    readScreen(infile, num_screens, screen_info);
    
    std::vector<int> init_types;
    ReadInitTypes(infile, init_types);

    double mu_0 = ReadMu0(infile);

    double charge = ReadSpeciesCharge(infile);

    //char dipole_field_type = ReadDipoleMagnetFieldType(infile);
    
///////////////////
    //Define Particle Beam:
    //std::cout << "Charge is " << charge << '\n';
    double mass       {1.0};    //hard-coded, no input from user
    double energy0 =         beam_info[1][0];   //Normalized Energy = gamma
    double energy_spread =   beam_spread_info[1][0];
    //double int_y = sqrt(-(1.0 - energy0*energy0)/(energy0*energy0))*energy0;
    //std::cout << std::setprecision(15);
    //del_time = 2*M_PI*energy0*100;
    //std::cout << int_y << std::endl;
    //std::cerr << beam_info[0][0] << std::endl;
    ThreeVec initial_position(beam_info[0][0], beam_info[0][1], beam_info[0][2]);
    ThreeVec initial_position_spread(beam_spread_info[0][0], beam_spread_info[0][1], beam_spread_info[0][2]);
    ThreeVec initial_angular_direction(beam_info[2][0], beam_info[2][1], beam_info[2][2]);
    ThreeVec initial_angular_spread(beam_spread_info[2][0], beam_spread_info[2][1], beam_spread_info[2][2]);
    //std::cerr << initial_angular_spread << "\n";


    Beam electron_beam(num_par, charge, mass, energy0, energy_spread, initial_position, initial_position_spread, initial_angular_direction, initial_angular_spread, static_cast<Beam::PositionInitializationTypes>(init_types[0]), static_cast<Beam::EnergyInitializationTypes>(init_types[1]), static_cast<Beam::DivergenceInitializationTypes>(init_types[2]));
    //std::cout << "Charge is " << charge << '\n';
    num_par = electron_beam.get_num_particles();

/////////////////////
    //Define Magnets:
    Magnet magnet[num_magnets];
    for(int ii=0; ii<num_magnets; ++ii) {
        
        magnet[ii].set_B0(0, 0.0);
        magnet[ii].set_B0(1, 0.0);
        magnet[ii].set_B0(2, 0.0);
        magnet[ii].set_axis_of_magnetization( magnet_axis_info[ii] );
        magnet[ii].set_height_of_dipole_block( PmagDim[ii] );
        
        if(magnet_axis_info[ii] == 'x') {
            magnet[ii].set_B0(0, magnet_info[2][ii][0]);
        }
        else if(magnet_axis_info[ii] == 'y') {
            magnet[ii].set_B0(1, magnet_info[2][ii][0]);
        }
        else {
            magnet[ii].set_B0(2, magnet_info[2][ii][0]);  // magnet_axis_info == 'z'
        }
        
        magnet[ii].set_pos(0, magnet_info[1][ii][0]); 
        magnet[ii].set_pos(1, magnet_info[1][ii][1]); 
        magnet[ii].set_pos(2, magnet_info[1][ii][2]);
        magnet[ii].set_length(magnet_info[0][ii][1]);
        magnet[ii].set_width(magnet_info[0][ii][0]);
        magnet[ii].set_height(magnet_info[0][ii][2]);
        magnet[ii].set_outfile(outfile_magnets);
        magnet[ii].set_type(magnet_type[ii]);
        if(magnet_type[ii] == 'q' || magnet_type[ii] == 'h')
        {
            magnet[ii].set_Br(magnet_info[2][ii][0]);
        }
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
            screen[ii].set_index(ii);
            screen[ii].set_pos(0, screen_info[1][ii][0]);
            screen[ii].set_pos(1, screen_info[1][ii][1]);
            screen[ii].set_pos(2, screen_info[1][ii][2]);
            screen[ii].set_length(screen_info[0][ii][0]);
            screen[ii].set_height(screen_info[0][ii][1]);
            screen[ii].set_angle_about_z(screen_info[2][ii][0], 'r');
            screen[ii].set_angle_about_y(screen_info[2][ii][1], 'r');
            screen[ii].set_angle_about_x(screen_info[2][ii][2], 'r');
            screen[ii].set_outfile(outfile_screens);
            screen[ii].set_outfile_particle_on_screen(outfile_particle_on_screen);
            outfile_screen_single (screen[ii], ii);
        }   //<---END OF SCREEN 'FOR' LOOP//
    }
    outfile_part_on_screen_first_line(screen[0]);
    
    double particle_time_limit = (2*M_PI*energy0)*10.0;
    bool intersected_screen_while_in_magnet;
    //MAIN LOOP FOR STEPPING PARTICLES THROUGH SYSTEM:
    omp_set_dynamic(0);     // Explicitly disable dynamic teams
    omp_set_num_threads(4); // Use 4 threads for all consecutive parallel regions
    #pragma omp parallel for ordered
    for(int ii{0}; ii < num_par; ++ii)
    {
        if(ii == 0){
        std::cout << "Number of parallel threads created: " << omp_get_num_threads() << '\n';
        }
        int magnet_counter = 0;
        int screen_counter = 0;
        time               = 0.0;
    
        //std::cout << "particle number " << (ii+1) << std::endl;

        if(ii>0)
            {
                
                electron_beam.next_particle();
            }
        
        Particle electron = electron_beam.get_particle();
        //std::cerr << "Energy after setting new electron is " << electron.get_energy() << std::endl;
        //std::cout << "Charge is " << electron.get_charge() << '\n';
        //std::cerr << "X-Position is " << electron.get_pos(0) << std::endl;
        //std::cerr << electron.get_pos(1) << std::endl;
        //std::cerr << "Energy is " << electron.get_energy() << std::endl;
        //std::cerr << "p is " << electron.get_p() << std::endl;
        electron.set_outfiles(outfile_time, outfile_xpos, outfile_ypos, outfile_zpos, outfile_px, outfile_py, outfile_pz, outfile_energy);

        electron.set_time(time);

        outfile_del_t << del_time << '\n';
        #pragma omp ordered
        {
            outfile_part_write(electron);
        }
        intersected_screen_while_in_magnet = move_through_magnets_general(magnet, num_magnets, electron, time, del_time, mu_0, particle_time_limit, screen, num_screens);
        //if(dipole_field_type=='d')
        //{
        //    intersected_screen_while_in_magnet = move_through_magnets_dipole(magnet, num_magnets, electron, time, del_time, mu_0, particle_time_limit, screen, num_screens);
        //}
        //else if(dipole_field_type=='u')
        //{
        //    intersected_screen_while_in_magnet = move_through_magnets_uniform(magnet, num_magnets, electron, time, del_time, mu_0, particle_time_limit, screen, num_screens);
        //}
        if(!intersected_screen_while_in_magnet || intersected_screen_while_in_magnet) 
        {
            move_to_screens(screen, num_screens, electron, ii);
        }
        #pragma omp ordered
        {
            outfile_part_newline(electron);
        }
        if(time_step_test)
            {half_time_step(del_time);}
    }   //<-END OF PARTICLE STEPPING 'FOR' LOOP

    outfile_readme.            close();
    //outfile_grid. close();
    outfile_time.              close();
    outfile_xpos.              close();
    outfile_ypos.              close();
    outfile_zpos.              close();
    outfile_px.                close();
    outfile_py.                close();
    outfile_pz.                close();
    outfile_magnets.           close();
    outfile_energy.            close();
    outfile_particle_on_screen.close();
    //outfile_array_size. close();
    outfile_del_t.             close();
    infile.close();

    return 0;
}
