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


//double M_c = 299792458;                         //SPEED OF LIGHT


int main(int argc, char *argv[])
{
    const int num_par {5};
    double time       {0.0};                            //Define a variable time in which will be stepped over
    double del_time   {0.01};                           //Define a time step

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
	std::ofstream outfile_dump;
    std::ifstream infile;


    //Define Magnets:
    const int num_magnets           = 2;
    Magnet magnet[num_magnets];
    double length_bmap[num_magnets] = {0.5, 1.0};
    double width_bmap[num_magnets]  = {1.0, 1.0};
    ThreeVec mag_B0[num_magnets];
        for(int ii{0}; ii<num_magnets; ii++)
        {
            switch(ii)
            {
            case 0:
                mag_B0[ii].setX(0.0); mag_B0[ii].setY(0.0); mag_B0[ii].setZ(1.0);
                break;
            
            case 1:
                mag_B0[ii].setX(0.0); mag_B0[ii].setY(0.0); mag_B0[ii].setZ(1.0);
                break;
            }
        }

    ThreeVec mag_pos[num_magnets];
        for(int ii{0}; ii<num_magnets; ii++)
        {
            switch(ii)
            {
            case 0:
                mag_pos[ii].setX(0.0); mag_pos[ii].setY(0.0); mag_pos[ii].setZ(0.0);
                break;
            
            case 1:
                mag_pos[ii].setX(1.0); mag_pos[ii].setY(0.0); mag_pos[ii].setZ(0.0);
                break;
            }
        }

    for(int ii{0}; ii < num_magnets; ii++)
    {
        magnet[ii].set_length(length_bmap[ii]);
        magnet[ii].set_width(width_bmap[ii]);
        magnet[ii].set_B0(mag_B0[ii]);
        magnet[ii].set_pos(mag_pos[ii]);
        magnet[ii].set_outfile(outfile_magnets);
        outfile_uniform_magnet(magnet[ii], ii);
    }

////////////////////
    //Define Screens:
    const int num_screens = 1;
    Screen screens[num_screens];
    double scr_angle[num_screens]  = {45.0};
    double scr_length[num_screens] = {1.0};
    double scr_height[num_screens] = {1.0};
    ThreeVec scr_position;
    for(int ii{0}; ii < num_screens; ii++)
    {
        switch(ii)
        {
            case 0:

                break;
            case 1:

                break;
        }
    }



///////////////////////
    //Define Particles:
     
    const double length_before{0.0};


    ThreeVec B0(0.0,0.0,1.0);

    ThreeVec energy0(1000.0,0.0,0.0); //Kinetic energy (in MeV)
    ThreeVec radius_energy0(100.0,0.0,0.0);

    //ThreeVec v0(1.0,0.0,0.0);                   //INITIAL CENTRAL VELOCITY OF UNIFORM BEAM
    //ThreeVec radius_v0(0.9,0.0,0.0);           //RADIUS OF INT VELOCITIES IN PHASE SPACE
    
    ThreeVec r0(0.0,0.0,0.0);                   //INITIAL CENTRAL POSTIION OF // //
    ThreeVec radius_r0(0.0,0.0,0.0);            //RADIUS OF INT POSITIIONS IN PHASE SPACE
    int qe = -1;                                //CHARGE OF PARTICLE SPECIES (normalizd to charge of proton)

    ThreeVec energy;

    //Particle electron(r0, v0, qe, time, outfile_time, outfile_xpos, outfile_ypos, outfile_zpos, outfile_px, outfile_py, outfile_pz, outfile_energy);
    Particle electron(r0, qe, energy0, time, outfile_time, outfile_xpos, outfile_ypos, outfile_zpos, outfile_px, outfile_py, outfile_pz, outfile_energy);

    
    double initial_x;
    double initial_y;
    double initial_z;
    double initial_enx;
    double initial_eny;
    double initial_enz;

    //enum InitializationTypes
    //{
    //    INITIALIZE_GAUSSIAN,
    //    INITIALIZE_UNIFORM_DIST,
    //};
    Particle::InitializationTypes initialize = electron.INITIALIZE_GAUSSIAN;
    int posx_counter = 0;
    int posy_counter = 0;
    int posz_counter = 0;
    int velx_counter = 0;
    int vely_counter = 0;
    int velz_counter = 0;



    for(int ii{0}; ii < num_par; ii++)
    {
        int magnet_counter = 0;
        time = 0.0;
        std::cerr << "particle number " << (ii+1) << '\n';
        
        switch (initialize)
        {
        case electron.INITIALIZE_GAUSSIAN:
            initial_x = 0.0;
            initial_y = gaussian()*(radius_r0.getY() ) + (r0.getY() );
            initial_z = gaussian()*(radius_r0.getZ() ) + (r0.getZ() );
            initial_enx = gaussian()*(radius_energy0.getX() ) +(energy0.getX() );
            initial_eny = gaussian()*(radius_energy0.getY() ) +(energy0.getY() );
            initial_enz = gaussian()*(radius_energy0.getZ() ) +(energy0.getZ() );
            break;
        
        case electron.INITIALIZE_UNIFORM_DIST:
            initial_x = 0.0;
            initial_y = uniform_dist_single(num_par, r0.getY(), radius_r0.getY(), posy_counter);
            initial_z = uniform_dist_single(num_par, r0.getZ(), radius_r0.getZ(), posz_counter);
            initial_enx = uniform_dist_single(num_par, energy0.getX(), radius_energy0.getX(), velx_counter);
            initial_eny = uniform_dist_single(num_par, energy0.getY(), radius_energy0.getY(), vely_counter);
            initial_enz = uniform_dist_single(num_par, energy0.getZ(), radius_energy0.getZ(), velz_counter);
            break;
        }
        
        ThreeVec r_int(initial_x, initial_y, initial_z);
        ThreeVec energy_int(initial_enx, initial_eny, initial_enz);
        electron.set_pos(r_int);
        electron.set_energy(energy_int);
        electron.set_time(time);

        outfile_part_writeAndComma(electron);
                            /*
                            *    Determine the vn and xn at time step time = 0-del_time using 
                            *    backward difference method method: x^n = x^n+1 - del_time*(y^n+1)
                            *    where y is vxB and x is v.
                            */

        

        //step_through_magnet(electron, vn_plus, vn_minus, B0, rn_plus, rn_minus, time, del_time, width_bmap[magnet_counter], length_bmap[magnet_counter]);
        //step_through_magnet_mag_leap(electron, magnet[magnet_counter], vn_plus, vn_minus, rn_plus, rn_minus, time, del_time);
        step_through_magnet_mag_boris(electron,magnet[magnet_counter],time,del_time);
        
        while(++magnet_counter < num_magnets)
        {   //Loop through magnets after first magnet
            //magnet_counter++;

            double dist_x_to_mag = magnet[(magnet_counter)].get_pos(0) - magnet[magnet_counter-1].get_pos(0);

            if( ( (dist_x_to_mag > 0) && ((electron.get_vel(0)) > 0) ) || ( (dist_x_to_mag < 0) && ((electron.get_vel(0)) < 0) ) )
            { //First check that electron is moving toward next magnet
                double time_btwn_mags = 0.0;
                time_btwn_mags = (((magnet[magnet_counter].get_pos(0))) - (electron.get_pos(0)))/(electron.get_vel(0)); 
                double y_at_time{0.0};
                y_at_time = (electron.get_pos(1)) + ((electron.get_vel(1))*(time_btwn_mags));
                double z_at_time{0.0};
                z_at_time = (electron.get_pos(2)) + ((electron.get_vel(2))*(time_btwn_mags));
                time = time + time_btwn_mags;

                if( ((y_at_time >= (magnet[magnet_counter].get_pos(1) - (magnet[magnet_counter].get_width()/2.0) ) ) && (y_at_time <= (magnet[magnet_counter].get_pos(1) + (magnet[magnet_counter].get_width()/2.0) ))))
                {   //Then check that particle ends up in magnet region
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
                    //std::cerr << "After: " << electron.get_pos() << '\n' << '~' << '\n';
                    
                }
            }
        }   //<- end of magnetic while loop

        outfile_part_newline(electron);
    }   //<- end of particle for loop


    outfile_readme.close();
    //outfile_grid.close();
    outfile_time.close();
    outfile_xpos.close();
    outfile_ypos.close();
    outfile_zpos.close();
    outfile_px.close();
    outfile_py.close();
    outfile_pz.close();
    outfile_magnets.close();
    outfile_energy.close();
    //infile.close();

    return 0;
}
