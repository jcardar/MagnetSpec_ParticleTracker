/*Leap Frog Numerical Method Code
 *for tracking a particle in static B and possible E field
 *Author: Jason Cardarelli
 *NERS 574: Computational Plasma Physics
 *Prof. Alexander Thomas
 */
#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include "threevector.h"
#include "threematrix.h"

using namespace std;

int main()
{

//ifstream file( "yourfilename.txt" )    //The Format to import data from data file
ifstream bmap( "bmap_in.txt" );   //7x21 x vs y matrix
ofstream outfile;
outfile.open("data.txt");
ThreeVec B0(0.0,0.0,1.0);
//ThreeVec E0(0.0,0.1,0.0);      //UNCOMMENT IF E-FIELD IS DESIRED
ThreeVec vn(1.0,0.0,0.0);
ThreeVec rn(0.0,0.0,0.0);


//for(int b_xn{0};b_xn<7;b_xn++)
//{
//    for(int b_yn{0};b_yn<21;b_yn++)
//    {
//
//    }
//}


int length_bmap{0};
int height_bmap{0};
int ii_byte{0};

while(bmap)
{
    char current;
    bmap.seekg(ii_byte);
    bmap >> current;
    if (current == 32)
        {length_bmap++;}
    else if(current == 92)
        {
            ii_byte++;
            bmap.seekg(ii_byte);
            bmap>> current;
            if(current == 110)
            {
                height_bmap++;
                length_bmap = 0;
            }
        }
}


while(bmap)
{

}


double time;                          //Define a variable time in which will be stepped over
double del_t = 0.05;                  //Define a time step
const double t_max = (4*(3.141592));  //Define a maximum time to run to


outfile << "t\t" << "vx\t" << "vy\t" << "vz\t" << "x\t" << "y\t" << "z\t" << endl;


/*
 *    Determine the vn and xn at time step time = 0-del_t using backward difference method
 *    method: x^n = x^n+1 - del_t*(y^n+1) where y is vxB and x is v.
 */

ThreeVec vn_minus = vn - ((vn^B0)*del_t);
ThreeVec vn_plus;

ThreeVec rn_minus = rn - (vn*del_t);
ThreeVec rn_plus;

for (time = 0; time <= t_max; time+=del_t)
    {
        outfile << time << "\t" << vn.get(0) << "\t" << vn.get(1) << "\t" << vn.get(2) << "\t"
        << rn.get(0) << "\t" << rn.get(1) << "\t" << rn.get(2) << endl;

        vn_plus = vn_minus + (vn^B0)*2*del_t;    //leapfrog algarithm vn

        rn_plus = rn_minus + vn*2*del_t;         //leapfrog algarithm rn


        rn_minus = rn;                          //iterate xn_minus up a timestep
        rn = rn_plus;                           //do the same to xn

        vn_minus = vn;                          //iterate vn_minus up a timestep
        vn = vn_plus;                           //do the same to vn
    }



outfile.close();
return 0;
}
