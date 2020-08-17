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

int main() {
    std::ofstream outfile_B("grid_maker_B_data.txt");
    std::ofstream outfile_Z("grid_maker_Z_coord.txt");
    std::ofstream outfile_Y("grid_maker_Y_coord.txt");

    Magnet magnet;

    // hard code input values for magnet
    magnet.set_pos(0, 0.0); 
    magnet.set_pos(1, 0.0); 
    magnet.set_pos(2, 0.0);

    magnet.set_length(176.0037714340999);
    magnet.set_width(41.06754666795665);
    magnet.set_height(17.60037714340999);

    magnet.set_B0(0, 0.0);
    magnet.set_B0(1, 0.0);
    magnet.set_B0(2, 1.0);

    double mag_dim = 46.9343390490933;

    double mu_0 = 2.0775073042165752e-11;

    char axis = 'z';

    double magnetization = find_magnetization(magnet, mag_dim, mu_0, axis);
    
    ThreeVec grid_point;
    grid_point.setX( magnet.get_pos(0) + 0.5*magnet.get_length() );  // x is middle of space

    int step_num = 1000;
    double z_step = magnet.get_height()/step_num;
    double y_step = magnet.get_width()/step_num;

    double z_coord = 0.0;  // start at 0 distance along width and height
    double y_coord = 0.0;  // i.e. bottom right corner of yz plane between magnets when using the diagram from input deck
    
    for(int i=0; i<step_num; ++i) {
        double offset_z = magnet.get_pos(2) - 0.5*magnet.get_height() + z_coord;
        grid_point.setZ(offset_z);

        for(int j=0; j<step_num; ++j) {
            double offset_y = magnet.get_pos(1) - 0.5*magnet.get_width() + y_coord;
            grid_point.setY(offset_y);

            ThreeVec grid_B_vec = calc_grid_point_B(grid_point, magnet, mag_dim, mu_0, axis, magnetization);
            double B_mag = grid_B_vec.mag();

            outfile_B << B_mag << ' ';
            outfile_Z << z_coord << ' ';
            outfile_Y << y_coord << ' ';
            y_coord += y_step;
        }
        outfile_B << '\n';
        outfile_Z << '\n';
        outfile_Y << '\n';

        z_coord += z_step;
        y_coord = 0.0;
    }
    outfile_B.close();
    outfile_Z.close();
    outfile_Y.close();
}