#ifndef MAGNET_H
#define MAGNET_H

#include "threevector.h"
#include <fstream>

class Magnet
{
private:
    ThreeVec m_position;                //position of y-center, left-hand-side of magnet
    ThreeVec m_Bfield;
    double* m_Bfield_grid;
    int m_grid_num;
    double m_width;
    double m_length;
    double m_height;

public:
    Magnet() {};

    Magnet(ThreeVec pos, double length, double width, double height, ThreeVec B0, std::ofstream& out_magnet);     //uniform B

    Magnet(ThreeVec pos, double length, double width, double height, double* Bmap, std::ofstream& out_magnet);    //discrete B

    void set_pos(ThreeVec pos);

    void set_pos(int index, double value);

    void set_B0(ThreeVec B0);

    void set_B0(int index, double value);

    void set_length(double length);

    void set_width(double width);

    void set_height(double height);

    ThreeVec get_pos()
    {
        return m_position;
    }

    double get_pos(int index);

    ThreeVec get_B0();

    double get_B0(int index);

    double get_length();

    double get_width();

    double get_height();

    void set_outfile(std::ofstream& out_mag);

    std::ofstream *m_out_magnet;
};



#endif