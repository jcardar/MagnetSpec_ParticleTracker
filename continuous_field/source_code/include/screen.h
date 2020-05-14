#ifndef SCREEN_H
#define SCREEN_H

#include "threevector.h"

class Screen
{
private:
    ThreeVec m_position;    //POSITION OF LOW-ENERGY EDGE, CENTER OF HEIGHT, IN MAGNET COORDINATE SPACE
    double m_angle;     //degrees from horizontal
    double m_length;
    double m_height;
public:

    Screen(ThreeVec position, double angle_deg, double length, double height);
    Screen(ThreeVec position, double length, double height);
    Screen() {};

    void set_pos(ThreeVec pos);

    void set_pos(int index, double value);

    void set_angle(double angle);

    void set_length(double length);

    void set_height(double height);

    ThreeVec get_pos();

    double get_pos(int index);

    double get_angle(char type = 'd');

    double get_length();

    double get_height();

    void set_outfile(std::ofstream& out_screen);

    std::ofstream *m_out_screen;
};


#endif