#ifndef SCREEN_H
#define SCREEN_H

#include "threevector.h"
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

class Screen
{
private:
    ThreeVec m_position;    //POSITION OF LOW-ENERGY EDGE, CENTER OF HEIGHT, IN MAGNET COORDINATE SPACE
    double m_angle_about_x;     //radians about x-axis
    double m_angle_about_y;     //radians about z-axis
    double m_angle_about_z;     //radians about z-axis
    double m_length;
    double m_height;
public:

    Screen(ThreeVec position, double angle_rad_around_x, double angle_rad_around_y, double angle_rad_around_z, double length, double height);
    ///// If using this instantiation, use radians for angles.
    Screen(ThreeVec position, double length, double height);
    Screen() {};

    void set_pos(ThreeVec pos);

    void set_pos(int index, double value);

    void set_pos(double posx, double posy, double posz);

    void set_angle_about_x(double angle, char unit = 'r');

    void set_angle_about_y(double angle, char unit = 'r');

    void set_angle_about_z(double angle, char unit = 'r');

    void set_length(double length);

    void set_height(double height);

    ThreeVec get_pos();

    double get_pos(int index);

    double get_angle_about_x(char type = 'r');

    double get_angle_about_y(char type = 'r');

    double get_angle_about_z(char type = 'r');

    double get_length();

    double get_height();

    void set_outfile(std::ofstream& out_screen);

    std::ofstream *m_out_screen;
};


#endif