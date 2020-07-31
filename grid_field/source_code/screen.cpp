#include <iostream>
#include "screen.h"
#include "threevector.h"
#include <cmath>
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

Screen::Screen(ThreeVec position, double angle_rad_about_x, double angle_rad_about_y, double angle_rad_about_z, double length, double height)
        : m_position(position), m_angle_about_x(angle_rad_about_x), m_angle_about_y(angle_rad_about_y), m_angle_about_z(angle_rad_about_z), m_length(length), m_height(height)
    {
        if(angle_rad_about_x > M_PI)
        {
            double multiple = angle_rad_about_x/M_PI;
            m_angle_about_x = angle_rad_about_x - ( static_cast<int>(multiple)*M_PI );
        }
        if(angle_rad_about_y > M_PI)
        {
            double multiple = angle_rad_about_y/M_PI;
            m_angle_about_y = angle_rad_about_y - ( static_cast<int>(multiple)*M_PI );
        }
        if(angle_rad_about_z > M_PI)
        {
            double multiple = angle_rad_about_z/M_PI;
            m_angle_about_z = angle_rad_about_z - ( static_cast<int>(multiple)*M_PI );
        }
    }

Screen::Screen(ThreeVec position, double length, double height)
        : m_position(position), m_angle_about_x(0.0), m_angle_about_y(0.0), m_angle_about_z(0.0), m_length(length), m_height(height)
    {
    }

void Screen::set_pos(ThreeVec pos)
    {
        m_position = pos;
    }

void Screen::set_pos(int index, double value)
{
    m_position.set(index, value);
}

void Screen::set_pos(double posx, double posy, double posz)
{
    m_position.set(0, posx);
    m_position.set(1, posy);
    m_position.set(2, posz);
}

void Screen::set_angle_about_x(double angle_deg, char unit)
{
    if(unit == 'd')
    {
        if(angle_deg > 180)
        {
            double multiple  = angle_deg/180;
            m_angle_about_x  = angle_deg - ( static_cast<int>(multiple)*180 );
            m_angle_about_x  = m_angle_about_x * (M_PI / 180.0);
        }
        else
        {
            m_angle_about_x  = angle_deg * (M_PI / 180.0);
        }   
    } 
    else //Assume radians if not 'd'
    {
        if(angle_deg > M_PI)
        {
            double multiple  = angle_deg/M_PI;
            m_angle_about_x  = angle_deg - ( static_cast<int>(multiple)*M_PI );
        }
        else
        {
            m_angle_about_x = angle_deg;
        } 
    }   
}

void Screen::set_angle_about_y(double angle, char unit)
{
    if(unit == 'd')
    {
        if(angle > 180)
        {
            double multiple = angle/180;
            m_angle_about_y = angle - ( static_cast<int>(multiple)*180 );
            m_angle_about_y = m_angle_about_y * (M_PI / 180.0);
        }
        else
        {
            m_angle_about_y = angle * (M_PI / 180.0);
        }   
    } 
    else //Assume radians if not 'd'
    {
        if(angle > M_PI)
        {
            double multiple  = angle/M_PI;
            m_angle_about_y  = angle - ( static_cast<int>(multiple)*M_PI );
        }
        else
        {
            m_angle_about_y = angle;
        } 
    }  
}

void Screen::set_angle_about_z(double angle_deg, char unit)
    {
        if(unit == 'd')
        {
            if(angle_deg > 180)
            {
                double multiple = angle_deg / 180;
                m_angle_about_z = angle_deg - ( static_cast<int>(multiple)*180 );
                m_angle_about_z = m_angle_about_z * (M_PI / 180.0);
            }
            else
            {
                m_angle_about_z = angle_deg * (M_PI / 180.0);
            }   
        } 
        else //Assume radians if not 'd'
        {
            if(angle_deg > M_PI)
            {
                double multiple = angle_deg/M_PI;
                m_angle_about_z = angle_deg - ( static_cast<int>(multiple)*M_PI );
            }
            else
            {
                m_angle_about_z = angle_deg;
            } 
        }
    }

void Screen::set_length(double length)
    {
        m_length = length;
    }

void Screen::set_height(double height)
    {
        m_height = height;
    }

ThreeVec Screen::get_pos()
    {
        return m_position;
    }

double Screen::get_pos(int index)
    {
        return m_position.get(index);
    }

double Screen::get_angle_about_x(char type)
{
    if(type == 'd')
    {
        return m_angle_about_x*(180.0/M_PI);
    }
    else if(type == 'r')
    {
        return m_angle_about_x;
    }
    else
    {
        std::cout << "Incorrect Screen angle_x type in get_angle_about_x(); radians used by default";
        return m_angle_about_x;
    }
    
}

double Screen::get_angle_about_y(char type)
{
    if(type == 'd')
    {
        return m_angle_about_y*(180.0/M_PI);
    }
    else if(type == 'r')
    {
        return m_angle_about_y;
    }
    else
    {
        std::cout << "Incorrect Screen angle_y type in get_angle_about_y(); radians used by default";
        return m_angle_about_y;
    }
    
}

double Screen::get_angle_about_z(char type)
{
    if(type == 'd')
    {
        return m_angle_about_z*(180.0/M_PI);
    }
    else if(type == 'r')
    {
        return m_angle_about_z;
    }
    else
    {
        std::cout << "Incorrect Screen angle_z type in get_angle_about_z(); radians used by default";
        return m_angle_about_z;
    }
    
}

double Screen::get_length()
    {
        return m_length;
    }

double Screen::get_height()
    {
        return m_height;
    }

void Screen::set_outfile(std::ofstream& out_screen)
    {
        m_out_screen = &out_screen;
    }