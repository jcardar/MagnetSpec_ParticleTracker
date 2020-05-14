#include <iostream>
#include "screen.h"
#include "threevector.h"
#include <cmath>

Screen::Screen(ThreeVec position, double angle_deg, double length, double height)
        : m_position(position), m_angle(angle_deg), m_length(length), m_height(height)
    {
        if(angle_deg > 180)
        {
            double multiple = angle_deg/180;
            m_angle         = angle_deg - ( static_cast<int>(multiple)*180 );
        }
    }

Screen::Screen(ThreeVec position, double length, double height)
        : m_position(position), m_angle(0.0), m_length(length), m_height(height)
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

void Screen::set_angle(double angle_deg)
    {
        if(angle_deg > 180)
        {
            double multiple = angle_deg/180;
            m_angle         = angle_deg - ( static_cast<int>(multiple)*180 );
        }
        else
        {
            m_angle = angle_deg;
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

double Screen::get_angle(char type)
    {
        if(type == 'd')
        {
            return m_angle;
        }
        else if(type == 'r')
        {
            return m_angle*(M_PI/180);
        }
        else
        {
            std::cout << "Incorrect Screen angle type; degrees used by default";
            return m_angle;
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