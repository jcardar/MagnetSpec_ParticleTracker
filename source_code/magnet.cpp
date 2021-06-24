#include "magnet.h"
#include <iostream>
#include "threevector.h"



Magnet::Magnet(ThreeVec pos, double length, double width, double height, double height_of_dipole_block, char type, double remanence, ThreeVec B0, std::ofstream& out_magnet)
        : m_position(pos), m_length(length), m_width(width), m_height(height), m_height_of_dipole_block(height_of_dipole_block), m_type(type), m_remanence(remanence), m_Bfield(B0), m_out_magnet(&out_magnet)
{

}

Magnet::Magnet(ThreeVec pos, double length, double width, double height, double height_of_dipole_block, char type, double remanence, double* Bmap, std::ofstream& out_magnet)
    : m_position(pos), m_length(length), m_width(width), m_height(height), m_height_of_dipole_block(height_of_dipole_block), m_type(type), m_remanence(remanence), m_Bfield_grid(Bmap), m_out_magnet(&out_magnet)
{

}

void Magnet::set_type(char type)
{
    m_type = type;
}

void Magnet::set_pos(ThreeVec pos)
{
    m_position = pos;
}

void Magnet::set_pos(int index, double value)
{
    m_position.set(index, value);
}

void Magnet::set_length(double length)
{
    m_length = length;
}

void Magnet::set_width(double width)
{
    m_width = width;
}

void Magnet::set_height(double height)
{
    m_height = height;
}

void Magnet::set_height_of_dipole_block(double height_of_block)
{
    m_height_of_dipole_block = height_of_block;
}

void Magnet::set_axis_of_magnetization(char axis_of_magnetization)
{
    m_axis_of_magnetization = axis_of_magnetization;
}

void Magnet::set_B0(ThreeVec B0)
{
    m_Bfield = B0;
}

void Magnet::set_B0(int index, double value)
{
    m_Bfield.set(index, value);
}

void Magnet::set_outfile(std::ofstream& out_mag)
{
    m_out_magnet = &out_mag;
}

double Magnet::get_pos(int index)
{
    return m_position.get(index);
}

ThreeVec Magnet::get_B0()
{
    return m_Bfield;
}

double Magnet::get_B0(int index)
{
    return m_Bfield.get(index);
}

double Magnet::get_length()
{
    return m_length;
}

double Magnet::get_width()
{
    return m_width;
}

double Magnet::get_height()
{
    return m_height;
}

double Magnet::get_height_of_dipole_block()
{
    return m_height_of_dipole_block;
}

char Magnet::get_axis_of_magnetization()
{
    return m_axis_of_magnetization;
}

char Magnet::get_type()
{
    return m_type;
}

double Magnet::get_Br()
{
    return m_remanence;
}