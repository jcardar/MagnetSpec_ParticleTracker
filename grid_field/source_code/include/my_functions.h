#ifndef MY_FUNCTIONS_H
#define MY_FUNCTIONS_H

#include <fstream>
#include <vector>
#include "threevector.h"
#include "threematrix.h"
#include "particle.h"
#include "magnet.h"
#include "screen.h"
#include "math.h"
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

void outfile_tab(double& time, std::ofstream& out_time, ThreeVec pos, std::ofstream& out_xpos, std::ofstream& out_ypos, std::ofstream& out_zpos, ThreeVec vel, std::ofstream& out_vx, std::ofstream& out_vy, std::ofstream& out_vz);

void outfile_part_writeAndTab(Particle& particle);

void outfile_part_writeAndComma(Particle& particle);

void outfile_part_write(Particle& particle);

void outfile_newline(std::ofstream& out_time, std::ofstream& out_xpos, std::ofstream& out_ypos, std::ofstream& out_zpos, std::ofstream& out_vx, std::ofstream& out_vy, std::ofstream& out_vz);

void outfile_part_newline(Particle& particle);

void outfile_part_comma(Particle& particle);

void outfile_uniform_magnet(Magnet& magnet, int counter);

void outfile_screen_single(Screen& screen, int counter);

void uniform_en_dist(double& initial_x, double& initial_y, double& initial_z, double& initial_enx, double& initial_eny, double& initial_enz, double length_before, int* posy_counter, int* posz_counter, ThreeVec r0, ThreeVec radius_r0, ThreeVec energy0, int num_par);

void uniform_pos_dist(double& initial_x, double& initial_y, double& initial_z, double& initial_enx, double& initial_eny, double& initial_enz, double length_before, int* posx_counter, int* posy_counter, int* posz_counter, ThreeVec p0, ThreeVec radius_p0, ThreeVec r0, int num_par, bool pointSource, bool dist_x, bool dist_y, bool dist_z);

double gaussian();

double uniform_dist_single(const int num_par_t, double vel0_t, double radius_v0_t, int &counter_t);

void step_through_magnet_leap(Particle &electron, ThreeVec vn_plus, ThreeVec vn_minus, ThreeVec B0, ThreeVec rn_plus, ThreeVec rn_minus, double& time, const double& del_time, double& width_bmap, double& length_bmap);

void step_through_magnet_mag_leap(Particle &electron, Magnet &magnet, ThreeVec vn_plus, ThreeVec vn_minus, ThreeVec rn_plus, ThreeVec rn_minus, double& time, const double& del_time);

void step_through_magnet_mag_leap(Particle *electron, Magnet &magnet, double& time, const double& del_time);

void boris(Particle &electron_t, Magnet &magnet_t, const double del_t, double mu_0);

void step_through_magnet_mag_boris(Particle &electron, Magnet &magnet, double& time, const double &del_time, double mu_0, double time_out = (2*M_PI*1000));

void first_half_position_step(Particle &electron_t, const double del_t);

bool inside_of_mag(Magnet magnet_t, Particle particle_t);

double time_to_magnet_boundary(Magnet magnet_t, Particle particle_t);

bool intersect_mag(Magnet magnet_t, Particle particle_t);

double dist_to_mag(Magnet magnet_t, Particle particle_t);

void move_particle_to_magnet(Magnet magnet_t, Particle &particle_t);

void move_through_magnets(Magnet magnet_t[], int num_mags, Particle &particle_t, double &time, double del_time, double mu_0, double time_limit);

void half_time_step(double &time_step);

void move_to_screens(Screen screen_t[], int num_screen, Particle particle_t);

void readUnits(std::ifstream &input_stream, std::vector<std::string> &desired_units);

int readNumOf(std::ifstream &input_stream);

void readMagnet(std::ifstream &input_stream, int &magNum, std::vector<std::vector<std::vector<double>>> &magInfo);

void readPermanentMagDim(std::ifstream &input_stream, int magNum, std::vector<double> &PmagDim);

void readMagAxis(std::ifstream &input_stream, int magNum, std::vector<char> &axisInfo);

void readBeam(std::ifstream &input_stream, std::vector<std::vector<double>> &beamInfo);

void readSpread(std::ifstream &input_stream, std::vector<std::vector<double>> &spreadInfo);

void readScreen(std::ifstream &input_stream, int &screenNum, std::vector<std::vector<std::vector<double>>> &screenInfo);

void ReadInitTypes(std::ifstream &input_stream, std::vector<int> &init_types);

double ReadMu0(std::ifstream &input_stream);

double find_magnetization(Magnet &magnet, double mag_dim, double mu_0, char axis);

void calc_grid_B_comps(double factor, double a, double b, double c, double x, double y, double z, double &temp_B1, double &temp_B2, double &temp_B3);

bool B_within_margin(double B_center_val, double B1, double B2, double B3);

ThreeVec calc_grid_point_B(ThreeVec &grid_point, Magnet &magnet, double mag_dim, double mu_0, char axis, double magnetization);

#endif