#import matplotlib as mp
#import matplotlib.pyplot as plt
#import matplotlib.patches as patches
#import mpl_toolkits.mplot3d.art3d as art3d
#import pandas as pd
#import csv
#import numpy as np
#import itertools
#import math
#import scipy.constants as const
"""
with open("../XPOS.csv") as csvfile:
    csv_reader = csv.reader(csvfile, delimiter=',')
    line_count = 0
    posx = list(csv_reader)
    posx = [[float(y) for y in x] for x in posx]

with open("../YPOS.csv") as csvfile:
    csv_reader = csv.reader(csvfile, delimiter=',',quoting=csv.QUOTE_NONNUMERIC)
    line_count = 0
    posy = list(csv_reader)
    posy = [[float(y) for y in x] for x in posy]

with open("../ZPOS.csv") as csvfile:
    csv_reader = csv.reader(csvfile, delimiter=',')
    line_count = 0
    posz = list(csv_reader)
#     posz = [[float(y) for y in x[0:-1]] for x in posz]
    posz = [[float(y) for y in x] for x in posz]

with open("../MOMENTUM_X.csv") as csvfile:
    csv_reader = csv.reader(csvfile, delimiter=',')
    line_count = 0
    px = list(csv_reader)
    px = [[float(y) for y in x] for x in px]

with open("../MOMENTUM_Y.csv") as csvfile:
    csv_reader = csv.reader(csvfile, delimiter=',');
    line_count = 0;
    py = list(csv_reader);
    py = [[float(y) for y in x] for x in py]

with open("../MOMENTUM_Z.csv") as csvfile:
    csv_reader = csv.reader(csvfile, delimiter=',');
    line_count = 0;
    pz = list(csv_reader);
    pz = [[float(y) for y in x] for x in pz]

with open("../TIME.csv") as csvfile:
    csv_reader = csv.reader(csvfile, delimiter=',');
    line_count = 0;
    time = list(csv_reader);
    time = [[float(y) for y in x] for x in time]
#     time = [[float(y) for y in x] for x in time]

magnet = pd.read_csv("../MAGNETS.csv")

screen = pd.read_csv("../SCREENS.csv")

del_time = pd.read_csv("../DEL_T.csv", dtype=float, header = -1)
"""
import csv
import sys
import numpy as np

def import_relevant_data():
    
    with open("../data/ENERGY.csv") as csvfile:
        csv_reader = csv.reader(csvfile, delimiter=',')
        line_count = 0
        energy = list(csv_reader)
        energy = [[float(y) for y in x] for x in energy]
    energy = [[energy[jj][ii]*0.511+0.511 for ii in range(len(energy[jj]))] for jj in range(len(energy))]

    with open("../data/PARTICLE_ON_SCREENS.csv") as csvfile:
        csv_reader = csv.reader(csvfile, delimiter=',')
        line_count = 0
        next(csv_reader)
        part_on_screen = list(csv_reader)
        part_on_screen = [[float(y) for y in x] for x in part_on_screen]
        part_on_screen_screen_index = [int(x[0]) for x in part_on_screen]
        part_on_screen_part_index = [int(x[1]) for x in part_on_screen]
        part_on_screen_pos_x = [x[2] for x in part_on_screen]
        part_on_screen_pos_y = [x[3] for x in part_on_screen]
        part_on_screen_pos_z = [x[4] for x in part_on_screen]

    with open("../data/XPOS.csv") as csvfile:
        csv_reader = csv.reader(csvfile, delimiter=',')
        line_count = 0
        posx = list(csv_reader)
        posx = [[float(y) for y in x] for x in posx]

    return energy, part_on_screen_screen_index, part_on_screen_part_index, part_on_screen_pos_x, part_on_screen_pos_y, part_on_screen_pos_z, posx

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def check_energy_range_captured(part_on_screen_part_index, energy, energy_range, posx):
    """
    Can work with total energy range of energy_range = [min_energy, max_energy], where 90% of particles in the system should be captured, or with
    subgroups of energy ranges energy_range = [[min_energy0, max_energy0],[min_energy1, max_energy1],[min_energy2, max_energy2],...]
    Subgroups meant to ensure no holes in energy at desired ranges.
    """
    if energy_range.shape == (2,):
        if (min([energy[part_on_screen_part_index[ii]][0] for ii in range(len(part_on_screen_part_index))]))>=energy_range[0] and (max([energy[part_on_screen_part_index[ii]][0] for ii in range(len(part_on_screen_part_index))]))<=energy_range[1]:
            captured = True
        else:
            captured = False
            return captured
        if len(part_on_screen_part_index) <= (0.90*len(posx)):
            captured = True
        else:
            captured = False
            return captured
    else:
        if (min([energy[part_on_screen_part_index[ii]][0] for ii in range(len(part_on_screen_part_index))]))>=min(map(min,energy_range)) and (max([energy[part_on_screen_part_index[ii]][0] for ii in range(len(part_on_screen_part_index))]))<=max(map(max,energy_range)):
            captured = True
        else:
            return False
        number_of_ranges = energy_range.shape[0]
        for jj in range(number_of_ranges):
            #energy ranges between the ranges supplied have no holes in them
            nearest_to_min_range = find_nearest([energy[part_on_screen_part_index[ii]][0] for ii in range(len(part_on_screen_part_index))], energy_range[jj][0])
            nearest_to_max_range = find_nearest([energy[part_on_screen_part_index[ii]][0] for ii in range(len(part_on_screen_part_index))], energy_range[jj][1])
                
    return captured



def energy_resolution(energy_range,normalizing_fom, isfirst):
    energy, part_on_screen_screen_index, part_on_screen_part_index, part_on_screen_pos_x, part_on_screen_pos_y, part_on_screen_pos_z, posx = import_relevant_data()
    energy_range_captured = check_energy_range_captured(part_on_screen_part_index, energy, energy_range, posx)
    number_of_screens = int(max(part_on_screen_screen_index)+1)
    if energy_range_captured == False and isfirst == False:
        energy_resolution = sys.float_info.max
        return energy_resolution
    elif isfirst == True:
        energy_resolution_sum = np.array([])
        for jj in range(number_of_screens):
            long_cord = [np.sqrt(part_on_screen_pos_x[ii]**2 + part_on_screen_pos_y[ii]**2) for ii in range(len(part_on_screen_pos_x)) if part_on_screen_screen_index[ii] == jj];
            indicies_kept = [ii for ii in range(len(part_on_screen_pos_x)) if part_on_screen_screen_index[ii] == jj]
            indicies_of_div_7_0 = [ii for ii in range(len(indicies_kept)) if part_on_screen_part_index[indicies_kept[ii]]%7==0]
            indicies_of_div_7_1 = [ii for ii in range(len(indicies_kept)) if part_on_screen_part_index[indicies_kept[ii]]%7==1]
            indicies_of_div_7_2 = [ii for ii in range(len(indicies_kept)) if part_on_screen_part_index[indicies_kept[ii]]%7==2]
            energies = [energy[int(part_on_screen_part_index[int(indicies_kept[ii])])][0] for ii in range(len(long_cord))]
            dE_dx    = [abs(energies[indicies_of_div_7_0[ii+1]]-energies[indicies_of_div_7_0[ii-1]])/abs(long_cord[indicies_of_div_7_0[ii+1]]-long_cord[indicies_of_div_7_0[ii-1]]) for ii in (range(1, int(len(indicies_of_div_7_0))-1, 1))]
            dE_dx_times_E = np.array(dE_dx)*np.array([energies[indicies_of_div_7_0[ii]] for ii in range(1, len(indicies_of_div_7_0)-1, 1)])
            np.append(energy_resolution_sum, dE_dx_times_E)
        energy_res_fom = np.sum(energy_resolution_sum)
        return energy_res_fom
    elif energy_range_captured == True and isfirst == False:
        energy_resolution_sum = np.array([])
        for jj in range(number_of_screens):
            long_cord = [np.sqrt(part_on_screen_pos_x[ii]**2 + part_on_screen_pos_y[ii]**2) for ii in range(len(part_on_screen_pos_x)) if part_on_screen_screen_index[ii] == jj];
            indicies_kept = [ii for ii in range(len(part_on_screen_pos_x)) if part_on_screen_screen_index[ii] == jj]
            indicies_of_div_7_0 = [ii for ii in range(len(indicies_kept)) if part_on_screen_part_index[indicies_kept[ii]]%7==0]
            indicies_of_div_7_1 = [ii for ii in range(len(indicies_kept)) if part_on_screen_part_index[indicies_kept[ii]]%7==1]
            indicies_of_div_7_2 = [ii for ii in range(len(indicies_kept)) if part_on_screen_part_index[indicies_kept[ii]]%7==2]
            energies = [energy[int(part_on_screen_part_index[int(indicies_kept[ii])])][0] for ii in range(len(long_cord))]
            dE_dx    = [abs(energies[indicies_of_div_7_0[ii+1]]-energies[indicies_of_div_7_0[ii-1]])/abs(long_cord[indicies_of_div_7_0[ii+1]]-long_cord[indicies_of_div_7_0[ii-1]]) for ii in (range(1, int(len(indicies_of_div_7_0))-1, 1))]
            dE_dx_times_E = np.array(dE_dx)*np.array([energies[indicies_of_div_7_0[ii]] for ii in range(1, len(indicies_of_div_7_0)-1, 1)])
            np.append(energy_resolution_sum, dE_dx_times_E)
        energy_res_fom = np.sum(energy_resolution_sum)
        energy_res_fom = energy_res_fom/normalizing_fom
        return energy_res_fom
    else:
        print("Something went wrong calculating energy resolution fom!")
