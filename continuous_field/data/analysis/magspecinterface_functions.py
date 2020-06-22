import ipywidgets as widgets
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Define units

units_length = widgets.Dropdown(
    options=['mm', 'cm', 'm'],
    value='cm',
    description='Length Unit:',
    disabled=False,
)

units_energy = widgets.Dropdown(
    options=['eV', 'MeV', 'GeV'],
    value='MeV',
    description='Energy Unit:',
    disabled=False,
)

units_angles = widgets.Dropdown(
    options=['mrad', 'Radians', 'Degrees'],
    value='mrad',
    description='Angle Unit:',
    disabled=False,
)

units_magnetic_field = widgets.Dropdown(
    options=['T', 'Gauss'],
    value='T',
    description='Field Unit:',
    disabled=False,
)

# Define Coordinate System attributes

global_max_x = widgets.FloatText(
    value=100,
    description='global x max',
    disabled=False
)

global_min_x = widgets.FloatText(
    value=0,
    description='global x min',
    disabled=False
)

global_max_y = widgets.FloatText(
    value=100,
    description='global y max',
    disabled=False
)

global_min_y = widgets.FloatText(
    value=0,
    description='global y min',
    disabled=False
)

global_max_z = widgets.FloatText(
    value=100,
    description='global z max',
    disabled=False
)

global_min_z = widgets.FloatText(
    value=0,
    description='global z min',
    disabled=False
)

global_bounds = [global_max_x, global_min_x, global_max_y, global_min_y, global_max_z, global_min_z]

# Define Magnet attributes

number_of_magnets = widgets.BoundedIntText(
    value=2,
    min=1,
    step=1,
    description='# of Magnets',
    disabled=False
)

def dynamicFloatValue_Magnet_Dimensions(num_of_magnets):
    listOfWidgets = []
    for i in range(num_of_magnets):
        widget1 = widgets.BoundedFloatText(
            value=0,
            min=0,
            description=f'width {i+1}',
        )
        widget2 = widgets.BoundedFloatText(
            value=0,
            min=0,
            description=f'length {i+1}',
        )
        widget3 = widgets.BoundedFloatText(
            value=0,
            min=0,
            description=f'height {i+1}',
        )
        listOfWidgets.append(widget1)
        listOfWidgets.append(widget2)
        listOfWidgets.append(widget3)
    return listOfWidgets

def dynamicFloatValue_Magnet_Position(num_of_magnets, magnet_dimensions, global_bounds):
    listOfWidgets = []
    j = 0
    for i in range(num_of_magnets):
        widget1 = widgets.BoundedFloatText(
            value=0,
            max=(global_bounds[0].value - magnet_dimensions[j+1].value),
            min=(global_bounds[1].value ),
            description=f'x pos {i+1}',
        )
        widget2 = widgets.BoundedFloatText(
            value=0,
            max=(global_bounds[2].value - 0.5*magnet_dimensions[j].value),
            min=(global_bounds[3].value + 0.5*magnet_dimensions[j].value),
            description=f'y pos {i+1}',
        )
        widget3 = widgets.BoundedFloatText(
            value=0,
            max=(global_bounds[4].value - 0.5*magnet_dimensions[j+2].value),
            min=(global_bounds[5].value + 0.5*magnet_dimensions[j+2].value),
            description=f'z pos {i+1}',
        )
        listOfWidgets.append(widget1)
        listOfWidgets.append(widget2)
        listOfWidgets.append(widget3)
        j += 3
    return listOfWidgets

def dynamicFloatValue_Magnetic_Field_Comps(num_of_magnets):
    listOfWidgets = []
    for i in range(num_of_magnets):
        widget1 = widgets.FloatText(
            value=0,
            description=f'x comp {i+1}',
        )
        widget2 = widgets.FloatText(
            value=0,
            description=f'y comp {i+1}',
        )
        widget3 = widgets.FloatText(
            value=0,
            description=f'z comp {i+1}',
        )
        listOfWidgets.append(widget1)
        listOfWidgets.append(widget2)
        listOfWidgets.append(widget3)
    return listOfWidgets

# Define Beam attributes

number_of_particles = widgets.BoundedIntText(
    value=1,
    min=1,
    step=1,
    description='# of Particles',
    disabled=False
)

def dynamicFloatValue_Beam_Start_Position(global_bounds):
    listOfWidgets = []
    widget1 = widgets.BoundedFloatText(
        value=0,
        max=global_bounds[0].value,
        min=global_bounds[1].value,
        description='x start pos',
    )
    widget2 = widgets.BoundedFloatText(
        value=0,
        max=global_bounds[2].value,
        min=global_bounds[3].value,
        description='y start pos',
    )
    widget3 = widgets.BoundedFloatText(
        value=0,
        max=global_bounds[4].value,
        min=global_bounds[5].value,
        description='z start pos',
    )
    listOfWidgets.extend([widget1, widget2, widget3])
    return listOfWidgets

beam_energy = widgets.BoundedFloatText(
    value=1,
    min=0,
    description='beam energy'
)

def dynamicFloatValue_Beam_Direction():
    listOfWidgets = []
    widget1 = widgets.FloatText(
        value=0,
        description='x angle',
    )
    widget2 = widgets.FloatText(
        value=0,
        description='y angle',
    )
    widget3 = widgets.FloatText(
        value=0,
        description='z angle',
    )
    listOfWidgets.extend([widget1, widget2, widget3])
    return listOfWidgets

# Define Beam Spread attributes

def dynamicFloatValue_Beam_Position_Spread():
    listOfWidgets = []
    widget1 = widgets.FloatText(
        value=0,
        description='x pos spread',
    )
    widget2 = widgets.FloatText(
        value=0,
        description='y pos spread',
    )
    widget3 = widgets.FloatText(
        value=0,
        description='z pos spread',
    )
    listOfWidgets.extend([widget1, widget2, widget3])
    return listOfWidgets

beam_energy_spread = widgets.BoundedFloatText(
    value=0,
    min=0,
    description='nrg spread'
)

def dynamicFloatValue_Beam_Divergence_Spread():
    listOfWidgets = []
    widget1 = widgets.FloatText(
        value=0,
        description='x divergence',
    )
    widget2 = widgets.FloatText(
        value=0,
        description='y divergence',
    )
    widget3 = widgets.FloatText(
        value=0,
        description='z divergence',
    )
    listOfWidgets.extend([widget1, widget2, widget3])
    return listOfWidgets

# Define Screen attriutes

number_of_screens = widgets.BoundedIntText(
    value=2,
    min=0,
    step=1,
    description='# of Screens',
    disabled=False
)

def dynamicFloatValue_Screen_Dimensions(num_of_screens):
    listOfWidgets = []
    for i in range(num_of_screens):
        widget1 = widgets.BoundedFloatText(
            value=0,
            min=0,
            description=f'length {i+1}',
        )
        widget2 = widgets.BoundedFloatText(
            value=0,
            min=0,
            description=f'height {i+1}',
        )
        listOfWidgets.append(widget1)
        listOfWidgets.append(widget2)
    return listOfWidgets

def dynamicFloatValue_Screen_Position(num_of_screens, global_bounds):
    listOfWidgets = []
    for i in range(num_of_screens):
        widget1 = widgets.BoundedFloatText(
            value=0,
            max=global_bounds[0].value,
            min=global_bounds[1].value,
            description=f'x pos {i+1}',
        )
        widget2 = widgets.BoundedFloatText(
            value=0,
            max=global_bounds[2].value,
            min=global_bounds[3].value,
            description=f'y pos {i+1}',
        )
        widget3 = widgets.BoundedFloatText(
            value=0,
            max=global_bounds[4].value,
            min=global_bounds[5].value,
            description=f'z pos {i+1}',
        )
        listOfWidgets.append(widget1)
        listOfWidgets.append(widget2)
        listOfWidgets.append(widget3)
    return listOfWidgets

def dynamicFloatValue_Screen_Angles(num_of_screens):
    listOfWidgets = []
    for i in range(num_of_screens):
        widget1 = widgets.FloatText(
            value=0,
            description=f'yaw angle {i+1}',
        )
        widget2 = widgets.FloatText(
            value=0,
            description=f'pitch angle {i+1}',
        )
        widget3 = widgets.FloatText(
            value=0,
            description=f'roll angle {i+1}',
        )
        listOfWidgets.append(widget1)
        listOfWidgets.append(widget2)
        listOfWidgets.append(widget3)
    return listOfWidgets

# Functions to aid in normalization/converting values

def averageB0(num_of_magnets, field_comps):
    all_fields_total = 0
    i = 0
    for j in range(num_of_magnets):
        one_field_total = 0
        
        for k in range(3):
            one_field_total += math.pow(field_comps[i].value,2)
            i += 1
        
        all_fields_total += math.sqrt(one_field_total)
    
    aveB_0 = all_fields_total / num_of_magnets
    
    return aveB_0

def convertAngles(units, beam_direction, divergence_spread, screen_angles):
    converted_beam_direction = []
    converted_divergence_spread = []
    converted_screen_angles = []

    angle_multiplier = 1 # Radians
    if units[4] == 'mrad':
        angle_multiplier = math.pow(10,-3)
    elif units[4] == 'Degrees':
        angle_multiplier = math.pi/180
    
    for i in range(3):
        converted_beam_direction.append( beam_direction[i].value * angle_multiplier )
    for i in range(3):
        converted_divergence_spread.append( divergence_spread[i].value * angle_multiplier )
    for i in range(len(screen_angles)):
        converted_screen_angles.append( screen_angles[i].value * angle_multiplier )
        
    return converted_beam_direction, converted_divergence_spread, converted_screen_angles

def normalizeValues(units, num_mag, mag_dim, mag_pos, fld_comps, beam_pos, beam_energy, scrn_dim, scrn_pos):
    norm_mag_dim = []
    norm_mag_pos = []
    norm_fld_comps = []
    norm_beam_pos = []
    norm_scrn_dim = []
    norm_scrn_pos = []
    
    aveB_0 = averageB0(num_mag, fld_comps)
    omega_div_c = (1.602177 * math.pow(10,-19) * aveB_0) / (9.109384 * math.pow(10,-31) * 2.997925 * math.pow(10,8))
    
    # distances
    distance_multiplier = 1 # m
    if units[0] == 'mm':
        distance_multiplier = math.pow(10,-3)
    elif units[0] == 'cm':
        distance_multiplier = math.pow(10,-2)
        
    for i in range(len(mag_dim)):
        norm_mag_dim.append( mag_dim[i].value * distance_multiplier * omega_div_c )
    for i in range(len(mag_pos)):
        norm_mag_pos.append( mag_pos[i].value * distance_multiplier * omega_div_c )
    for i in range(len(beam_pos)):
        norm_beam_pos.append( beam_pos[i].value * distance_multiplier * omega_div_c )
    for i in range(len(scrn_dim)):
        norm_scrn_dim.append( scrn_dim[i].value * distance_multiplier * omega_div_c )
    for i in range(len(scrn_pos)):
        norm_scrn_pos.append( scrn_pos[i].value * distance_multiplier * omega_div_c )
    
    # magnetic field
    magnetic_multiplier = 1 # Tesla
    if units[6] == 'Gauss':
        magnetic_multiplier = math.pow(10,-4)
    
    for i in range(len(fld_comps)):
        norm_fld_comps.append( (fld_comps[i].value * magnetic_multiplier) / aveB_0 )
    
    # energy
    rest_energy = 0.511 # MeV
    if units[2] == 'eV':
        rest_energy = 0.511 * math.pow(10,6)
    elif units[2] == 'GeV':
        rest_energy = 0.511 * math.pow(10,-3)
    
    norm_beam_energy = (rest_energy + beam_energy)/rest_energy
    
    return norm_mag_dim, norm_magn_pos, norm_fld_comps, norm_beam_pos, norm_beam_energy, norm_scrn_dim, norm_scrn_pos

# Functions to aid in plotting

def showDiagram(num_mag, mag_dim, mag_pos, beam_pos, beam_dir, num_scrn, scrn_dim, scrn_pos, scrn_angl):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.autoscale_view()
    
    # plot magnet(s)
    i = 0
    for k in range(num_mag):
        mx1 = [mag_pos[i].value, mag_pos[i].value + mag_dim[i+1].value]
        mx2 = [mag_pos[i].value, mag_pos[i].value]
        mx3 = [mag_pos[i].value + mag_dim[i+1].value, mag_pos[i].value + mag_dim[i+1].value]
        my1 = [mag_pos[i+1].value - (0.5 * mag_dim[i].value), mag_pos[i+1].value + (0.5 * mag_dim[i].value)]
        my2 = [mag_pos[i+1].value + (0.5 * mag_dim[i].value), mag_pos[i+1].value + (0.5 * mag_dim[i].value)]
        my3 = [mag_pos[i+1].value - (0.5 * mag_dim[i].value), mag_pos[i+1].value - (0.5 * mag_dim[i].value)]
        mz1 = [mag_pos[i+2].value - (0.5 * mag_dim[i+2].value), mag_pos[i+2].value + (0.5 * mag_dim[i+2].value)]
        mz2 = [mag_pos[i+2].value + (0.5 * mag_dim[i+2].value), mag_pos[i+2].value + (0.5 * mag_dim[i+2].value)]
        mz3 = [mag_pos[i+2].value - (0.5 * mag_dim[i+2].value), mag_pos[i+2].value - (0.5 * mag_dim[i+2].value)]
    
        ax.plot(mx2, my1, mz2, c='b', alpha=0.5)
        ax.plot(mx2, my1, mz3, c='b', alpha=0.5)
        ax.plot(mx2, my2, mz1, c='b', alpha=0.5)
        ax.plot(mx2, my3, mz1, c='b', alpha=0.5)
        ax.plot(mx3, my1, mz2, c='b', alpha=0.5)
        ax.plot(mx3, my1, mz3, c='b', alpha=0.5)
        ax.plot(mx3, my2, mz1, c='b', alpha=0.5)
        ax.plot(mx3, my3, mz1, c='b', alpha=0.5)
        ax.plot(mx1, my2, mz2, c='b', alpha=0.5)
        ax.plot(mx1, my3, mz2, c='b', alpha=0.5)
        ax.plot(mx1, my2, mz3, c='b', alpha=0.5)
        ax.plot(mx1, my3, mz3, c='b', alpha=0.5)
        i += 3
    
    # plot beam
    length = mag_pos[0].value - beam_pos[0].value
    startx = beam_pos[0].value
    starty = beam_pos[1].value
    startz = beam_pos[2].value
    endx = length * math.cos(beam_dir[0])   # using direction cosines
    endy = length * math.cos(beam_dir[1])
    endz = length * math.cos(beam_dir[2])
    
    ax.plot([startx, endx], [starty, endy], [startz, endz], c='c', alpha=0.5)
    
    # plot screen(s)
    j = 0
    for k in range(num_scrn):
        yaw = scrn_angl[j].value
        pitch = scrn_angl[j+1].value
        roll = scrn_angl[j+2].value
        

# Functions to aid in output

def createList(values_list):
    newlist = []
    for i in range(len(values_list)):
        newlist.append(f'{values_list[i]}')
        newlist.append(' ')
    return newlist

def createOutput(num_mag, norm_mag_dim, norm_mag_pos, norm_fld_comps, num_particles, norm_beam_pos, norm_beam_energy, converted_beam_dir,
                 pos_sprd, energy_sprd, div_sprd, num_scrn, norm_scrn_dim, norm_scrn_pos, norm_scrn_angl):
    mag_num = [f'{num_mag}', ' ']
    mag_dim = createList(norm_mag_dim)
    mag_pos = createList(norm_mag_pos)
    mag_field = createList(norm_fld_comps)
    mag_info = mag_num + mag_dim + mag_pos + mag_field

    particle_num = [f'{num_particles}', ' ']
    beam_pos = createList(norm_beam_pos)
    beam_nrg = [f'{norm_beam_energy}', ' ']
    beam_dir = createList(norm_beam_dir)
    beam_info = particle_num + beam_pos + beam_nrg + beam_dir

    spread_pos = createList(pos_sprd)
    spread_nrg = [f'{energy_sprd}', ' ']
    spread_div = createList(div_sprd)
    spread_info = spread_pos + spread_nrg + spread_div

    screen_num = [f'{num_scrn}', ' ']
    screen_dim = createList(norm_scrn_dim)
    screen_pos = createList(norm_scrn_pos)
    screen_angles = createList(norm_scrn_angl)
    screen_info = screen_num + screen_dim + screen_pos + screen_angles
    
    return mag_info, beam_info, spread_info, screen_info

def writeOutput(units, num_mag, mag_dim, mag_pos, fld_comps, num_particles, beam_pos, beam_energy, converted_beam_dir, pos_sprd,
                energy_sprd, converted_div_spread, num_scrn, scrn_dim, scrn_pos, converted_scrn_angl):
    
    outfile = open('input_deck.txt', 'w')
    
    norm_mag_dim, norm_magt_pos, norm_fld_comps, norm_beam_pos,\
    norm_beam_energy, norm_scrn_dim, norm_scrn_pos = normalizeValues(units, num_mag, mag_dim, mag_pos, fld_comps, beam_pos, beam_energy, 
                                                                     scrn_dim, scrn_pos)
    
    mag_info, beam_info, spread_info, screen_info = createOutput(num_mag, norm_mag_dim, norm_mag_pos, norm_fld_comps, num_particles,
                                                                 norm_beam_pos, norm_beam_energy, converted_beam_dir, num_scrns, pos_sprd,
                                                                 energy_sprd, converted_div_sprd, num_scrn, norm_scrn_dim, norm_scrn_pos,
                                                                 converted_scrn_angl)
    outfile.writelines(units)
    outfile.write('\n')
    outfile.writelines(mag_info)
    outfile.write('\n')
    outfile.writelines(beam_info)
    outfile.write('\n')
    outfile.writelines(spread_info)
    outfile.write('\n')
    outfile.writelines(screen_info)
    
    outfile.close()