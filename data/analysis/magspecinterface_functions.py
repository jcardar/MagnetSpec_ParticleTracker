import ipywidgets as widgets
import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Define units

units_length = widgets.Dropdown(
    options=['mm', 'cm', 'm'],
    value='cm',
    description='Length Unit',
    disabled=False,
)

units_energy = widgets.Dropdown(
    options=['eV', 'MeV', 'GeV'],
    value='MeV',
    description='Energy Unit',
    disabled=False,
)

units_angles = widgets.Dropdown(
    options=['mrad', 'Radians', 'Degrees'],
    value='mrad',
    description='Angle Unit',
    disabled=False,
)

units_magnetic_field = widgets.Dropdown(
    options=['Tesla', 'Gauss'],
    value='Tesla',
    description='Field Unit',
    disabled=False,
)

# Define initialization types

init_position = widgets.Dropdown(
    options=[('Gaussian',0) , ('Uniform',1) , ('Scan',2)],
    value=2,
    description='Position',
    disabled=False,
)

init_energy = widgets.Dropdown(
    options=[('Gaussian',0) , ('Uniform',1) , ('Log',2)],
    value=1,
    description='Energy',
    disabled=False,
)

init_divergence = widgets.Dropdown(
    options=[('Gaussian',0) , ('Uniform',1) , ('Scan',2)],
    value=2,
    description='Divergence',
    disabled=False,
)

# Define Coordinate System attributes

global_max_x = widgets.FloatText(
    value=100,
    description='global x max',
    disabled=False
)

global_min_x = widgets.FloatText(
    value=-100,
    description='global x min',
    disabled=False
)

global_max_y = widgets.FloatText(
    value=100,
    description='global y max',
    disabled=False
)

global_min_y = widgets.FloatText(
    value=-100,
    description='global y min',
    disabled=False
)

global_max_z = widgets.FloatText(
    value=100,
    description='global z max',
    disabled=False
)

global_min_z = widgets.FloatText(
    value=-100,
    description='global z min',
    disabled=False
)

global_bounds = [global_max_x, global_min_x, global_max_y, global_min_y, global_max_z, global_min_z]

# Define Magnet attributes

number_of_magnets = widgets.BoundedIntText(
    value=1,
    min=1,
    step=1,
    description='# of Magnets',
    disabled=False
)

def dynamicFloatValue_Magnet_Dimensions(num_of_magnets, global_bounds):
    listOfWidgets = []
    for i in range(num_of_magnets):
        widget1 = widgets.BoundedFloatText(
            value=0,
            max=global_bounds[2].value,
            min=0,
            description=f'width {i+1}',
        )
        widget2 = widgets.BoundedFloatText(
            value=0,
            max=global_bounds[0].value,
            min=0,
            description=f'length {i+1}',
        )
        widget3 = widgets.BoundedFloatText(
            value=0,
            max=global_bounds[4].value,
            min=0,
            description=f'height {i+1}',
        )
        listOfWidgets.extend([widget1, widget2, widget3])
    return listOfWidgets

def dynamicFloatValue_Permanent_Magnet_Dimension(num_of_magnets):
    listOfWidgets = []
    for i in range(num_of_magnets):
        widget = widgets.FloatText(
            value=0,
            max=999999,
            min=0,
            description=f'dimension {i+1}',
        )
        listOfWidgets.append(widget)
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
        listOfWidgets.extend([widget1, widget2, widget3])
        j += 3
    return listOfWidgets

def dynamicFloatValue_Magnetic_Field_Value(num_of_magnets):
    listOfWidgets = []
    for i in range(num_of_magnets):
        widget = widgets.FloatText(
            value=1,
            description=f'B field {i+1}',
        )
        listOfWidgets.append(widget)
    return listOfWidgets

def dynamicFloatValue_Magnetic_Field_Axis(num_of_magnets):
    listOfWidgets = []
    for i in range(num_of_magnets):
        widget = widgets.Dropdown(
            options=[('x-axis','x') ,('y-axis','y') , ('z-axis','z')],
            value='z',
            description='along the',
        )
        listOfWidgets.append(widget)
    return listOfWidgets

# Define Beam attributes

number_of_particles = widgets.BoundedIntText(
    value=1,
    max=999999,
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
    max=9999999,
    description='beam energy'
)

def dynamicFloatValue_Beam_Direction(angle_unit):
    listOfWidgets = []
    init_val = math.pi*500 # mrad
    if angle_unit == 'Radians':
        init_val = math.pi/2
    elif angle_unit == 'Degrees':
        init_val = 90
    
    widget1 = widgets.FloatText(
        value=0,
        description='x angle',
    )
    widget2 = widgets.FloatText(
        value=init_val,
        description='y angle',
    )
    widget3 = widgets.FloatText(
        value=init_val,
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
    max=9999999,
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
    value=1,
    min=0,
    step=1,
    description='# of Screens',
    disabled=False
)

def dynamicFloatValue_Screen_Dimensions(num_of_screens, global_bounds):
    listOfWidgets = []
    for i in range(num_of_screens):
        widget1 = widgets.BoundedFloatText(
            value=0,
            max=global_bounds[0].value,
            min=0,
            description=f'length {i+1}',
        )
        widget2 = widgets.BoundedFloatText(
            value=0,
            max=global_bounds[4].value,
            min=0,
            description=f'height {i+1}',
        )
        listOfWidgets.extend([widget1, widget2])
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
        listOfWidgets.extend([widget1, widget2, widget3])
    return listOfWidgets

def dynamicFloatValue_Screen_Position(angle_unit, num_of_screens, global_bounds, scrn_dim, scrn_angl):
    listOfWidgets = []
    j = 0
    k = 0
    angle_multiplier = 1 # Radians
    if angle_unit == 'mrad':
        angle_multiplier = math.pow(10,-3)
    elif angle_unit == 'Degrees':
        angle_multiplier = math.pi/180
        
    for i in range(num_of_screens):
        alpha = scrn_angl[j].value * angle_multiplier
        beta = scrn_angl[j+1].value * angle_multiplier
        gamma = scrn_angl[j+2].value * angle_multiplier
        
        c12x = abs(scrn_dim[k+1].value*0.5*(np.cos(alpha)*np.sin(beta)*np.cos(gamma) + np.sin(alpha)*np.sin(gamma)))
        c3x =  scrn_dim[k].value*np.cos(alpha)*np.cos(beta) + scrn_dim[k+1].value*0.5*(np.cos(alpha)*np.sin(beta)*np.cos(gamma) +\
np.sin(alpha)*np.sin(gamma))
        c4x =  scrn_dim[k].value*np.cos(alpha)*np.cos(beta) - scrn_dim[k+1].value*0.5*(np.cos(alpha)*np.sin(beta)*np.cos(gamma) +\
np.sin(alpha)*np.sin(gamma))
        c12y = abs(scrn_dim[k+1].value*0.5*(np.sin(alpha)*np.sin(beta)*np.cos(gamma) - np.cos(alpha)*np.sin(gamma)))
        c3y =  (scrn_dim[k].value*np.sin(alpha)*np.cos(beta)) + scrn_dim[k+1].value*0.5*(np.sin(alpha)*np.sin(beta)*np.cos(gamma) -\
np.cos(alpha)*np.sin(gamma))
        c4y =  (scrn_dim[k].value*np.sin(alpha)*np.cos(beta)) - scrn_dim[k+1].value*0.5*(np.sin(alpha)*np.sin(beta)*np.cos(gamma) -\
np.cos(alpha)*np.sin(gamma))
        c12z = abs(scrn_dim[k+1].value*0.5*np.cos(beta)*np.cos(gamma))
        c3z =  -1*(scrn_dim[k].value*np.sin(beta) + scrn_dim[k+1].value*0.5*np.cos(beta)*np.cos(gamma))
        c4z =  -1*(scrn_dim[k].value*np.sin(beta) - scrn_dim[k+1].value*0.5*np.cos(beta)*np.cos(gamma))
        
        xmax = max(c12x, c3x, c4x)
        xmin = min(c12x, c3x, c4x)
        ymax = max(c12y, c3y, c4y)
        ymin = min(c12y, c3y, c4y)
        zmax = max(c12z, c3z, c4z)
        zmin = min(c12z, c3z, c4z)
        widget1 = widgets.BoundedFloatText(
            value=0,
            max=global_bounds[0].value - xmax,
            min=global_bounds[1].value + xmin,
            description=f'x pos {i+1}',
        )
        widget2 = widgets.BoundedFloatText(
            value=0,
            max=global_bounds[2].value - ymax,
            min=global_bounds[3].value + ymin,
            description=f'y pos {i+1}',
        )
        widget3 = widgets.BoundedFloatText(
            value=0,
            max=global_bounds[4].value - zmax,
            min=global_bounds[5].value + zmin,
            description=f'z pos {i+1}',
        )
        listOfWidgets.extend([widget1, widget2, widget3])
        j += 3
        k += 2
    return listOfWidgets

# Functions to aid in standardizing values

def averageB0(num_of_magnets, field_values):
    all_fields_total = 0.0
    for i in range(num_of_magnets):
        all_fields_total += field_values[i].value
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

def normalizeValues(units, num_mag, mag_dim, Pmag_dim, mag_pos, fld_vals, beam_pos, beam_energy, pos_sprd, energy_sprd, scrn_dim, 
                    scrn_pos):
    n_mag_dim = []
    n_Pmag_dim = []
    n_mag_pos = []
    n_fld_vals = []
    n_beam_pos = []
    n_pos_sprd = []
    n_scrn_dim = []
    n_scrn_pos = []
    
    aveB_0 = averageB0(num_mag, fld_vals)
    q_e = 1.602177 * math.pow(10,-19)
    m_e = 9.109384 * math.pow(10,-31)
    c = 2.997925 * math.pow(10,8)
    omega_div_c = (q_e * aveB_0) / (m_e * c)
    
    # distances
    distance_multiplier = 1 # m
    if units[0] == 'mm':
        distance_multiplier = math.pow(10,-3)
    elif units[0] == 'cm':
        distance_multiplier = math.pow(10,-2)
        
    for i in range(len(mag_dim)):
        n_mag_dim.append( mag_dim[i].value * distance_multiplier * omega_div_c )
    for i in range(len(Pmag_dim)):
        n_Pmag_dim.append( Pmag_dim[i].value * distance_multiplier * omega_div_c )
    for i in range(len(mag_pos)):
        n_mag_pos.append( mag_pos[i].value * distance_multiplier * omega_div_c )
    for i in range(3):
        n_beam_pos.append( beam_pos[i].value * distance_multiplier * omega_div_c )
    for i in range(3):
        n_pos_sprd.append( pos_sprd[i].value * distance_multiplier * omega_div_c )
    for i in range(len(scrn_dim)):
        n_scrn_dim.append( scrn_dim[i].value * distance_multiplier * omega_div_c )
    for i in range(len(scrn_pos)):
        n_scrn_pos.append( scrn_pos[i].value * distance_multiplier * omega_div_c )
    
    # magnetic field
    magnetic_multiplier = 1 # Tesla
    if units[6] == 'Gauss':
        magnetic_multiplier = math.pow(10,-4)
    
    for i in range(len(fld_vals)):
        n_fld_vals.append( (fld_vals[i].value * magnetic_multiplier) / aveB_0 )
    
    # energy
    rest_energy = 0.511 # MeV
    if units[2] == 'eV':
        rest_energy = 0.511 * math.pow(10,6)
    elif units[2] == 'GeV':
        rest_energy = 0.511 * math.pow(10,-3)
    
    n_beam_energy = (rest_energy + beam_energy.value)/rest_energy
    n_energy_sprd = (rest_energy + energy_sprd.value)/rest_energy
    
    return n_mag_dim, n_Pmag_dim, n_mag_pos, n_fld_vals, n_beam_pos, n_beam_energy, n_pos_sprd, n_energy_sprd, n_scrn_dim, n_scrn_pos

def normalizeMu0(num_mag, fld_vals):
    aveB_0 = averageB0(num_mag, fld_vals)
    q_e = 1.602177 * math.pow(10,-19)
    m_e = 9.109384 * math.pow(10,-31)
    c = 2.997925 * math.pow(10,8)
    omega_div_c = (q_e * aveB_0) / (m_e * c)
    
    mu_0 = 4 * math.pi * math.pow(10,-7) * ((omega_div_c * q_e * q_e) / m_e)
    
    return mu_0

# Functions to aid in plotting

def checkObjBounds(xyz_ranges, xmax, xmin, ymax, ymin, zmax, zmin, color, bound_check_bool):
    """
    xyz_ranges -> list of previously calculated ?min ?max values for non-intersecting objects [min, max]
    ?min/max is the bounds for the current object being looked at
    """
    #screen_bool argument
    ## xx = 0.0, yy = 0.0, zz = 0.0
    j = 0
    for i in range(len(xyz_ranges)//6):
        if (xmax>=xyz_ranges[j] and xmax<=xyz_ranges[j+1]) or (xmin>=xyz_ranges[j] and xmin<=xyz_ranges[j+1]) or \
        (xyz_ranges[j]>=xmin and xyz_ranges[j]<=xmax and xyz_ranges[j+1]>=xmin and xyz_ranges[j+1]<=xmax):
            if (ymax>=xyz_ranges[j+2] and ymax<=xyz_ranges[j+3]) or (ymin>=xyz_ranges[j+2] and ymin<=xyz_ranges[j+3]) or \
            (xyz_ranges[j+2]>=ymin and xyz_ranges[j+2]<=ymax and xyz_ranges[j+3]>=ymin and xyz_ranges[j+3]<=ymax):
                if (zmax>=xyz_ranges[j+4] and zmax<=xyz_ranges[j+5]) or (zmin>=xyz_ranges[j+4] and zmin<=xyz_ranges[j+5]) or \
                (xyz_ranges[j+4]>=zmin and xyz_ranges[j+4]<=zmax and xyz_ranges[j+5]>=zmin and xyz_ranges[j+5]<=zmax):
                    #if screen_bool == True:
                        #if xx (... check collision bounds like above conditional statements) and xx != np.nan:
                            #if yy (...): 
                        #check the point 
                    #else:
                    color = 'r'
                    bound_check_bool = True
                    return color, bound_check_bool
        j += 6
    return color, bound_check_bool

# Functions to aid in output

def createList(values_list):
    newlist = []
    for i in range(len(values_list)):
        newlist.append(f'{values_list[i]}')
        newlist.append(' ')
    return newlist

def createAxesList(axes_list):
    newlist = []
    for i in range(len(axes_list)):
        newlist.append(axes_list[i].value)
        newlist.append(' ')
    return newlist

def createOutput(num_mag, mag_dim, Pmag_dim, mag_pos, fld_vals, fld_axs, num_particles, beam_pos, beam_energy, beam_dir, pos_sprd, 
                 energy_sprd, div_sprd, num_scrn, scrn_dim, scrn_pos, scrn_angl):
    mag_num = [f'{num_mag}', ' ']
    dim_mag = createList(mag_dim)
    pos_mag = createList(mag_pos)
    mag_field = createList(fld_vals)
    pMag_dim = createList(Pmag_dim)
    field_axes = createAxesList(fld_axs)
    mag_info = mag_num + dim_mag + pos_mag + mag_field + pMag_dim + field_axes

    particle_num = [f'{num_particles}', ' ']
    pos_beam = createList(beam_pos)
    beam_nrg = [f'{beam_energy}', ' ']
    dir_beam = createList(beam_dir)
    beam_info = particle_num + pos_beam + beam_nrg + dir_beam

    spread_pos = createList(pos_sprd)
    spread_nrg = [f'{energy_sprd}', ' ']
    spread_div = createList(div_sprd)
    spread_info = spread_pos + spread_nrg + spread_div

    screen_num = [f'{num_scrn}', ' ']
    screen_dim = createList(scrn_dim)
    screen_pos = createList(scrn_pos)
    screen_angles = createList(scrn_angl)
    screen_info = screen_num + screen_dim + screen_pos + screen_angles
    
    return mag_info, beam_info, spread_info, screen_info

def writeOutput(units, num_mag, mag_dim, Pmag_dim, mag_pos, fld_vals, fld_axs, num_particles, beam_pos, beam_energy, beam_dir, pos_sprd, 
                energy_sprd, div_sprd, num_scrn, scrn_dim, scrn_pos, scrn_angl, init_types):
    
    outfile = open('input_deck.txt', 'w')
    
    n_mag_dim, n_Pmag_dim, n_mag_pos, n_fld_vals, n_beam_pos, n_beam_energy, \
    n_pos_sprd, n_energy_sprd, n_scrn_dim, n_scrn_pos = normalizeValues(units, num_mag, mag_dim, Pmag_dim, mag_pos, fld_vals, beam_pos,
                                                                        beam_energy, pos_sprd, energy_sprd, scrn_dim, scrn_pos)
    mu_0 = normalizeMu0(num_mag, fld_vals)
    
    mag_info, beam_info, spread_info, screen_info = createOutput(num_mag, n_mag_dim, n_Pmag_dim, n_mag_pos, n_fld_vals, fld_axs, 
                                                                 num_particles, n_beam_pos, n_beam_energy, beam_dir, n_pos_sprd, 
                                                                 n_energy_sprd, div_sprd, num_scrn, n_scrn_dim, n_scrn_pos, scrn_angl)
    outfile.writelines(units)
    outfile.write('\n')
    outfile.writelines(mag_info)
    outfile.write('\n')
    outfile.writelines(beam_info)
    outfile.write('\n')
    outfile.writelines(spread_info)
    outfile.write('\n')
    outfile.writelines(screen_info)
    outfile.write('\n')
    outfile.writelines(init_types)
    outfile.write('\n')
    outfile.write(f'{mu_0}')
    
    outfile.close()

# Plots and outputs

def DisplayAndOutput(global_bounds, units, num_mag, mag_dim, Pmag_dim, mag_pos, fld_vals, fld_axs, num_particles, beam_pos, beam_energy,
                     beam_dir, pos_sprd, energy_sprd, div_spread, num_scrn, scrn_dim, scrn_pos, scrn_angl, init_types):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlabel(f'x position ({units[0]})')
    ax.set_ylabel(f'y position ({units[0]})')
    ax.set_zlabel(f'z position ({units[0]})')
    ax.set_title('User Defined Magnetic Spectrometer')
    ax.autoscale_view()
    
    # plot global boundaries
    xmax = global_bounds[0].value
    xmin = global_bounds[1].value
    ymax = global_bounds[2].value
    ymin = global_bounds[3].value
    zmax = global_bounds[4].value
    zmin = global_bounds[5].value
    
    ax.plot([xmin, xmin], [ymin, ymax], [zmax, zmax], linestyle=':', c='k', alpha=0.5)
    ax.plot([xmin, xmin], [ymin, ymax], [zmin, zmin], linestyle=':', c='k', alpha=0.5)
    ax.plot([xmin, xmin], [ymax, ymax], [zmin, zmax], linestyle=':', c='k', alpha=0.5)
    ax.plot([xmin, xmin], [ymin, ymin], [zmin, zmax], linestyle=':', c='k', alpha=0.5)
    ax.plot([xmax, xmax], [ymin, ymax], [zmax, zmax], linestyle=':', c='k', alpha=0.5)
    ax.plot([xmax, xmax], [ymin, ymax], [zmin, zmin], linestyle=':', c='k', alpha=0.5)
    ax.plot([xmax, xmax], [ymax, ymax], [zmin, zmax], linestyle=':', c='k', alpha=0.5)
    ax.plot([xmax, xmax], [ymin, ymin], [zmin, zmax], linestyle=':', c='k', alpha=0.5)
    ax.plot([xmin, xmax], [ymax, ymax], [zmax, zmax], linestyle=':', c='k', alpha=0.5)
    ax.plot([xmin, xmax], [ymin, ymin], [zmax, zmax], linestyle=':', c='k', alpha=0.5)
    ax.plot([xmin, xmax], [ymax, ymax], [zmin, zmin], linestyle=':', c='k', alpha=0.5)
    ax.plot([xmin, xmax], [ymin, ymin], [zmin, zmin], linestyle=':', c='k', alpha=0.5)
    
    bound_check_bool = False # check for overlapping boundaries
    
    # plot magnet(s)
    mag_xyz_ranges = [] # order of min-max of x, y, z
    i = 0
    for l in range(num_mag):
        mxmax = mag_pos[i].value + mag_dim[i+1].value
        mxmin = mag_pos[i].value
        mymax = mag_pos[i+1].value + (0.5 * mag_dim[i].value)
        mymin = mag_pos[i+1].value - (0.5 * mag_dim[i].value)
        mzmax = mag_pos[i+2].value + (0.5 * mag_dim[i+2].value)
        mzmin = mag_pos[i+2].value - (0.5 * mag_dim[i+2].value)
       
        clr = 'g'
        if l > 0:
            clr, bound_check_bool = checkObjBounds(mag_xyz_ranges, mxmax, mxmin, mymax, mymin, mzmax, mzmin, clr, bound_check_bool)
       
        ax.plot([mxmin, mxmin], [mymin, mymax], [mzmax, mzmax], c=clr, alpha=0.5)
        ax.plot([mxmin, mxmin], [mymin, mymax], [mzmin, mzmin], c=clr, alpha=0.5)
        ax.plot([mxmin, mxmin], [mymax, mymax], [mzmin, mzmax], c=clr, alpha=0.5)
        ax.plot([mxmin, mxmin], [mymin, mymin], [mzmin, mzmax], c=clr, alpha=0.5)
        ax.plot([mxmax, mxmax], [mymin, mymax], [mzmax, mzmax], c=clr, alpha=0.5)
        ax.plot([mxmax, mxmax], [mymin, mymax], [mzmin, mzmin], c=clr, alpha=0.5)
        ax.plot([mxmax, mxmax], [mymax, mymax], [mzmin, mzmax], c=clr, alpha=0.5)
        ax.plot([mxmax, mxmax], [mymin, mymin], [mzmin, mzmax], c=clr, alpha=0.5)
        ax.plot([mxmin, mxmax], [mymax, mymax], [mzmax, mzmax], c=clr, alpha=0.5)
        ax.plot([mxmin, mxmax], [mymin, mymin], [mzmax, mzmax], c=clr, alpha=0.5)
        ax.plot([mxmin, mxmax], [mymax, mymax], [mzmin, mzmin], c=clr, alpha=0.5)
        ax.plot([mxmin, mxmax], [mymin, mymin], [mzmin, mzmin], c=clr, alpha=0.5)
        
        mag_xyz_ranges.extend([mxmin, mxmax, mymin, mymax, mzmin, mzmax])
        i += 3
    
    # plot beam
    length = mag_pos[0].value - beam_pos[0].value
    startx = beam_pos[0].value
    starty = beam_pos[1].value
    startz = beam_pos[2].value
    endx = length * math.cos(beam_dir[0]) + startx   # using direction cosines
    endy = length * math.cos(beam_dir[1]) + starty
    endz = length * math.cos(beam_dir[2]) + startz
    
    ax.plot([startx, endx], [starty, endy], [startz, endz], linestyle='--', c='m', alpha=0.7)
    
    # plot screen(s)
    scrn_xyz_ranges = [] # order of min-max of x, y, z
    j = 0
    k = 0
    for l in range(num_scrn):
        alpha = scrn_angl[j]
        beta = scrn_angl[j+1]
        gamma = scrn_angl[j+2]
        
        plane_point = np.array([scrn_pos[j].value, scrn_pos[j+1].value, scrn_pos[j+2].value])
        normal_vec  = np.array([np.sin(alpha)*(np.cos(beta)**2)*np.cos(gamma) + np.sin(beta)*(np.sin(alpha)*np.sin(beta)*np.cos(gamma) -\
np.cos(alpha)*np.sin(gamma)), -(np.cos(alpha)*(np.cos(beta)**2)*np.cos(gamma) + np.sin(beta)*\
(np.cos(alpha)*np.sin(beta)*np.cos(gamma) + np.sin(alpha)*np.sin(gamma))), np.cos(alpha)*np.cos(beta)*\
(np.sin(alpha)*np.sin(beta)*np.cos(gamma) - np.cos(alpha)*np.sin(gamma)) - np.sin(alpha)*np.cos(beta)*\
(np.cos(alpha)*np.sin(beta)*np.cos(gamma) + np.sin(alpha)*np.sin(gamma))])
        d_plane     = plane_point.dot(normal_vec)

        corner1_x = scrn_pos[j].value + scrn_dim[k+1].value*0.5*(np.cos(alpha)*np.sin(beta)*np.cos(gamma) + np.sin(alpha)*np.sin(gamma))
        corner2_x = scrn_pos[j].value - scrn_dim[k+1].value*0.5*(np.cos(alpha)*np.sin(beta)*np.cos(gamma) + np.sin(alpha)*np.sin(gamma))
        corner3_x = scrn_pos[j].value + scrn_dim[k].value*np.cos(alpha)*np.cos(beta) + scrn_dim[k+1].value*0.5*\
(np.cos(alpha)*np.sin(beta)*np.cos(gamma) + np.sin(alpha)*np.sin(gamma))
        corner4_x = scrn_pos[j].value + scrn_dim[k].value*np.cos(alpha)*np.cos(beta) - scrn_dim[k+1].value*0.5*\
(np.cos(alpha)*np.sin(beta)*np.cos(gamma) + np.sin(alpha)*np.sin(gamma))

        corner1_y = scrn_pos[j+1].value + scrn_dim[k+1].value*0.5*(np.sin(alpha)*np.sin(beta)*np.cos(gamma) - np.cos(alpha)*np.sin(gamma))
        corner2_y = scrn_pos[j+1].value - scrn_dim[k+1].value*0.5*(np.sin(alpha)*np.sin(beta)*np.cos(gamma) - np.cos(alpha)*np.sin(gamma))
        corner3_y = scrn_pos[j+1].value + (scrn_dim[k].value*np.sin(alpha)*np.cos(beta)) + scrn_dim[k+1].value*0.5*\
(np.sin(alpha)*np.sin(beta)*np.cos(gamma) - np.cos(alpha)*np.sin(gamma))
        corner4_y = scrn_pos[j+1].value + (scrn_dim[k].value*np.sin(alpha)*np.cos(beta)) - scrn_dim[k+1].value*0.5*\
(np.sin(alpha)*np.sin(beta)*np.cos(gamma) - np.cos(alpha)*np.sin(gamma))

        corner1_z = scrn_pos[j+2].value + scrn_dim[k+1].value*0.5*np.cos(beta)*np.cos(gamma)
        corner2_z = scrn_pos[j+2].value - scrn_dim[k+1].value*0.5*np.cos(beta)*np.cos(gamma)
        corner3_z = scrn_pos[j+2].value - scrn_dim[k].value*np.sin(beta) + scrn_dim[k+1].value*0.5*np.cos(beta)*np.cos(gamma)
        corner4_z = scrn_pos[j+2].value - scrn_dim[k].value*np.sin(beta) - scrn_dim[k+1].value*0.5*np.cos(beta)*np.cos(gamma)
        
        corner1  = np.array([corner1_x, corner1_y, corner1_z])
        corner2  = np.array([corner2_x, corner2_y, corner2_z])
        corner3  = np.array([corner3_x, corner3_y, corner3_z])
        corner4  = np.array([corner4_x, corner4_y, corner4_z])
        
        sxmax = max(corner1_x, corner2_x, corner3_x, corner4_x)
        sxmin = min(corner1_x, corner2_x, corner3_x, corner4_x)
        symax = max(corner1_y, corner2_y, corner3_y, corner4_y)
        symin = min(corner1_y, corner2_y, corner3_y, corner4_y)
        szmax = max(corner1_z, corner2_z, corner3_z, corner4_z)
        szmin = min(corner1_z, corner2_z, corner3_z, corner4_z)
        
        clr = 'b'
        
        
        if normal_vec[0] != 0:
            yy     = np.linspace(symin, symax, num=200)
            zz     = np.linspace(szmin, szmax, num=200)
            yy, zz = np.meshgrid(yy,zz)
            xx     = ( -normal_vec[1]*yy - normal_vec[2]*zz + d_plane ) / normal_vec[0]
        elif normal_vec[1] != 0:
            zz     = np.linspace(szmin, szmax, num=200)
            xx     = np.linspace(sxmin, sxmax, num=200)
            xx, zz = np.meshgrid(xx,zz)
            yy     = ( -normal_vec[0]*xx - normal_vec[2]*zz + d_plane ) / normal_vec[1]
        elif normal_vec[2] != 0:
            xx = np.linspace(sxmin, sxmax, num=200)
            yy = np.linspace(symin, symax, num=200)
            xx, yy = np.meshgrid(xx,yy)
            zz = ( -normal_vec[0]*xx - normal_vec[1]*yy + d_plane ) / normal_vec[2]
            
        area0 = scrn_dim[k+1].value*scrn_dim[k].value
        tmin = [[(-((np.array(corner1 - np.array([xx[jj][kk], yy[jj][kk], zz[jj][kk]])) @ np.array(corner2 - corner1))/\
abs(np.array(corner2 - corner1)@np.array(corner2 - corner1)))) for jj in range(len(xx[0]))] for kk in range(len(xx))]
        closest_point = [[[corner1_x + (corner2_x - corner1_x)*tmin[jj][kk], corner1_y + (corner2_y - corner1_y)*tmin[jj][kk], corner1_z +\
(corner2_z - corner1_z)*tmin[jj][kk]] for jj in range(len(xx[0]))] for kk in range(len(xx))]
        dist_to_intersect = [[np.sqrt(((closest_point[jj][kk][0]-xx[jj][kk])**2) + ((closest_point[jj][kk][1]-yy[jj][kk])**2) + \
((closest_point[jj][kk][2]-zz[jj][kk])**2)) for jj in range(len(xx[0]))] for kk in range(len(xx))];
        area12 = np.array([[0.5*scrn_dim[k+1].value*dist_to_intersect[jj][kk] for jj in range(len(xx[0]))] for kk in range(len(xx))])

        tmin = [[(-((np.array(corner3 - np.array([xx[jj][kk], yy[jj][kk], zz[jj][kk]])) @ np.array(corner1 - corner3))/\
abs(np.array(corner1 - corner3)@np.array(corner1 - corner3)))) for jj in range(len(xx[0]))] for kk in range(len(xx))]
        closest_point = [[[corner3_x + (corner1_x - corner3_x)*tmin[jj][kk], corner3_y + (corner1_y - corner3_y)*tmin[jj][kk], corner3_z +\
(corner1_z - corner3_z)*tmin[jj][kk]] for jj in range(len(xx[0]))] for kk in range(len(xx))]
        dist_to_intersect = [[np.sqrt(((closest_point[jj][kk][0]-xx[jj][kk])**2) + ((closest_point[jj][kk][1]-yy[jj][kk])**2) + \
((closest_point[jj][kk][2]-zz[jj][kk])**2)) for jj in range(len(xx[0]))] for kk in range(len(xx))]
        area13 = np.array([[0.5*scrn_dim[k].value*dist_to_intersect[jj][kk] for jj in range(len(xx[0]))] for kk in range(len(xx))])

        tmin = [[(-((np.array(corner4 - np.array([xx[jj][kk], yy[jj][kk], zz[jj][kk]])) @ np.array(corner2 - corner4))/\
abs(np.array(corner2 - corner4)@np.array(corner2 - corner4)))) for jj in range(len(xx[0]))] for kk in range(len(xx))]
        closest_point = [[[corner4_x + (corner2_x - corner4_x)*tmin[jj][kk], corner4_y + (corner2_y - corner4_y)*tmin[jj][kk], corner4_z +\
(corner2_z - corner4_z)*tmin[jj][kk]] for jj in range(len(xx[0]))] for kk in range(len(xx))]
        dist_to_intersect = [[np.sqrt(((closest_point[jj][kk][0]-xx[jj][kk])**2) + ((closest_point[jj][kk][1]-yy[jj][kk])**2) + \
((closest_point[jj][kk][2]-zz[jj][kk])**2)) for jj in range(len(xx[0]))] for kk in range(len(xx))]
        area24 = np.array([[0.5*scrn_dim[k].value*dist_to_intersect[jj][kk] for jj in range(len(xx[0]))] for kk in range(len(xx))])

        tmin = [[(-((np.array(corner3 - np.array([xx[jj][kk], yy[jj][kk], zz[jj][kk]])) @ np.array(corner4 - corner3))/\
abs(np.array(corner4 - corner3)@np.array(corner4 - corner3)))) for jj in range(len(xx[0]))] for kk in range(len(xx))]
        closest_point = [[[corner3_x + (corner4_x - corner3_x)*tmin[jj][kk], corner3_y + (corner4_y - corner3_y)*tmin[jj][kk], corner3_z +\
(corner4_z - corner3_z)*tmin[jj][kk]] for jj in range(len(xx[0]))] for kk in range(len(xx))]
        dist_to_intersect = [[np.sqrt(((closest_point[jj][kk][0]-xx[jj][kk])**2) + ((closest_point[jj][kk][1]-yy[jj][kk])**2) + \
((closest_point[jj][kk][2]-zz[jj][kk])**2)) for jj in range(len(xx[0]))] for kk in range(len(xx))];
        area34 = np.array([[0.5*scrn_dim[k+1].value*dist_to_intersect[jj][kk] for jj in range(len(xx[0]))] for kk in range(len(xx))])
        if normal_vec[0] == min(normal_vec):
            xx[(area12+ area13+ area24+ area34) > (area0)+(area0*0.001)] = np.nan
        elif normal_vec[1] == min(normal_vec):
            yy[(area12+ area13+ area24+ area34) > (area0)+(area0*0.001)] = np.nan
        elif normal_vec[2] == min(normal_vec):
            zz[(area12+ area13+ area24+ area34) > (area0)+(area0*0.001)] = np.nan
            
        clr, bound_check_bool = checkObjBounds(mag_xyz_ranges, sxmax, sxmin, symax, symin, szmax, szmin, clr, bound_check_bool)
        if l > 0:
            clr, bound_check_bool = checkObjBounds(scrn_xyz_ranges, sxmax, sxmin, symax, symin, szmax, szmin, clr, bound_check_bool)
            
        ax.plot_surface(xx, yy, zz, alpha=0.4, color=clr)
        ax.scatter(scrn_pos[j].value,scrn_pos[j+1].value,scrn_pos[j+2].value, color='c')
        ax.scatter([corner1_x, corner2_x, corner3_x, corner4_x], [corner1_y, corner2_y, corner3_y, corner4_y], 
                   [corner1_z, corner2_z, corner3_z, corner4_z], color='b', alpha=0.5)
        scrn_xyz_ranges.extend([sxmin, sxmax, symin, symax, szmin, szmax])
        j += 3
        k += 2
    
    # output data
    dscrpt = 'Save Inputs'
    icn = 'check'
    if bound_check_bool == True:
        dscrpt = 'Objects Overlap'
        icn = 'exclamation'
    button = widgets.Button(description=dscrpt, icon=icn, disabled=bound_check_bool, layout=widgets.Layout(width='50%', height='80px'))
    display(button)
    
    def on_button_clicked(b):
        writeOutput(units, num_mag, mag_dim, Pmag_dim, mag_pos, fld_vals, fld_axs, num_particles, beam_pos, beam_energy, beam_dir, 
                    pos_sprd, energy_sprd, div_spread, num_scrn, scrn_dim, scrn_pos, scrn_angl, init_types)
        print('Inputs saved and exported!')
    button.on_click(on_button_clicked)

# Genetic algorithm

def createBoundsList(bounds_list, length_unit, num_mag, field_values):
    newlist = []
    
    aveB_0 = averageB0(num_mag, field_values)
    q_e = 1.602177 * math.pow(10,-19)
    m_e = 9.109384 * math.pow(10,-31)
    c = 2.997925 * math.pow(10,8)
    omega_div_c = (q_e * aveB_0) / (m_e * c)
    
    distance_multiplier = 1 # m
    if length_unit == 'mm':
        distance_multiplier = math.pow(10,-3)
    elif length_unit == 'cm':
        distance_multiplier = math.pow(10,-2)
    
    for i in range(len(bounds_list)):
        newlist.append(f'{(bounds_list[i].value)* distance_multiplier * omega_div_c }')
        newlist.append(' ')
        
    return newlist

def bool_output(check_bool):
    bool_num = 0
    if check_bool:
        bool_num = -1
    return bool_num

def genetic_algorithm_setup(global_bounds, num_mag, num_screen, length_unit, field_values):
    button_label = widgets.Label(value='Do you want to use the genetic algorithm? ')
    button = widgets.Button(description='Yes', icon='check', disabled=False)
    display(button_label, button)
    
    check_label = widgets.Label(value='The genetic algorithm has access to mutate:')
    mag_pos_label = widgets.Label(value='Magnet position')
    screen_pos_label = widgets.Label(value='Screen position')
    screen_angles_label = widgets.Label(value='Screen angles')
    
    mag_pos_checks = []
    mag_pos_steps = []
    
    for i in range(num_mag):
        check1 = widgets.Checkbox(value=False, description=f'x pos {i+1}')
        check2 = widgets.Checkbox(value=False, description=f'y pos {i+1}')
        check3 = widgets.Checkbox(value=False, description=f'z pos {i+1}')
        mag_pos_checks.extend([check1, check2, check3])
        
        step1 = widgets.BoundedFloatText(value=0.0, min=0.0, max=100.0, description='step size')
        step2 = widgets.BoundedFloatText(value=0.0, min=0.0, max=100.0, description='step size')
        step3 = widgets.BoundedFloatText(value=0.0, min=0.0, max=100.0, description='step size')
        mag_pos_steps.extend([step1, step2, step3])
    
    mag_pos_left_box = widgets.VBox(mag_pos_checks)
    mag_pos_right_box = widgets.VBox(mag_pos_steps)
    mag_pos_box = widgets.HBox([mag_pos_left_box, mag_pos_right_box])
    
    screen_pos_checks = []
    screen_pos_steps = []
    
    for i in range(num_screen):
        check1 = widgets.Checkbox(value=False, description=f'x pos {i+1}')
        check2 = widgets.Checkbox(value=False, description=f'y pos {i+1}')
        check3 = widgets.Checkbox(value=False, description=f'z pos {i+1}')
        screen_pos_checks.extend([check1, check2, check3])
        
        step1 = widgets.BoundedFloatText(value=0.0, min=0.0, max=100.0, description='step size')
        step2 = widgets.BoundedFloatText(value=0.0, min=0.0, max=100.0, description='step size')
        step3 = widgets.BoundedFloatText(value=0.0, min=0.0, max=100.0, description='step size')
        screen_pos_steps.extend([step1, step2, step3])
    
    screen_pos_left_box = widgets.VBox(screen_pos_checks)
    screen_pos_right_box = widgets.VBox(screen_pos_steps)
    screen_pos_box = widgets.HBox([screen_pos_left_box, screen_pos_right_box])
    
    screen_angles_checks = []
    screen_angles_steps = []
    
    for i in range(num_screen):
        check1 = widgets.Checkbox(value=False, description=f'yaw {i+1}')
        check2 = widgets.Checkbox(value=False, description=f'pitch {i+1}')
        check3 = widgets.Checkbox(value=False, description=f'roll {i+1}')
        screen_angles_checks.extend([check1, check2, check3])
        
        step1 = widgets.BoundedFloatText(value=0.0, min=0.0, max=100.0, description='step size')
        step2 = widgets.BoundedFloatText(value=0.0, min=0.0, max=100.0, description='step size')
        step3 = widgets.BoundedFloatText(value=0.0, min=0.0, max=100.0, description='step size')
        screen_angles_steps.extend([step1, step2, step3])
    
    screen_angles_left_box = widgets.VBox(screen_angles_checks)
    screen_angles_right_box = widgets.VBox(screen_angles_steps)
    screen_angles_box = widgets.HBox([screen_angles_left_box, screen_angles_right_box])
    
    save_label = widgets.Label(value='Click below to save and output selection')
    save_button = widgets.Button(description='Save Checked Boxes', icon='check', disabled=False, layout=widgets.Layout(width='50%', 
                                                                                                                       height='80px'))
    
    def on_button_clicked(b):
        display(check_label, mag_pos_label, mag_pos_box, screen_pos_label, screen_pos_box, screen_angles_label, screen_angles_box, 
                save_label, save_button)
        button.disabled = True
    button.on_click(on_button_clicked)
    
    def on_save_button_clicked(b):
        outfile = open('algorithm_access.txt', 'w')
        
        aveB_0 = averageB0(num_mag, field_values)
        q_e = 1.602177 * math.pow(10,-19)
        m_e = 9.109384 * math.pow(10,-31)
        c = 2.997925 * math.pow(10,8)
        omega_div_c = (q_e * aveB_0) / (m_e * c)

        distance_multiplier = 1 # m
        if length_unit == 'mm':
            distance_multiplier = math.pow(10,-3)
        elif length_unit == 'cm':
            distance_multiplier = math.pow(10,-2)
        
        bounds_info = createBoundsList(global_bounds, length_unit, num_mag, field_values)
        outfile.writelines(bounds_info)
        outfile.write('\n')
        
        mag_pos_index = 0
        for i in range(num_mag):
            for j in range(3):
                bool_num = bool_output(mag_pos_checks[mag_pos_index+j].value)
                outfile.write(f'{bool_num} ')
                if bool_num == -1:
                    outfile.write(f'{mag_pos_steps[mag_pos_index+j].value*distance_multiplier*omega_div_c} ')
            mag_pos_index += 3
        outfile.write('\n')
        
        screen_pos_index = 0
        for i in range(num_screen):
            for j in range(3):
                bool_num = bool_output(screen_pos_checks[screen_pos_index+j].value)
                outfile.write(f'{bool_num} ')
                if bool_num == -1:
                    outfile.write(f'{screen_pos_steps[screen_pos_index+j].value*distance_multiplier*omega_div_c} ')
            screen_pos_index += 3
        outfile.write('\n')
        
        screen_angles_index = 0
        for i in range(num_screen):
            for j in range(3):
                bool_num = bool_output(screen_angles_checks[screen_angles_index+j].value)
                outfile.write(f'{bool_num} ')
                if bool_num == -1:
                    outfile.write(f'{screen_angles_steps[screen_angles_index+j].value} ')
            screen_angles_index += 3
        
        outfile.close
        
        print('Choices saved and exported!')
    save_button.on_click(on_save_button_clicked)
