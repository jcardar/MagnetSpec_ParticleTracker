# Define global widgets here
import ipywidgets as widgets

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
        widgetx = widgets.BoundedFloatText(
            value=0,
            min=0,
            description=f'width {i+1}',
        )
        widgety = widgets.BoundedFloatText(
            value=0,
            min=0,
            description=f'length {i+1}',
        )
        widgetz = widgets.BoundedFloatText(
            value=0,
            min=0,
            description=f'height {i+1}',
        )
        listOfWidgets.append(widgetx)
        listOfWidgets.append(widgety)
        listOfWidgets.append(widgetz)
    return listOfWidgets

def dynamicFloatValue_Magnet_Position(num_of_magnets):
    listOfWidgets = []
    for i in range(num_of_magnets):
        widgetx = widgets.BoundedFloatText(
            value=0,
            min=0,
            description=f'x pos {i+1}',
        )
        widgety = widgets.FloatText(
            value=0,
            description=f'y pos {i+1}',
        )
        widgetz = widgets.FloatText(
            value=0,
            description=f'z pos {i+1}',
        )
        listOfWidgets.append(widgetx)
        listOfWidgets.append(widgety)
        listOfWidgets.append(widgetz)
    return listOfWidgets

def dynamicFloatValue_Magnetic_Field_Strength(num_of_magnets):
    listOfWidgets = []
    for i in range(num_of_magnets):
        widgetx = widgets.FloatText(
            value=0,
            description=f'x comp {i+1}',
        )
        widgety = widgets.FloatText(
            value=0,
            description=f'y comp {i+1}',
        )
        widgetz = widgets.FloatText(
            value=0,
            description=f'z comp {i+1}',
        )
        listOfWidgets.append(widgetx)
        listOfWidgets.append(widgety)
        listOfWidgets.append(widgetz)
    return listOfWidgets

# Define Beam attributes
number_of_particles = widgets.BoundedIntText(
    value=1,
    min=1,
    step=1,
    description='# of Particles',
    disabled=False
)

def dynamicFloatValue_Beam_Start_Position():
    listOfWidgets = []
    widgetx = widgets.FloatText(
        value=0,
        description='x start pos',
    )
    widgety = widgets.FloatText(
        value=0,
        description='y start pos',
    )
    widgetz = widgets.FloatText(
        value=0,
        description='z start pos',
    )
    listOfWidgets.extend([widgetx, widgety, widgetz])
    return listOfWidgets

beam_energy = widgets.BoundedFloatText(
    value=1,
    min=0,
    description='beam energy'
)

def dynamicFloatValue_Beam_Direction():
    listOfWidgets = []
    widgetx = widgets.FloatText(
        value=0,
        description='x angle',
    )
    widgety = widgets.FloatText(
        value=0,
        description='y angle',
    )
    widgetz = widgets.FloatText(
        value=0,
        description='z angle',
    )
    listOfWidgets.extend([widgetx, widgety, widgetz])
    return listOfWidgets

# Define Beam spread
def dynamicFloatValue_Beam_Position_Spread():
    listOfWidgets = []
    widgetx = widgets.FloatText(
        value=0,
        description='x pos spread',
    )
    widgety = widgets.FloatText(
        value=0,
        description='y pos spread',
    )
    widgetz = widgets.FloatText(
        value=0,
        description='z pos spread',
    )
    listOfWidgets.extend([widgetx, widgety, widgetz])
    return listOfWidgets

beam_energy_spread = widgets.BoundedFloatText(
    value=1,
    min=0,
    description='nrg spread'
)

beam_divergence_spread = widgets.BoundedFloatText(
    value=1,
    min=0,
    description='div spread'
)

# Define Screen attriutes
number_of_screens = widgets.BoundedIntText(
    value=2,
    min=1,
    step=1,
    description='# of Screens',
    disabled=False
)

def dynamicFloatValue_Screen_Dimensions(num_of_screens):
    listOfWidgets = []
    for i in range(num_of_screens):
        widgety = widgets.BoundedFloatText(
            value=0,
            min=0,
            description=f'length {i+1}',
        )
        widgetz = widgets.BoundedFloatText(
            value=0,
            min=0,
            description=f'height {i+1}',
        )
        listOfWidgets.append(widgety)
        listOfWidgets.append(widgetz)
    return listOfWidgets

def dynamicFloatValue_Screen_Position(num_of_screens):
    listOfWidgets = []
    for i in range(num_of_screens):
        widgetx = widgets.FloatText(
            value=0,
            description=f'x pos {i+1}',
        )
        widgety = widgets.FloatText(
            value=0,
            description=f'y pos {i+1}',
        )
        widgetz = widgets.FloatText(
            value=0,
            description=f'z pos {i+1}',
        )
        listOfWidgets.append(widgetx)
        listOfWidgets.append(widgety)
        listOfWidgets.append(widgetz)
    return listOfWidgets

def dynamicFloatValue_Screen_Angles(num_of_screens):
    listOfWidgets = []
    for i in range(num_of_screens):
        widgetz = widgets.FloatText(
            value=0,
            description=f'yaw angle {i+1}',
        )
        widgety = widgets.FloatText(
            value=0,
            description=f'pitch angle {i+1}',
        )
        widgetx = widgets.FloatText(
            value=0,
            description=f'roll angle {i+1}',
        )
        listOfWidgets.append(widgetz)
        listOfWidgets.append(widgety)
        listOfWidgets.append(widgetx)
    return listOfWidgets

# Functions to aid in output
def convertValues(units, beam_angles, screen_angles, beam_nrg):
    if units[4] == 'mrad':
        for i in range(len(beam_angles)):
            beam_angles[i].value /= 1000
        for j in range(len(screen_angles)):
            screen_angles[j].value /= 1000
    elif units[4] == 'Degrees':
        for i in range(len(beam_angles)):
            beam_angles[i].value *= (3.14159/180)
        for j in range(len(screen_angles)):
            screen_angles[j].value *= (3.14159/180)
    
    rest_nrg = 0.511 # MeV
    if units[2] == 'eV':
        rest_nrg = 0.511 * (10**6)
    elif units[2] == 'GeV':
        rest_nrg = 0.511 * (10**-3)
    beam_nrg = (rest_nrg + beam_nrg)/rest_nrg
    
    return beam_angles, screen_angles, beam_nrg

def createList(widgets_list):
    newlist = []
    for i in range(len(widgets_list)):
        newlist.append(f'{widgets_list[i].value}')
        newlist.append(' ')
    return newlist

def createOutput(magNum, magDim, magPos, magFld, partNum, beamPos, beamNrg, beamDir, sprdPos, sprdNrg, sprdDiv, scrnNum, scrnDim, scrnPos, scrnAngl):
    mag_num = [f'{magNum.value}', ' ']
    mag_dim = createList(magDim)
    mag_pos = createList(magPos)
    mag_field = createList(magFld)
    mag_info = mag_num + mag_dim + mag_pos + mag_field

    particle_num = [f'{partNum.value}', ' ']
    beam_pos = createList(beamPos)
    beam_nrg = [f'{beamNrg.value}', ' ']
    beam_dir = createList(beamDir)
    beam_info = particle_num + beam_pos + beam_nrg + beam_dir

    spread_pos = createList(sprdPos)
    spread_nrg = [f'{sprdNrg.value}', ' ']
    spread_div = [f'{sprdDiv.value}', ' ']
    spread_info = spread_pos + spread_nrg + spread_div

    screen_num = [f'{scrnNum.value}', ' ']
    screen_dim = createList(scrnDim)
    screen_pos = createList(scrnPos)
    screen_angles = createList(scrnAngl)
    screen_info = screen_num + screen_dim + screen_pos + screen_angles
    
    return mag_info, beam_info, spread_info, screen_info

def writeOutput(unitLENval, unitNRGval, unitANGLval, unitFLDval, MAGnum, MAGdim, MAGpos, MAGfld, PARTnum, BEAMpos, BEAMnrg, BEAMdir, SPRDpos, SPRDnrg, SPRDdiv, SCRNnum, SCRNdim, SCRNpos, SCRNangl):
    
    outfile = open('input_deck.txt', 'w')
    
    unit_info = [units_length.value, ' ', units_energy.value, ' ', units_angles.value, ' ', units_magnetic_field.value]
    
    BEAMdir, SCRNangl, BEAMnrg.value = convertValues(unit_info, BEAMdir, SCRNangl, BEAMnrg.value)
    
    mag_info, beam_info, spread_info, screen_info = createOutput(MAGnum, MAGdim, MAGpos, MAGfld, PARTnum, BEAMpos, BEAMnrg, BEAMdir,
                                                                 SPRDpos, SPRDnrg, SPRDdiv, SCRNnum, SCRNdim, SCRNpos, SCRNangl)
    
    outfile.writelines(unit_info)
    outfile.write('\n')
    outfile.writelines(mag_info)
    outfile.write('\n')
    outfile.writelines(beam_info)
    outfile.write('\n')
    outfile.writelines(spread_info)
    outfile.write('\n')
    outfile.writelines(screen_info)
    
    outfile.close()