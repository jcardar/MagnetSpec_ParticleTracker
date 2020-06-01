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

# Define Magnet attributes
number_of_magnets = widgets.BoundedIntText(
    value=2,
    min=0,
    step=1,
    description='# of Magnets',
    disabled=False
)

def dynamicFloatValue_Magnet_Position(num_of_magnets):
    listOfWidgets = []
    for i in range(num_of_magnets):
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
    min=0,
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
    min=0,
    step=1,
    description='# of Screens',
    disabled=False
)

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

def dynamicFloatValue_Screen_Angles(num_of_screens):
    listOfWidgets = []
    for i in range(num_of_screens):
        widgetz = widgets.FloatText(
            value=0,
            description=f'z angle {i+1}',
        )
        widgetx = widgets.FloatText(
            value=0,
            description=f'x angle {i+1}',
        )
        listOfWidgets.append(widgetz)
        listOfWidgets.append(widgetx)
    return listOfWidgets