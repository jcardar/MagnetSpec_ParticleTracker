#Define global widgets here
import ipywidgets as widgets

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

number_of_magnets = widgets.BoundedIntText(
    value=2,
    min=0,
    step=1,
    description='# of Magnets',
    disabled=False
)


def dynamicFloatValue_Magnet_Position(num_of_widgets):
    listOfWidgets = [];
    for i in range(num_of_widgets):
        widgetx = widgets.FloatText(
            value=0,
            description=f'x of Mag {i+1}',
        )
        widgety = widgets.FloatText(
            value=0,
            description=f'y of Mag {i+1}',
        )
        widgetz = widgets.FloatText(
            value=0,
            description=f'z of Mag {i+1}',
        )
        listOfWidgets.append(widgetx)
        listOfWidgets.append(widgety)
        listOfWidgets.append(widgetz)
    return listOfWidgets