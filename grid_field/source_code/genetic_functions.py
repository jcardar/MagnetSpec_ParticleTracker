import os

def read_global_bounds(access_infile):
    bounds_strings = access_infile.readline().split()
    x_max = float(bounds_strings[0])
    x_min = float(bounds_strings[1])
    y_max = float(bounds_strings[2])
    y_min = float(bounds_strings[3])
    z_max = float(bounds_strings[4])
    z_min = float(bounds_strings[5])

    global_bounds = (x_max, x_min, y_max, y_min, z_max, z_min)
    return global_bounds

def read_access_and_mutation_sizes(access_infile):
    mutation_sizes = []
    mag_access = []
    screen_access = []

    mag_strings = access_infile.readline().split()
    #print("Mag strings are ") 
    #print(mag_strings)
    #mag_strings.pop()
    i = 0
    #print(len(mag_strings))
    while i < len(mag_strings):
        #print(i)
        mag_access_value = int(mag_strings[i])
        i += 1
        mag_access.append(mag_access_value)
        if mag_access_value == -1:
            mag_mutation_value = float(mag_strings[i])
            mutation_sizes.append(mag_mutation_value)
            i += 1
        
    
    screen_pos_strings = access_infile.readline().split()
    #print(screen_pos_strings)
    #screen_pos_strings.pop()
    j = 0
    while j < len(screen_pos_strings):
        #print(j)
        screen_access_value = int(screen_pos_strings[j])
        screen_access.append(screen_access_value)
        if screen_access_value == -1:
            j += 1
            screen_mutation_value = float(screen_pos_strings[j])
            mutation_sizes.append(screen_mutation_value)
        j += 1

    screen_angles_strings = access_infile.readline().split()
    #screen_angles_strings.pop()
    k = 0
    while k < len(screen_angles_strings):
        screen_access_value = int(screen_angles_strings[k])
        screen_access.append(screen_access_value)
        if screen_access_value == -1:
            k += 1
            screen_mutation_value = float(screen_angles_strings[k])
            mutation_sizes.append(screen_mutation_value)
        k += 1
    #print(screen_access)
    return mutation_sizes, mag_access, screen_access

def read_starting_points(mag_access, screen_access):
    infile = open(os.path.join(os.path.dirname(__file__),os.pardir,'data','analysis','input_deck.txt'), 'r')
    outfile = open(os.path.join(os.path.dirname(__file__),os.pardir,'data','analysis','original_input_deck.txt'), 'w')
    starting_points = []

    first_line = infile.readline()
    outfile.write(first_line)

    second_line = infile.readline().split()
    num_mag = int(second_line[0])
    i = 0
    offset_2l = num_mag*3 + 1
    #print(len(mag_access))
    for ii in range(num_mag*3):
        if mag_access[i] == -1:
            starting_point_value = float(second_line[i+offset_2l])
            starting_points.append(starting_point_value)
        i += 1
    for jj in range(len(second_line)-1):
        outfile.write(second_line[jj])
        outfile.write(' ')
    outfile.write(second_line[-1])

    third_line = infile.readline()
    outfile.write(third_line)

    fourth_line = infile.readline()
    outfile.write(fourth_line)

    fifth_line = infile.readline().split()
    #print(fifth_line)
    num_screen = int(fifth_line[0])
    j = 0
    offset_5l = num_screen*2 + 1
    #print(len(screen_access))
    #print(screen_access)
    for ii in range(num_screen*6):
     #   print(j)
        if screen_access[j] == -1:
            starting_point_value = float(fifth_line[j+offset_5l])
            starting_points.append(starting_point_value)
        j += 1
    for jj in range(len(fifth_line)-1):
        outfile.write(fifth_line[jj])
        outfile.write(' ')
    outfile.write(fifth_line[-1])

    sixth_line = infile.readline()
    outfile.write(sixth_line)

    seventh_line = infile.readline()
    outfile.write(seventh_line)

    infile.close()
    outfile.close()
    return starting_points

def edit_input_deck(mutated_values, mag_access, screen_access):
    #If this is ever moved from "grid_field/source_code/.", test whether __file__ refer to this file or file that imports this file.
    infile = open(os.path.join(os.path.dirname(__file__),os.pardir,'data','analysis','input_deck.txt'), 'r')
    outfile = open(os.path.join(os.path.dirname(__file__),os.pardir,'data','analysis','temp_input_deck.txt'), 'w')

    mutate_value_count = 0

    first_line = infile.readline()
    outfile.write(first_line)

    second_line = infile.readline().split()
    num_mag = int(second_line[0])
    i = 0
    offset_2l = num_mag*3 + 1
    for ii in range(num_mag*3):
        if mag_access[i] == -1:
            new_value = mutated_values[mutate_value_count]
            second_line[i+offset_2l] = f'{new_value}'
            mutate_value_count += 1
        i += 1
    for jj in range(len(second_line)-1):
        outfile.write(second_line[jj])
        outfile.write(' ')
    outfile.write(second_line[-1])

    third_line = infile.readline()
    outfile.write(third_line)

    fourth_line = infile.readline()
    outfile.write(fourth_line)

    fifth_line = infile.readline().split()
    num_screen = int(fifth_line[0])
    j = 0
    offset_5l = num_screen*2 + 1
    for ii in range(num_screen*6):
        if screen_access[j] == -1:
            new_value = mutated_values[mutate_value_count]
            fifth_line[j+offset_5l] = f'{new_value}'
            mutate_value_count += 1
        j += 1
    for jj in range(len(fifth_line)-1):
        outfile.write(fifth_line[jj])
        outfile.write(' ')
    outfile.write(fifth_line[-1])

    sixth_line = infile.readline()
    outfile.write(sixth_line)

    seventh_line = infile.readline()
    outfile.write(seventh_line)

    infile.close()
    outfile.close()

    os.replace(os.path.join(os.path.dirname(__file__),os.pardir,'data','analysis','temp_input_deck.txt'),os.path.join(os.path.dirname(__file__),os.pardir,'data','analysis','input_deck.txt'))

def read_bounds(mag_access, screen_access, global_bounds):
    infile = open(os.path.join(os.path.dirname(__file__),os.pardir,'data','analysis','input_deck.txt'), 'r')
    lbound = []
    ubound = [] 
    first_line = infile.readline()
    #outfile.write(first_line)

    second_line = infile.readline().split()
    num_mag = int(second_line[0])
    mag_dims = [float(second_line[ii]) for ii in range(1,3*num_mag,1)]#width,length,height
    i = 0
    offset_2l = num_mag*3 + 1
    #print(global_bounds)
    #print(global_bounds[0])
    for ii in range(num_mag*3):
        if mag_access[i] == -1 and i%3 == 0:
            #print("This is true")
            print(global_bounds[0])
            lbound.append(global_bounds[1])
            ubound.append(global_bounds[0] - mag_dims[i])
            print(global_bounds[0]-mag_dims[i])
            #starting_points.append(starting_point_value)
        if mag_access[i] == -1 and i%3 == 1:
            lbound.append(global_bounds[3] + 0.5*mag_dims[i])
            ubound.append(global_bounds[2] - 0.5*mag_dims[i])
        if mag_access[i] == -1 and i%3 == 2:
            lbound.append(global_bounds[5] + 0.5*mag_dims[i])
            ubound.append(global_bounds[4] - 0.5*mag_dims[i])
        i += 1
    #for jj in range(len(second_line)-1):
        #outfile.write(second_line[jj])
        #outfile.write(' ')
    #outfile.write(second_line[-1])

    third_line = infile.readline()
    #outfile.write(third_line)

    fourth_line = infile.readline()
    #outfile.write(fourth_line)

    fifth_line = infile.readline().split()
    num_screen = int(fifth_line[0])
    screen_dims = [float(fifth_line[jj]) for jj in range(1, 2*num_screen, 1)] #[length, height]
    j = 0
    screen_count = 0
    offset_5l = num_screen*2 + 1
    for ii in range(num_screen*6):
        if screen_access[i] == -1 and i%6 == 0:
            lbound.append(global_bounds[1] + screen_dims[0 + 2*screen_count])
            ubound.append(global_bounds[0] - screen_dims[0 + 2*screen_count])
            #starting_points.append(starting_point_value)
        if screen_access[i] == -1 and i%6 == 1:
            lbound.append(global_bounds[3])
            ubound.append(global_bounds[2] - screen_dims[0 + 2*screen_count])
        if screen_access[i] == -1 and i%6 == 2:
            lbound.append(global_bounds[5] + screen_dims[1 + 2*screen_count])
            ubound.append(global_bounds[4] - screen_dims[1 + 2*screen_count])
        if screen_access[i] == -1 and i%6 == 3:
            lbound.append(0)
            ubound.append(3.14159)
            #starting_points.append(starting_point_value)
        if screen_access[i] == -1 and i%6 == 4:
            lbound.append(0)
            ubound.append(3.14159)
        if screen_access[i] == -1 and i%6 == 5:
            lbound.append(0)
            ubound.append(3.14159)
        j += 1
        if j==5:
            screen_count += 1
    #for jj in range(len(fifth_line)-1):
        #outfile.write(fifth_line[jj])
        #outfile.write(' ')
    #outfile.write(fifth_line[-1])

    sixth_line = infile.readline()
    #outfile.write(sixth_line)

    seventh_line = infile.readline()
    #outfile.write(seventh_line)

    infile.close()
    #outfile.close()
    print(lbound)
    print(ubound)
    return lbound, ubound