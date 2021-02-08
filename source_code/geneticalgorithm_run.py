
from opt import OptimizationProblem, Goal
from opt.ga import DifferentialGaDriver, BasicGaDriver
import os
import time
import datetime
import numpy as np
import numpy
import matplotlib.pyplot as plt
import copy
import math
from spectrometer_energy_resolution import *
from genetic_functions import *
from python_run_cpp_code_mac import *
#from phaseCal_agrt import *

#Recommended to use differential driver: changes mutation size to correspond to how close it is to optimized value

####z0 = 0.25
####s0 = 3.5e-6
####x, y, phi0, I0 = shadowgraphy(zeff = z0, sigma = s0, Lx = 10*13e-6, EkeV = [5], Nx = 2048, Ny = 512)

save_all_mag_spec_images = False

access_file = open(os.path.join(os.path.dirname(__file__),os.pardir,'data','analysis','algorithm_access.txt'), 'r')
global_bounds = read_global_bounds(access_file)
mutation_size, mag_access, screen_access = read_access_and_mutation_sizes(access_file)
starting_point = read_starting_points(mag_access, screen_access)
#print(starting_point)
#print(len(starting_point))
#print(mutation_size)
#print(len(mutation_size))
num_parent   = 20
num_children = 5
num_elitism  = 5
isfirst = True
run_spectrometer_code()
#import_relevant_data()
normalizing_fom = 0
energy_range = np.array([500.51121400000005, 2501.53407])
#normalizing_fom = energy_weighted_resolution(energy_range, normalizing_fom, isfirst)
normalizing_fom = energy_and_divergence_resolution(energy_range, normalizing_fom, isfirst)
print(f"Normalizing fom is {normalizing_fom}")
#isfirst = False
best_fom = 1.5
#Change these:
#starting_point = [0.5, 1.5e-6]
#Mutation size only two parameters - first changes first starting point, second changes second starting point
#mutation_size = (0.01, 0.05e-6)

iter_num = 100

print("Done with pre- main()")
def main():
    
    # bound limit, zeff [0.05, 7], sigma [1e-6, 5e-6]
    print("main(). About to read bounds")
    lbound, ubound = read_bounds(mag_access, screen_access, global_bounds)
    print("Bounds read. About to set beset_fom = 1.5")
    best_fom = 1.5
    #print(len(lbound))
    #lbound = (0.05, 1e-6) 
    #ubound = (1.0 , 7e-6)

    population_id = 0

    def validate(params):
        return True
    
    # Fom
    def synthetic_diagnostic(params):
        #energy_resolution(part_on_screen_part_index = params[0], energy = params[1], energy_range = params[2], posx = params[3], normalizing_fom = params[4], isfirst = isfirst)
        global isfirst
        if isfirst == True:
            print("FIRST RUN")
            isfirst = False
            fom = energy_and_divergence_resolution(energy_range, normalizing_fom, isfirst)
            return fom

        edit_input_deck(params, mag_access, screen_access)
        run_spectrometer_code()
        #fom = energy_weighted_resolution(energy_range, normalizing_fom, isfirst)
        fom = energy_and_divergence_resolution(energy_range, normalizing_fom, isfirst)

        
        return fom
        ##x, y, phi, I = shadowgraphy(zeff = params[0], sigma = params[1], Lx = 10*13e-6, EkeV = [5], Nx = 2048, Ny = 512)
        ##fom = np.sum(abs(I - I0)**0.5) #/ np.sum(abs(I0)**0.5)
    
    def evaluate(params):
        #FOM here, synthetic diagnostic
        global best_fom
        time.sleep(0.001)

        nonlocal population_id
        population_id +=1
        if population_id % num_parent == 1:
            population_id = 1

        fom = synthetic_diagnostic(params)
        #print(f"Best fom is {best_fom}")
        if fom < best_fom:
            writeout_best_fom_input()
            best_fom = fom

        print('Iter No.: %d, Child: %d, Fom: %.1f'%
                (ga.generation_number, population_id, fom))
        
        if save_raw_option.lower() == 'y':
            pass
#            np.savetxt(rawpath+'raw_iter%d_people%d.txt'%(ga.generation_number, population_id), np.asarray(im))
        return fom

    ## TODO, detect keyboard
#     print("\nYou are about to run the GAD code. \n\nYou can stop running by pressing 'q' at any time.") 
    print("\nYou are about to run the GAD code. ") 

    while True:
        save_raw_option = input('\nDo you want to save all raw data (y/n)?  ')
        if save_raw_option.lower() == 'n':
            break
        elif save_raw_option.lower() == 'y':
            run_num = input('\nPlease input a run number or a file name, so I can make a new folder for you: ')
            date =  datetime.datetime.now().strftime('%Y-%m-%d-')
            rawpath = './save/%srun-%s-RAW/'%(date, run_num)
            if not os.path.exists(rawpath):
                os.makedirs(rawpath)
                break
            else:
                print('\nPath exists! Try again!')
                continue

    while True:
        save_raw_option = input('\nDo you want to save figure-of-merit progression (y/n)?  ')
        if save_raw_option.lower() == 'n':
            break
        elif save_raw_option.lower() == 'y':
            run_num = input('\nPlease input a run number or a file name, so I can make a new folder for you: ')
            date =  datetime.datetime.now().strftime('%Y-%m-%d-')
            progpath = './save/%srun-%s-progression/'%(date, run_num)
            if not os.path.exists(progpath):
                os.makedirs(progpath)
                break
            else:
                print('\nPath exists! Try again!')
                continue
            
    goal = int(input('\nWhat is the goal of your optimization  problem?\nTips: Answer "1" for minimization or "0" for maximization : '))
    if goal not in [0, 1]:
        goal = int(input('\nAnswer "0" or "1", please try again:'))
    if goal == 0:
        goal = Goal.maximize
    else:
        goal = Goal.minimize
    
    
#CORE OF THE CODE:
    problem = OptimizationProblem((len(lbound),), evaluate, goal,
                   lbound=lbound, ubound=ubound, validator=validate)
    #print(problem)
                
    

    diffdriver = input('\nDo you want to use Differential Driver? \n y/n?')

    if diffdriver.lower() == 'n':
        print("setting diff driver settings, setting basic problem including mutation size")
        ga = BasicGaDriver(problem , mutation_size)
        print("about to set num_parents")
        ga['population'] = num_parent
        print("about to set num_children")
        ga['selection']  = num_children
        print("about to set num_elitism")
        ga['elitism']    = num_elitism
        print("about to set mutation_size for diff dirver")
        ga['mutation_size'] = 1.0
        ga['mutation_probability'] = 1.0
        print("about to populate diff dirver")
        print(f"starting point is {[starting_point[ii] for ii in range(len(starting_point))]}")
        ga.populate(tuple([starting_point[jj] for jj in range(len(starting_point))]), scale = 1)
        print("diff driver populated")
        
    elif diffdriver.lower() == 'y':
        ga = DifferentialGaDriver(problem)
        ga['population'] = num_parent
        ga['selection']  = num_children
        ga['elitism']    = num_elitism
        ga['mutation_size'] = 2.0
        ga['mutation_probability'] = 1.0
        ga.populate(tuple([starting_point[jj] for jj in range(len(starting_point))]), scale = mutation_size)  # for DifferentialGaDriver
    
    print("Figures for iteration tracking about to create")
    fig = plt.figure(figsize = (12,8))
    plt.subplots_adjust(wspace = 0.45 , hspace = 0.35, left =0.1 , right = 0.95, top = 0.95, bottom = 0.07)
    plt.ion()
    plt.show()
    plt.tight_layout()
    ax1 = plt.subplot(231)
    ax2 = plt.subplot(232)
    ax3 = plt.subplot(233)

    ax4 = plt.subplot(223)
    ax5 = plt.subplot(224)

    # ax1.set_title('Figure of merit progression')
    ax3.set_xlabel('Iterations')
    ax2.set_xlabel('Iterations')
    ax1.set_xlabel('Iterations')

    ax1.set_ylabel('Figures of merit')
    ax2.set_ylabel('Gene[0] (Magnet y-position, normalized)')
    ax3.set_ylabel('Gene[1] (Screen x-position, normalized)')
    

    #ax2.axhline(z0, ls = '--', color = 'gray', lw = 2)
    #ax3.axhline(s0*1e6, ls = '--', color = 'gray', lw = 2)

    ax4.set_title('Gene[2] (Screen y-position, normalized)')
    ax5.set_title('Gene[3] (Screen yaw-angle, radians)')

    fig1, ax = plt.subplots(figsize = (10,10))
    fig1.set_tight_layout(True)
    print("Figures created. Lists for GA values being instantiated")
    FOM =  []
    GENE = []
    FOM_ALL = []
    GENE_ALL =[]
    print("Going into iteration loop.")


    i = 0
    while i <= iter_num:
        print("running through iterations. about to ga.generate()")
        #if i == 1:
        #    isfirst = False
    # while True:

        ga.generate()
        # time.sleep(3)
        print("fom_generate done")
        fom_all  = copy.copy(ga.fitness) 
        fom      = copy.copy(ga.selection_fitness)
        gene_all = np.array(copy.copy(ga.population)).ravel()
        gene     = np.array(copy.copy(ga.selection )).ravel()

        plt.figure(fig.number)
        ax1.plot([ga.generation_number] * ga['population'], fom_all, '+', color = 'gray')
        ax1.plot([ga.generation_number] * ga['selection'],  fom,     '*', color = 'r')
        if max(fom_all)==3.14159:
            ax1.set_ylim([(min(fom)-0.1),1.4])
        else:
            ax1.set_ylim([(min(fom)-0.1),1.2])
        ax2.plot([ga.generation_number] * ga['population'], gene_all.reshape((ga['population'], len(lbound)))[:,0], '+', color = 'gray')
        ax2.plot([ga.generation_number] * ga['selection'],  gene.reshape(    (ga['selection'],len(lbound))  )[:,0], 'r*')
        ax3.plot([ga.generation_number] * ga['population'], gene_all.reshape((ga['population'], len(lbound)))[:,1], '+', color = 'gray')
        ax3.plot([ga.generation_number] * ga['selection'],  gene.reshape(    (ga['selection'],len(lbound))  )[:,1], 'r*')
        ax4.plot([ga.generation_number] * ga['population'], gene_all.reshape((ga['population'], len(lbound)))[:,2], '+', color = 'gray')
        ax4.plot([ga.generation_number] * ga['selection'],  gene.reshape(    (ga['selection'],len(lbound))  )[:,2], 'r*')
        ax5.plot([ga.generation_number] * ga['population'], gene_all.reshape((ga['population'], len(lbound)))[:,3], '+', color = 'gray')
        ax5.plot([ga.generation_number] * ga['selection'],  gene.reshape(    (ga['selection'],len(lbound))  )[:,3], 'r*')


        #####_, _, _, I_best = shadowgraphy(zeff = gene[0], sigma = gene[1], Lx = 10*13e-6, EkeV = [5], Nx =  2048, Ny = 512)
        #E_fom_best = energy_resolution()
        #####p4 = ax4.pcolormesh(x*1e6, y*1e6, I0, vmin = 0, vmax = 0.5, rasterized = True)
        #####p5 = ax5.pcolormesh(x*1e6, y*1e6, I_best, vmin = 0, vmax = 0.5, rasterized = True)
        # plt.colorbar(p4, ax = ax4)
        # plt.colorbar(p5, ax = ax5)

        plt.savefig(progpath+'prog_iter%d.png'%(ga.generation_number), dpi = 300)

        plt.figure(fig1.number)
        ax.plot([ga.generation_number] * ga['population'], fom_all, '+', color = 'gray')
        ax.plot([ga.generation_number] * ga['selection'],  fom,     '*', color = 'r')
        ax.set_title(f"Figure of Merit Progression Up To Iteration {ga.generation_number}")
        ax.set_ylabel("Figure of Merit")
        ax.set_xlabel("Iteration")
        if max(fom_all)==3.14159:
            ax.set_ylim([(min(fom)-0.1),1.4])
        else:
            ax.set_ylim([(min(fom)-0.1),1.2])
        plt.savefig(progpath+'prog_iter%d_JUST_FOM.png'%(ga.generation_number), dpi = 300)

        plt.draw()
        FOM.append(fom)
        GENE.append(gene)
        FOM_ALL.append(fom_all)
        GENE_ALL.append(gene_all)

        #plt.pause(0.001)


        i+=1

    #     # Quit loop when focus on terminal
    #     if msvcrt.kbhit():
    #         keyboard_input = msvcrt.getwche()
    #         if keyboard_input.lower() == 'q':
    #             print('Quiting...')
    #             break

    #     # Quit loop when focus on figure
    #     def on_key(event):
    #         global mpl_key_event
    #         mpl_key_event = event.key

    #     fig.canvas.mpl_connect( 'key_press_event' , on_key )
    #     if mpl_key_event:
    #         if mpl_key_event.lower() == 'q':
    #             print('Quiting...')
    #             break

    while True:
        save_option = input('Do you want to save data (including best/original input decks) (y/n)?  ')
        if save_option.lower() == 'n':
            break
        elif save_option.lower() == 'y':
            run_num = input('Please input a run number or a file name, so I can make a new folder for you: ')
            date =  datetime.datetime.now().strftime('%Y-%m-%d-')
            newpath = './save/%srun-%s/'%(date, run_num)
            from shutil import copy2
            if not os.path.exists(newpath):
                os.makedirs(newpath)
                np.savetxt(newpath + 'fom.txt', np.array(FOM), header = 'Run %s, fom'%run_num )
                np.savetxt(newpath + 'gene.txt', np.array(GENE), header = 'Run %s, gene\n %i of genes'%(run_num, len(starting_point)))
                np.savetxt(newpath + 'gene_all.txt', np.array(GENE_ALL), header = 'Run %s, gene all\n %i of genes'%(run_num, len(starting_point)))
                np.savetxt(newpath + 'fom_all.txt', np.array(FOM_ALL), header = 'Run %s, fom all'%run_num)
                plt.savefig(newpath + 'fom.pdf', dpi = 300, bbox_inches = 'tight')
                #nput_files = ['../data/analysis/original_input_deck.txt', '../data/analysis/best_input_deck.txt']
                best_deck = os.path.join(os.path.dirname(__file__),os.pardir,'data','analysis','best_input_deck.txt')
                first_deck = os.path.join(os.path.dirname(__file__),os.pardir,'data','analysis','original_input_deck.txt')
                copy2(best_deck, newpath)
                copy2(first_deck, newpath)
                break
            else:
                print('Path exists! Try again!')
                continue
        else:
            print("I don't understand your command, try again? ")

if __name__ == '__main__':
    main()

