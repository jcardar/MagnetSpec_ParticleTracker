
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
from phaseCal_agrt import *

#Recommended to use differential driver: changes mutation size to correspond to how close it is to optimized value

####z0 = 0.25
####s0 = 3.5e-6
####x, y, phi0, I0 = shadowgraphy(zeff = z0, sigma = s0, Lx = 10*13e-6, EkeV = [5], Nx = 2048, Ny = 512)


num_parent   = 20
num_children = 5
num_elitism  = 5

#Change these:
starting_point = [0.5, 1.5e-6]
#Mutation size only two parameters - first changes first starting point, second changes second starting point
mutation_size = (0.01, 0.05e-6)

iter_num = 40


def main():
    
    # bound limit, zeff [0.05, 7], sigma [1e-6, 5e-6]
    lbound = (0.05, 1e-6) 
    ubound = (1.0 , 7e-6)

    population_id = 0

    def validate(params):
        return True
    
    # Fom
    def synthetic_diagnostic(params):

        ##x, y, phi, I = shadowgraphy(zeff = params[0], sigma = params[1], Lx = 10*13e-6, EkeV = [5], Nx = 2048, Ny = 512)
        ##fom = np.sum(abs(I - I0)**0.5) #/ np.sum(abs(I0)**0.5)

        return fom, I
    
    def evaluate(params):
        #FOM here, synthetic diagnostic
        time.sleep(0.001)

        nonlocal population_id
        population_id +=1
        if population_id % num_parent == 1:
            population_id = 1

        fom, I = synthetic_diagnostic(params)

        print('Iter No.: %d, Child: %d, Fom: %.1f, Genes: %.3f, %.2f'%
                (ga.generation_number, population_id, fom, params[0], params[1]*1e6))

        if save_raw_option.lower() == 'y':
            pass
#             np.savetxt(rawpath+'raw_iter%d_people%d.txt'%(ga.generation_number, population_id), np.asarray(im))
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
                
    

    diffdriver = input('\nDo you want to use Differential Driver? \n y/n?')

    if diffdriver.lower() == 'n':
        
        ga = BasicGaDriver(problem , mutation_size )
        ga['population'] = num_parent
        ga['selection']  = num_children
        ga['elitism']    = num_elitism
        ga['mutation_size'] = 1.0
        ga['mutation_probability'] = 1.0
        ga.populate((starting_point[0],starting_point[1]), scale = 1)
        
    elif diffdriver.lower() == 'y':
        ga = DifferentialGaDriver(problem)
        ga['population'] = num_parent
        ga['selection']  = num_children
        ga['elitism']    = num_elitism
        ga['mutation_size'] = 2.0
        ga['mutation_probability'] = 1.0
        ga.populate((starting_point[0],starting_point[1]), scale = mutation_size)  # for DifferentialGaDriver
    
    
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
    ax2.set_ylabel('Zeff')
    ax3.set_ylabel('Sigma (1e-6)')

    ax2.axhline(z0, ls = '--', color = 'gray', lw = 2)
    ax3.axhline(s0*1e6, ls = '--', color = 'gray', lw = 2)

    ax4.set_title('Goal')
    ax5.set_title('Current best')

    FOM =  []
    GENE = []
    FOM_ALL = []
    GENE_ALL =[]


    i = 0
    while i <= iter_num:
    # while True:

        ga.generate()
        # time.sleep(3)
        fom_all  = copy.copy(ga.fitness) 
        fom      = copy.copy(ga.selection_fitness)
        gene_all = np.array(copy.copy(ga.population)).ravel()
        gene     = np.array(copy.copy(ga.selection )).ravel()

        ax1.plot([ga.generation_number] * ga['population'], fom_all, '+', color = 'gray')
        ax1.plot([ga.generation_number] * ga['selection'],  fom,     '*', color = 'r')
        ax2.plot([ga.generation_number] * ga['population'], gene_all.reshape((ga['population'], len(lbound)))[:,0], '+', color = 'gray')
        ax2.plot([ga.generation_number] * ga['selection'],  gene.reshape(    (ga['selection'],len(lbound))  )[:,0], 'r*')
        ax3.plot([ga.generation_number] * ga['population'], 1e6*gene_all.reshape((ga['population'], len(lbound)))[:,1], '+', color = 'gray')
        ax3.plot([ga.generation_number] * ga['selection'],  1e6*gene.reshape(    (ga['selection'],len(lbound))  )[:,1], 'r*')
        
        _, _, _, I_best = shadowgraphy(zeff = gene[0], sigma = gene[1], Lx = 10*13e-6, EkeV = [5], Nx =  2048, Ny = 512)

        p4 = ax4.pcolormesh(x*1e6, y*1e6, I0, vmin = 0, vmax = 0.5, rasterized = True)
        p5 = ax5.pcolormesh(x*1e6, y*1e6, I_best, vmin = 0, vmax = 0.5, rasterized = True)
        # plt.colorbar(p4, ax = ax4)
        # plt.colorbar(p5, ax = ax5)

        plt.savefig(progpath+'prog_iter%d.png'%(ga.generation_number), dpi = 300)



        plt.draw()
        FOM.append(fom)
        GENE.append(gene)
        FOM_ALL.append(fom_all)
        GENE_ALL.append(gene_all)

        plt.pause(0.001)


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
        save_option = input('Do you want to save data (y/n)?  ')
        if save_option.lower() == 'n':
            break
        elif save_option.lower() == 'y':
            run_num = input('Please input a run number or a file name, so I can make a new folder for you: ')
            date =  datetime.datetime.now().strftime('%Y-%m-%d-')
            newpath = './save/%srun-%s/'%(date, run_num)
            if not os.path.exists(newpath):
                os.makedirs(newpath)
                np.savetxt(newpath + 'fom.txt', np.array(FOM), header = 'Run %s, fom'%run_num )
                np.savetxt(newpath + 'gene.txt', np.array(GENE), header = 'Run %s, gene'%run_num)
                np.savetxt(newpath + 'gene_all.txt', np.array(GENE_ALL), header = 'Run %s, gene all'%run_num)
                np.savetxt(newpath + 'fom_all.txt', np.array(FOM_ALL), header = 'Run %s, fom all'%run_num)
                plt.savefig(newpath + 'fom.pdf', dpi = 300, bbox_inches = 'tight')
                break
            else:
                print('Path exists! Try again!')
                continue
        else:
            print("I don't understand your command, sir, try again? ")

if __name__ == '__main__':
    main()

