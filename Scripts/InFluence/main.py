import argparse
from InFluenceSimulation import InFluenceSimulation
from ImageObjectConstructor import image_pixels

from matplotlib import pyplot as plt
import numpy as np

from filepaths import *
from InFluenceSimulation import *

import itertools
import sys
import time


def twinkling_stars_spinner():
    emojis = ['✨', '  ']  # Customize the emojis as you like
    spinner_characters = itertools.cycle(emojis)
    while True:
        sys.stdout.write(next(spinner_characters))  # Print the next emoji
        sys.stdout.flush()
        time.sleep(0.25)  # Adjust the delay (in seconds) to control the flickering speed
        sys.stdout.write('\b\b')  # Move the cursor back to the beginning of the line


def RunSimulation(simulation_type):

    print('_____________________________________________________________')

    print(' ✨ ✨ ✨ ✨ ✨ ✨ ✨ ✨ ✨ ✨ ✨ ✨ ✨ ✨ ✨ ✨ ✨ ✨ ✨ ✨ ')

    print('                   *  InFluence  *  ')

    print('                                                         ')
    
    print('InFluence simulation will run on ' + simulation_type + '...')
    
    # Start the twinkling stars spinner
    #spinner_process = multiprocessing.Process(target=twinkling_stars_spinner)
    #spinner_process.start()
    image_object = image_pixels(filepaths.InputImageLocation)
    influence_simulator = InFluenceSimulation(image_object)

    print('Dose is:     ', params.dose )
    print('Simulation starting...')
    
    if simulation_type =='GPU':
        SimulationResults = influence_simulator.Run_GPU_KERNEL_MCScatteringSimulation()
    if simulation_type == "PartialVec":
        SimulationResults = influence_simulator.RunPartiallyVectorizedMCScatteringSimulation()
    else:
        SimulationResults = influence_simulator.RunMCSimulation()

    # Stop the twinkling stars spinner
    #spinner_process.terminate()

    
    plt.imshow(np.uint32(SimulationResults))
    plt.axis('off')
    plt.show() 

    output_file = filepaths.SaveLocation
    plt.savefig(output_file)

    print('    ...Simulation complete! ✨')

    print('                                                         ')

    print('You can view the output image at: ' + filepaths.SaveLocation )

    print('                                                         ')

    print(' ✨ ✨ ✨ ✨ ✨ ✨ ✨ ✨ ✨ ✨ ✨ ✨ ✨ ✨ ✨ ✨ ✨ ✨ ✨ ✨ ')
    print('_____________________________________________________________')


def profiled_run():
    RunSimulation('CPU')

def main():
    start_time = time.time()
    
    RunSimulation('CPU')
    
    end_time = time.time()
    elapsed_time = end_time - start_time
    
    print(f"Simulation on 'CPU' mode took {elapsed_time:.2f} seconds.")





if __name__ == "__main__":
    main()


