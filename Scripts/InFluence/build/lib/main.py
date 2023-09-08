
import itertools
import sys
import time

from matplotlib import pyplot as plt
import numpy as np

from filepaths import Filepaths as filepaths
from image_preprocessing.image_object_constructor import image_pixels

from simulation.simulation_interface import InFluenceSimulation 
from parameters import SimulationParameters


def RunSimulation(simulation_type, gui=False):
    
    params = SimulationParameters()
    params.refresh_config()
    image_object = image_pixels(filepaths.InputImageLocation, params)
    influence_simulator = InFluenceSimulation(image_object, params)

    if simulation_type =='GPU':
        SimulationResults = influence_simulator.Run_GPU_KERNEL_MCScatteringSimulation()
    elif simulation_type == "CPP":
        SimulationResults = influence_simulator.RunCPPMCScatteringSimulation()
    elif simulation_type == "PySparse":
        SimulationResults = influence_simulator.RunMCSimulation_Py("PySparse")
    elif simulation_type == "PyDense":
        SimulationResults = influence_simulator.RunMCSimulation_Py("PyDense")
    elif simulation_type == "Plot":
        SimulationResults = influence_simulator.RunMCSimulation_Py("Plot")


    if gui==True:
        print('Simulation complete.')
        return SimulationResults
    
    else:
        
                
        plt.imshow(np.uint32(SimulationResults))
        plt.axis('off')
        plt.show() 
        

        output_file = filepaths.SaveLocation
        plt.savefig(output_file)
        
        # Save as an image
        plt.savefig(output_file + '.png')

        # Save as a .npy file
        np.save(output_file + '.npy', SimulationResults)

        return
        

import time  # Import the time module

def main():
    # Record the start time
    start_time = time.time()

    # Run the simulation
    RunSimulation('PySparse')

    # Record the end time
    end_time = time.time()

    # Calculate the elapsed time
    elapsed_time = end_time - start_time
    

    print(f"Simulation completed in {elapsed_time:.2f} seconds.")  # Print the elapsed time


if __name__ == "__main__":
    try:
        main()
    except ValueError as e:
        print(e)  # Print the error message if an exception is raised
        sys.exit(1)  # Exit the program with an error code (1 indicates an error)
