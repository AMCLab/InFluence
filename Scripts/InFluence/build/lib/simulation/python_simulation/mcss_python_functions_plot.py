 
import math
import random
import numpy as np
from tqdm import tqdm
from multiprocessing import Pool, cpu_count
from joblib import Parallel, delayed
from numba.typed import Dict
from numba.types import UniTuple, int64, float64
import numba
 

from simulation.python_simulation.helper_functions import (
    custom_round_down, evaluate_alpha, evaluate_cross_section_opt,
    evaluate_path_length, evaluate_step, evaluate_phi,
    evaluate_pho, evaluate_direction_cosine_a, evaluate_direction_cosine_b,
    evaluate_direction_cosine_c, evaluate_energy_loss_rate, initialise_postions
)
from simulation.python_simulation.electron_object_constructor import ElectronEvolver

 
from numba import njit, prange, jit, float64, int32, typed, typeof
from numba.experimental import jitclass
from numba.typed import List




def MonteCarloScatteringSimulationPlot_function(
    pixels, 
    E_i, 
    dE_threshold, 
    MinimumEnergy, 
    Density, 
    t_counting, 
    N, 
    ProbeDiameter, 
    image_shape_0, 
    image_shape_1):
   

    #counters
    number_transmitted = 0
    number_backscattered = 0
    number_stopped = 0
    vector_coordinates_x = List([np.float64(0)])
    vector_coordinates_y = List([np.float64(0)])
    vector_coordinates_z = List([np.float64(0)])
    points_per_trajectory = List([np.float64(0)])
    
    ModulatedImage_flat = np.zeros(image_shape_0 * image_shape_1, dtype=np.float64)
    local_math_floor = math.floor
    local_numpy_floor = np.floor

    ModulatedImage_total = np.zeros(image_shape_0 * image_shape_1, dtype=np.float64)
    
    counter = 0

    for i in range(len(pixels)):
        
        counter +=1
        print(counter)


        pixel = pixels[i]
        for j in range(pixel.electron_count):
            
            point_counter = 1
           
            pixel_x_dimension = pixel.x_dimension
            pixel_y_dimension = pixel.y_dimension
            pixel_z_dimension = pixel.z_dimension
            pixel_i = pixel.i
            pixel_j = pixel.j
 
            electron = ElectronEvolver(E_i, ProbeDiameter, Density)
            
            vector_coordinates_x.append(electron.x0)
            vector_coordinates_y.append(electron.y0)
            vector_coordinates_z.append(electron.z0)
 
            ElectronHoleChargeCounter_index = List([np.float64(0)])
            ElectronHoleChargeCounter_values = List([np.float64(0)])
   
 
            while electron.condition:  # Monte Carlo Loop until backscatter, trasmission or stopping
                
                # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- #
 
                if electron.E <= MinimumEnergy:                                             #if electron stops in the material
                    number_stopped = number_stopped + 1
                    electron.condition = False
                   
                if electron.z0 < 0.01:                                                                   #if electron backscatters
                    number_backscattered = number_backscattered + 1
                    electron.condition = False                                                              #not sure how to deal with electron scattering outside of material
                   
                if electron.z0 > pixel_z_dimension:                                                        # if electron penetrates material #change to have de/ds as the counting
                                                                                       
                    number_transmitted = number_transmitted + 1
                    electron.condition = False            
                                                              
            
                electron.update_state()
               

                number_stopped=0
                number_backscattered =0
                number_transmitted=0
                
                vector_coordinates_x.append(electron.x0)
                vector_coordinates_y.append(electron.y0)
                vector_coordinates_z.append(electron.z0)
                point_counter += 1
 
            points_per_trajectory.append(point_counter)
            
               
    return np.reshape(ModulatedImage_flat, (image_shape_0, image_shape_1)), vector_coordinates_x, vector_coordinates_y, vector_coordinates_z, points_per_trajectory




 