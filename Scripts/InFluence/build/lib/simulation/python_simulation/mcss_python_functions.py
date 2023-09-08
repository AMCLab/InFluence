
import math
import random
import numpy as np
from tqdm import tqdm
from multiprocessing import Pool, cpu_count
from joblib import Parallel, delayed
from numba.typed import Dict
from numba.types import UniTuple, int64, float64


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
from threading import Lock
from scipy.sparse import csr_matrix, lil_matrix



@jit(nopython=True, parallel=True)
def MonteCarloScatteringSimulation_function(
    pixels, 
    E_i, 
    dE_threshold, 
    MinimumEnergy, 
    Density, 
    t_counting, 
    N, 
    ProbeDiameter, 
    image_shape_0, 
    image_shape_1
    ):

    
    ModulatedImage = np.zeros((image_shape_0, image_shape_1), dtype=np.float64)
    local_random_uniform = random.uniform

    print(len(pixels))

    for i in prange(len(pixels)):
        
        pixel = pixels[i]
        for j in prange(pixel.electron_count):
            
            pixel_x_dimension = pixel.x_dimension
            pixel_y_dimension = pixel.y_dimension
            pixel_z_dimension = pixel.z_dimension
            pixel_i = pixel.i
            pixel_j = pixel.j

            electron = ElectronEvolver(
                E_i, ProbeDiameter, Density)
            
            ElectronHoleChargeCounter = np.zeros((image_shape_0, image_shape_1), dtype = np.float64)
    

            while electron.condition:  
                                                               
            
                electron.update_state()
                
                electron_x0 = electron.x0
                electron_y0 = electron.y0
                electron_z0 = electron.z0
                number_stopped=0
                number_backscattered =0
                number_transmitted=0

                # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- #

                if electron.E <= MinimumEnergy:                                             #if electron stops in the material
                    number_stopped = number_stopped + 1
                    electron.condition = False
                    
                if electron_z0 < 0.01:                                                                   #if electron backscatters
                    number_backscattered = number_backscattered + 1
                    electron.condition = False                                                             
                    
                if electron_z0 > pixel_z_dimension:                                                        # if electron penetrates material #change to have de/ds as the counting
                                                                                        
                    number_transmitted = number_transmitted + 1
                    electron.condition = False             

                # ------------------------------------------------------------------------------ Electron Deposits Sufficient Energy ------------------------------------------------------------------------------ #    
                            
                if -1*electron.dE >= dE_threshold:                                                            
                    new_eh_pairs = math.floor(-1*electron.dE/dE_threshold)
                    
                # ------------------------------------------------------------------------------ Electron Stays Within Pixel Boundary ------------------------------------------------------------------------------ #    
                    
                    if (electron_x0 <= 1*pixel_x_dimension) and (electron_x0 >= -1*pixel_x_dimension) and (electron_y0 <= 1*pixel_y_dimension)and (electron_y0 >=-1*pixel_y_dimension):                                        #electron stays within pixel region

                        ElectronHoleChargeCounter[pixel_i, pixel_j] += new_eh_pairs 

    
                # ------------------------------------------------------------------------------ Electron Moves Positively in X and Y ------------------------------------------------------------------------------ #
                    
                    elif (electron_x0 > 1*pixel_x_dimension) and (electron_y0 > 1*pixel_y_dimension):                   
                        translation_x = custom_round_down(x = electron_x0/(pixel_x_dimension))
                        translation_y = custom_round_down(x = electron_y0/(2*pixel_y_dimension))
                        
                        if pixel_i + translation_x  <= image_shape_0 - 1 and pixel_j + translation_y <= image_shape_0 - 1:
                                ElectronHoleChargeCounter[pixel_i + translation_x, pixel_j + translation_y] += new_eh_pairs 

                # ------------------------------------------------------------------------------ Electron Moves Negatively in X and Y ------------------------------------------------------------------------------ #

                    elif (electron_x0 < -1*pixel_x_dimension) and (electron_y0 < -1*pixel_y_dimension):                 #electron moves negatively in x
                            # and y
                        translation_x = custom_round_down(x=electron_x0 / (2 * pixel_x_dimension))
                        translation_y = custom_round_down(x=electron_y0 / (2 * pixel_y_dimension))
                        if pixel_i - translation_x >= 0 and pixel_j - translation_y >= 0:
                            ElectronHoleChargeCounter[pixel_i - translation_x, pixel_j - translation_y] += new_eh_pairs 

                # ----------------------------------------------------------------------- Electron Moves Positively in X and Negatively in Y ------------------------------------------------------------------------ #


                    elif electron_x0 > (1*pixel_x_dimension) and electron_y0 < (-1*pixel_y_dimension):  
                        
                        translation_x = custom_round_down(x=electron_x0 / (2 * pixel_x_dimension))
                        translation_y = custom_round_down(x=electron_y0 / (2 * pixel_y_dimension))
                        if pixel_j - translation_y >= 0 and pixel_i + translation_x  <= image_shape_0 - 1:

                            ElectronHoleChargeCounter[pixel_i + translation_x, pixel_j - translation_y] += new_eh_pairs 


                # ----------------------------------------------------------------------- Electron Moves Negatively in X and Positively in Y ------------------------------------------------------------------------ #                         


                    elif electron_x0 < (-1*pixel_x_dimension) and electron_y0 > (1*pixel_y_dimension):  

                        translation_x = custom_round_down(x=electron_x0 / (2 * pixel_x_dimension))
                        translation_y = custom_round_down(x=electron_y0 / (2 * pixel_y_dimension))
                        if pixel_i - translation_x >= 0 and pixel_j + translation_y <= image_shape_0 - 1:
                            ElectronHoleChargeCounter[pixel_i - translation_x, pixel_j + translation_y] += new_eh_pairs 


                # ------------------------------------------------------------------------------ Electron Moves Positively in X Direction ----------------------------------------------------------------------------- # 


                    elif electron_x0 > (1*pixel_x_dimension):                                               
                        translation_x = custom_round_down(x=electron_x0 / (2 * pixel_x_dimension))
                        if pixel_i + \
                            translation_x <= image_shape_0 - 1:

                            ElectronHoleChargeCounter[pixel_i + translation_x, pixel_j] += new_eh_pairs 


                # ------------------------------------------------------------------------------ Electron Moves in Negative in X Direction ----------------------------------------------------------------------------- # 


                    elif electron_x0 < (-1*pixel_x_dimension):                                          
                        translation_x = custom_round_down(x=electron_x0 / (2 * pixel_x_dimension))
                        if pixel_i - translation_x >= 0:

                            ElectronHoleChargeCounter[pixel_i - translation_x, pixel_j] += new_eh_pairs 


                # --------------------------------------------------------------------------------- Electron Moves in Positive Y Direction ------------------------------------------------------------------------------ # 
               

                    elif electron_y0 > (1*pixel_y_dimension):                      
                        translation_y = custom_round_down(x=electron_y0 / (2 * pixel_y_dimension))
                        if pixel_j + \
                            translation_y <= image_shape_0 - 1:

                            ElectronHoleChargeCounter[pixel_i, pixel_j + translation_y] += new_eh_pairs 

                # ---------------------------------------------------------------------------------- Electron Moves in Negative Y Direction ------------------------------------------------------------------------------- # 

                    elif electron_y0 < (-1*pixel_x_dimension): #electron moves in negative y direction
                        translation_y = custom_round_down(x=electron_y0 / (2 * pixel_y_dimension))
                        
                        if pixel_j - translation_y >= 0:

                            ElectronHoleChargeCounter[pixel_i, pixel_j - translation_y] += new_eh_pairs 

                # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- #


            ElectronHoleChargeCounter = np.floor(dE_threshold * ElectronHoleChargeCounter / t_counting)

            ModulatedImage += ElectronHoleChargeCounter 
    
    return ModulatedImage



