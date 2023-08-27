
import math
import random
import numpy as np
from tqdm import tqdm
from multiprocessing import Pool, cpu_count
from joblib import Parallel, delayed
from numba.typed import Dict
from numba.types import UniTuple, int64, float64

from parameters import params
from HelperFunctions import (
    custom_round_down, evaluate_alpha, evaluate_cross_section_opt, 
    evaluate_path_length, evaluate_step, evaluate_phi, 
    evaluate_pho, evaluate_direction_cosine_a, evaluate_direction_cosine_b, 
    evaluate_direction_cosine_c, evaluate_energy_loss_rate, initialise_postions
)
from ElectronObjectConstructor import ElectronEvolver
from PlotTrajectoryFunctions import PlotTrajectories
from PixelObjectConstructor import pixel_constructor

from numba import njit, prange, jit, float64, int32, typed, typeof
from numba.experimental import jitclass
from numba.typed import List
from threading import Lock
from scipy.sparse import csr_matrix, lil_matrix




@jit(nopython=True)
def MonteCarloScatteringSimulationPlot_function(pixels, E_i, dE_threshold, MinimumEnergy, ProtonNum, AtomicMass, Density, t_counting, N, ProbeDiameter, AlphaMultiplier, CrossSectionNumorator, CrossSectionLogArgMultiplier, CrossSectionDenominatorA, CrossSectionDenominatorB, PathLengthMultiplier, EnergyLossMultiplierA, EnergyLossMultiplierB, new_image_MCS, image_shape_0, image_shape_1):
    
    ModulatedImage = np.zeros((image_shape_0, image_shape_1), dtype=np.float64)
    local_random_uniform = random.uniform

            # Counters
    number_transmitted = 0
    number_backscattered = 0
    number_stopped = 0
    vector_coordinates_x = List([np.float64(0)])
    vector_coordinates_y = List([np.float64(0)])
    vector_coordinates_z = List([np.float64(0)])
    points_per_trajectory = List([np.float64(0)])


    for i in prange(len(pixels)):


        pixel = pixels[i]
        for j in prange(pixel.electron_count):

            point_counter = 1
            
            pixel_x_dimension = pixel.x_dimension
            pixel_y_dimension = pixel.y_dimension
            pixel_z_dimension = pixel.z_dimension
            pixel_i = pixel.i
            pixel_j = pixel.j

            electron = ElectronEvolver(E_i, pixel_x_dimension, Density)
            ElectronHoleChargeCounter = np.zeros((image_shape_0, image_shape_1), dtype = np.float64)
    

            while electron.condition:  # Monte Carlo Loop until backscatter, trasmission or stopping
                                                               
            
                electron.update_state
                
                electron_x0 = electron.x0
                electron_y0 = electron.y0
                electron_z0 = electron.z0



                vector_coordinates_x.append(electron_x0)
                vector_coordinates_y.append(electron_y0)
                vector_coordinates_z.append(electron_z0)
                point_counter += 1
    

                number_stopped=0
                number_backscattered =0
                number_transmitted=0

                # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- #

                if electron.E <= MinimumEnergy:                                             #if electron stops in the material
                    number_stopped = number_stopped + 1
                    electron.condition = False
                    
                if electron_z0 < 0.01:                                                                   #if electron backscatters
                    number_backscattered = number_backscattered + 1
                    electron.condition = False                                                              #not sure how to deal with electron scattering outside of material
                    
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

            points_per_trajectory.append(point_counter)
            ElectronHoleChargeCounter = np.floor(dE_threshold * ElectronHoleChargeCounter / t_counting)
            new_image_MCS += ElectronHoleChargeCounter 

    
    
    return new_image_MCS, vector_coordinates_x, vector_coordinates_y, vector_coordinates_z, points_per_trajectory



