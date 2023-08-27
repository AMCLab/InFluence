 
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
from PixelObjectConstructor import pixel_constructor
 
from numba import njit, prange, jit, float64, int32, typed, typeof
from numba.experimental import jitclass
from numba.typed import List
from threading import Lock
from scipy.sparse import csr_matrix, lil_matrix




#@jit(nopython=True)

 

@jit(nopython=True, parallel = True)
def MonteCarloScatteringSimulationSparse_function(pixels, E_i, dE_threshold, MinimumEnergy, ProtonNum, AtomicMass, Density, t_counting, N, ProbeDiameter, AlphaMultiplier, CrossSectionNumorator, CrossSectionLogArgMultiplier, CrossSectionDenominatorA, CrossSectionDenominatorB, PathLengthMultiplier, EnergyLossMultiplierA, EnergyLossMultiplierB, image_shape_0, image_shape_1):
   
    print(len(pixels))

    ModulatedImage_flat = np.zeros(image_shape_0 * image_shape_1, dtype=np.float64)
    local_random_uniform = random.uniform
    local_math_floor = math.floor
    local_numpy_floor = np.floor

    ModulatedImage_total = np.zeros(image_shape_0 * image_shape_1, dtype=np.float64)

    for i in prange(len(pixels)):


       
        pixel = pixels[i]
        for j in prange(pixel.electron_count):
           
            pixel_x_dimension = pixel.x_dimension
            pixel_y_dimension = pixel.y_dimension
            pixel_z_dimension = pixel.z_dimension
            pixel_i = pixel.i
            pixel_j = pixel.j
 
            electron = ElectronEvolver(E_i, pixel_x_dimension, Density)
 
            ElectronHoleChargeCounter_index = List([np.float64(0)])
            ElectronHoleChargeCounter_values = List([np.float64(0)])
   
 
            while electron.condition:  # Monte Carlo Loop until backscatter, trasmission or stopping
                                                              
            
                electron.update_state()
               
                electron_x0 = electron.x0
                electron_y0 = electron.y0
                electron_z0 = electron.z0
                electron_dE = electron.dE
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
                            
                if -1*electron_dE >= dE_threshold:                                                           
                    new_eh_pairs = local_math_floor(-1*electron_dE/dE_threshold)
                   
                # ------------------------------------------------------------------------------ Electron Stays Within Pixel Boundary ------------------------------------------------------------------------------ #   
                    
                    if (electron_x0 <= 1*pixel_x_dimension) and (electron_x0 >= -1*pixel_x_dimension) and (electron_y0 <= 1*pixel_y_dimension)and (electron_y0 >=-1*pixel_y_dimension):                                        #electron stays within pixel region
 
                        pixel_1D_index = pixel_i*image_shape_1+pixel_j
 
                        if pixel_1D_index in ElectronHoleChargeCounter_index:
 
                            i = ElectronHoleChargeCounter_index.index(pixel_1D_index)
                            ElectronHoleChargeCounter_values[i] += new_eh_pairs
                       
                        else:
 
                            ElectronHoleChargeCounter_index.append(pixel_1D_index)
                            ElectronHoleChargeCounter_values.append(new_eh_pairs)
   
                # ------------------------------------------------------------------------------ Electron Moves Positively in X and Y ------------------------------------------------------------------------------ #
                    
                    elif (electron_x0 > 1*pixel_x_dimension) and (electron_y0 > 1*pixel_y_dimension):                  
                        translation_x = custom_round_down(x = electron_x0/(pixel_x_dimension))
                        translation_y = custom_round_down(x = electron_y0/(2*pixel_y_dimension))
                       
                        if pixel_i + translation_x  <= image_shape_0 - 1 and pixel_j + translation_y <= image_shape_0 - 1:
 
                            pixel_1D_index = (pixel_i - translation_x) * image_shape_1 + (pixel_j - translation_y)
 
                            if pixel_1D_index in ElectronHoleChargeCounter_index:
 
                                i = ElectronHoleChargeCounter_index.index(pixel_1D_index)
                                ElectronHoleChargeCounter_values[i] += new_eh_pairs
                           
                            else:
 
                                ElectronHoleChargeCounter_index.append(pixel_1D_index)
                                ElectronHoleChargeCounter_values.append(new_eh_pairs)
 
                # ------------------------------------------------------------------------------ Electron Moves Negatively in X and Y ------------------------------------------------------------------------------ #
 
                    elif (electron_x0 < -1*pixel_x_dimension) and (electron_y0 < -1*pixel_y_dimension):                 #electron moves negatively in x
                            # and y
                        translation_x = custom_round_down(x=electron_x0 / (2 * pixel_x_dimension))
                        translation_y = custom_round_down(x=electron_y0 / (2 * pixel_y_dimension))
                        if pixel_i - translation_x >= 0 and pixel_j - translation_y >= 0:
 
                            pixel_1D_index =  (pixel_i - translation_x) * image_shape_1 + (pixel_j - translation_y)
 
                            if pixel_1D_index in ElectronHoleChargeCounter_index:
 
                                i = ElectronHoleChargeCounter_index.index(pixel_1D_index)
                                ElectronHoleChargeCounter_values[i] += new_eh_pairs
                           
                            else:
 
                                ElectronHoleChargeCounter_index.append(pixel_1D_index)
                                ElectronHoleChargeCounter_values.append(new_eh_pairs)
 
                # ----------------------------------------------------------------------- Electron Moves Positively in X and Negatively in Y ------------------------------------------------------------------------ #
 
                    elif electron_x0 > (1*pixel_x_dimension) and electron_y0 < (-1*pixel_y_dimension): 
                        
                        translation_x = custom_round_down(x=electron_x0 / (2 * pixel_x_dimension))
                        translation_y = custom_round_down(x=electron_y0 / (2 * pixel_y_dimension))
                        if pixel_j - translation_y >= 0 and pixel_i + translation_x  <= image_shape_0 - 1:
 
                            pixel_1D_index = (pixel_i + translation_x) * image_shape_1 + (pixel_j - translation_y)
 
                            if pixel_1D_index in ElectronHoleChargeCounter_index:
 
                                i = ElectronHoleChargeCounter_index.index(pixel_1D_index)
                                ElectronHoleChargeCounter_values[i] += new_eh_pairs
                           
                            else:
 
                                ElectronHoleChargeCounter_index.append(pixel_1D_index)
                                ElectronHoleChargeCounter_values.append(new_eh_pairs)
 
                # ----------------------------------------------------------------------- Electron Moves Negatively in X and Positively in Y ------------------------------------------------------------------------ #                        
 

                    elif electron_x0 < (-1*pixel_x_dimension) and electron_y0 > (1*pixel_y_dimension): 
 
                        translation_x = custom_round_down(x=electron_x0 / (2 * pixel_x_dimension))
                        translation_y = custom_round_down(x=electron_y0 / (2 * pixel_y_dimension))
                        if pixel_i - translation_x >= 0 and pixel_j + translation_y <= image_shape_0 - 1:
 
                            pixel_1D_index = (pixel_i - translation_x) * image_shape_1 + (pixel_j + translation_y)
 
                            if pixel_1D_index in ElectronHoleChargeCounter_index:
 
                                i = ElectronHoleChargeCounter_index.index(pixel_1D_index)
                                ElectronHoleChargeCounter_values[i] += new_eh_pairs
                           
                            else:
 
                                ElectronHoleChargeCounter_index.append(pixel_1D_index)
                                ElectronHoleChargeCounter_values.append(new_eh_pairs)                          
 
                # ------------------------------------------------------------------------------ Electron Moves Positively in X Direction ----------------------------------------------------------------------------- #
 
                    elif electron_x0 > (1*pixel_x_dimension):                                              
                        translation_x = custom_round_down(x=electron_x0 / (2 * pixel_x_dimension))
                        if pixel_i + \
                            translation_x <= image_shape_0 - 1:
 
                            pixel_1D_index = (pixel_i + translation_x) * image_shape_1 + pixel_j
 
                            if pixel_1D_index in ElectronHoleChargeCounter_index:
 
                                i = ElectronHoleChargeCounter_index.index(pixel_1D_index)
                                ElectronHoleChargeCounter_values[i] += new_eh_pairs
                            
                            else:
 
                                ElectronHoleChargeCounter_index.append(pixel_1D_index)
                                ElectronHoleChargeCounter_values.append(new_eh_pairs) 
 
                # ------------------------------------------------------------------------------ Electron Moves in Negative in X Direction ----------------------------------------------------------------------------- #
 
                    elif electron_x0 < (-1*pixel_x_dimension):                                          
                        translation_x = custom_round_down(x=electron_x0 / (2 * pixel_x_dimension))
                        if pixel_i - translation_x >= 0:
 
                            pixel_1D_index = (pixel_i - translation_x) * image_shape_1 + pixel_j
                           
                            if pixel_1D_index in ElectronHoleChargeCounter_index:
 
                                i = ElectronHoleChargeCounter_index.index(pixel_1D_index)
                                ElectronHoleChargeCounter_values[i] += new_eh_pairs
                           
                            else:
 
                                ElectronHoleChargeCounter_index.append(pixel_1D_index)
                                ElectronHoleChargeCounter_values.append(new_eh_pairs) 
 
                # --------------------------------------------------------------------------------- Electron Moves in Positive Y Direction ------------------------------------------------------------------------------ #
               
 
                    elif electron_y0 > (1*pixel_y_dimension):                     
                        translation_y = custom_round_down(x=electron_y0 / (2 * pixel_y_dimension))
                        if pixel_j + \
                            translation_y <= image_shape_0 - 1:
 
                            pixel_1D_index = pixel_i * image_shape_1 + (pixel_j + translation_y)
                            if pixel_1D_index in ElectronHoleChargeCounter_index:
 
                                i = ElectronHoleChargeCounter_index.index(pixel_1D_index)
                                ElectronHoleChargeCounter_values[i] += new_eh_pairs
                           
                            else:
 
                                ElectronHoleChargeCounter_index.append(pixel_1D_index)
                                ElectronHoleChargeCounter_values.append(new_eh_pairs) 
 
                # ---------------------------------------------------------------------------------- Electron Moves in Negative Y Direction ------------------------------------------------------------------------------- #
 
                    elif electron_y0 < (-1*pixel_x_dimension): #electron moves in negative y direction
                        translation_y = custom_round_down(x=electron_y0 / (2 * pixel_y_dimension))
                       
                        if pixel_j - translation_y >= 0:
 
                            pixel_1D_index = pixel_i * image_shape_1 + (pixel_j - translation_y)
 
                            if pixel_1D_index in ElectronHoleChargeCounter_index:
 
                                i = ElectronHoleChargeCounter_index.index(pixel_1D_index)
                                ElectronHoleChargeCounter_values[i] += new_eh_pairs
                           
                            else:
 
                                ElectronHoleChargeCounter_index.append(pixel_1D_index)
                                ElectronHoleChargeCounter_values.append(new_eh_pairs) 
 

                # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- #
 
            for k in prange(len(ElectronHoleChargeCounter_index)):
                idx = int(ElectronHoleChargeCounter_index[k])
                value = ElectronHoleChargeCounter_values[k]
                ModulatedImage_flat[idx] +=  local_numpy_floor(dE_threshold * value / t_counting)
            
               
    return ModulatedImage_flat




 