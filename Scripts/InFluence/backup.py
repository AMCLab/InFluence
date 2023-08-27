import math
import numpy as np
import random
from tqdm import tqdm
from parameters import params
from ElectronObjectConstructor import ElectronObject
from HelperFunctions import custom_round_down, CalculateConstants, evaluate_alpha, evaluate_cross_section_opt, evaluate_path_length, evaluate_step, evaluate_phi, evaluate_pho, evaluate_direction_cosine_a, evaluate_direction_cosine_b, evaluate_direction_cosine_c, evaluate_energy_loss_rate, initialise_postions
import multiprocessing
import numba 
from numba import njit, prange, jit
from numba import int32, float32    # import the types
from numba.experimental import jitclass

import multiprocessing
from multiprocessing import Pool



@njit
def _MCScatteringSimulation_SingleElectron(pixel, image_object, E_i, dE_threshold, MinimumEnergy, ProtonNum, AtomicMass, Density, t_counting, N, ProbeDiameter, AlphaMultiplier, CrossSectionNumorator, CrossSectionLogArgMultiplier, CrossSectionDenominatorA, CrossSectionDenominatorB, PathLengthMultiplier, EnergyLossMultiplierA, EnergyLossMultiplierB, new_image_MCS):
    electron = ElectronObject(image_object, E_i, dE_threshold, MinimumEnergy, ProtonNum, AtomicMass, Density, t_counting, N, ProbeDiameter)

    ElectronHoleChargeCounter = np.zeros((image_object.perfect_image.shape[0], image_object.perfect_image.shape[1]), dtype = np.float64)

    eh_charge_counter = np.zeros(image_object.perfect_image.shape, dtype=np.float64)


    while electron.condition:                                                                   # Monte Carlo Loop until backscatter, trasmission or stopping
            


        RND_phi = random.uniform(a=0, b=1)                                                                 # generate random number for the phi angle
        RND_step = random.uniform(a = 0.000001, b = 0.999999)
        RND_pho = random.uniform(a = 0, b = 1)
        alpha = evaluate_alpha(electron.E, AlphaMultiplier)                                                            # calc screening constant, function of previous energy
        cross = evaluate_cross_section_opt(electron.E, CrossSectionLogArgMultiplier, CrossSectionNumorator, CrossSectionDenominatorA, CrossSectionDenominatorB)                                             # calc cross section
        path_length = evaluate_path_length(electron.CrossSection, PathLengthMultiplier)                                      # calc mean free path length
        step = evaluate_step(path_length, RND_step)                             # calculate step of this iteration
        dE = step * Density * evaluate_energy_loss_rate(electron.E, EnergyLossMultiplierA, EnergyLossMultiplierB)
        electron.E = electron.E + dE                                                                              # calc new energy #separate out dE/dS
        phi = evaluate_phi(RND_phi, alpha=alpha)                                            # calc scattering angle
        psi = evaluate_pho(RND_pho, math.pi)                                                         # calc other scattering angle
        ca = evaluate_direction_cosine_a(phi, psi, electron.cosineX, electron.cosineY, electron.cosineZ)                                  # calc direction cosines
        cb = evaluate_direction_cosine_b(phi, psi, electron.cosineX, electron.cosineY, electron.cosineZ)
        cc = evaluate_direction_cosine_c(phi, psi, electron.cosineZ)
        electron.x0 = electron.x0 + step * ca                                                                     # find and reset to new positions
        electron.y0 = electron.y0 + step * cb
        electron.z0 = electron.z0 + step * cc
        electron.cosineX = ca                                                                                 # reset direction cosines
        electron.cosineY = cb
        electron.cosineZ = cc

        number_stopped=0

        if electron.E <= MinimumEnergy:                                             #if electron stops in the material
            number_stopped = number_stopped + 1
            electron.condition = False
            
        if electron.z0 < 0.01:                                                                   #if electron backscatters
            number_backscattered = number_backscattered + 1
            electron.condition = False                                                              #not sure how to deal with electron scattering outside of material
            
        if electron.z0 > pixel.x_dimension:                                                        # if electron penetrates material #change to have de/ds as the counting
                                                                                
            number_transmitted = number_transmitted + 1
            electron.condition = False
                    
        if -1*dE >= dE_threshold:                                                            #if electron deposits sufficient energyo
            new_eh_pairs = math.floor(-1*dE/dE_threshold)
            number_eh_pairs = int(number_eh_pairs + new_eh_pairs)

            if (electron.x0 <= 1*pixel.x_dimension) and (electron.x0 >= -1*pixel.x_dimension) and (electron.y0 <= 1*pixel.y_dimension)and (electron.y0 >=-1*pixel.y_dimension):                                        #electron stays within pixel region

                eh_charge_counter[pixel.i, pixel.j] = new_eh_pairs + eh_charge_counter[pixel.i, pixel.j]

            elif (electron.x0 > 1*pixel.x_dimension) and (electron.y0 > 1*pixel.y_dimension):                   #electron moves positively in x
                    # and y
                translation_x = custom_round_down(x = electron.x0/(pixel.x_dimension))
                translation_y = custom_round_down(x = electron.y0/(2*pixel.y_dimension))
                if pixel.i + translation_x  <= image_object.perfect_image.shape[0] - 1 and pixel.j + \
                        translation_y <= image_object.perfect_image.shape[0] - 1:
                        eh_charge_counter[pixel.i + translation_x, pixel.j + translation_y] = new_eh_pairs + eh_charge_counter[pixel.i +translation_x, pixel.j+translation_y]

            elif (electron.x0 < -1*pixel.x_dimension) and (electron.y0 < -1*pixel.y_dimension):                 #electron moves negatively in x
                    # and y
                translation_x = custom_round_down(x=electron.x0 / (2 * pixel.x_dimension))
                translation_y = custom_round_down(x=electron.y0 / (2 * pixel.y_dimension))
                if pixel.i - translation_x >= 0 and pixel.j - translation_y >= 0:
                    eh_charge_counter[pixel.i - translation_x, pixel.j - translation_y] = new_eh_pairs + eh_charge_counter[pixel.i -translation_x, pixel.j -translation_y]


            elif electron.x0 > (1*pixel.x_dimension) and electron.y0 < (-1*pixel.y_dimension):  # electron moves psoitively in
                    # x and negatively in y
                translation_x = custom_round_down(x=electron.x0 / (2 * pixel.x_dimension))
                translation_y = custom_round_down(x=electron.y0 / (2 * pixel.y_dimension))
                if pixel.j - translation_y >= 0 and pixel.i + translation_x  <= \
                    image_object.perfect_image.shape[0] - 1:
                    eh_charge_counter[pixel.i + translation_x, pixel.j - translation_y] = new_eh_pairs + eh_charge_counter[pixel.i +translation_x,pixel.j -translation_y]


            elif electron.x0 < (-1*pixel.x_dimension) and electron.y0 > (1*pixel.y_dimension):  # electron moves negatively in
                            # x and positively in y
                translation_x = custom_round_down(x=electron.x0 / (2 * pixel.x_dimension))
                translation_y = custom_round_down(x=electron.y0 / (2 * pixel.y_dimension))
                if pixel.i - translation_x >= 0 and pixel.j + \
                    translation_y <= image_object.perfect_image.shape[0] - 1:
                    eh_charge_counter[pixel.i - translation_x, pixel.j + translation_y] = new_eh_pairs + eh_charge_counter[pixel.i -translation_x,pixel.j +translation_y]


            elif electron.x0 > (1*pixel.x_dimension):                                               #electron moves in positive x direction
                translation_x = custom_round_down(x=electron.x0 / (2 * pixel.x_dimension))
                if pixel.i + \
                    translation_x <= image_object.perfect_image.shape[0] - 1:

                    eh_charge_counter[pixel.i + translation_x, pixel.j] = new_eh_pairs + eh_charge_counter[pixel.i +translation_x, pixel.j]

            elif electron.x0 < (-1*pixel.x_dimension):                                          #electron moves in negative x direction
                translation_x = custom_round_down(x=electron.x0 / (2 * pixel.x_dimension))
                if pixel.i - translation_x >= 0:

                    eh_charge_counter[pixel.i - translation_x, pixel.j] = new_eh_pairs + eh_charge_counter[pixel.i -translation_x,pixel.j]

            elif electron.y0 > (1*pixel.y_dimension):                                       #electron moves in positive y direction
                translation_y = custom_round_down(x=electron.y0 / (2 * pixel.y_dimension))
                if pixel.j + \
                    translation_y <= image_object.perfect_image.shape[0] - 1:

                    eh_charge_counter[pixel.i, pixel.j + translation_y] = new_eh_pairs + eh_charge_counter[pixel.i, pixel.j +translation_y]

            elif electron.y0 < (-1*pixel.x_dimension): #electron moves in negative y direction
                translation_y = custom_round_down(x=electron.y0 / (2 * pixel.y_dimension))

                if pixel.j - translation_y >= 0:

                    eh_charge_counter[pixel.i, pixel.j - translation_y] = new_eh_pairs + eh_charge_counter[pixel.i, pixel.j -translation_y]

    eh_charge_counter = (np.floor(dE_threshold*eh_charge_counter/t_counting))
    new_image_MCS += eh_charge_counter



class MonteCarloScatteringSimulation:
    
    # Constants
    AvagadrosConstant = params.N  # Assumes 'params' is globally available
    pi = math.pi

    # Shared counters - shared memory (if they are truly shared across instances)
    number_transmitted = 0
    number_backscattered = 0
    number_stopped = 0
    number_eh_pairs = 0
    new_image_MCS = 0
    COUNTER = 0

    def __init__(self, image_object, E_i, dE_threshold, MinimumEnergy, ProtonNum, AtomicMass, Density, t_counting, N):
        self.image_object = image_object
        self.E_i = E_i
        self.dE_threshold = dE_threshold
        self.MinimumEnergy = MinimumEnergy
        self.ProtonNum = ProtonNum
        self.AtomicMass = AtomicMass
        self.Density = Density
        self.t_counting = t_counting
        self.N = N
        self.ProbeDiameter = image_object.pixel_dimensions[0]
        
        (self.AlphaMultiplier, self.CrossSectionNumorator, self.CrossSectionLogArgMultiplier, 
         self.CrossSectionDenominatorA, self.CrossSectionDenominatorB, 
         self.PathLengthMultiplier, self.EnergyLossMultiplierA, 
         self.EnergyLossMultiplierB) = CalculateConstants(params)
        
        self.new_image_MCS = np.zeros(self.image_object.perfect_image.shape, dtype=np.float64)
        
        
                                                                                    
    
    def MCScatteringSimulation_SingleElectron(self, pixel):
        
        electron = ElectronObject(self.image_object, self.E_i, self.dE_threshold, self.MinimumEnergy, self.ProtonNum, self.AtomicMass, self.Density, self.t_counting, self.N, self.ProbeDiameter)

        ElectronHoleChargeCounter = np.zeros((self.image_object.perfect_image.shape[0], self.image_object.perfect_image.shape[1]), dtype = np.float64)

        eh_charge_counter = np.zeros(self.image_object.perfect_image.shape, dtype=np.float64)


        while electron.condition:                                                                   # Monte Carlo Loop until backscatter, trasmission or stopping
            


            RND_phi = random.uniform(a=0, b=1)                                                                 # generate random number for the phi angle
            RND_step = random.uniform(a = 0.000001, b = 0.999999)
            RND_pho = random.uniform(a = 0, b = 1)
            alpha = evaluate_alpha(electron.E, self.AlphaMultiplier)                                                            # calc screening constant, function of previous energy
            cross = evaluate_cross_section_opt(electron.E, self.CrossSectionLogArgMultiplier, self.CrossSectionNumorator, self.CrossSectionDenominatorA, self.CrossSectionDenominatorB)                                             # calc cross section
            path_length = evaluate_path_length(electron.CrossSection, self.PathLengthMultiplier)                                      # calc mean free path length
            step = evaluate_step(path_length, RND_step)                             # calculate step of this iteration
            dE = step * self.Density * evaluate_energy_loss_rate(electron.E, self.EnergyLossMultiplierA, self.EnergyLossMultiplierB)
            electron.E = electron.E + dE                                                                              # calc new energy #separate out dE/dS
            phi = evaluate_phi(RND_phi, alpha=alpha)                                            # calc scattering angle
            psi = evaluate_pho(RND_pho, self.pi)                                                         # calc other scattering angle
            ca = evaluate_direction_cosine_a(phi, psi, electron.cosineX, electron.cosineY, electron.cosineZ)                                  # calc direction cosines
            cb = evaluate_direction_cosine_b(phi, psi, electron.cosineX, electron.cosineY, electron.cosineZ)
            cc = evaluate_direction_cosine_c(phi, psi, electron.cosineZ)
            electron.x0 = electron.x0 + step * ca                                                                     # find and reset to new positions
            electron.y0 = electron.y0 + step * cb
            electron.z0 = electron.z0 + step * cc
            electron.cosineX = ca                                                                                 # reset direction cosines
            electron.cosineY = cb
            electron.cosineZ = cc

            if electron.E <= self.MinimumEnergy:                                             #if electron stops in the material
                self.number_stopped = self.number_stopped + 1
                electron.condition = False
            
            if electron.z0 < 0.01:                                                                   #if electron backscatters
                self.number_backscattered = self.number_backscattered + 1
                electron.condition = False                                                              #not sure how to deal with electron scattering outside of material
            
            if electron.z0 > pixel.x_dimension:                                                        # if electron penetrates material #change to have de/ds as the counting
                                                                                
                self.number_transmitted = self.number_transmitted + 1
                electron.condition = False
                    
            if -1*dE >= self.dE_threshold:                                                            #if electron deposits sufficient energyo
                self.new_eh_pairs = math.floor(-1*dE/self.dE_threshold)
                self.number_eh_pairs = int(self.number_eh_pairs + self.new_eh_pairs)

                if (electron.x0 <= 1*pixel.x_dimension) and (electron.x0 >= -1*pixel.x_dimension) and (electron.y0 <= 1*pixel.y_dimension)and (electron.y0 >=-1*pixel.y_dimension):                                        #electron stays within pixel region

                    eh_charge_counter[pixel.i, pixel.j] = self.new_eh_pairs + eh_charge_counter[pixel.i, pixel.j]

                elif (electron.x0 > 1*pixel.x_dimension) and (electron.y0 > 1*pixel.y_dimension):                   #electron moves positively in x
                    # and y
                    translation_x = custom_round_down(x = electron.x0/(pixel.x_dimension))
                    translation_y = custom_round_down(x = electron.y0/(2*pixel.y_dimension))
                    if pixel.i + translation_x  <= self.image_object.perfect_image.shape[0] - 1 and pixel.j + \
                            translation_y <= self.image_object.perfect_image.shape[0] - 1:
                            eh_charge_counter[pixel.i + translation_x, pixel.j + translation_y] = self.new_eh_pairs + eh_charge_counter[pixel.i +translation_x, pixel.j+translation_y]

                elif (electron.x0 < -1*pixel.x_dimension) and (electron.y0 < -1*pixel.y_dimension):                 #electron moves negatively in x
                    # and y
                    translation_x = custom_round_down(x=electron.x0 / (2 * pixel.x_dimension))
                    translation_y = custom_round_down(x=electron.y0 / (2 * pixel.y_dimension))
                    if pixel.i - translation_x >= 0 and pixel.j - translation_y >= 0:
                        eh_charge_counter[pixel.i - translation_x, pixel.j - translation_y] = self.new_eh_pairs + eh_charge_counter[pixel.i -translation_x, pixel.j -translation_y]


                elif electron.x0 > (1*pixel.x_dimension) and electron.y0 < (-1*pixel.y_dimension):  # electron moves psoitively in
                    # x and negatively in y
                    translation_x = custom_round_down(x=electron.x0 / (2 * pixel.x_dimension))
                    translation_y = custom_round_down(x=electron.y0 / (2 * pixel.y_dimension))
                    if pixel.j - translation_y >= 0 and pixel.i + translation_x  <= \
                        self.image_object.perfect_image.shape[0] - 1:
                        eh_charge_counter[pixel.i + translation_x, pixel.j - translation_y] = self.new_eh_pairs + eh_charge_counter[pixel.i +translation_x,pixel.j -translation_y]


                elif electron.x0 < (-1*pixel.x_dimension) and electron.y0 > (1*pixel.y_dimension):  # electron moves negatively in
                            # x and positively in y
                    translation_x = custom_round_down(x=electron.x0 / (2 * pixel.x_dimension))
                    translation_y = custom_round_down(x=electron.y0 / (2 * pixel.y_dimension))
                    if pixel.i - translation_x >= 0 and pixel.j + \
                        translation_y <= self.image_object.perfect_image.shape[0] - 1:
                        eh_charge_counter[pixel.i - translation_x, pixel.j + translation_y] = self.new_eh_pairs + eh_charge_counter[pixel.i -translation_x,pixel.j +translation_y]


                elif electron.x0 > (1*pixel.x_dimension):                                               #electron moves in positive x direction
                    translation_x = custom_round_down(x=electron.x0 / (2 * pixel.x_dimension))
                    if pixel.i + \
                        translation_x <= self.image_object.perfect_image.shape[0] - 1:

                        eh_charge_counter[pixel.i + translation_x, pixel.j] = self.new_eh_pairs + eh_charge_counter[pixel.i +translation_x, pixel.j]

                elif electron.x0 < (-1*pixel.x_dimension):                                          #electron moves in negative x direction
                    translation_x = custom_round_down(x=electron.x0 / (2 * pixel.x_dimension))
                    if pixel.i - translation_x >= 0:

                        eh_charge_counter[pixel.i - translation_x, pixel.j] = self.new_eh_pairs + eh_charge_counter[pixel.i -translation_x,pixel.j]

                elif electron.y0 > (1*pixel.y_dimension):                                       #electron moves in positive y direction
                    translation_y = custom_round_down(x=electron.y0 / (2 * pixel.y_dimension))
                    if pixel.j + \
                        translation_y <= self.image_object.perfect_image.shape[0] - 1:

                        eh_charge_counter[pixel.i, pixel.j + translation_y] = self.new_eh_pairs + eh_charge_counter[pixel.i, pixel.j +translation_y]

                elif electron.y0 < (-1*pixel.x_dimension): #electron moves in negative y direction
                    translation_y = custom_round_down(x=electron.y0 / (2 * pixel.y_dimension))

                    if pixel.j - translation_y >= 0:

                        eh_charge_counter[pixel.i, pixel.j - translation_y] = self.new_eh_pairs + eh_charge_counter[pixel.i, pixel.j -translation_y]

        eh_charge_counter = (np.floor(self.dE_threshold*eh_charge_counter/self.t_counting))
        self.new_image_MCS += eh_charge_counter




    def RunMCScatteringSimulation_SinglePixel(self,pixel):
        pixel_count = 0
        for i in prange(1, pixel.electron_count+1):
            self.MCScatteringSimulation_SingleElectron(pixel)
            print('1')


    def RunMCScatteringSimulation(self, num_threads=None):

        num_threads = multiprocessing.cpu_count()
        print(f'Number of threads: {num_threads}')

        for pixel in self.image_object.pixels:
            self.RunMCScatteringSimulation_SinglePixel(pixel)        
        return self.new_image_MCS
    


