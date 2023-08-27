import math
import numpy as np
import random
from tqdm import tqdm
from parameters import params
from ElectronObjectConstructor import ElectronObject
from Scripts.InFluence.HelperFunctions import custom_round_down, CalculateConstants, evaluate_alpha, evaluate_cross_section_opt, evaluate_path_length, evaluate_step, evaluate_phi, evaluate_pho, evaluate_direction_cosine_a, evaluate_direction_cosine_b, evaluate_direction_cosine_c, evaluate_energy_loss_rate, initialise_postions
from Scripts.InFluence.PartiallyVectorizedMCScatteringSimulation import *


class MonteCarloScatteringSimulation:
    
    # input paramaters

    image_object = None
    E_i = None
    dE_threshold = None
    MinimumEnergy = None
    ProtonNum = None
    AtomicMass = None
    Density = None
    t_counting = None
    N = None
    AvagadrosConstant = params.N                                            #This will always be the same so opted not to pass as an argument to the class from the Influence class

    # constants (defined out of loop improve performence)

    AlphaMultiplier = None
    CrossSectionNumorator = None
    CrossSectionLogArgMultiplier = None
    CrossSectionDenominatorA = None
    CrossSectionDenominatorB = None
    PathLengthMultiplier = None
    pi = math.pi
    EnergyLossMultiplierA = None
    EnergyLossMultiplierB = None
    ProbeDiameter = None

    # derived scattering parameters

    alpha = None 
    CrossSection = None 
    PathLength = None 
    RND_step =  None 
    step = None 
    InitialPosition = None 
    cosineX = None 
    cosineY = None 
    cosineZ = None 
    x0 = None
    y0 = None
    z0 = None
    condition = None                                                                            #What is this
    E = None

    # global counters - shared memory

    number_transmitted = 0
    number_backscattered = 0
    number_stopped = 0
    number_eh_pairs = 0
    new_image_MCS = 0
    number_eh_pairs = 0
    eh_charge_counter = None
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


        self.AlphaMultiplier, self.CrossSectionNumorator, self.CrossSectionLogArgMultiplier, self.CrossSectionDenominatorA, self.CrossSectionDenominatorB, self.PathLengthMultiplier, self.EnergyLossMultiplierA, self.EnergyLossMultiplierB = CalculateConstants(params)

        self.new_image_MCS = np.zeros((self.image_object.perfect_image.shape[0], self.image_object.perfect_image.shape[1]), dtype=np.float64)

        self.loop_counter = 0

        self.eh_charge_counter = np.zeros((self.image_object.perfect_image.shape[0], self.image_object.perfect_image.shape[1]), dtype = np.float64)
                                                                                    
       

        

    def RunMCScatteringSimulation(self):

        print('Simulation starting ...')

        pixels_array = np.array([
        (pixel.electron_count,                  # Column 0: Number of electrons in the pixel
         pixel.i,                               # Column 1: Pixel index i
         pixel.j,                               # Column 2: Pixel index j
         pixel.x_dimension,                     # Column 3: Dimension along the x-axis                  
         pixel.y_dimension,                     # Column 4: Dimension along the y-axis
         pixel.z_dimension)                     # Column 5: Dimension along the z-axis
        
        for pixel in self.image_object.pixels], dtype=[('electron_count', int), ('i', int), ('j', int), ('x_dimension', float), ('y_dimension', float), ('z_dimension', float)])

        
        print('Pixel data prepared ...')

        print('Beginning main simulation loop...')

        for pixel in pixels_array:


            RunMCScatteringSimulation_Vectorized(pixel, self.E_i, self.ProbeDiameter)
        
        print('Simulation complete.')








                

        
        

