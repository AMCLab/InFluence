

from parameters import params

from MonteCarloScatteringSimulation_Python import MonteCarloScatteringSimulation
from ImageObjectConstructor import *
from HelperFunctions import ConstantsArray
from numba.typed import List







class InFluenceSimulation:

    E_i = params.E_i                                              #whats this - rename
    dE_threshold = params.dE_threshold                                           #whast this -rename
    MinimumEnergy = params.minimum_energy           
    ProtonNum = params.ProtonNum
    AtomicMass = params.AtomicMass                                           #check this is atomic mass not weight / put them in th 
    Density = params.Density
    t_counting = params.t_counting
    N = params.N
    image_object = None
    modulated_image = None
    image_object = None 
    ProbeDiameter = None

    def __init__(self, image_object):
        self.image_object = image_object
        self.image_object.generate_scaled_perfect_image()
        self.image_object.generate_perfect_image_pixel_objects()
        self.ProbeDiameter = self.image_object.pixel_dimensions[0]

    def RunMCSimulation(self):

        pixels_typed = List()
        for pixel in self.image_object.pixels:
            pixels_typed.append(pixel)

        model = MonteCarloScatteringSimulation(
            pixels_typed, 
            self.E_i, 
            self.dE_threshold, 
            self.MinimumEnergy, 
            self.ProtonNum, 
            self.AtomicMass, 
            self.Density, 
            self.t_counting,
            self.N,
            self.image_object.perfect_image.shape[0],
            self.image_object.perfect_image.shape[1],
            params.N,
            ConstantsArray

        )

       
        return  model.RunMCScatteringSimulation_Sparse()
    
        
    
    def RunCPPMCScatteringSimulation(self):   
                                 #Requires CPU

        from MonteCarloScatteringSimulation_CPP import CPP_MCScatteringSimulation
        model = CPP_MCScatteringSimulation(
            self.image_object.pixels,
            self.E_i, 
            self.ProbeDiameter, 
            self.MinimumEnergy, 
            self.image_object.pixel_dimensions, 
            self.dE_threshold, 
            self.image_object.perfect_image.shape[0], 
            self.image_object.perfect_image.shape[1], 
            self.Density, 
            self.t_counting

        )
        return model
    
    def Run_GPU_KERNEL_MCScatteringSimulation(self):

        from MonteCarloScatteringSimualtion_GPU import GPU_KERNEL_MCScatteringSimulation

        model = GPU_KERNEL_MCScatteringSimulation(
            self.image_object.pixels,
            self.E_i, 
            self.ProbeDiameter, 
            self.MinimumEnergy, 
            self.image_object.pixel_dimensions, 
            self.dE_threshold, 
            self.image_object.perfect_image.shape[0], 
            self.image_object.perfect_image.shape[1], 
            self.Density, 
            self.t_counting

        )
        return model



    def PlotTrajectories():

        pass

    def SaveModulatedImage():

        pass







