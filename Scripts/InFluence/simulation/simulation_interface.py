from numba.typed import List
from simulation.python_simulation.python_simulation_interface import MonteCarloScatteringSimulation



class InFluenceSimulation:

    def __init__(self, image_object, params):
        self.image_object = image_object
        self.image_object.generate_scaled_perfect_image()
        self.image_object.generate_perfect_image_pixel_objects()
        self.ProbeDiameter = self.image_object.pixel_dimensions[0]
        self.E_i = params.E_i 
        self.dE_threshold = params.dE_threshold                                       
        self.MinimumEnergy = params.minimum_energy           
        self.ProtonNum = params.ProtonNum
        self.AtomicMass = params.AtomicMass                                         
        self.Density = params.Density
        self.t_counting = params.t_counting
        self.N = params.N
        image_object = None
        image_object = None
        self.params = params

    def RunMCSimulation(self):


        simulation = MonteCarloScatteringSimulation(
            self.image_object, self.params
        )
        return  simulation.RunMCScatteringSimulation_Sparse()
    
        
    
    def RunCPPMCScatteringSimulation(self):   

        from cpp_simulation import CPP_MCScatteringSimulation
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
        from cuda_simulation import gpu_simulation_interface

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


    def SaveModulatedImage():

        pass







