from simulation.python_simulation.mcss_python_functions import MonteCarloScatteringSimulation_function
from simulation.python_simulation.mcss_python_functions_sparse import MonteCarloScatteringSimulationSparse_function
from simulation.python_simulation.mcss_python_functions_plot import MonteCarloScatteringSimulationPlot_function
from simulation.python_simulation.plot_trajectory_functions import PlotTrajectories


from math import pi
import numpy as np

from numba.typed import Dict
from numba.types import UniTuple, int64, float64

from numba.typed import List


class MonteCarloScatteringSimulation:

    def __init__(self, image_object, params):


        pixels_typed = List()
        for pixel in image_object.pixels:
            pixels_typed.append(pixel)
        
        self.pixels = pixels_typed
        self.image_object = image_object

        self.E_i = params.E_i
        self.dE_threshold = params.dE_threshold
        self.MinimumEnergy = params.minimum_energy
        self.ProtonNum = params.ProtonNum
        self.AtomicMass = params.AtomicMass
        self.Density = params.Density
        self.t_counting = params.t_counting
        self.N = params.N
        self.ProbeDiameter = image_object.pixel_dimensions[0]
        self.image_shape_0 = image_object.perfect_image.shape[0] 
        self.image_shape_1 = image_object.perfect_image.shape[1] 
        self.pi = pi
        self.params = params

        
       
       
        self.new_image_MCS = np.zeros((self.image_shape_0, self.image_shape_1), dtype=np.float64)
    

    def RunMCScatteringSimulation(self):
        SimulationResult = MonteCarloScatteringSimulation_function(self.pixels, 
                                                                        self.E_i, 
                                                                        self.dE_threshold, 
                                                                        self.MinimumEnergy, 
                                                                        self.Density, 
                                                                        self.t_counting, 
                                                                        self.N, 
                                                                        self.ProbeDiameter, 
                                                                        self.image_shape_0, 
                                                                        self.image_shape_1)
        
        return SimulationResult
    
    def RunMCScatteringSimulation_Sparse(self):
        SimulationResult = MonteCarloScatteringSimulationSparse_function(self.pixels, 
                                                                        self.E_i, 
                                                                        self.dE_threshold, 
                                                                        self.MinimumEnergy, 
                                                                        self.Density, 
                                                                        self.t_counting, 
                                                                        self.N, 
                                                                        self.ProbeDiameter, 
                                                                        self.image_shape_0, 
                                                                        self.image_shape_1)
        
        
        return np.reshape(SimulationResult, (self.image_shape_0, self.image_shape_1))
    
    
    def RunMCScatteringSimulation_Plot(self):
        if self.params.dose > 5000:
        # Raise an exception and exit the program with an error code
            error_message = "Dose is too high. Please enter an electron dose lower than 5000."
            raise ValueError(error_message)

        SimulationResult, vector_coordinates_x, vector_coordinates_y, vector_coordinates_z, points_per_trajectory = MonteCarloScatteringSimulationPlot_function(self.pixels, 
                                                                                                                                                                self.E_i, 
                                                                                                                                                                self.dE_threshold, 
                                                                                                                                                                self.MinimumEnergy, 
                                                                                                                                                                self.Density, 
                                                                                                                                                                self.t_counting, 
                                                                                                                                                                self.N, 
                                                                                                                                                                self.ProbeDiameter, 
                                                                                                                                                                self.image_shape_0, 
                                                                                                                                                                self.image_shape_1)
        PlotTrajectories(vector_coordinates_x, vector_coordinates_y, vector_coordinates_z, points_per_trajectory)
        return SimulationResult








    

    


