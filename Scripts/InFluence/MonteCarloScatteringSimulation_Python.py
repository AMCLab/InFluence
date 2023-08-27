from MCSS_PythonFunctions import MonteCarloScatteringSimulation_function
from MCSS_PythonFunctions_Sparse import MonteCarloScatteringSimulationSparse_function
from MCSS_PythonFunctions_Plot import MonteCarloScatteringSimulationPlot_function
from PlotTrajectoryFunctions import PlotTrajectories
from math import pi
import numpy as np

from numba.typed import Dict
from numba.types import UniTuple, int64, float64


class MonteCarloScatteringSimulation:

    def __init__(self, pixels, E_i, dE_threshold, MinimumEnergy, ProtonNum, AtomicMass, Density, t_counting, N, image_shape_0, image_shape_1, AvagadrosConstant, ConstantsArray):
       
 
        self.pixels = pixels
        self.E_i = E_i
        self.dE_threshold = dE_threshold
        self.MinimumEnergy = MinimumEnergy
        self.ProtonNum = ProtonNum
        self.AtomicMass = AtomicMass
        self.Density = Density
        self.t_counting = t_counting
        self.N = N
        self.ProbeDiameter = image_shape_0
        self.image_shape_0 = image_shape_0
        self.image_shape_1 = image_shape_1
        self.AvagadrosConstant = AvagadrosConstant  
        self.pi = pi
       
        self.AlphaMultiplier = ConstantsArray[0]
        self.CrossSectionNumorator = ConstantsArray[1]
        self.CrossSectionLogArgMultiplier = ConstantsArray[2]
        self.CrossSectionDenominatorA = ConstantsArray[3]
        self.CrossSectionDenominatorB = ConstantsArray[4]
        self.PathLengthMultiplier = ConstantsArray[5]
        self.EnergyLossMultiplierA = ConstantsArray[6]
        self.EnergyLossMultiplierB = ConstantsArray[7]
       
        self.new_image_MCS = np.zeros((image_shape_0, image_shape_1), dtype=np.float64)
    

    def RunMCScatteringSimulation(self):
        SimulationResult = MonteCarloScatteringSimulation_function(self.pixels, self.E_i, self.dE_threshold, self.MinimumEnergy, self.ProtonNum, self.AtomicMass, self.Density, self.t_counting, self.N, self.ProbeDiameter, self.AlphaMultiplier, self.CrossSectionNumorator, self.CrossSectionLogArgMultiplier, self.CrossSectionDenominatorA, self.CrossSectionDenominatorB, self.PathLengthMultiplier, self.EnergyLossMultiplierA, self.EnergyLossMultiplierB, self.new_image_MCS, self.image_shape_0, self.image_shape_1)
        return SimulationResult
    
    def RunMCScatteringSimulation_Sparse(self):
        SimulationResult = MonteCarloScatteringSimulationSparse_function(self.pixels, self.E_i, self.dE_threshold, self.MinimumEnergy, self.ProtonNum, self.AtomicMass, self.Density, self.t_counting, self.N, self.ProbeDiameter, self.AlphaMultiplier, self.CrossSectionNumorator, self.CrossSectionLogArgMultiplier, self.CrossSectionDenominatorA, self.CrossSectionDenominatorB, self.PathLengthMultiplier, self.EnergyLossMultiplierA, self.EnergyLossMultiplierB, self.image_shape_0, self.image_shape_1)
        return np.reshape(SimulationResult, (self.image_shape_0, self.image_shape_1))
    
    
    def RunMCScatteringSimulation_Plot(self):
        SimulationResult, vector_coordinates_x, vector_coordinates_y, vector_coordinates_z, points_per_trajectory = MonteCarloScatteringSimulationPlot_function(self.pixels, self.E_i, self.dE_threshold, self.MinimumEnergy, self.ProtonNum, self.AtomicMass, self.Density, self.t_counting, self.N, self.ProbeDiameter, self.AlphaMultiplier, self.CrossSectionNumorator, self.CrossSectionLogArgMultiplier, self.CrossSectionDenominatorA, self.CrossSectionDenominatorB, self.PathLengthMultiplier, self.EnergyLossMultiplierA, self.EnergyLossMultiplierB, self.new_image_MCS, self.image_shape_0, self.image_shape_1)
        #PlotTrajectories(vector_coordinates_x, vector_coordinates_y, vector_coordinates_z, points_per_trajectory)
        return SimulationResult








    

    


