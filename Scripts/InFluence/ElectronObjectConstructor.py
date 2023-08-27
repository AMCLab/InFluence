# electron.py

from numba import float64, boolean
from numba.experimental import jitclass

import math
import numpy as np
import random

from HelperFunctions import custom_round_down, CalculateConstants, evaluate_alpha, evaluate_cross_section_opt, evaluate_path_length, evaluate_step, evaluate_phi, evaluate_pho, evaluate_direction_cosine_a, evaluate_direction_cosine_b, evaluate_direction_cosine_c, evaluate_energy_loss_rate, initialise_postions, AlphaMultiplier, CrossSectionNumorator, CrossSectionLogArgMultiplier, CrossSectionDenominatorA, CrossSectionDenominatorB, PathLengthMultiplier, EnergyLossMultiplierA, EnergyLossMultiplierB

local_random_uniform = random.uniform
spec = {
    'ProbeDiameter': float64,
    'E_i': float64,
    'E': float64,
    'alpha': float64,
    'CrossSection': float64,
    'PathLength': float64,
    'RND_step': float64,
    'step': float64,
    'InitialPosition': float64[:],
    'cosineX': float64,
    'cosineY': float64,
    'cosineZ': float64,
    'z0': float64,
    'y0': float64,
    'x0': float64,
    'condition': boolean,
    'dE': float64,
    'Density': float64
}

@jitclass(spec)
class ElectronEvolver:
    
    def __init__(self, E_i, ProbeDiameter, Density):
        self.ProbeDiameter = ProbeDiameter
        self.E_i = E_i
        self.E = E_i
        self.alpha = evaluate_alpha(self.E_i, AlphaMultiplier)
        self.CrossSection = evaluate_cross_section_opt(self.E_i, CrossSectionLogArgMultiplier, CrossSectionNumorator, CrossSectionDenominatorA, CrossSectionDenominatorB)
        self.PathLength = evaluate_path_length(self.CrossSection, PathLengthMultiplier)
        self.RND_step = np.random.uniform(0.000001, 0.999999)
        self.step = evaluate_step(self.PathLength, self.RND_step)
        self.InitialPosition = initialise_postions(self.step, self.ProbeDiameter)         #d = probe diameter
        self.cosineX = self.InitialPosition[0]
        self.cosineY = self.InitialPosition[1]
        self.cosineZ = self.InitialPosition[2]
        self.z0 = self.InitialPosition[3]
        self.y0 = self.InitialPosition[4]
        self.x0 = self.InitialPosition[5]
        self.condition = True
        self.dE = 0
        self.Density = Density

    def update_state(self):
        """Updates the state of the electron based on input parameters and methods."""
    
        RND_phi = local_random_uniform(0, 1)
        RND_step = local_random_uniform(0.000001, 0.999999)
        RND_pho = local_random_uniform(0, 1)
    
        alpha = evaluate_alpha(self.E, AlphaMultiplier)
        cross = evaluate_cross_section_opt(self.E, CrossSectionLogArgMultiplier, CrossSectionNumorator, CrossSectionDenominatorA, CrossSectionDenominatorB)
        path_length = evaluate_path_length(self.CrossSection, PathLengthMultiplier)
        step = evaluate_step(path_length, RND_step)
        self.dE = step * self.Density * evaluate_energy_loss_rate(self.E, EnergyLossMultiplierA, EnergyLossMultiplierB)
        self.E += self.dE
        phi = evaluate_phi(RND_phi, alpha=alpha)
        psi = evaluate_pho(RND_pho, math.pi)
    
        ca = evaluate_direction_cosine_a(phi, psi, self.cosineX, self.cosineY, self.cosineZ)
        cb = evaluate_direction_cosine_b(phi, psi, self.cosineX, self.cosineY, self.cosineZ)
        cc = evaluate_direction_cosine_c(phi, psi, self.cosineZ)
    
        self.x0 += step * ca
        self.y0 += step * cb
        self.z0 += step * cc
    
        self.cosineX = ca
        self.cosineY = cb
        self.cosineZ = cc
    
        


        







