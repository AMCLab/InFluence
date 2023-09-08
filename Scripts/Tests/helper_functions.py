import numpy as np
import random
import os
from parameters import SimulationParameters

params = SimulationParameters()
params.refresh_config()

local_log10 = np.log10
local_log = np.log
local_arccos = np.arccos
local_sin = np.sin
local_cos = np.cos
local_sqrt = np.sqrt
local_exp = np.exp
local_ceil = np.ceil
local_abs = np.abs


def custom_round_down(x, decimals=0):
    multiplier = 10 ** decimals
    rounded_value = local_ceil(local_abs(x) * multiplier - 0.5)
    rounded_value /= multiplier
    return int(rounded_value)

def CalculateConstants(ProgramParameters):                                                           
        
    AlphaMultiplier = (3.4*10**-3)*(ProgramParameters.ProtonNum**0.67)                          
    CrossSectionNumorator = (4.7 * 10 ** -18) * (ProgramParameters.ProtonNum ** 1.33 + 0.032 * ProgramParameters.ProtonNum  ** 2)
    CrossSectionLogArgMultiplier = 8 * ProgramParameters.ProtonNum  ** -1.33
    CrossSectionDenominatorA = 0.0155 * ProgramParameters.ProtonNum ** 1.33
    CrossSectionDenominatorB = 0.02 * ProgramParameters.ProtonNum  ** 0.5
    PathLengthMultiplier = ProgramParameters.AtomicMass/(ProgramParameters.N * ProgramParameters.Density)
    EnergyLossMultiplierA = -78500*ProgramParameters.ProtonNum/ProgramParameters.AtomicMass
    EnergyLossMultiplierB = (9.76*ProgramParameters.ProtonNum + 58.5/ProgramParameters.ProtonNum**0.19)*10**-3

    return AlphaMultiplier, CrossSectionNumorator, CrossSectionLogArgMultiplier, CrossSectionDenominatorA, CrossSectionDenominatorB, PathLengthMultiplier, EnergyLossMultiplierA, EnergyLossMultiplierB

AlphaMultiplier, CrossSectionNumorator, CrossSectionLogArgMultiplier, CrossSectionDenominatorA, CrossSectionDenominatorB, PathLengthMultiplier, EnergyLossMultiplierA, EnergyLossMultiplierB = CalculateConstants(params)


def evaluate_alpha(E, AlphaMultiplier):
    return AlphaMultiplier / E                             
        

def evaluate_cross_section_opt(E, CrossSectionLogArgMultiplier, CrossSectionNumorator, CrossSectionDenominatorA, CrossSectionDenominatorB):
    LogArg = local_log10(E * CrossSectionLogArgMultiplier)
    CrossSection = CrossSectionNumorator / ((E + (E ** 0.5) * CrossSectionDenominatorA) * (1 - local_exp(-1 * LogArg ** 2) * CrossSectionDenominatorB))
    return CrossSection


def evaluate_path_length(CrossSection, PathLengthMultiplier):
    PathLength = PathLengthMultiplier * (1/CrossSection)
    return PathLength


def evaluate_step(path_length, RandomStep):
    return -path_length*local_log(RandomStep)


def evaluate_phi(RandomNum, alpha):
    cosPhi = 1 - (2*alpha*RandomNum**2)/(1+alpha-RandomNum)
    return local_arccos(cosPhi)        


def evaluate_pho(RandomNum, pi):
    return 2*pi*RandomNum
   

def evaluate_direction_cosine_a(phi, psi, cosineX, cosineY, cosineZ):
    alpha = local_sin(psi) * local_sin(phi)
    beta = local_sin(phi) * local_cos(psi)
    gamma = local_cos(phi)
    cos_1 = cosineZ
    sin_1 = local_sqrt(1 - cosineZ**2)
    cos_2 = cosineY / sin_1
    sin_2 = cosineX / sin_1
    return alpha * cos_2 + sin_2 * (beta * cos_1 + gamma * sin_1)

def evaluate_direction_cosine_b(phi, psi, cosineX, cosineY, cosineZ):
    alpha = local_sin(phi) * local_sin(psi)
    beta = local_sin(phi) * local_cos(psi)
    gamma = local_cos(phi)
    cos_1 = cosineZ
    sin_1 = local_sqrt(1 - cosineZ ** 2)
    cos_2 = cosineY / sin_1
    sin_2 = cosineX / sin_1
    return -alpha * sin_2 + cos_2 * (beta * cos_1 + gamma * sin_1)


def evaluate_direction_cosine_c(phi, psi, cosineZ):
    beta = local_sin(phi) * local_cos(psi)
    gamma = local_cos(phi)
    cos_1 = cosineZ
    sin_1 = local_sqrt(1 - cosineZ ** 2)
    return -beta*sin_1 + gamma*cos_1


def evaluate_energy_loss_rate(E, EnergyLossMultiplierA, EnergyLossMultiplier):
    return EnergyLossMultiplierA*(1/E)*local_log(1.166*(E+0.85*EnergyLossMultiplier)/EnergyLossMultiplier)


def initialise_postions(step, ProbeDiameter):
    z0 = 10**-2 + step                              
    x0 = ProbeDiameter * random.uniform(-1, 1)
    y0 = ProbeDiameter * random.uniform(-1,  1)
    vector_length = local_sqrt(x0**2 + y0**2 + z0**2)
    cosineX = x0/vector_length
    cosineY = y0/vector_length
    cosineZ = z0/vector_length
    return np.array([cosineX, cosineY, cosineZ, z0, y0, x0, vector_length])


