import sys
sys.path.append("../../")
from parameters import SimulationParameters

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

params = SimulationParameters()
params.refresh_config()
AlphaMultiplier, CrossSectionNumorator, CrossSectionLogArgMultiplier, CrossSectionDenominatorA, CrossSectionDenominatorB, PathLengthMultiplier, EnergyLossMultiplierA, EnergyLossMultiplierB = CalculateConstants(params)