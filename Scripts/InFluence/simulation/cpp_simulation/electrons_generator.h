#ifndef ELECTRONS_GENERATOR_H
#define ELECTRONS_GENERATOR_H

#include <vector>
#include <random>
#include <cmath>
#include <random>

std::vector<std::vector<double>> generate_electrons(int electron_count, double E_i, double ProbeDiameter, double AlphaMultiplier,
    double CrossSectionNumorator,
    double CrossSectionLogArgMultiplier,
    double CrossSectionDenominatorA,
    double CrossSectionDenominatorB,
    double PathLengthMultiplier,
    double EnergyLossMultiplierA,
    double EnergyLossMultiplierB);

#endif // ELECTRONS_GENERATOR_H
