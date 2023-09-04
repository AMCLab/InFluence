
#ifndef MCLOOP_H
#define MCLOOP_H

#include <vector>

std::vector<std::vector<double>> RunMCScatteringSimulation(const std::vector<std::vector<double>>& pixels, double E_i, double ProbeDiameter, double MinimumEnergy,
                              double dE_threshold, int perfect_image_0, 
                              int perfect_image_1, double Density, double t_counting);

#endif MCLOOP_H