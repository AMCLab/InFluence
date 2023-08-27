#ifndef HELPER_FUNCTIONS_H
#define HELPER_FUNCTIONS_H

#include <vector>
#include <cmath> 


int custom_round_down(double x, int decimals = 0);
std::vector<int> custom_round_down_vector(const std::vector<double>& input_vec, int decimals = 0);

double evaluate_alpha(double E, double AlphaMultiplier);
std::vector<double> evaluate_alpha_vector(const std::vector<double>& energies, double AlphaMultiplier);

double evaluate_cross_section_opt(double E, double CrossSectionLogArgMultiplier, double CrossSectionNumorator, double CrossSectionDenominatorA, double CrossSectionDenominatorB);
std::vector<double> evaluate_cross_section_opt_vector(const std::vector<double>& energies, double CrossSectionLogArgMultiplier, double CrossSectionNumorator, double CrossSectionDenominatorA, double CrossSectionDenominatorB);

double evaluate_path_length(double CrossSection, double PathLengthMultiplier);
std::vector<double> evaluate_path_length_vector(const std::vector<double>& cross_sections, double PathLengthMultiplier);

double evaluate_step(double path_length, double RandomStep);
std::vector<double> evaluate_step(const std::vector<double>& path_lengths, const std::vector<double>& RandomSteps);


double evaluate_phi(double RandomNum, double alpha);
std::vector<double> evaluate_phi(const std::vector<double>& RandomNums, double alpha);

double evaluate_pho(double RandomNum, double pi);
std::vector<double> evaluate_pho(const std::vector<double>& RandomNums, double pi);

double evaluate_direction_cosine_a(double phi, double psi, double cosineY, double cosineZ);


std::vector<double> evaluate_direction_cosine_a(const std::vector<double>& phis,
                                                const std::vector<double>& psis,
                                                const std::vector<double>& cosineYs,
                                                const std::vector<double>& cosineZs);

double evaluate_direction_cosine_b(double phi, double psi, double cosineX, double cosineY, double cosineZ);
std::vector<double> evaluate_direction_cosine_b(const std::vector<double>& phis,
                                                const std::vector<double>& psis,
                                                const std::vector<double>& cosineXs,
                                                const std::vector<double>& cosineYs,
                                                const std::vector<double>& cosineZs);


double evaluate_direction_cosine_c(double phi, double psi, double cosineZ);
std::vector<double> evaluate_direction_cosine_c(const std::vector<double>& phis,
                                                const std::vector<double>& psis,
                                                const std::vector<double>& cosineZs);


double evaluate_energy_loss_rate(double E, double EnergyLossMultiplierA, double EnergyLossMultiplier);
std::vector<double> evaluate_energy_loss_rate(const std::vector<double>& energies,
                                              double EnergyLossMultiplierA,
                                              double EnergyLossMultiplier);

double get_random_uniform(double a, double b);

std::vector<double> initialize_positions(double step, double ProbeDiameter);
std::vector<std::vector<double>> initialize_positions(const std::vector<double>& steps, double ProbeDiameter);



#endif // HELPER_FUNCTIONS_H



