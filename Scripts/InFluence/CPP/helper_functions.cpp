#include "helper_functions.h"
#include <cmath>
#include <random>




int custom_round_down(double x, int decimals) {
    double multiplier = pow(10, decimals);
    double rounded_value = ceil(fabs(x) * multiplier - 0.5);
    rounded_value /= multiplier;
    return static_cast<int>(rounded_value);
}

std::vector<int> custom_round_down_vector(const std::vector<double>& input_vec, int decimals) {
    std::vector<int> result_vec;

    for (const double& value : input_vec) {
        int result = custom_round_down(value, decimals);
        result_vec.push_back(result);
    }

    return result_vec;
}

double evaluate_alpha(double E, double AlphaMultiplier) {
    return AlphaMultiplier / E;
}

std::vector<double> evaluate_alpha_vector(const std::vector<double>& energies, double AlphaMultiplier) {
    std::vector<double> alpha_values;

    for (const double& E : energies) {
        double alpha = evaluate_alpha(E, AlphaMultiplier);
        alpha_values.push_back(alpha);
    }

    return alpha_values;
}

double evaluate_cross_section_opt(double E, double CrossSectionLogArgMultiplier, double CrossSectionNumorator, double CrossSectionDenominatorA, double CrossSectionDenominatorB) {
    double LogArg = log10(E * CrossSectionLogArgMultiplier);
    double CrossSection = CrossSectionNumorator / ((E + sqrt(E) * CrossSectionDenominatorA) * (1 - exp(-1 * LogArg * LogArg) * CrossSectionDenominatorB));
    return CrossSection;
}

std::vector<double> evaluate_cross_section_opt_vector(const std::vector<double>& energies, double CrossSectionLogArgMultiplier, double CrossSectionNumorator, double CrossSectionDenominatorA, double CrossSectionDenominatorB) {
    std::vector<double> cross_section_values;

    for (const double& E : energies) {
        double cross_section = evaluate_cross_section_opt(E, CrossSectionLogArgMultiplier, CrossSectionNumorator, CrossSectionDenominatorA, CrossSectionDenominatorB);
        cross_section_values.push_back(cross_section);
    }

    return cross_section_values;
}

double evaluate_path_length(double CrossSection, double PathLengthMultiplier) {
    double PathLength = PathLengthMultiplier * (1 / CrossSection);
    return PathLength;
}

std::vector<double> evaluate_path_length_vector(const std::vector<double>& cross_sections, double PathLengthMultiplier) {
    std::vector<double> path_length_values;

    for (const double& cross_section : cross_sections) {
        double path_length = evaluate_path_length(cross_section, PathLengthMultiplier);
        path_length_values.push_back(path_length);
    }

    return path_length_values;
}


double evaluate_step(double path_length, double RandomStep) {
    return -path_length * std::log(RandomStep);
}


std::vector<double> evaluate_step(const std::vector<double>& path_lengths, const std::vector<double>& RandomSteps) {
    std::vector<double> results;
    results.reserve(path_lengths.size());

    for (size_t i = 0; i < path_lengths.size(); ++i) {
        results.push_back(-path_lengths[i] * std::log(RandomSteps[i]));
    }

    return results;
}

double evaluate_phi(double RandomNum, double alpha) {
    double cosPhi = 1.0 - (2.0 * alpha * std::pow(RandomNum, 2)) / (1.0 + alpha - RandomNum);
    return std::acos(cosPhi);
}

std::vector<double> evaluate_phi(const std::vector<double>& RandomNums, double alpha) {
    std::vector<double> cosPhis;
    cosPhis.reserve(RandomNums.size());

    for (double RandomNum : RandomNums) {
        double cosPhi = 1.0 - (2.0 * alpha * std::pow(RandomNum, 2)) / (1.0 + alpha - RandomNum);
        cosPhis.push_back(cosPhi);
    }

    std::vector<double> phis;
    phis.reserve(cosPhis.size());

    for (double cosPhi : cosPhis) {
        phis.push_back(std::acos(cosPhi));
    }

    return phis;
}

double evaluate_pho(double RandomNum, double pi) {
    return 2.0 * pi * RandomNum;
}

std::vector<double> evaluate_pho(const std::vector<double>& RandomNums, double pi) {
    std::vector<double> phos;
    phos.reserve(RandomNums.size());

    for (double RandomNum : RandomNums) {
        phos.push_back(2.0 * pi * RandomNum);
    }

    return phos;
}

double evaluate_direction_cosine_a(double phi, double psi, double cosineY, double cosineZ) {
    double alpha = std::sin(psi) * std::sin(phi);
    double beta = std::sin(phi) * std::cos(psi);
    double gamma = std::cos(phi);
    double cos_1 = cosineZ;
    double sin_1 = std::sqrt(1.0 - cosineZ * cosineZ);
    double cos_2 = cosineY / sin_1;
    double sin_2 = cosineY / sin_1;
    return alpha * cos_2 + sin_2 * (beta * cos_1 + gamma * sin_1);
}

std::vector<double> evaluate_direction_cosine_a(const std::vector<double>& phis,
                                                const std::vector<double>& psis,
                                                const std::vector<double>& cosineYs,
                                                const std::vector<double>& cosineZs) {
    std::vector<double> results;
    results.reserve(phis.size());

    for (size_t i = 0; i < phis.size(); ++i) {
        double phi = phis[i];
        double psi = psis[i];
        double cosineY = cosineYs[i];
        double cosineZ = cosineZs[i];

        double alpha = std::sin(psi) * std::sin(phi);
        double beta = std::sin(phi) * std::cos(psi);
        double gamma = std::cos(phi);
        double sin_1 = std::sqrt(1.0 - cosineZ * cosineZ);
        double cos_2 = cosineY / sin_1;
        double sin_2 = cosineY / sin_1;

        results.push_back(alpha * cos_2 + sin_2 * (beta * cosineZ + gamma * sin_1));
    }

    return results;
}

double evaluate_direction_cosine_b(double phi, double psi, double cosineX, double cosineY, double cosineZ) {
    double alpha = std::sin(phi) * std::sin(psi);
    double beta = std::sin(phi) * std::cos(psi);
    double gamma = std::cos(phi);
    double sin_1 = std::sqrt(1.0 - cosineZ * cosineZ);
    double cos_2 = cosineY / sin_1;
    double sin_2 = cosineX / sin_1;

    return -alpha * sin_2 + cos_2 * (beta * cosineZ + gamma * sin_1);
}

std::vector<double> evaluate_direction_cosine_b(const std::vector<double>& phis,
                                                const std::vector<double>& psis,
                                                const std::vector<double>& cosineXs,
                                                const std::vector<double>& cosineYs,
                                                const std::vector<double>& cosineZs) {
    std::vector<double> results;
    results.reserve(phis.size());

    for (size_t i = 0; i < phis.size(); ++i) {
        double alpha = std::sin(phis[i]) * std::sin(psis[i]);
        double beta = std::sin(phis[i]) * std::cos(psis[i]);
        double gamma = std::cos(phis[i]);
        double sin_1 = std::sqrt(1.0 - cosineZs[i] * cosineZs[i]);
        double cos_2 = cosineYs[i] / sin_1;
        double sin_2 = cosineXs[i] / sin_1;

        results.push_back(-alpha * sin_2 + cos_2 * (beta * cosineZs[i] + gamma * sin_1));
    }

    return results;
}

double evaluate_direction_cosine_c(double phi, double psi, double cosineZ) {
    double beta = std::sin(phi) * std::cos(psi);
    double gamma = std::cos(phi);
    double sin_1 = std::sqrt(1.0 - cosineZ * cosineZ);
    double cos_1 = cosineZ;

    return -beta * sin_1 + gamma * cos_1;
}


std::vector<double> evaluate_direction_cosine_c(const std::vector<double>& phis,
                                                const std::vector<double>& psis,
                                                const std::vector<double>& cosineZs) {
    std::vector<double> results;
    results.reserve(phis.size());

    for (size_t i = 0; i < phis.size(); ++i) {
        double beta = std::sin(phis[i]) * std::cos(psis[i]);
        double gamma = std::cos(phis[i]);
        double sin_1 = std::sqrt(1.0 - cosineZs[i] * cosineZs[i]);
        double cos_1 = cosineZs[i];

        results.push_back(-beta * sin_1 + gamma * cos_1);
    }

    return results;
}


double evaluate_energy_loss_rate(double E, double EnergyLossMultiplierA, double EnergyLossMultiplier) {
    return EnergyLossMultiplierA * (1.0 / E) * std::log(1.166 * (E + 0.85 * EnergyLossMultiplier) / EnergyLossMultiplier);
}

std::vector<double> evaluate_energy_loss_rate(const std::vector<double>& energies,
                                              double EnergyLossMultiplierA,
                                              double EnergyLossMultiplier) {
    std::vector<double> results;
    results.reserve(energies.size());

    for (size_t i = 0; i < energies.size(); ++i) {
        results.push_back(evaluate_energy_loss_rate(energies[i], EnergyLossMultiplierA, EnergyLossMultiplier));
    }

    return results;
}

double get_random_uniform(double a, double b) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(a, b);
    return dis(gen);
}

std::vector<double> initialize_positions(double step, double ProbeDiameter) {
    double z0 = 1e-2 + step;
    double x0 = ProbeDiameter * get_random_uniform(-1.0, 1.0);
    double y0 = ProbeDiameter * get_random_uniform(-1.0, 1.0);
    double vector_length = std::sqrt(x0 * x0 + y0 * y0 + z0 * z0);
    double cosineX = x0 / vector_length;
    double cosineY = y0 / vector_length;
    double cosineZ = z0 / vector_length;

    return {cosineX, cosineY, cosineZ, z0, y0, x0, vector_length};
}

std::vector<std::vector<double>> initialize_positions(const std::vector<double>& steps, double ProbeDiameter) {
    std::vector<std::vector<double>> results;
    results.reserve(steps.size());

    for (double step : steps) {
        results.push_back(initialize_positions(step, ProbeDiameter));
    }

    return results;
}



