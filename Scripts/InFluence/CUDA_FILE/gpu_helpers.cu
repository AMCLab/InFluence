#include "gpu_helpers.cuh"
#include <cmath>


__device__ int custom_round_down(double x, int decimals) {
    double multiplier = pow(10.0, static_cast<double>(decimals));
    double rounded_value = ceil(fabs(x) * multiplier - 0.5);
    rounded_value /= multiplier;
    return static_cast<int>(rounded_value);
}

__device__ double evaluate_alpha(double E, double AlphaMultiplier) {
    return AlphaMultiplier / E;
}

__device__ double evaluate_path_length(double CrossSection, double PathLengthMultiplier) {
    double PathLength = PathLengthMultiplier * (1.0 / CrossSection);
    return PathLength;
}

__device__ float evaluate_cross_section_opt(float E, float CrossSectionLogArgMultiplier, float CrossSectionNumorator, float CrossSectionDenominatorA, float CrossSectionDenominatorB) {
    float LogArg = log10f(E * CrossSectionLogArgMultiplier);
    float CrossSection = CrossSectionNumorator / ((E + sqrtf(E) * CrossSectionDenominatorA) * (1.0f - expf(-1.0f * LogArg * LogArg) * CrossSectionDenominatorB));
    return CrossSection;
}

__device__ double evaluate_step(double path_length, double RandomStep) {
    return -path_length * logf(RandomStep);
}

__device__ double evaluate_phi(double RandomNum, double alpha) {
    double cosPhi = 1.0 - (2.0 * alpha * pow(RandomNum, 2)) / (1.0 + alpha - RandomNum);
    return acosf(cosPhi);
}

__device__ double evaluate_pho(double RandomNum, double pi) {
    return 2.0 * pi * RandomNum;
}


__device__ double evaluate_direction_cosine_a(double phi, double psi, double cosineY, double cosineZ) {
    float alpha = sinf(psi) * sinf(phi);
    float beta = sinf(phi) * cosf(psi);
    float gamma = cosf(phi);
    float cos_1 = static_cast<float>(cosineZ);
    float sin_1 = sqrtf(1.0f - cos_1 * cos_1);
    float cos_2 = static_cast<float>(cosineY) / sin_1;
    float sin_2 = sqrtf(1.0f - cos_2 * cos_2);
    return static_cast<double>(alpha * cos_2 + sin_2 * (beta * cos_1 + gamma * sin_1));
}


__device__ double evaluate_direction_cosine_c(double phi, double psi, double cosineZ) {
    float beta = sinf(static_cast<float>(phi)) * cosf(static_cast<float>(psi));
    float gamma = cosf(static_cast<float>(phi));
    float sin_1 = sqrtf(1.0f - static_cast<float>(cosineZ) * static_cast<float>(cosineZ));
    float cos_1 = static_cast<float>(cosineZ);

    return static_cast<double>(-beta * sin_1 + gamma * cos_1);
}

__device__ double evaluate_energy_loss_rate(double E, double EnergyLossMultiplierA, double EnergyLossMultiplier) {
    float log_arg = logf(1.166f * (E + 0.85f * static_cast<float>(EnergyLossMultiplier)) / static_cast<float>(EnergyLossMultiplier));
    return static_cast<double>(EnergyLossMultiplierA * (1.0f / static_cast<float>(E)) * log_arg);
}

__device__ double get_random_uniform(double a, double b) {
    curandState state;
    curand_init(clock64(), 0, 0, &state); // Seed the random number generator

    double random_value = curand_uniform_double(&state); // Generate a random double between 0 and 1
    return a + random_value * (b - a); // Scale to [a, b]
}

__device__ void initialize_positions(double step, double ProbeDiameter, double* position_data) {
    double z0 = 1e-2 + step;
    double x0 = ProbeDiameter * get_random_uniform(-1.0, 1.0);
    double y0 = ProbeDiameter * get_random_uniform(-1.0, 1.0);
    double vector_length = sqrtf(x0 * x0 + y0 * y0 + z0 * z0);
    double cosineX = x0 / vector_length;
    double cosineY = y0 / vector_length;
    double cosineZ = z0 / vector_length;

    // Store the results in the custom data structure (array)
    position_data[0] = cosineX;
    position_data[1] = cosineY;
    position_data[2] = cosineZ;
    position_data[3] = z0;
    position_data[4] = y0;
    position_data[5] = x0;
    position_data[6] = vector_length;
}