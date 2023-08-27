#include "electrons_generator.h"
#include "helper_functions.h" // Include any necessary custom functions
#include "constants.h"

std::vector<std::vector<double>> generate_electrons(int electron_count, double E_i, double ProbeDiameter) {
    std::vector<std::vector<double>> electrons(electron_count, std::vector<double>(16));

    // Pre-calculate some common values
    double alpha = evaluate_alpha(E_i, AlphaMultiplier);
    double crossSection = evaluate_cross_section_opt(E_i, CrossSectionLogArgMultiplier, CrossSectionNumorator, CrossSectionDenominatorA, CrossSectionDenominatorB);
    double pathLength = evaluate_path_length(crossSection, PathLengthMultiplier);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(0.000001, 0.999999);

    for (int i = 0; i < electron_count; ++i) {
        electrons[i][0] = E_i; // E_i
        electrons[i][1] = E_i; // E
        electrons[i][2] = alpha; // Alpha
        electrons[i][3] = crossSection; // CrossSection
        electrons[i][4] = pathLength; // PathLength
        electrons[i][5] = dis(gen); // RND_step
        electrons[i][6] = evaluate_step(pathLength, electrons[i][5]); // step
        electrons[i][7] = ProbeDiameter; // ProbeDiameter

        std::vector<double> initialPositions = initialize_positions(electrons[i][6], ProbeDiameter); // InitialPosition
        electrons[i][9] = initialPositions[0]; // cosineX
        electrons[i][10] = initialPositions[1]; // cosineY
        electrons[i][11] = initialPositions[2]; // cosineZ
        electrons[i][12] = initialPositions[3]; // z0
        electrons[i][13] = initialPositions[4]; // y0
        electrons[i][14] = initialPositions[5]; // x0
        electrons[i][15] = 0; // dE
    }

    return electrons;
}
