#include <iostream>
#include <vector>

#include "helper_functions.h"
#include "electrons_generator.h"
#include "MCloop.h"

int main() {
    // Create the sample 2D vector "Pixels"


    double x_dimension = 1.0;
    double y_dimension = 2.0;
    double z_dimension = 3.0;

    // Number of rows in the array
    int num_rows = 10000;

    // Initialize the vector of vectors
    std::vector<std::vector<double>> Pixels;

    for (int i = 0; i < num_rows; ++i) {
        int count = i % 4 + 1;
        int x = i % 10 + 1;
        int y = i % 10 + 1;
        std::vector<double> row = {static_cast<double>(count), static_cast<double>(x), static_cast<double>(y), x_dimension, y_dimension, z_dimension};
        Pixels.push_back(row);
    }




   double E_i = 200; // Replace with the desired value of E_i
    double ProbeDiameter = 0.1; // Replace with the desired value of ProbeDiameter
    double MinimumEnergy = 0.05; // Replace with the desired value of MinimumEnergy


    double dE_threshold = 0.00362; // Replace with the desired value of dE_threshold
    int perfect_image_0 = 10000; // Replace with the desired value of perfect_image_0
    int perfect_image_1 = 10000; // Replace with the desired value of perfect_image_1
    double Density = 2.33; // Replace with the desired value of Density
    double t_counting = 40; // Replace with the desired value of t_counting

    // Call the function for each row of Pixels
    RunMCScatteringSimulation(Pixels, E_i, ProbeDiameter, MinimumEnergy, dE_threshold, perfect_image_0, perfect_image_1, Density, t_counting);

    std::cout<<"Simulation complete."<<std::endl;



    return 0;

}