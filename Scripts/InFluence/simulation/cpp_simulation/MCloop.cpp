#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <iostream>

#include "MCloop.h"
#include "helper_functions.h"
#include <Eigen/Sparse>
#include <omp.h>

std::vector<std::vector<double>> RunMCScatteringSimulation(const std::vector<std::vector<double>>& pixels, double E_i, double ProbeDiameter, double MinimumEnergy,
                              double dE_threshold, int perfect_image_0, int perfect_image_1, double Density, double t_counting, double AlphaMultiplier,
    double CrossSectionNumorator,
    double CrossSectionLogArgMultiplier,
    double CrossSectionDenominatorA,
    double CrossSectionDenominatorB,
    double PathLengthMultiplier,
    double EnergyLossMultiplierA,
    double EnergyLossMultiplierB) {

    // Initialize 2D arrays as vectors of vectors

    Eigen::SparseMatrix<double> eh_charge_counter(perfect_image_0, perfect_image_1);
    std::vector<std::vector<double>> new_image_MCS(perfect_image_0, std::vector<double>(perfect_image_1, 0.0));

    int number_stopped = 0;
    int number_backscattered = 0;
    int number_transmitted = 0;
    int number_eh_pairs = 0;

    int pixel_counter = 0;
    int electron_per_pixel = 0;
    


    #pragma omp parallel
    for (const std::vector<double>& pixel : pixels) {
        int count = static_cast<int>(pixel[0]);
        double i_coordinate = pixel[1];
        double j_coordinate = pixel[2];
        double x_dimension = pixel[3];
        double y_dimension = pixel[4];
        double z_dimension = pixel[5];



        if (count == 0) {
        
        int do_nothing = 0;
        }

        else{



        for (int i = 0; i < count; ++i) {


            for (int k = 0; k < eh_charge_counter.outerSize(); ++k) {
                for (Eigen::SparseMatrix<double>::InnerIterator it(eh_charge_counter, k); it; ++it) {
                    it.valueRef() = 0.0;  // Set the non-zero element to zero
                }
            }
            
            //std::cout<<i<<std::endl;
            
            // Initiate electron conditions
            //std::vector<std::vector<double>> eh_charge_counter(perfect_image_0, std::vector<double>(perfect_image_1, 0.0));
            double alpha = evaluate_alpha(E_i, AlphaMultiplier);
            double CrossSection = evaluate_cross_section_opt(E_i, CrossSectionLogArgMultiplier, CrossSectionNumorator, CrossSectionDenominatorA, CrossSectionDenominatorB);
            double PathLength = evaluate_path_length(CrossSection, PathLengthMultiplier);
            double RND_step = get_random_uniform(0.000001, 0.999999);
            double step = evaluate_step(PathLength, RND_step);
            std::vector<double> ip = initialize_positions(step, ProbeDiameter);
            double cosineX = ip[0];
            double cosineY = ip[1];
            double cosineZ = ip[2];
            double z0 = ip[3];
            double y0 = ip[4];
            double x0 = ip[5];
            double E = E_i;
            bool condition = true;
            
            pixel_counter +=1;
            



            while (condition == true) {
                

                // Generate random numbers
                double RND_phi = get_random_uniform(0, 1);
                double RND_step = get_random_uniform(0.000001, 0.999999);
                double RND_pho = get_random_uniform(0, 1);

                // Calculate alpha, cross section, path length, and step
                double alpha = evaluate_alpha(E_i, AlphaMultiplier);
                double CrossSection = evaluate_cross_section_opt(E, CrossSectionLogArgMultiplier, CrossSectionNumorator, CrossSectionDenominatorA, CrossSectionDenominatorB);
                double PathLength = evaluate_path_length(CrossSection, PathLengthMultiplier);
                double step = evaluate_step(PathLength, RND_step);

                // Calculate energy loss and update energy
                double dE = step * Density * evaluate_energy_loss_rate(E, EnergyLossMultiplierA, EnergyLossMultiplierB);
                E = E + dE;

                // Calculate scattering angles and direction cosines
                double phi = evaluate_phi(RND_phi, alpha);
                double psi = evaluate_pho(RND_pho, M_PI); 

                double ca = evaluate_direction_cosine_a(phi, psi, cosineY, cosineZ);
                double cb = evaluate_direction_cosine_b(phi, psi, cosineX, cosineY, cosineZ);
                double cc = evaluate_direction_cosine_c(phi, psi, cosineZ);

                // Update positions
                x0 = x0 + step * ca;
                y0 = y0 + step * cb;
                z0 = z0 + step * cc;

                // Update direction cosines
                cosineX = ca;
                cosineY = cb;
                cosineZ = cc;


                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                //      Scenario 1 : Electron stops in material.                                                                 //
                //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                if (E <= MinimumEnergy) {
                    number_stopped++;
                    condition = false;          
                }

                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                //      Scenario 2 : Electron backscatters.                                                                     //
                //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                if (z0 < std::pow(10, -2)) {
                    number_backscattered++;
                    condition = false; 
                }

                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                //      Scenario 3 : Electron penetrates the material.                                                          //
                //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                if (z0 > z_dimension) {
                    number_transmitted++;
                    condition = false;
                }

                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                //      Scenario 4 : Electron deposits sufficient energy.                                                        //
                //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                if (-dE >= dE_threshold) {
                    int new_eh_pairs = std::floor(-dE / dE_threshold);
                    number_eh_pairs += new_eh_pairs;

                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    //      Scenario 5 : Electron deposits sufficient energy and stays within pixel boundary.               //
                    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                    if ((x0 <= x_dimension) && (x0 >= -x_dimension) && (y0 <= y_dimension) && (y0 >= -y_dimension)) {
                        #pragma omp critical
                        eh_charge_counter.coeffRef(i_coordinate, j_coordinate) += new_eh_pairs;
                        
                    }

                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    //      Scenario 6 : Electron deposits sufficient energy and moves positively in x and y.                       //
                    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                    else if ((x0 > x_dimension) && (y0 > y_dimension)) {
                        int translation_x = custom_round_down(x0 / (2 * x_dimension));
                        int translation_y = custom_round_down(y0 / (2 * y_dimension));
                        if ((i_coordinate + translation_x <= perfect_image_0 - 1) && (j_coordinate + translation_y <= perfect_image_1 - 1)) {
                                #pragma omp critical                           
                                eh_charge_counter.coeffRef(i_coordinate + translation_x, j_coordinate + translation_y) += new_eh_pairs;
                            }
                        }

                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    //      Scenario 7 : Electron deposits sufficient energy and moves negatively in x and y.                       //
                    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                    else if ((x0 < -x_dimension) && (y0 < -y_dimension)) {
                            int translation_x = custom_round_down(x0 / (2 * x_dimension));
                            int translation_y = custom_round_down(y0 / (2 * y_dimension));
                            if ((i_coordinate - translation_x >= 0) && (j_coordinate - translation_y >= 0)) {

                                    #pragma omp critical
                                    eh_charge_counter.coeffRef(i_coordinate - translation_x, j_coordinate - translation_y) += new_eh_pairs;
                                }
                            }
                        
                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    //      Scenario 8 : Electron deposits sufficient energy and moves positively in x only.                        //
                    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                    else if (x0 > x_dimension) {
                        int translation_x = custom_round_down(x0 / (2 * x_dimension));
                        if (i_coordinate + translation_x <= perfect_image_0 - 1) {

                                #pragma omp critical
                                eh_charge_counter.coeffRef(i_coordinate + translation_x, j_coordinate) += new_eh_pairs;
                            }
                        }

                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    //      Scenario 9 : Electron deposits sufficient energy and moves negatively in x only.                        //
                    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                    else if (x0 < -x_dimension) {
                        int translation_x = custom_round_down(x0 / (2 * x_dimension));
                        if (i_coordinate - translation_x >= 0) {
                           
                                #pragma omp critical
                                eh_charge_counter.coeffRef(i_coordinate - translation_x, j_coordinate) += new_eh_pairs;
                            
                        }
                    }
                    
                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    //      Scenario 10 : Electron deposits sufficient energy and moves positively in y only.                       //
                    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                    else if (y0 > y_dimension) {
                        int translation_y = custom_round_down(y0 / (2 * y_dimension));
                        if (j_coordinate + translation_y <= perfect_image_0 - 1) {

                                #pragma omp critical
                                eh_charge_counter.coeffRef(i_coordinate, j_coordinate + translation_y) += new_eh_pairs;
                            }
                        }
                    
                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    //      Scenario 11 : Electron deposits sufficient energy and moves negatively in y only.                       //
                    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                    else if (y0 < -y_dimension) {
                        int translation_y = custom_round_down(y0 / (2 * y_dimension));
                        if (j_coordinate - translation_y >= 0) {

                                #pragma omp critical
                                eh_charge_counter.coeffRef(i_coordinate, j_coordinate - translation_y) += new_eh_pairs;
                            
                        }
                    }



//---------------------------------------------------------------------------------------------------------------------------------------------------------------

                
                }               // Sufficient energy loop
            }                   // While loop 
                                // Add updated eh_charge_counter to new_image_MCS
            
            for (int k = 0; k < eh_charge_counter.outerSize(); ++k) {
                for (Eigen::SparseMatrix<double>::InnerIterator it(eh_charge_counter, k); it; ++it) {
                    new_image_MCS[it.row()][it.col()] += it.value();
                    it.valueRef() = 0.0; // Reset the non-zero elements for the next iteration
                }
            }
        


        }    //electron level loop
        }

    }   //pixels

return new_image_MCS;
}   
