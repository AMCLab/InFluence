#include <curand_kernel.h>
#include <cuda_runtime.h>




__device__ float getRand(float a, float b, curandState* state)
{
    return a + (b - a) * curand_uniform(state);
}


extern "C" {


__device__ int custom_round_down(double x, int decimals=0) {
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
__device__ double evaluate_direction_cosine_a(double phi, double psi, double cosineX, double cosineY, double cosineZ) {
    float alpha = sinf(psi) * sinf(phi);
    float beta = sinf(phi) * cosf(psi);
    float gamma = cosf(phi);
    float cos_1 = static_cast<float>(cosineZ);
    float sin_1 = sqrtf(1.0f - cos_1 * cos_1);
    float cos_2 = static_cast<float>(cosineY) / sin_1;
    float sin_2 = static_cast<float>(cosineX) / sin_1;  // Corrected line
    return static_cast<double>(alpha * cos_2 + sin_2 * (beta * cos_1 + gamma * sin_1));
}

__device__ double evaluate_direction_cosine_b(double phi, double psi, double cosineX, double cosineY, double cosineZ) {
    float alpha = sinf(static_cast<float>(phi)) * sinf(static_cast<float>(psi));
    float beta = sinf(static_cast<float>(phi)) * cosf(static_cast<float>(psi));
    float gamma = cosf(static_cast<float>(phi));
    float sin_1 = sqrtf(1.0f - static_cast<float>(cosineZ) * static_cast<float>(cosineZ));
    float cos_2 = static_cast<float>(cosineY) / sin_1;
    float sin_2 = static_cast<float>(cosineX) / sin_1;

    return static_cast<double>(-alpha * sin_2 + cos_2 * (beta * static_cast<float>(cosineZ) + gamma * sin_1));
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


__device__ void initialize_positions(double step, double ProbeDiameter, double* position_data, curandState* state) {
    double z0 = 1e-2 + step;
    double x0 = ProbeDiameter * getRand(-1.0, 1.0, state);
    double y0 = ProbeDiameter * getRand(-1.0, 1.0, state);
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



__device__ void atomicAdd_double(double* address, double value) {
        unsigned long long int* address_as_ull = (unsigned long long int*)address;
        unsigned long long int old = *address_as_ull, assumed;

        do {
            assumed = old;
            old = atomicCAS(address_as_ull, assumed,
                            __double_as_longlong(value +
                            __longlong_as_double(assumed)));
        } while (assumed != old);
    }


__global__ void MCScatteringSimulationKernel(
    const double* pixels, 
    int numPixels, 
    double E_i, 
    double ProbeDiameter,
    double MinimumEnergy, 
    double dE_threshold, 
    int perfect_image_0,
    int perfect_image_1, 
    double Density, 
    double t_counting,
    double AlphaMultiplier, 
    double CrossSectionNumorator, 
    double CrossSectionLogArgMultiplier, 
    double CrossSectionDenominatorA, 
    double CrossSectionDenominatorB, 
    double PathLengthMultiplier, 
    double EnergyLossMultiplierA, 
    double EnergyLossMultiplierB,
    double* new_image_MCS) 
                                                
    {
    
    int threadIndex_x = blockIdx.x * blockDim.x + threadIdx.x;
    int threadIndex_y = blockIdx.y * blockDim.y + threadIdx.y;   

    

    // Convert 2D thread indices to 1D index
    int threadIndex = threadIndex_y * perfect_image_1 + threadIndex_x;

    curandState state;
    curand_init(1234, threadIndex, 0, &state); 
    

    // Check if the thread index is within the valid range of pixels
    if (threadIndex < numPixels) {
        // Get pixel data for this thread
        int pixelOffset = threadIndex * 6; // Each pixel has 6 values: count, i_coordinate, j_coordinate, x_dimension, y_dimension, z_dimension
        int count = static_cast<int>(pixels[pixelOffset]);
        double i_coordinate = pixels[pixelOffset + 1];
        double j_coordinate = pixels[pixelOffset + 2];
        double x_dimension = pixels[pixelOffset + 3];
        double y_dimension = pixels[pixelOffset + 4];
        double z_dimension = pixels[pixelOffset + 5];


        // Initialize 2D arrays as device arrays (raw pointers)
        int pixelIndex = static_cast<int>(i_coordinate) * perfect_image_0 + static_cast<int>(j_coordinate);
        //eh_charge_counter[pixelIndex] = 0.0;
        //new_image_MCS[pixelIndex] = 0.0;

        int number_transmitted = 0;
        int number_eh_pairs = 0;
        int number_stopped  = 0;
        int number_backscattered = 0;


        if (count == 0) {
            // Do nothing for this pixel

        } else {
            // Loop 'count' times for this pixel



            for (int k = 0; k < count; ++k) {

                const int max_nnz = 50;  // Maximum number of non-zero values
                double values[max_nnz];
                int indices[max_nnz];
                int nnz = 0;  // This will keep track of the number of non-zero values currently stored
                int localCount = 0;
                // Initialize variables for electron conditions
                double alpha = evaluate_alpha(E_i, AlphaMultiplier);
                double CrossSection = evaluate_cross_section_opt(E_i, CrossSectionLogArgMultiplier, CrossSectionNumorator, CrossSectionDenominatorA, CrossSectionDenominatorB);
                double PathLength = evaluate_path_length(CrossSection, PathLengthMultiplier);
                double RND_step = getRand(0.000001, 0.999999, &state);
                double step = evaluate_step(PathLength, RND_step);  

                // Initialize position data
                double position_data[7];
                initialize_positions(step, ProbeDiameter, position_data, &state);
                double cosineX = position_data[0];
                double cosineY = position_data[1];
                double cosineZ = position_data[2];
                double z0 = position_data[3];
                double y0 = position_data[4];
                double x0 = position_data[5];
                double E = E_i;
                bool condition = true;


                while (condition == true) {
                
                
                // Generate random numbers
                double RND_phi = getRand(0, 1, &state);
                double RND_step = getRand(0.000001, 0.999999, &state);
                double RND_pho = getRand(0, 1, &state);

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

                double ca = evaluate_direction_cosine_a(phi, psi, cosineX, cosineY, cosineZ);
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

                    if ((x0 <= x_dimension) && (x0 >= -x_dimension) &&
                        (y0 <= y_dimension) && (y0 >= -y_dimension)) {
                        int pixel_1D_index = static_cast<int>(i_coordinate) * perfect_image_1 + static_cast<int>(j_coordinate);
                        
                        bool found = false;
                        for (int i = 0; i < nnz; ++i) {
                            if (indices[i] == pixel_1D_index) {
                                values[i] += static_cast<double>(new_eh_pairs);
                                found = true;
                                break;
                            }
                        }
                        
                        if (!found && nnz < max_nnz) {
                            indices[nnz] = pixel_1D_index;
                            values[nnz] = static_cast<double>(new_eh_pairs);
                            nnz++;
                        }

                    }

                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    //      Scenario 6 : Electron deposits sufficient energy and moves positively in x and y.                       //
                    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                    else if ((x0 > x_dimension) && (y0 > y_dimension)) {
                        
                        // Electron moves positively in x and y

                        

                        // Calculate the translations
                        int translation_x = custom_round_down(x0 / (2 * x_dimension));
                        int translation_y = custom_round_down(y0 / (2 * y_dimension));

                        

                        if ((i_coordinate + translation_x <= perfect_image_0 - 1) && (j_coordinate + translation_y <= perfect_image_1 - 1)) {
                            // Update the value in the device array using an atomic add
                            
                            int pixel_1D_index = (i_coordinate + translation_x) * perfect_image_1 + (j_coordinate + translation_y);
                            

                            bool found = false;
                            for (int i = 0; i < nnz; ++i) {
                                if (indices[i] == pixel_1D_index) {
                                    values[i] += static_cast<double>(new_eh_pairs);
                                    found = true;
                                    break;
                                }
                            }
                            
                            if (!found && nnz < max_nnz) {
                                indices[nnz] = pixel_1D_index;
                                values[nnz] = static_cast<double>(new_eh_pairs);
                                nnz++;
                            }

                        }
                    }                 

                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    //      Scenario 7 : Electron deposits sufficient energy and moves negatively in x and y.                       //
                    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                    else if ((x0 < -x_dimension) && (y0 < -y_dimension)) {


                        // Calculate the translations
                        int translation_x = custom_round_down(x0 / (2 * x_dimension));
                        int translation_y = custom_round_down(y0 / (2 * y_dimension));

                        if ((i_coordinate - translation_x >= 0) && (j_coordinate - translation_y >= 0)) {

                        
                        int pixel_1D_index = (i_coordinate - translation_x) * perfect_image_1 + (j_coordinate - translation_y);

                            bool found = false;
                            for (int i = 0; i < nnz; ++i) {
                                if (indices[i] == pixel_1D_index) {
                                    values[i] += static_cast<double>(new_eh_pairs);
                                    found = true;
                                    break;
                                }
                            }
                            
                            if (!found && nnz < max_nnz) {
                                indices[nnz] = pixel_1D_index;
                                values[nnz] = static_cast<double>(new_eh_pairs);
                                nnz++;
                            }

                        }
                    }

                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    //      Scenario 8 : Electron deposits sufficient energy and moves positively in x only.                        //
                    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                    else if (x0 > x_dimension) {
                        int translation_x = custom_round_down(x0 / (2 * x_dimension));
                        if (i_coordinate + translation_x <= perfect_image_0 - 1) {
                        // Update the value in the device array using an atomic add
                        
                        int pixel_1D_index = (i_coordinate + translation_x) * perfect_image_1 + j_coordinate;
                        
                            bool found = false;
                            for (int i = 0; i < nnz; ++i) {
                                if (indices[i] == pixel_1D_index) {
                                    values[i] += static_cast<double>(new_eh_pairs);
                                    found = true;
                                    break;
                                }
                            }
                            
                            if (!found && nnz < max_nnz) {
                                indices[nnz] = pixel_1D_index;
                                values[nnz] = static_cast<double>(new_eh_pairs);
                                nnz++;
                            }

                        }
                    }

                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    //      Scenario 9 : Electron deposits sufficient energy and moves negatively in x only.                        //
                    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                    else if (x0 < -x_dimension) {
                        int translation_x = custom_round_down(x0 / (2 * x_dimension));
                        if (i_coordinate - translation_x >= 0) {
                        // Update the value in the device array using an atomic add
                        
                        int pixel_1D_index = (i_coordinate - translation_x) * perfect_image_1 + j_coordinate;

                            bool found = false;
                            for (int i = 0; i < nnz; ++i) {
                                if (indices[i] == pixel_1D_index) {
                                    values[i] += static_cast<double>(new_eh_pairs);
                                    found = true;
                                    break;
                                }
                            }
                            
                            if (!found && nnz < max_nnz) {
                                indices[nnz] = pixel_1D_index;
                                values[nnz] = static_cast<double>(new_eh_pairs);
                                nnz++;
                            }

                        }
                        }

                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    //      Scenario 10 : Electron deposits sufficient energy and moves positively in y only.                       //
                    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                    else if (y0 > y_dimension) {
                        int translation_y = custom_round_down(y0 / (2 * y_dimension));
                        if (j_coordinate + translation_y <= perfect_image_1 - 1) {
                        // Update the value in the device array using an atomic add
                        
                        int pixel_1D_index = i_coordinate * perfect_image_1 + (j_coordinate + translation_y);
                        

                        bool found = false;
                        for (int i = 0; i < nnz; ++i) {
                            if (indices[i] == pixel_1D_index) {
                                values[i] += static_cast<double>(new_eh_pairs);
                                found = true;
                                break;
                            }
                        }
                            
                            if (!found && nnz < max_nnz) {
                                indices[nnz] = pixel_1D_index;
                                values[nnz] = static_cast<double>(new_eh_pairs);
                                nnz++;
                            }

                        }    
                           
                        }
                    

                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    //      Scenario 11 : Electron deposits sufficient energy and moves negatively in y only.                       //
                    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                    else if (y0 < -y_dimension) {
                        int translation_y = custom_round_down(y0 / (2 * y_dimension));
                        if (j_coordinate - translation_y >= 0) {
                        // Update the value in the device array using an atomic add
                        
                        int pixel_1D_index = i_coordinate * perfect_image_1 + (j_coordinate - translation_y);
                        

                        bool found = false;
                        for (int i = 0; i < nnz; ++i) {
                            if (indices[i] == pixel_1D_index) {
                                values[i] += static_cast<double>(new_eh_pairs);
                                found = true;
                                break;
                            }
                        }
                            
                            if (!found && nnz < max_nnz) {
                                indices[nnz] = pixel_1D_index;
                                values[nnz] = static_cast<double>(new_eh_pairs);
                                nnz++;
                            }


                        }
                    }
                
                }  // Sufficient energy loop
            } // wile loop

            for (int i = 0; i < nnz; ++i) { // nnz is the number of non-zero entries 
                int globalIndex = indices[i];
                values[i] = floor(dE_threshold * values[i] / t_counting);
                atomicAdd_double(&(new_image_MCS[globalIndex]), values[i]);
                
            }
     

            } //electron level
        
        } //else

    } //pixel_lop


}

} // c linkage


