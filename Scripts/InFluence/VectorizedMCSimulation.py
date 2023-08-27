import math
import numpy as np
import random
from tqdm import tqdm
from parameters import params
from ElectronObjectConstructor import ElectronObject
from HelperFunctions import custom_round_down, CalculateConstants, evaluate_alpha, evaluate_cross_section_opt, evaluate_path_length, evaluate_step, evaluate_phi, evaluate_pho, evaluate_direction_cosine_a, evaluate_direction_cosine_b, evaluate_direction_cosine_c, evaluate_energy_loss_rate, initialise_postions
from ElectronArrayConstructor import generate_electrons
import numba 



def Single_Pixel(pixel, AlphaMultiplier, CrossSectionNumorator, CrossSectionLogArgMultiplier, CrossSectionDenominatorA, CrossSectionDenominatorB, PathLengthMultiplier, EnergyLossMultiplierA, EnergyLossMultiplierB, E_i, ProbeDiameter, MinimumEnergy, pixel_dimensions, dE_threshold, perfect_image_0, perfect_image_1, Density, t_counting):
   

        def generate_electrons_from_pixels(pixels, E_i, ProbeDiameter):
            electron_data = []

            for pixel in pixels:
                electron_count = int(pixel[0])  # Extract electron count from the first column
                electrons = generate_electrons(electron_count, E_i, ProbeDiameter)
                pixel_data = np.repeat([pixel], electron_count, axis=0)
                electrons_with_pixels = np.hstack((electrons, pixel_data))
                electron_data.append(electrons_with_pixels)
            result = np.concatenate(electron_data, axis=0)
            return result

        pixel = pixel[pixel[:, 0] != 0]
        

        electrons = generate_electrons_from_pixels(pixel, E_i, ProbeDiameter)

         #electron array index        
            # 0 : E_i
            # 1 : E
            # 2 : Alpha
            # 3 : CrossSection
            # 4 : PathLength
            # 5 : RND_step
            # 6 : step
            # 7 : ProbeDiameter
            # 8 : InitialPosition
            # 9 : cosineX
            # 10 : cosineY
            # 11 : cosineZ
            # 12 : z0
            # 13 : y0
            # 14 : x0
            # 15 : dE
            # 16 : count        pixel[:, 0]      electrons[:,16]
            # 17 : i            pixel[:, 1]      electrons[:,17] 
            # 18 : j            pixel[:, 2]      electrons[:,18]
            # 19 : Dimension along the x-axis   pixel[:, 3]   electrons[:,19]
            # 20 : Dimension along the y-axis   pixel[:, 4]   electrons[:,20]
            # 21 : Dimension along the z-axis   pixel[:, 5]   electrons[:,21]

        
        eh_charge_counter = np.zeros((perfect_image_0, perfect_image_1), dtype = np.float64)
        ElectronHoleChargeCounter = np.zeros((perfect_image_0, perfect_image_1), dtype = np.float64)

        new_image_MCS = np.zeros((perfect_image_0, perfect_image_1), dtype=np.float64)

    
        number_stopped = 0
        number_backscattered = 0
        number_transmitted = 0
        number_eh_pairs = 0
        
        pi = math.pi

        
        loop_count = 0
        while electrons.size> 0:  
            loop_count += 1

            print(loop_count)
            

                                                                                                                                                                     # Monte Carlo Loop until backscatter, trasmission or stopping
                    
            RND_phi = random.uniform(a=0, b=1)                                                                                                                                                                   # generate random number for the phi angle
            RND_step = random.uniform(a = 0.000001, b = 0.999999)
            RND_pho = random.uniform(a = 0, b = 1)

            alpha = evaluate_alpha(electrons[:,1], AlphaMultiplier)                                                                                                                                   # calc screening constant, function of previous energy
            cross = evaluate_cross_section_opt(electrons[:,1], CrossSectionLogArgMultiplier, CrossSectionNumorator, CrossSectionDenominatorA, CrossSectionDenominatorB)                                             # calc cross section
            path_length = evaluate_path_length(electrons[:,3], PathLengthMultiplier)                                                                                                   # calc mean free path length
            step = evaluate_step(path_length, RND_step)                                                                                                                                                  # calculate step of this iteration
            electrons[:,15] = step * Density * evaluate_energy_loss_rate(electrons[:,1], EnergyLossMultiplierA, EnergyLossMultiplierB)
            electrons[:,1] = electrons[:,1] + electrons[:,15]                                                                                                                                                     # calc new energy #separate out dE/dS
            phi = evaluate_phi(RND_phi, alpha=alpha)                                                                                                                                                     # calc scattering angle
            psi = evaluate_pho(RND_pho, pi)                                                                                                                                                              # calc other scattering angle
            ca = evaluate_direction_cosine_a(phi, psi, electrons[:,9], electrons[:,10], electrons[:,11])                                                                             # calc direction cosines
            cb = evaluate_direction_cosine_b(phi, psi, electrons[:,9], electrons[:,10], electrons[:,11])
            cc = evaluate_direction_cosine_c(phi, psi, electrons[:,11])
            electrons[:,14] = electrons[:,14] + step * ca                                                                                                                                                 # find and reset to new positions
            electrons[:,13] = electrons[:,13] + step * cb
            electrons[:,12] = electrons[:,12] + step * cc
            electrons[:,9] = ca                                                                                                                                                                        # reset direction cosines
            electrons[:,10] = cb
            electrons[:,11] = cc


            #  Scenario 1: Electron stops in material
            
            ElectronStopsInMaterial_mask = electrons[:,1] <= MinimumEnergy    # Make a mask
            number_stopped += np.sum(ElectronStopsInMaterial_mask)              # If true for any rows, add the number of rows to counter
            electrons = electrons[~ElectronStopsInMaterial_mask]                # Remove the electrons that stop in the material


            # Scenario 2: Electron backscatters
            BackscatteredElectrons_mask = electrons[:,12] < 0.01
            number_backscattered += np.sum(BackscatteredElectrons_mask)
            electrons = electrons[~BackscatteredElectrons_mask]

            # Scenario 3: Electron penetrates material
            PenetratingElectrons_mask = electrons[:,12]> pixel_dimensions[2]
            number_transmitted += np.sum(PenetratingElectrons_mask)
            electrons = electrons[~PenetratingElectrons_mask]

            # Scenario 4: Electron deposits sufficient energy
            SufficientEnergy_mask = -1 * electrons[:,15] >= dE_threshold
            new_eh_pairs = np.floor(-1 * electrons[:,15][SufficientEnergy_mask] / dE_threshold).astype(int)
            number_eh_pairs += np.sum(new_eh_pairs)

            # Select those electrons with sufficient energy 
            electrons_SufficientEnergy = electrons[SufficientEnergy_mask]


            #______________________________________________________________________________________________________________________

            # Scenario 5: Electron has sufficient energy and stays within pixel boundary

            # Create mask where electrons stay within the pixel region
            within_pixel_mask = (
            (electrons_SufficientEnergy[:, 14] <= 1 * electrons_SufficientEnergy[:, 19]) &
            (electrons_SufficientEnergy[:, 14] >= -1 * electrons_SufficientEnergy[:, 19]) &
            (electrons_SufficientEnergy[:, 13] <= 1 * electrons_SufficientEnergy[:, 20]) &
            (electrons_SufficientEnergy[:, 13] >= -1 * electrons_SufficientEnergy[:, 20])
            )

            if not np.any(within_pixel_mask):
                pass
            else:
                eh_charge_counter[electrons_SufficientEnergy[:, 17][within_pixel_mask].astype(int), electrons_SufficientEnergy[:, 18][within_pixel_mask].astype(int)] += np.sum(new_eh_pairs[within_pixel_mask])
            #______________________________________________________________________________________________________________________

            # Scenario 5: Electron has sufficient energy and moves positively in x and y

            # Create mask where electrons move positively in x and y
            positive_move_mask = (
            (electrons_SufficientEnergy[:, 14] > 1 * electrons_SufficientEnergy[:, 19]) &
            (electrons_SufficientEnergy[:, 13] > 1 * electrons_SufficientEnergy[:, 20])
            )
            
            if not np.any(positive_move_mask):
                pass
            else:
                # Round the distance the electron moved in terms of pixel units
                translation_x = custom_round_down(electrons_SufficientEnergy[:, 14][positive_move_mask] / (2 * electrons_SufficientEnergy[:, 19][positive_move_mask]))
                translation_y = custom_round_down(electrons_SufficientEnergy[:, 13][positive_move_mask] / (2 * electrons_SufficientEnergy[:, 20][positive_move_mask]))

                # Convert to integer indices
                translation_x = translation_x.astype(int)
                translation_y = translation_y.astype(int)

                # Create mask where updated coordinates are within the image dimensions
                within_image_mask = (
                (electrons_SufficientEnergy[:, 17][positive_move_mask] + translation_x <= perfect_image_0 - 1) &
                (electrons_SufficientEnergy[:, 18][positive_move_mask] + translation_y <= perfect_image_1 - 1)
                )

                # Update the electron-hole charge counter for the pixel
                for x, y in zip(translation_x[within_image_mask], translation_y[within_image_mask]):
                    eh_charge_counter[electrons_SufficientEnergy[:, 17][positive_move_mask][within_image_mask].astype(int) + x, electrons_SufficientEnergy[:, 18][positive_move_mask][within_image_mask].astype(int) + y] += new_eh_pairs[positive_move_mask][within_image_mask]


            
            #______________________________________________________________________________________________________________________#

            #    Scenario 6: Electron has sufficient energy and moves negatively in x and y
                    
            # Create mask where electrons move negatively in x and y
            negative_move_mask = (
                (electrons_SufficientEnergy[:, 14] < -1 * electrons_SufficientEnergy[:, 19]) &
                (electrons_SufficientEnergy[:, 13] < -1 * electrons_SufficientEnergy[:, 20])
            )

            if not np.any(negative_move_mask):
                pass
            else:

                # Round the distance the electron moved in terms of pixel units
                translation_x = custom_round_down(electrons_SufficientEnergy[:, 14][negative_move_mask] / (2 * electrons_SufficientEnergy[:, 19][negative_move_mask]))
                translation_y = custom_round_down(electrons_SufficientEnergy[:, 13][negative_move_mask] / (2 * electrons_SufficientEnergy[:, 20][negative_move_mask]))

                # Convert to integer indices
                translation_x = translation_x.astype(int)
                translation_y = translation_y.astype(int)

                # Create mask where updated coordinates are within the image dimensions
                within_image_mask = (
                    (electrons_SufficientEnergy[:, 17][negative_move_mask] - translation_x >= 0) &
                    (electrons_SufficientEnergy[:, 18][negative_move_mask] - translation_y >= 0)
                )


                # Update the electron-hole charge counter for the pixel
                for x, y in zip(translation_x[within_image_mask], translation_y[within_image_mask]):
                    eh_charge_counter[electrons_SufficientEnergy[:, 17][negative_move_mask][within_image_mask].astype(int) - x, electrons_SufficientEnergy[:, 18][negative_move_mask][within_image_mask].astype(int) - y] += new_eh_pairs[negative_move_mask][within_image_mask]

            #______________________________________________________________________________________________________________________#
            #    Scenario 7: Electron has sufficient energy and moves positively in x and negatively in y                           
              
            # Create mask where electrons move positively in x and negatively in y
            positive_negative_move_mask = (
            (electrons_SufficientEnergy[:, 14] > 1 * electrons_SufficientEnergy[:, 19]) &
            (electrons_SufficientEnergy[:, 13] < -1 * electrons_SufficientEnergy[:, 20])
            )

            if not np.any(positive_negative_move_mask):
                pass
            else:
                # Round the distance the electron moved in terms of pixel units
                translation_x = custom_round_down(electrons_SufficientEnergy[:, 14][positive_negative_move_mask] / (2 * electrons_SufficientEnergy[:, 19][positive_negative_move_mask]))
                translation_y = custom_round_down(electrons_SufficientEnergy[:, 13][positive_negative_move_mask] / (2 * electrons_SufficientEnergy[:, 20][positive_negative_move_mask]))

                # Convert to integer indices
                translation_x = translation_x.astype(int)
                translation_y = translation_y.astype(int)

                # Create mask where updated coordinates are within the image dimensions
                within_image_mask = (
                (electrons_SufficientEnergy[:, 17][positive_negative_move_mask] + translation_x <= perfect_image_0 - 1) &
                (electrons_SufficientEnergy[:, 18][positive_negative_move_mask] - translation_y >= 0)
                )

                # Update the electron-hole charge counter for the pixel
                for x, y in zip(translation_x[within_image_mask], translation_y[within_image_mask]):
                    eh_charge_counter[electrons_SufficientEnergy[:, 17][positive_negative_move_mask].astype(int) + x, electrons_SufficientEnergy[:, 18][positive_negative_move_mask].astype(int) - y] += new_eh_pairs[positive_negative_move_mask][within_image_mask]

            #______________________________________________________________________________________________________________________#
            #    Scenario 8: Electron has sufficient energy and moves negatively in x and positively in y
            
            # Create mask where electrons move negatively in x and positively in y
            neg_pos_move_mask = (
            (electrons_SufficientEnergy[:, 14] < -1 * electrons_SufficientEnergy[:, 19]) &
            (electrons_SufficientEnergy[:, 13] > 1 * electrons_SufficientEnergy[:, 20])
            )

            if not np.any(neg_pos_move_mask):
                pass
            else:

                # Round the distance the electron moved in terms of pixel units
                translation_x = np.floor(electrons_SufficientEnergy[:, 14][neg_pos_move_mask] / (2 * electrons_SufficientEnergy[:, 19][neg_pos_move_mask]))
                translation_y = np.floor(electrons_SufficientEnergy[:, 13][neg_pos_move_mask] / (2 * electrons_SufficientEnergy[:, 20][neg_pos_move_mask]))


                # Convert to integer indices
                translation_x = translation_x.astype(int)
                translation_y = translation_y.astype(int)

                # Create mask where updated coordinates are within the image dimensions
                within_image_mask = (
                (electrons_SufficientEnergy[:, 17][neg_pos_move_mask] - translation_x >= 0) &
                (electrons_SufficientEnergy[:, 18][neg_pos_move_mask] + translation_y <= perfect_image_0 - 1)
                )

                # Update the electron-hole charge counter for the pixel
                for x, y in zip(translation_x[within_image_mask], translation_y[within_image_mask]):
                    eh_charge_counter[electrons_SufficientEnergy[:, 17][neg_pos_move_mask].astype(int) - x, electrons_SufficientEnergy[:, 18][neg_pos_move_mask].astype(int) + y] += new_eh_pairs[neg_pos_move_mask][within_image_mask]

            #______________________________________________________________________________________________________________________#
            #   Scenario 9: Electron has sufficient energy and moves positively in x only

            # Create mask where electrons move positively in x only
            positive_x_move_mask = electrons_SufficientEnergy[:, 14] > 1 * electrons_SufficientEnergy[:, 19]

            if not np.any(positive_x_move_mask):
                pass
            else:

                # Round the distance the electron moved in terms of pixel units
                translation_x = np.floor(electrons_SufficientEnergy[:, 14][positive_x_move_mask] / (2 * electrons_SufficientEnergy[:, 19][positive_x_move_mask]))

                # Convert to integer indices
                translation_x = translation_x.astype(int)

                # Create mask where updated coordinates are within the image dimensions
                within_image_mask = (electrons_SufficientEnergy[:, 17][positive_x_move_mask] + translation_x <= perfect_image_0 - 1)

                # Update the electron-hole charge counter for the pixel
                for x in translation_x[within_image_mask]:
                    eh_charge_counter[electrons_SufficientEnergy[:, 17][positive_x_move_mask].astype(int) + translation_x[within_image_mask], electrons_SufficientEnergy[:, 18][positive_x_move_mask].astype(int)] += new_eh_pairs[positive_x_move_mask][within_image_mask]


            #______________________________________________________________________________________________________________________#
            #   Scenario 10: Electron has sufficient energy and moves negatively in x only

            # Create mask where electrons move negatively in x only
            negative_x_move_mask = electrons_SufficientEnergy[:,14] < -1*electrons_SufficientEnergy[:,19]

            if not np.any(negative_x_move_mask):
                pass
            else:

                # Round the distance the electron moved in terms of pixel units
                translation_x = np.floor(-1 * electrons_SufficientEnergy[:,14][negative_x_move_mask] / (2*electrons_SufficientEnergy[:,19][negative_x_move_mask]))

                # Convert to integer indices
                translation_x = translation_x.astype(int)

                # Create mask where updated coordinates are within the image dimensions
                within_image_mask = (electrons_SufficientEnergy[:, 17][negative_x_move_mask] - translation_x >= 0)


                # Update the electron-hole charge counter for the pixel
                for x in translation_x[within_image_mask]:
                    eh_charge_counter[electrons_SufficientEnergy[:, 17][negative_x_move_mask].astype(int) - translation_x[within_image_mask], electrons_SufficientEnergy[:, 18][negative_x_move_mask].astype(int)] += new_eh_pairs[negative_x_move_mask][within_image_mask]


            #______________________________________________________________________________________________________________________#
            #   Scenario 11: Electron has sufficient energy and moves positively in y only
             
            # Create mask where electrons move positively in y only
            positive_y_move_mask = electrons_SufficientEnergy[:,13] > 1*electrons_SufficientEnergy[:,20]

            if not np.any(positive_y_move_mask):
                pass
            else:

                # Round the distance the electron moved in terms of pixel units
                translation_y = np.floor(electrons_SufficientEnergy[:,13][positive_y_move_mask] / (2*electrons_SufficientEnergy[:, 20][positive_y_move_mask]))

                # Convert to integer indices
                translation_y = translation_y.astype(int)

                # Create mask where updated coordinates are within the image dimensions
                within_image_mask = (electrons_SufficientEnergy[:,18][positive_y_move_mask] + translation_y <= perfect_image_0 - 1)

                # Update the electron-hole charge counter for the pixel
                for y in translation_y[within_image_mask]:
                    eh_charge_counter[electrons_SufficientEnergy[:, 17][positive_y_move_mask].astype(int), electrons_SufficientEnergy[:, 18][positive_y_move_mask].astype(int) + translation_y[within_image_mask]] += new_eh_pairs[positive_y_move_mask][within_image_mask]
                

            #______________________________________________________________________________________________________________________#
            #    Scenario 12: Electron has sufficient energy and moves negatively in y only
            

            # Create mask where electrons move negatively in y only
            negative_y_move_mask = electrons_SufficientEnergy[:,13] < -1*electrons_SufficientEnergy[:,20]

            if not np.any(negative_y_move_mask):
                pass
            else:

                # Round the distance the electron moved in terms of pixel units
                translation_y = np.floor(electrons_SufficientEnergy[:,13][negative_y_move_mask] / -(2*electrons_SufficientEnergy[:,20][negative_y_move_mask]))

                # Convert to integer indices
                translation_y = translation_y.astype(int)

                # Create mask where updated coordinates are within the image dimensions
                within_image_mask = (electrons_SufficientEnergy[:,18][negative_y_move_mask] - translation_y >= 0)


                # Update the electron-hole charge counter for the pixel
                for y in translation_y[within_image_mask]:
                    eh_charge_counter[electrons_SufficientEnergy[:,17][negative_y_move_mask].astype(int), electrons_SufficientEnergy[:,18][negative_y_move_mask].astype(int) - y] += new_eh_pairs[negative_y_move_mask][within_image_mask]



        eh_charge_counter = (np.floor(dE_threshold*eh_charge_counter/t_counting))
        new_image_MCS += eh_charge_counter  

        nonzero_indices = np.nonzero(new_image_MCS)
        nonzero_elements = new_image_MCS[nonzero_indices]

        print(nonzero_elements)

        return new_image_MCS

                

            
def VectorizedMCScatteringSimulation(pixels, E_i, ProbeDiameter, MinimumEnergy, pixel_dimensions, dE_threshold, perfect_image_0, perfect_image_1, Density, t_counting):

    
    # Create a 2D array from the 'pixels' list
    pixels_array = np.vstack((
        [pixel.electron_count for pixel in pixels],      # Column 0: Number of electrons in the pixel
        [pixel.i for pixel in pixels],                   # Column 1: Pixel index i
        [pixel.j for pixel in pixels],                   # Column 2: Pixel index j
        [pixel.x_dimension for pixel in pixels],         # Column 3: Dimension along the x-axis      
        [pixel.y_dimension for pixel in pixels],         # Column 4: Dimension along the y-axis
        [pixel.z_dimension for pixel in pixels]          # Column 5: Dimension along the z-axis
    )).T  # Transpose to get the desired shape (rows represent pixels, columns represent attributes)

    # Convert the data type of the first three columns to int and the last three columns to float
    

    print('Pixel data prepared ...')

    print(pixels_array[1])

    AlphaMultiplier, CrossSectionNumorator, CrossSectionLogArgMultiplier, CrossSectionDenominatorA, CrossSectionDenominatorB, PathLengthMultiplier, EnergyLossMultiplierA, EnergyLossMultiplierB = CalculateConstants(params)

    pi = math.pi

    print('Beginning main simulation loop...')


    result = Single_Pixel(pixels_array, AlphaMultiplier, CrossSectionNumorator, CrossSectionLogArgMultiplier, CrossSectionDenominatorA, CrossSectionDenominatorB, PathLengthMultiplier, EnergyLossMultiplierA, EnergyLossMultiplierB, E_i, ProbeDiameter, MinimumEnergy, pixel_dimensions, dE_threshold, perfect_image_0, perfect_image_1, Density, t_counting)
                    

    return result
            



 #electron array index        
            # 0 : E_i
            # 1 : E
            # 2 : Alpha
            # 3 : CrossSection
            # 4 : PathLength
            # 5 : RND_step
            # 6 : step
            # 7 : ProbeDiameter
            # 8 : InitialPosition
            # 9 : cosineX
            # 10 : cosineY
            # 11 : cosineZ
            # 12 : z0
            # 13 : y0
            # 14 : x0