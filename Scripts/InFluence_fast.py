import time
from params import *
from common_functions import *
import simplejson
import matplotlib.pyplot as plt
from numba_progress import ProgressBar
from numba.typed import List
import warnings

#warnings.filterwarnings("ignore")

#Main function
@jit(nopython = True, parallel = True)
def MCS_ChargeCounting(minimum_energy, t_counting, dE_threshold, pixel_dimensions, pixel_information, distrib, perfect_image):
    # counters
    # load coordinates from file
    number_transmitted = 0
    number_backscattered = 0
    number_stopped = 0
    number_eh_pairs = 0
    new_image_MCS = np.zeros((perfect_image.shape[0], perfect_image.shape[1]), dtype=np.float64)

    for pixel_counter in prange(len(pixel_information)):
            number_of_electrons = int(pixel_information[pixel_counter][0])
            i_coordinate = int(pixel_information[pixel_counter][1])
            j_coordinate = int(pixel_information[pixel_counter][2])

            for _ in range(1, number_of_electrons + 1):
                eh_charge_counter = np.zeros((perfect_image.shape[0], perfect_image.shape[1]), dtype=np.float64)
                #print(eh_charge_counter.sum())
                # initial conditions
                random_index = np.random.randint(len(distrib))
                random_traj_from_list = distrib[random_index]
                condition = True
                q_counter = 0
                while condition: #Monte Carlo Loop until backscatter, trasmission or stopping

                    z0 = random_traj_from_list[2][q_counter]
                    y0 = random_traj_from_list[1][q_counter]
                    x0 = random_traj_from_list[0][q_counter]
                    E = random_traj_from_list[3][q_counter]
                    dE = random_traj_from_list[3][q_counter+1] - random_traj_from_list[3][q_counter]

                    if E <= minimum_energy: #if electron stops in the material
                        number_stopped = number_stopped + 1
                        condition = False
                    if z0 < 10**-2: #if electron backscatters
                        number_backscattered = number_backscattered + 1
                        condition = False #not sure how to deal with electron scattering outside of material
                    if z0 > pixel_dimensions[2]:  # if electron penetrates material #change to have de/ds as the counting
                        # thin
                        number_transmitted = number_transmitted + 1
                        condition = False
                    if -1*dE >= dE_threshold: #if electron deposits sufficient energyo
                        new_eh_pairs = math.floor(-1*dE/dE_threshold)
                        number_eh_pairs = int(number_eh_pairs + new_eh_pairs)

                        if (x0 <= 1*pixel_dimensions[0]) and (x0 >= -1*pixel_dimensions[0]) and (y0 <= 1*pixel_dimensions[
                            1])and (y0 >=-1*pixel_dimensions[1]): #electron stays within pixel region

                            eh_charge_counter[i_coordinate, j_coordinate] = eh_charge_counter[i_coordinate, j_coordinate]\
                                                                            + new_eh_pairs

                        elif (x0 > 1*pixel_dimensions[0]) and (y0 > 1*pixel_dimensions[1]): #electron moves positively in x
                            # and y
                            translation_x = round_half_down(x = x0/(2*pixel_dimensions[0]))
                            translation_y = round_half_down(x =y0/(2*pixel_dimensions[1]))
                            if i_coordinate + translation_x  <= perfect_image.shape[0] - 1 and j_coordinate + \
                                    translation_y <= perfect_image.shape[0] - 1:
                                    eh_charge_counter[i_coordinate + translation_x, j_coordinate + translation_y] = eh_charge_counter[
                                                                                                                    i_coordinate +
                                                                                                                    translation_x,
                                                                                                                    j_coordinate +
                                                                                                                    translation_y] + new_eh_pairs


                        elif (x0 < -1*pixel_dimensions[0]) and (y0 < -1*pixel_dimensions[1]): #electron moves negatively in x
                            # and y
                            translation_x = round_half_down(x=x0 / (2 * pixel_dimensions[0]))
                            translation_y = round_half_down(x=y0 / (2 * pixel_dimensions[1]))
                            if i_coordinate - translation_x >= 0 and j_coordinate - translation_y >= 0:
                                    eh_charge_counter[i_coordinate - translation_x, j_coordinate - translation_y] = eh_charge_counter[
                                                                                                                   i_coordinate -
                                                                                                                   translation_x,
                                                                                                                   j_coordinate -
                                                                                                                   translation_y] + new_eh_pairs



                        elif x0 > (1*pixel_dimensions[0]) and y0 < (-1*pixel_dimensions[1]):  # electron moves psoitively in
                            # x and negatively in y
                            translation_x = round_half_down(x=x0 / (2 * pixel_dimensions[0]))
                            translation_y = round_half_down(x=y0 / (2 * pixel_dimensions[1]))
                            if j_coordinate - translation_y >= 0 and i_coordinate + translation_x  <= \
                                    perfect_image.shape[0] - 1:
                                    eh_charge_counter[i_coordinate + translation_x, j_coordinate - translation_y] = eh_charge_counter[
                                                                                                                    i_coordinate +
                                                                                                                    translation_x,
                                                                                                                    j_coordinate -
                                                                                                                    translation_y] + new_eh_pairs

                        elif x0 < (-1*pixel_dimensions[0]) and y0 > (1*pixel_dimensions[1]):  # electron moves negatively in
                            # x and positively in y
                            translation_x = round_half_down(x=x0 / (2 * pixel_dimensions[0]))
                            translation_y = round_half_down(x=y0 / (2 * pixel_dimensions[1]))
                            if i_coordinate - translation_x >= 0 and j_coordinate + \
                                    translation_y <= perfect_image.shape[0] - 1:
                                    eh_charge_counter[i_coordinate - translation_x, j_coordinate + translation_y] = eh_charge_counter[
                                                                                                                    i_coordinate -
                                                                                                                    translation_x,
                                                                                                                    j_coordinate +
                                                                                                                    translation_y] + new_eh_pairs

                        elif x0 > (1*pixel_dimensions[0]): #electron moves in positive x direction
                            translation_x = round_half_down(x=x0 / (2 * pixel_dimensions[0]))
                            if i_coordinate + \
                                    translation_x <= perfect_image.shape[0] - 1:

                                    eh_charge_counter[i_coordinate + translation_x, j_coordinate] = eh_charge_counter[i_coordinate +
                                                                                                                      translation_x,
                                                                                                                      j_coordinate] + new_eh_pairs

                        elif x0 < (-1*pixel_dimensions[0]): #electron moves in negative x direction
                            translation_x = round_half_down(x=x0 / (2 * pixel_dimensions[0]))
                            if i_coordinate - translation_x >= 0:

                                    eh_charge_counter[i_coordinate - translation_x, j_coordinate] = eh_charge_counter[i_coordinate -
                                                                                                                      translation_x,
                                                                                                                      j_coordinate] + new_eh_pairs

                        elif y0 > (1*pixel_dimensions[1]): #electron moves in positive y direction
                            translation_y = round_half_down(x=y0 / (2 * pixel_dimensions[1]))
                            if j_coordinate + \
                                    translation_y <= perfect_image.shape[0] - 1:

                                    eh_charge_counter[i_coordinate, j_coordinate + translation_y] = eh_charge_counter[i_coordinate,
                                                                                                                      j_coordinate +
                                                                                                                      translation_y] + new_eh_pairs

                        elif y0 < (-1*pixel_dimensions[0]): #electron moves in negative y direction
                            translation_y = round_half_down(x=y0 / (2 * pixel_dimensions[1]))

                            if j_coordinate - translation_y >= 0:

                                    eh_charge_counter[i_coordinate, j_coordinate - translation_y] = eh_charge_counter[i_coordinate,
                                                                                                                      j_coordinate -
                                                                                                                      translation_y] + new_eh_pairs
                    q_counter +=1

                eh_charge_counter_norm = (np.floor(dE_threshold*eh_charge_counter/t_counting))

                #print('Charge counter image is', eh_charge_counter_norm.sum())
                new_image_MCS += eh_charge_counter_norm #+ new_image_MCS
    print('Final image is', new_image_MCS.sum())
    return new_image_MCS





if __name__ == '__main__':
    # image
    if path_to_image.endswith('.dat'):  # opens Dr Probe image
        perfect_image = opendatfile(image=path_to_image)
    elif path_to_image.endswith('.tif'):
        from PIL import Image
        perfect_image = Image.open(path_to_image)
        perfect_image = np.array(perfect_image)
    else:  # opens npy files
        perfect_image = np.load(path_to_image)  
     
    chosen_pixels = find_distribution(perfect_image, dose)
    perfect_image = num_samples * distribute_electrons(perfect_image, chosen_pixels)
    
    # output image
    with open(path_to_unmodulated_image, 'wb') as f:#############
        np.save(f, perfect_image/num_samples)

    pixel_information = get_pixel_info(image=perfect_image)
    pixel_information, _ = make_2D_array(pixel_information)



    print('\nInFluence is running...')

    start = time.time()
    with open(sample_distribution_file, 'r') as f:
        distribution_of_trajectories_arroflists = simplejson.load(f)
    distribution_of_trajectories = List([List(xi) for xi in distribution_of_trajectories_arroflists])
    print('Time to load distribution {}'.format(time.time() - start))
    modulated_image = np.uint32(MCS_ChargeCounting(minimum_energy, t_counting, dE_threshold,  pixel_dimensions,pixel_information, distrib = distribution_of_trajectories, perfect_image = perfect_image))
    print('\nInFluence has generated the modulated image, saved to: ', path_to_modulated_image, '\n\nRun time: ', time.time() - start)
    #output image
    with open(path_to_modulated_image, 'wb') as f:
        np.save(f, modulated_image/num_samples)


