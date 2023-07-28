import time

import numpy as np

from params import *
from common_functions import *

from tqdm import tqdm
import matplotlib.pyplot as plt
from numba_progress import ProgressBar


#Main function
#@jit(nopython = True, parallel =True)
def MCS_ChargeCounting(E_i, minimum_energy, t_counting, dE_threshold, N, Z, A, p, pixel_dimensions, pixel_information, perfect_image):
    # counters
    number_transmitted = 0
    number_backscattered = 0
    number_stopped = 0
    number_eh_pairs = 0
    new_image_MCS = np.zeros((perfect_image.shape[0], perfect_image.shape[1]), dtype=np.float64)
    for pixel_counter in range(len(pixel_information)):
            number_of_electrons = int(pixel_information[pixel_counter][0])
            i_coordinate = int(pixel_information[pixel_counter][1])
            j_coordinate = int(pixel_information[pixel_counter][2])
            for _ in range(1, number_of_electrons + 1):
                eh_charge_counter = np.zeros((perfect_image.shape[0], perfect_image.shape[1]), dtype = np.float64)
                # initial conditions
                alpha = evaluate_alpha(E_i, Z)
                cross = evaluate_cross_section(E_i, Z, alpha)
                path_length = evaluate_path_length(A, N, p, cross)
                RND_step = RND(a=0.000001, b=0.999999)
                step = evaluate_step(path_length = path_length, RND=RND_step)
                ip = initialise_postions(step=step, d = pixel_dimensions[0]) #d = probe diameter
                cx = ip[0]
                cy = ip[1]
                cz = ip[2]
                z0 = ip[3]
                y0 = ip[4]
                x0 = ip[5]
                condition = True
                E = E_i
                while condition: #Monte Carlo Loop until backscatter, trasmission or stopping
                    RND_phi = RND(a=0, b=1)  # generate random number for the phi angle
                    RND_step = RND(a = 0.000001, b = 0.999999)
                    RND_pho = RND(a = 0, b = 1)

                    alpha = evaluate_alpha(E, Z)  # calc screening constant, function of previous energy
                    cross = evaluate_cross_section(E, Z, alpha)  # calc cross section
                    path_length = evaluate_path_length(A, N, p, cross)  # calc mean free path length
                    step = evaluate_step(path_length=path_length, RND=RND_step)  # calculate step of this iteration
                    dE = step * p * evaluate_energy_loss_rate(E, Z, A)
                    E = E + dE  # calc new energy #separate out dE/dS
                    phi = evaluate_phi(RND=RND_phi, alpha=alpha)  # calc scattering angle
                    psi = evaluate_pho(RND=RND_pho)  # calc other scattering angle
                    ca = evaluate_direction_cosine_a(phi, psi, cx, cy, cz)  # calc direction cosines
                    cb = evaluate_direction_cosine_b(phi, psi, cx, cy, cz)
                    cc = evaluate_direction_cosine_c(phi, psi, cz)
                    x0 = x0 + step * ca  # find and reset to new positions
                    y0 = y0 + step * cb
                    z0 = z0 + step * cc
                    cx = ca  # reset direction cosines
                    cy = cb
                    cz = cc

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
###########################################################################################################################
                        if (x0 <= 1*pixel_dimensions[0]) and (x0 >= -1*pixel_dimensions[0]) and (y0 <= 1*pixel_dimensions[
                            1])and (y0 >=-1*pixel_dimensions[1]): #electron stays within pixel region

                            eh_charge_counter[i_coordinate, j_coordinate] = new_eh_pairs + eh_charge_counter[i_coordinate, j_coordinate]

                        elif (x0 > 1*pixel_dimensions[0]) and (y0 > 1*pixel_dimensions[1]): #electron moves positively in x
                            # and y
                            translation_x = round_half_down(x = x0/(2*pixel_dimensions[0]))
                            translation_y = round_half_down(x =y0/(2*pixel_dimensions[1]))
                            if i_coordinate + translation_x  <= perfect_image.shape[0] - 1 and j_coordinate + \
                                    translation_y <= perfect_image.shape[0] - 1:
                                    eh_charge_counter[i_coordinate + translation_x, j_coordinate + translation_y] = new_eh_pairs + eh_charge_counter[i_coordinate +translation_x,j_coordinate +translation_y]

###########################################################################################################################

                        elif (x0 < -1*pixel_dimensions[0]) and (y0 < -1*pixel_dimensions[1]): #electron moves negatively in x
                            # and y
                            translation_x = round_half_down(x=x0 / (2 * pixel_dimensions[0]))
                            translation_y = round_half_down(x=y0 / (2 * pixel_dimensions[1]))
                            if i_coordinate - translation_x >= 0 and j_coordinate - translation_y >= 0:
                                    eh_charge_counter[i_coordinate - translation_x, j_coordinate - translation_y] =new_eh_pairs + eh_charge_counter[i_coordinate -translation_x,j_coordinate -translation_y]



                        elif x0 > (1*pixel_dimensions[0]) and y0 < (-1*pixel_dimensions[1]):  # electron moves psoitively in
                            # x and negatively in y
                            translation_x = round_half_down(x=x0 / (2 * pixel_dimensions[0]))
                            translation_y = round_half_down(x=y0 / (2 * pixel_dimensions[1]))
                            if j_coordinate - translation_y >= 0 and i_coordinate + translation_x  <= \
                                    perfect_image.shape[0] - 1:
                                    eh_charge_counter[i_coordinate + translation_x, j_coordinate - translation_y] = new_eh_pairs + eh_charge_counter[i_coordinate +translation_x,j_coordinate -translation_y]

                        elif x0 < (-1*pixel_dimensions[0]) and y0 > (1*pixel_dimensions[1]):  # electron moves negatively in
                            # x and positively in y
                            translation_x = round_half_down(x=x0 / (2 * pixel_dimensions[0]))
                            translation_y = round_half_down(x=y0 / (2 * pixel_dimensions[1]))
                            if i_coordinate - translation_x >= 0 and j_coordinate + \
                                    translation_y <= perfect_image.shape[0] - 1:
                                    eh_charge_counter[i_coordinate - translation_x, j_coordinate + translation_y] = new_eh_pairs + eh_charge_counter[i_coordinate -translation_x,j_coordinate +translation_y]
#############################################################################################################
                        elif x0 > (1*pixel_dimensions[0]): #electron moves in positive x direction
                            translation_x = round_half_down(x=x0 / (2 * pixel_dimensions[0]))
                            if i_coordinate + \
                                    translation_x <= perfect_image.shape[0] - 1:

                                    eh_charge_counter[i_coordinate + translation_x, j_coordinate] = new_eh_pairs + eh_charge_counter[i_coordinate +translation_x, j_coordinate]

                        elif x0 < (-1*pixel_dimensions[0]): #electron moves in negative x direction
                            translation_x = round_half_down(x=x0 / (2 * pixel_dimensions[0]))
                            if i_coordinate - translation_x >= 0:

                                    eh_charge_counter[i_coordinate - translation_x, j_coordinate] = new_eh_pairs + eh_charge_counter[i_coordinate -translation_x,j_coordinate]

                        elif y0 > (1*pixel_dimensions[1]): #electron moves in positive y direction
                            translation_y = round_half_down(x=y0 / (2 * pixel_dimensions[1]))
                            if j_coordinate + \
                                    translation_y <= perfect_image.shape[0] - 1:

                                    eh_charge_counter[i_coordinate, j_coordinate + translation_y] = new_eh_pairs + eh_charge_counter[i_coordinate,j_coordinate +translation_y]

                        elif y0 < (-1*pixel_dimensions[0]): #electron moves in negative y direction
                            translation_y = round_half_down(x=y0 / (2 * pixel_dimensions[1]))

                            if j_coordinate - translation_y >= 0:

                                    eh_charge_counter[i_coordinate, j_coordinate - translation_y] = new_eh_pairs + eh_charge_counter[i_coordinate,j_coordinate -translation_y]

                eh_charge_counter = (np.floor(dE_threshold*eh_charge_counter/t_counting))
                new_image_MCS += eh_charge_counter #+ new_image_MCS
    return new_image_MCS


#image
if path_to_image.endswith('.dat'): #opens Dr Probe image
    perfect_image = opendatfile(image = path_to_image)
elif path_to_image.endswith('.tif'):
    from PIL import Image
    perfect_image = Image.open(path_to_image)
    perfect_image = np.array(perfect_image)
else: #opens most standard files
    perfect_image = np.load(path_to_image)#io.imread(path_to_image)##
#perfect_image = np.ones((256,256))
chosen_pixels = find_distribution(perfect_image, dose)

print(chosen_pixels)
perfect_image = num_samples*distribute_electrons(perfect_image, chosen_pixels)
#output image
with open(path_to_unmodulated_image, 'wb') as f:#############
    np.save(f, perfect_image/num_samples)


pixel_information = get_pixel_info(image = perfect_image)



pixel_information, _ = make_2D_array(pixel_information)

print('\nInFluence is running...')



start = time.time()
modulated_image = np.uint32(MCS_ChargeCounting(E_i, minimum_energy, t_counting, dE_threshold, N, Z, A, p, pixel_dimensions,pixel_information, perfect_image))
print('\nInFluence has generated the modulated image, saved to: ', path_to_modulated_image, '\n\nRun time: ', time.time() - start)

#output image
with open(path_to_modulated_image, 'wb') as f:
    np.save(f, modulated_image/num_samples)
 

plt.imshow(modulated_image/num_samples)  # 'gray' colormap for grayscale images, omit for RGB
plt.axis('off')  # Turn off axes and ticks for better visualization (optional)
plt.show()


output_file = 'legacy.png'
plt.savefig(output_file)


non_zero_values = perfect_image [perfect_image  != 0]

# Print the non-zero values
print(non_zero_values)

# Add code to perform other tasks as needed (e.g., plot trajectories, save modulated image)


print(E_i, minimum_energy, t_counting, dE_threshold, N, Z, A, p, pixel_dimensions)




print(f"E_i: {E_i}")
print(f"dE_threshold: {dE_threshold}")
print(f"MinimumEnergy: {minimum_energy}")
print(f"ProtonNum: {Z}")
print(f"AtomicMass: {A}")
print(f"Density: {p}")
print(f"t_counting: {t_counting}")
print(f"N: {N}")