import math
import random
import numpy as np
from numba import jit, prange, cuda
from numba.typed import List
from skimage import io


##IGNORE FOR NOW - KEEP AS HELPER FUNCTION
def opendatfile(image):
    datfile = np.fromfile(image,
                          dtype='float32')
    image_shape = int(math.sqrt(datfile.shape[0]))  # obtains the square dimensions for the image from the one dimensional input file
    image_shift = int(image_shape / 2)
    datfile = datfile.reshape(image_shape, image_shape)
    datfile = np.roll(datfile, shift=0, axis=1)
    datfile = np.roll(datfile, shift=0, axis=0)
    return datfile


#get information about the pixel: counts, i, j coords
#@jit(nopython = True)               
def get_pixel_info(image):
    pixel_info = []
    rows, cols = image.shape
    for i in prange(rows):              ### this loop can be removd
        for j in prange(cols):
            counts = image[i, j]
            pixel_info.append(np.array([counts, i, j]))
    return pixel_info

#A modified version of finding pixel information to scale the input image based on initial intensity.
#@jit(nopython = True)
def get_pixel_info_mod(image):
    pixel_info = []
    count_distribution = []
    rows, cols = image.shape
    for i in range(rows):
        for j in range(cols):
            counts = image[i, j]                #this loop can be removed
            pixel_info.append([i, j])
            count_distribution.append(counts)

        
    return pixel_info, count_distribution

#function to distribute a budget of electrons based on probability of a pixel being bright.
def find_distribution(image, dose):
    image = image.astype('float64')
    image = image/np.sum(image)
    p_coordinates, count_distribution = get_pixel_info_mod(image)
    p_coordinates1 = np.empty(len(p_coordinates), dtype=object) 
    
    p_coordinates1[:] = p_coordinates

    return np.random.choice(p_coordinates1, dose, p = count_distribution)


#@jit(nopython = True)
def distribute_electrons(image, chosen_pixels):

    scaled_perfect_image = np.zeros((image.shape[0], image.shape[1]))
    for coordinates in chosen_pixels:
        i, j = coordinates[0], coordinates[1]
        scaled_perfect_image[i, j] = scaled_perfect_image[i, j] + 1


    return scaled_perfect_image




def make_2D_array(lis):                                         #incorperate into method
    """Funciton to get 2D array from a list of lists
    """
    n = len(lis)
    lengths = np.array([len(x) for x in lis])
    max_len = np.max(lengths)
    arr = np.zeros((n, max_len))

    for i in range(n):
        arr[i, :lengths[i]] = lis[i]
    return arr, lengths


#@jit(nopython = True)
def round_half_down(x, decimals = 0):

    multiplier = 10**decimals
    return int(math.ceil(abs(x)*multiplier - 0.5) / multiplier)
#binning
def rebin(arr, new_shape):
    shape = (new_shape[0], arr.shape[0] // new_shape[0],
             new_shape[1], arr.shape[1] // new_shape[1])
    return arr.reshape(shape).mean(-1).mean(1)





#The next set of functions are for Monte Carlo scattering
#generate a random variable


#@jit(nopython = True)
def RND(a, b):                              #no point in being own function
    return random.uniform(a, b)
#evaluate cross section and mean free path
#@jit(nopython = True)
def evaluate_alpha(E, Z):
    return (3.4*10**-3)*(Z**0.67)/E



#@jit(nopython = True)
def evaluate_cross_section(E, Z, alpha):
    # return (((E+511)/(E+1024))**2)*(5.21*10**-21)*((Z/E)**2)*(4*math.pi)/(alpha*(1+alpha))
    u = math.log10(8 * E * Z ** -1.33)
    return (4.7 * 10 ** -18) * (Z ** 1.33 + 0.032 * Z ** 2) / (
                (E + (E ** 0.5) * 0.0155 * Z ** 1.33) * (1 - math.exp(-1*u**2) * 0.02 * Z ** 0.5))


#@jit(nopython = True)
def evaluate_path_length(A, N, p, cross):
    return A/(N*p*cross)




#evaluate step and angula changes
#@jit(nopython = True)
def evaluate_step(path_length, RND):
    return -path_length*np.log(RND)


#@jit(nopython = True)
def evaluate_phi(RND, alpha):
    cosphi = 1 - (2*alpha*RND**2)/(1+alpha-RND)
    return math.acos(cosphi)


#@jit(nopython = True)
def evaluate_pho(RND):
    return 2*math.pi*RND

##evaluate direction cosines
##evaluate direction cosines



#@jit(nopython = True)
def evaluate_direction_cosine_a(phi, psi, cx, cy, cz):
    alpha = math.sin(psi)*math.sin(phi)
    beta = math.sin(phi)*math.cos(psi)
    gamma = math.cos(phi)
    cos_1 = cz
    sin_1 = math.sqrt(1 -cz**2)
    cos_2 = cy/sin_1
    sin_2 = cx/sin_1

    return alpha*cos_2 + sin_2*(beta*cos_1  + gamma*sin_1)

#@jit(nopython = True)
def evaluate_direction_cosine_b(phi, psi, cx, cy, cz):
    alpha = math.sin(phi) * math.sin(psi)
    beta = math.sin(phi) * math.cos(psi)
    gamma = math.cos(phi)
    cos_1 = cz
    sin_1 = math.sqrt(1 - cz ** 2)
    cos_2 = cy / sin_1
    sin_2 = cx / sin_1
    return -alpha*sin_2 + cos_2*(beta*cos_1 + gamma*sin_1)

#@jit(nopython = True)
def evaluate_direction_cosine_c(phi, psi, cz):

    beta = math.sin(phi) * math.cos(psi)
    gamma = math.cos(phi)
    cos_1 = cz
    sin_1 = math.sqrt(1 - cz ** 2)
    return -beta*sin_1 + gamma*cos_1

#calculate the energy loss rate.
##@jit(nopython = True)
def evaluate_energy_loss_rate(E, Z, A):
    J = (9.76*Z + 58.5/Z**0.19)*10**-3
    return -78500*Z/(A*E)*np.log(1.166*(E+0.85*J)/J)

#evaluate initial positions
#@jit(nopython = True)
def initialise_postions(step, d):
    z0 = 10**-2 + step #10**-2 is to ensure all electrons are incident approximately normally on the surface

    x0 = d*RND(a = -1, b = 1)
    y0 = d*RND(a = -1, b = 1)
    vector_length = math.sqrt(x0**2 + y0**2 + z0**2)
    cx = x0/vector_length
    cy = y0/vector_length
    cz = z0/vector_length
    return np.array([cx, cy, cz, z0, y0, x0, vector_length])
