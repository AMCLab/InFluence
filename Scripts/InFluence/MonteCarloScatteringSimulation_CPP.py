from CPP import InFluence
from parameters import params
import numpy as np

def CPP_MCScatteringSimulation(pixels,
    E_i, 
    ProbeDiameter, 
    MinimumEnergy, 
    pixel_dimensions, 
    dE_threshold, 
    perfect_image_0, 
    perfect_image_1, 
    Density,
    t_counting
    ):

    print('Simulation starting ...')


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
    
    pixels_list = pixels_array.tolist()

    print('Pixel data prepared ...')

    
    print(len(pixels_list))
    print('Beginning main simulation loop...')

    import time
    
    start_time = time.time()
    

    FlattenedImage = InFluence.RunMCScatteringSimulation(pixels_list, E_i, ProbeDiameter, MinimumEnergy, dE_threshold,
    perfect_image_0, perfect_image_1, Density, t_counting)

    end_time = time.time()

    # Calculate the elapsed time
    duration = end_time - start_time
    print("Execution time:", duration, "seconds")
    
    OutputImage = np.array(FlattenedImage)
    
   
    print('Simulation complete.')

    return OutputImage