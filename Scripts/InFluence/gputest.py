

import numpy as np

first_column = np.random.choice(4, 9, replace=True)
other_columns = np.random.randint(0, 10, size=(9, 8))

pixels = np.column_stack((first_column, other_columns))

columns = 3
rows = 3
charge_arr = np.ones((rows, columns))


def gpu_t(pixels, charge_arr):

    for pixel in pixels:
        electron_count = pixel[0]

        for i in range(electron_count):
            
            print(charge_arr*electron_count)
            


gpu_t(pixels, charge_arr)