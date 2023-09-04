from numba import int32, float64
from numba.experimental import jitclass




# Define the specification for the pixel_constructor class
pixel_spec = [
    ('electron_count', int32),
    ('i', int32),
    ('j', int32),
    ('x_dimension', float64),
    ('y_dimension', float64),
    ('z_dimension', float64)
]

@jitclass(pixel_spec)
class pixel_constructor:
    
    def __init__(self, electron_count, i, j, x_dimension, y_dimension, z_dimension):
        self.electron_count = int(electron_count)
        self.i = int(i)
        self.j = int(j)
        self.x_dimension = x_dimension
        self.y_dimension = y_dimension
        self.z_dimension = z_dimension