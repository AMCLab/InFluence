
#This class finds and stores the pixel information from the input image 

#Approach is to get itto 

from PIL import Image
import numpy as np
from image_preprocessing.pixel_object_constructor import pixel_constructor



class image_pixels:
    
    

    def __init__(self, filepath, params):

        self.original_image = None
        self.perfect_image = None
        self.electron_dose = params.dose
        self.pixels = None                                
        self.pixel_dimensions = params.pixel_dimensions

        self.params = params


        if filepath.endswith('.dat'):
            self.original_image = opendatfile(image=filepath)
        elif filepath.endswith('.tif'):
            self.original_image = np.array(Image.open(filepath))
        else:
            self.original_image = np.load(filepath)
        

    def create_pixel_objects(self, image =  None):
       
        if image is None:
            image = self.original_image

        rows, cols = image.shape
        indices = np.indices((rows, cols))
        counts = image.flatten()
        i_coordinates = indices[0].flatten()                        #redesigned to remove loops
        j_coordinates = indices[1].flatten()
        pixel_info = np.column_stack((counts, i_coordinates, j_coordinates))
        pixel_info = pixel_info.tolist()

        self.pixels = [
            pixel_constructor(
                electron_count,
                i,
                j,
                self.pixel_dimensions[0],
                self.pixel_dimensions[1],
                self.pixel_dimensions[2]
            ) for electron_count, i, j in pixel_info
        ]

        return self.pixels
    
    def find_distribution(self):
        
        electron_counts = self.original_image.flatten()
        probabilities = electron_counts / np.sum(electron_counts)
        selected_indices = np.random.choice(len(self.pixels), size=self.params.dose, p=probabilities)
        selected_pixels = [self.pixels[i] for i in selected_indices]
        return selected_pixels
    

    def distribute_electrons(self, selected_pixels):
        scaled_perfect_image = np.zeros(self.original_image.shape)     
        for pixel in selected_pixels:                               # can't think how to avoid this loop
            i, j = pixel.i, pixel.j
            scaled_perfect_image[i, j] += 1
        return scaled_perfect_image                                


    def generate_scaled_perfect_image(self):                    #design choice: kept seperate so each function has 1 concern
        
        self.create_pixel_objects()                                         #no image variable passed to this will take the original image to make pixels
        selected_pixels = self.find_distribution()                         
        self.perfect_image = self.distribute_electrons(selected_pixels)  
        self.pixels = None                                                  #for memory purposes  - ensure pixels stored only for perfect image to avoid confusion         
    
    def generate_perfect_image_pixel_objects(self):
        self.pixels = self.create_pixel_objects(self.perfect_image)      # Replace the pixels variable with the scaled "perfect image" pixels


    def save_image(self):
        pass


    


