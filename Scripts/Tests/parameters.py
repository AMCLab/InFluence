import os
from pathlib import Path
import configparser


class SimulationParameters:

    # Define class variables with default values
    E_i = None
    minimum_energy = None
    t_counting = None
    N = None
    ProtonNum = None
    AtomicMass = None
    Density = None
    dE_threshold = None
    pixel_dimensions = None
    dose = None
    
    def __init__(self):
        self.config_path = self.get_config_path()
        self.config = self.load_config()
        self.refresh_config()
        

    def get_config_path(self):
        # Get the directory where this Python script resides
        current_dir = Path(os.path.dirname(os.path.abspath(__file__)))

        # Define the location of the config.ini relative to the script's location
        config_path = current_dir / "config.ini"
        return config_path

    def load_config(self):
        config = configparser.ConfigParser()
        config.read(self.config_path)
        return config

    def refresh_config(self):
        # Reload the configuration from the file
        self.config = self.load_config()

        # Update class variables using the get_ functions
        self.E_i = self.get_E_i()
        self.minimum_energy = self.get_minimum_energy()
        self.t_counting = self.get_t_counting()
        self.N = self.get_N()
        self.ProtonNum = self.get_ProtonNum()
        self.AtomicMass = self.get_AtomicMass()
        self.Density = self.get_Density()
        self.dE_threshold = self.get_dE_threshold()
        self.pixel_dimensions = self.get_pixel_dimensions()
        self.dose = self.get_electron_dose()

    def get_E_i(self):
        return self.config.getint('Energies', 'e_i')

    def get_minimum_energy(self):
        return self.config.getfloat('Energies', 'minimum_energy')

    def get_t_counting(self):
        return self.config.getint('Energies', 't_counting')

    def get_N(self):
        return self.config.getfloat('Material properties', 'n')

    def get_ProtonNum(self):
        return self.config.getint('Material properties', 'ProtonNum')

    def get_AtomicMass(self):
        return self.config.getfloat('Material properties', 'AtomicMass')

    def get_Density(self):
        return self.config.getfloat('Material properties', 'Density')

    def get_dE_threshold(self):
        return self.config.getfloat('Material properties', 'dE_threshold')

    def get_pixel_dimensions(self):
        pixel_dimensions_str = self.config.get('Material properties', 'pixel_dimensions')
        pixel_dimensions_list = [float(val) for val in pixel_dimensions_str.split(',')]
        return pixel_dimensions_list


    def get_electron_dose(self):
        return self.config.getint('Electrons', 'dose')
