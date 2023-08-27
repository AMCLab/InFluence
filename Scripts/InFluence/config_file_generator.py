import configparser
import scipy.constants as scc


# Create a ConfigParser object
config = configparser.ConfigParser()

# Set the values for the configuration options
config['Energies'] = {
    'E_i': '200',  # in keV
    'minimum_energy': '0.05',  # in keV
    't_counting': '40'  # in keV
}

config['Material properties'] = {
    'N': str(scc.physical_constants['Avogadro constant'][0]),
    'Z': '14',  # number of protons
    'A': '28.1',  # atomic weight in Daltons
    'p': '2.33',  # density of material in g/cm cubed
    'dE_threshold': '0.00362',  # eh pair generation of material in keV
    'pixel_dimensions': '[0.0275, 0.0275, 0.04]',  # [x, y, 10^-2 + z] in cm - change x, y, z only #correction required here
}

config['Choose electron dose'] = {
    'dose': '50'  # choose dose in number of electrons
}

# Write the configuration to a file
with open('config.ini', 'w') as configfile:
    config.write(configfile)
