import configparser

class SimulationParameters:
    
    config = configparser.ConfigParser()
    config.read('config.ini')

    E_i = config.getint('Energies', 'e_i')
    minimum_energy = config.getfloat('Energies', 'minimum_energy')
    t_counting = config.getint('Energies', 't_counting')
    N = config.getfloat('Material properties', 'n')
    ProtonNum = config.getint('Material properties', 'ProtonNum')
    AtomicMass = config.getfloat('Material properties', 'AtomicMass')
    Density = config.getfloat('Material properties', 'Density')
    dE_threshold = config.getfloat('Material properties', 'dE_threshold')
    pixel_dimensions = list(map(float, config.get('Material properties', 'pixel_dimensions').split(',')))
    dose = config.getint('Electrons', 'dose')

