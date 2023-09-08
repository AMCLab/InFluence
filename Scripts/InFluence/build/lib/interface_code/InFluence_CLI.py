import sys
import argparse
import configparser
from inquirer import prompt
import inquirer

sys.path.append('..')  # Adds the parent directory to the sys.path
from main import RunSimulation
sys.path.remove('..')
from gui import InFluenceGUI


class SimulationApp:

    def __init__(self, config_file='../config.ini'):
        self.config_file = config_file
        self.config = configparser.ConfigParser()
        self.config.read(self.config_file)
        self.parser = self.create_parser()


    def update_config(self,config_file, answers):
        config = configparser.ConfigParser()
        config.read(config_file)

        config['Energies']['e_i'] = str(answers['e_i'])
        config['Energies']['minimum_energy'] = str(answers['minimum_energy'])
        config['Energies']['t_counting'] = str(answers['t_counting'])
        config['Material properties']['n'] = str(answers['n'])
        config['Material properties']['ProtonNum'] = str(answers['ProtonNum'])
        config['Material properties']['AtomicMass'] = str(answers['AtomicMass'])
        config['Material properties']['Density'] = str(answers['Density'])
        config['Material properties']['de_threshold'] = str(answers['de_threshold'])
        config['Material properties']['pixel_dimensions'] = ", ".join(map(str.strip, answers['pixel_dimensions'].split(',')))
        config['Electrons']['dose'] = str(answers['dose'])
        config['Filepaths']['img_in'] = str(answers['InputImageLocation'])
        config['Filepaths']['img_out'] = str(answers['SaveImageLocation'])

        with open(config_file, 'w') as configfile:
            config.write(configfile)


    def create_parser(self):


        parser = argparse.ArgumentParser(description='InFluence, a package to simulate electron-detector interactions in TEM imaging.')

        parser.add_argument('--gui', action='store_true', help='Launch InFluence GUI.')

        parser.add_argument('--gpu', action='store_true', help='Run the simulation on the GPU using existing parameters in config file.')
        parser.add_argument('--pysparse', action='store_true', help='Run the simulation on the CPU with PySparse using existing parameters in config file.')
        parser.add_argument('--pydense', action='store_true', help='Run the simulation on the CPU with PyDense using existing parameters in config file.')
        parser.add_argument('--cpp', action='store_true', help='Run the simulation on the CPU with CPP using existing parameters in config file')
        parser.add_argument('--plot', action='store_true', help='Run the simulation on the CPU and plot electron trajectories using existing parameters in config file, the dose must not exceed 5,000 when using this option.') #add input check 
        parser.add_argument('--go', action='store_true', help='Use the existing parameters from config.ini and run simulation of default simulation type (PySparse)')
        parser.add_argument('--show_params', action='store_true', help='Display stored simulation parameters')      
        
        parser.add_argument('--update_dose', type=int, help='Set the electron dose value')
        parser.add_argument('--e_i', type=int, help='Set the e_i value')
        parser.add_argument('--minimum_energy', type=float, help='Set the minimum_energy')
        parser.add_argument('--t_counting', type=int, help='Set the t_counting')
        parser.add_argument('--n', type=float, help='Set the n value')
        parser.add_argument('--proton_num', type=int, help='Set the proton number of the detector material')
        parser.add_argument('--atomic_mass', type=float, help='Set the atomic mass of the detector material')
        parser.add_argument('--Density', type=float, help='Set the density of the detector material')
        parser.add_argument('--de_threshold', type=float, help='Set the detector threshold')
        parser.add_argument('--pixel_dimensions', type=str, help='Set the pixel dimensions')
        parser.add_argument('--InputImageLocation', type=str, help='Set the input image file path')
        parser.add_argument('--SaveImageLocation', type=str, help='Set the saved image file path')

        
        return parser


    def display_stored_params(self, config_file):
        config = configparser.ConfigParser()
        config.read(config_file)

        print("*" * 84)
        print("*" + " " * 34 + "\U0001F31F" + " Stored Parameters " + "\U0001F31F" + " " * 34 )
        for section in config.sections():
            print(f"* [{section}]")
            for option in config[section]:
                value = config.get(section, option)
                print(f"* {option} = {value}")
        print("*" + "-" * 82 )

    def print_output(self, simulation_type, dose):
        print("*" * 90)
        print("*" + " " * 34 + "\U0001F31F" + " InFluence " + "\U0001F31F" + " " * 34 )

        print()
        print("  InFluence simulation will run using ", simulation_type, " for ", dose, " electrons...")
        print()




    def run_app(self):

        config_file = '../config.ini'

        parser = self.create_parser()
        args = parser.parse_args()
        config = configparser.ConfigParser()
        config.read(config_file)

        if any(arg is not None for arg in [args.update_dose, args.e_i, args.minimum_energy, args.t_counting, args.n, args.proton_num, args.atomic_mass, args.Density, args.de_threshold, args.pixel_dimensions, args.InputImageLocation, args.SaveImageLocation]):
            if args.update_dose is not None:
                config['Electrons']['dose'] = str(args.update_dose)
            if args.e_i is not None:
                config['Energies']['e_i'] = str(args.e_i)
            if args.minimum_energy is not None:
                config['Energies']['minimum_energy'] = str(args.minimum_energy)
            if args.t_counting is not None:
                config['Energies']['t_counting'] = str(args.t_counting)
            if args.n is not None:
                config['Material properties']['n'] = str(args.n)
            if args.proton_num is not None:
                config['Material properties']['ProtonNum'] = str(args.ProtonNum)
            if args.atomic_mass is not None:
                config['Material properties']['AtomicMass'] = str(args.AtomicMass)
            if args.Density is not None:
                config['Material properties']['Density'] = str(args.Density)
            if args.de_threshold is not None:
                config['Material properties']['de_threshold'] = str(args.de_threshold)
            if args.pixel_dimensions is not None:
                config['Material properties']['pixel_dimensions'] = args.pixel_dimensions
            if args.InputImageLocation is not None:
                config['Filepaths']['img_in'] = args.InputImageLocation
            if args.SaveImageLocation is not None:
                config['Filepaths']['img_out'] = args.SaveImageLocation

            # Save updated config
            with open(config_file, 'w') as configfile:
                config.write(configfile)
            print("Simulation configuration updated.")

        if all(arg is None for arg in [args.update_dose, args.e_i, args.minimum_energy, args.t_counting, args.n, args.proton_num, args.atomic_mass, args.Density, args.de_threshold, args.pixel_dimensions, args.InputImageLocation, args.SaveImageLocation]):




            if len(sys.argv) >= 3:
                    print("Error: Too many arguments passed.")
                    sys.exit(1)
            
            if args.gui:
                gui = InFluenceGUI()
                gui.build_gui()
                

            else:

                if args.show_params:
                    print('woo')
                    self.display_stored_params(config_file)
                    return

                elif args.update_dose is not None:
                    # Update dose value in the configuration file
                    config = configparser.ConfigParser()
                    config.read(config_file)
                    config['Electrons']['dose'] = str(args.update_dose)

                    with open(config_file, 'w') as configfile:
                        config.write(configfile)
                    print(f"New stored dose value updated to: {args.update_dose}")

                elif len(sys.argv)==1:
                    # Read the existing values from the config file
                    config = configparser.ConfigParser()
                    config.read(config_file)
                    stored_ei = config.getint('Energies', 'e_i')
                    stored_min_energy = config.getfloat('Energies', 'minimum_energy')
                    stored_t_counting = config.getint('Energies', 't_counting')
                    stored_n = config.getfloat('Material properties', 'n')
                    stored_proton_num = config.getint('Material properties', 'ProtonNum')
                    stored_atomic_mass = config.getfloat('Material properties', 'AtomicMass')
                    stored_density = config.getfloat('Material properties', 'Density')
                    stored_de_threshold = config.getfloat('Material properties', 'de_threshold')
                    stored_pixel_dimensions = config.get('Material properties', 'pixel_dimensions')
                    stored_dose = config.getint('Electrons', 'dose')
                    stored_img_in = config.get('Filepaths', 'img_in')
                    stored_img_out = config.get('Filepaths', 'img_out')

                    # Create a dictionary to store the user's input
                    answers = {}

                    # Existing code for prompting user for input and updating the configuration
                    questions = [
                        inquirer.Text('e_i', message='Enter value for e_i:', default=stored_ei),
                        inquirer.Text('minimum_energy', message='Enter value for minimum_energy:', default=stored_min_energy),
                        inquirer.Text('t_counting', message='Enter value for t_counting:', default=stored_t_counting),
                        inquirer.Text('n', message='Enter value for n:', default=stored_n),
                        inquirer.Text('ProtonNum', message='Enter value for ProtonNum:', default=stored_proton_num),
                        inquirer.Text('AtomicMass', message='Enter value for AtomicMass:', default=stored_atomic_mass),
                        inquirer.Text('Density', message='Enter value for Density:', default=stored_density),
                        inquirer.Text('de_threshold', message='Enter value for de_threshold:', default=stored_de_threshold),
                        inquirer.Text('pixel_dimensions', message='Enter pixel_dimensions (comma-separated x,y,z):', default=stored_pixel_dimensions),
                        inquirer.Text('dose', message='Enter value for electron dose:', default=stored_dose),
                        inquirer.Text('InputImageLocation', message='Enter filepath for input image.', default=stored_img_in),
                        inquirer.Text('SaveImageLocation', message='Enter filepath for saved image.', default=stored_img_out),
                        inquirer.Text('simulation_type', message='Enter simulation type, GPU, CPP, PySparse, PyDense.', default='PySparse')  
                    ]

                    answers = inquirer.prompt(questions)
                    self.update_config(config_file, answers)
                    simulation_type = answers['simulation_type']
                    self.print_output('PySparse simulation on CPU', stored_dose)
                    RunSimulation(simulation_type)

                if len(sys.argv)==2:  

                    if args.gpu:
                        self.print_output('GPU', stored_dose)
                        RunSimulation('GPU')
                        return

                    elif args.pysparse:
                        self.print_output('PySparse simulation on CPU', stored_dose)
                        RunSimulation('PySparse')
                        return

                    if args.pysdense:
                        self.print_output('PyDense simulation on CPU', stored_dose)
                        RunSimulation('PyDense')
                        return

                    if args.cpp:
                        self.print_output('CPP simulation on CPU', stored_dose)
                        RunSimulation('CPP')
                        return

                    if args.pysparse:
                        self.print_output('using CPU and plot electron tragectories.', stored_dose)
                        RunSimulation('Plot')
                        return

            
                if args.go:

                    if len(sys.argv) == 2:
                        self.print_output('PySparse simulation on CPU', stored_dose)
                        RunSimulation('CPU')
                        return
                
                print("                     Simulation complete. Output saved to ", stored_img_out)
                print()
                print()
                print("*" * 90)
                
            return 

                 

if __name__ == "__main__":
    app = SimulationApp()

    app.run_app()



