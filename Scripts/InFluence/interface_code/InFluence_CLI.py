import argparse
import configparser
from inquirer import prompt
import inquirer
from main import RunSimulation

def update_config(config_file, answers):
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
    # Remove spaces from pixel_dimensions before updating the config
    config['Material properties']['pixel_dimensions'] = ", ".join(map(str.strip, answers['pixel_dimensions'].split(',')))
    config['Electrons']['dose'] = str(answers['dose'])
    config['Filepaths']['img_in'] = str(answers['InputImageLocation'])
    config['Filepaths']['img_out'] = str(answers['SaveImageLocation'])

    with open(config_file, 'w') as configfile:
        config.write(configfile)



def create_parser():
    parser = argparse.ArgumentParser(description='Your application description')
    parser.add_argument('--gpu', action='store_true', help='Run the simulation on the GPU')
    parser.add_argument('--stored_params', action='store_true', help='Use the existing parameters from config.ini')
    parser.add_argument('--input_img', nargs='?', type=str, help='File path of the input image')
    parser.add_argument('--save_img', type=str, help='Specify the file path to save the image')    
    parser.add_argument('--update_dose', type=int, help='Set the dose value')
    return parser


def display_stored_params(config_file):
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

def main():
    config_file = 'config.ini'
    parser = create_parser()
    args = parser.parse_args()
    config = configparser.ConfigParser()
    config.read(config_file)

    if args.stored_params:
        display_stored_params(config_file)
        return

    if args.update_dose is not None:
        # Update dose value in the configuration file
        config = configparser.ConfigParser()
        config.read(config_file)
        config['Electrons']['dose'] = str(args.update_dose)

        with open(config_file, 'w') as configfile:
            config.write(configfile)
        print(f"New stored dose value updated to: {args.update_dose}")

    elif not args.stored_params:
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
            inquirer.Text('SaveImageLocation', message='Enter filepath for saved image.', default=stored_img_out)
        ]

        answers = inquirer.prompt(questions)
        update_config(config_file, answers)
        
        # if argument gpu used then:

        if args.gpu:
            RunSimulation('GPU')
        


if __name__ == "__main__":
    main()






