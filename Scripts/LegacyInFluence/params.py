import scipy.constants as scc
from common_functions import *
import pathlib
import os

#Energies
E_i = 200 #in keV
minimum_energy = 0.05 #in keV
t_counting = 40 #in keV

#Material properties
N = scc.physical_constants['Avogadro constant'][0]
Z = 14 #number of protons
A = 28.1 #atomic weight in Daltons
p = 2.33 #density of material in g/cm cubed
dE_threshold = 0.00362 #eh pair generation of material in keV
pixel_dimensions = (np.around([27.5*10**-4, 27.5*10**-4, 10**-2 + 300*10**-4], decimals=6)) # [x, y, 10**-2 + z] in cm -change, x, y, z only

#Choose electron dose
dose = 1000000 #choose dose in number of electrons

#Output file type
ext = '.png'
ext_2 = '.svg'

#images
num_samples = 1 #relative precision = 1/root(n)

#unmodulated image
path_to_image = r"/workspaces/InFluence/ExampleData/slanted_edge.npy"

image_type = 'specimen' #knife_edge, white_noise or specimen

#No need to change these
save_file_folder = os.path.join(str(pathlib.Path(path_to_image).parent.parent), '{}_keV_HT_and_{}_keV_countingthreshold'.format(E_i, t_counting),str(dose) + '_electrons')
os.makedirs(save_file_folder, exist_ok=True)



path_to_unmodulated_image = os.path.join(save_file_folder, 'Unmodulated_{}_{}_electrons.npy'.format(image_type,dose))

path_to_modulated_image = os.path.join(save_file_folder, 'Modulated_{}_{}_electrons.npy'.format(image_type,dose))

path_to_parameter_file = os.path.join(save_file_folder, 'Modulated_{}_parameter_file_for_{}_electrons.txt'.format(image_type,dose))

path_to_trajectory_plot = os.path.join(save_file_folder,'Trajectory_plot_view_1' + ext)
path_to_trajectory_plot_svg = os.path.join(save_file_folder,'Trajectory_plot_view_1' + ext_2)
path_to_trajectory_plot_v2 = os.path.join(save_file_folder, 'Trajectory_plot_view_2' + ext)
path_to_trajectory_plot_v2_svg = os.path.join(save_file_folder,'Trajectory_plot_view_2' + ext_2)
path_to_trajectory_histogram_plot = os.path.join(save_file_folder, 'Trajectory_histogram_plot' + ext)
path_to_trajectory_histogram_plot_svg = os.path.join(save_file_folder,'Trajectory_histogram_plot' + ext_2)

#Do not change these
print('\nParameters used in simulation: \n \n', ' Incident energy (keV): ', E_i, '\n  Mininum energy (keV): ', minimum_energy, '\n  Counting threshold (keV): ', t_counting,
      '\n  Electron-hole pair energy threshold (keV): ', dE_threshold, '\n  Z: ', Z, '\n  Atomic weight (u): ', A, '\n  Density (g/cm\N{SUPERSCRIPT THREE}): ', p,
      '\n  Pixel dimensions (cm): ', pixel_dimensions, '\n  Fluence (e): ', dose, '\n  Number of samples: ', num_samples)

with open(path_to_parameter_file, "w+") as text_file:
    print('Parameters used in simulation: \n \n', ' Incident energy (keV): ', E_i, '\n  Mininum energy (keV): ', minimum_energy, '\n  Counting threshold (keV): ', t_counting,
      '\n  Electron-hole pair energy threshold (keV): ', dE_threshold, '\n  Z: ', Z, '\n  Atomic weight (u): ', A, '\n  Density (g/cm\N{SUPERSCRIPT THREE}): ', p,
      '\n  Pixel dimensions (cm): ', pixel_dimensions, '\n  Fluence (e): ', dose, '\n  Number of samples: ', num_samples, file=text_file)
print('\nParameters used in this experiment are saved to: ', path_to_parameter_file)
