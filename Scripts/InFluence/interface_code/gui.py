import tkinter as tk
from tkinter import Tk, Label
from tkinter import ttk
from tkinter import filedialog
from PIL import Image, ImageTk
from mendeleev import element
import platform
import psutil
from threading import Thread, Event
import GPUtil
import configparser
import sys
import matplotlib.pyplot as plt
import io
import importlib
import tkinter.messagebox as messagebox
import os
import tkinter as tk
from tkinter import Toplevel, Text, Scrollbar, Button

from PIL import Image, ImageDraw, ImageFont
from queue import Queue
import numpy as np

sys.path.append('..')
from main import RunSimulation
sys.path.remove('..')



class InFluenceGUI:

    def __init__(self):
        
        self.root = tk.Tk()
        self.root.title("InFluence")
        self.image_label = None  
        
        self.simulation_type = "PySparse" 
        
        self.stop_event = Event()
        self.simulation_thread = None
        self.queue = Queue()

        self.config_file = '../config.ini'
        self.config = configparser.ConfigParser()
        self.config.read(self.config_file)


    def start_simulation(self):
        print("Starting simulation...")
        self.create_simulating_display_image()
        if self.simulation_thread and self.simulation_thread.is_alive():
            print("A simulation is already running.")
            return
        self.stop_event.clear()

        self.simulation_thread = Thread(target=self.run_simulation)
        self.simulation_thread.start()

    def run_simulation(self):
        try:
            if not self.stop_event.is_set():

            
                self.simulation_result = RunSimulation(self.simulation_type, gui=True)
                self.display_simulation_results()
                self.queue.put("Done")
        except KeyboardInterrupt:
            print("Simulation stopped.")
        finally:
            self.reset_simulation()

            

    def stop_simulation(self):
        self.stop_event.set()
        print("Stopping simulation...")
        if self.simulation_thread:
            self.simulation_thread.join()  # Wait for the simulation thread to finish

    def reset_simulation(self):
        print("Resetting simulation...")
        self.stop_event.set()
        print('Simulation reset.')
        
        # Here, check for messages in the queue and update the GUI as needed

        
    def plot(self):
        self.simulation_type = 'Plot'
        self.start_simulation()
        messagebox.showinfo("Notification", "Plots will be saved to current directory.")
        


    def create_image_with_text(self, text):
        width, height = 128, 128
        influence_image = Image.new('RGB', (width, height), 'black')
        draw = ImageDraw.Draw(influence_image)
        color = '#E0EEEE'
        font = ImageFont.load_default()
        text_width, text_height = draw.textsize(text, font)
        text_x = (width - text_width) // 2
        text_y = (height - text_height) // 2
        draw.text((text_x, text_y), text, fill=color, font=font)
        resized_image = influence_image.resize((700,700), Image.ANTIALIAS)
        return ImageTk.PhotoImage(image=resized_image)

    def create_influence_display_image(self):
        influence_image_tk = self.create_image_with_text("InFluence")
        self.image_label = Label(self.frame, image=influence_image_tk, bg='#E0EEEE', padx=40, pady=40)
        self.image_label.grid(row=0, column=6, rowspan=15)
        self.image_label.image = influence_image_tk

    def create_simulating_display_image(self):
        simulating_image_tk = self.create_image_with_text("Simulating...")
        self.image_label.config(image=simulating_image_tk)
        self.image_label.image = simulating_image_tk
        self.frame.update_idletasks()
        
    def numpy_to_image(self):
        # Normalize the array to scale between 0-255
        normalized_array = ((self.simulation_result - self.simulation_result.min()) * (1/(self.simulation_result.max() - self.simulation_result.min()) * 255)).astype('uint8')
        image = Image.fromarray(normalized_array, 'L')  # 'L' indicates grayscale
        resized_image = image.resize((700, 700), Image.ANTIALIAS)
        return ImageTk.PhotoImage(image=resized_image)


    def display_simulation_results(self):
        simulation_results_tk = self.numpy_to_image()
        self.image_label.config(image=simulation_results_tk)
        self.image_label.image = simulation_results_tk
        self.frame.update_idletasks()
        
    def show_user_tutorial(self):
        new_window = Toplevel(self.root)
        new_window.title("User Tutorial")

        text_widget = Text(new_window, wrap='word')
        text_widget.pack(expand=1, fill='both')

        # Add a Scrollbar
        scrollbar = Scrollbar(new_window, command=text_widget.yview)
        scrollbar.pack(side='right', fill='y')
        text_widget['yscrollcommand'] = scrollbar.set

        # Read text from file and insert into the Text widget
        with open('gui_user_guide.txt', 'r') as f:
            content = f.read()
            text_widget.insert('1.0', content)

        # Add Close button
        close_button = Button(new_window, text='Close', command=new_window.destroy)
        close_button.pack()


    
    def choose_input_file(self):
        file_path = filedialog.askopenfilename()
        self.file_path_entry_in.delete(0, tk.END)  # Remove current text in entry
        self.file_path_entry_in.insert(0, file_path)  # Insert the file path

        
    def choose_ouput_file(self):
        file_path = filedialog.askopenfilename()
        self.file_path_entry_out.delete(0, tk.END)  # Remove current text in entry
        self.file_path_entry_out.insert(0, file_path)  # Insert the file path

            
    def set_file_paths(self):
        input_path = self.file_path_entry_in.get().strip()
        output_path = self.file_path_entry_out.get().strip()

        if not input_path:
            input_path = os.getcwd()  # Sets to current directory if entry is empty
        if not output_path:
            output_path = os.getcwd()  # Sets to current directory if entry is empty
            
        self.config['Filepaths']['img_out'] = str(output_path)
        self.config['Filepaths']['img_in'] = str(input_path)
        
        with open(self.config_file, 'w') as configfile:
            self.config.write(configfile)
        print(f"Set input path to: {input_path}")
        print(f"Set output path to: {output_path}")


        
    def update_electron_dose_scale(self, *args):
        self.electron_dose_entry.delete(0, tk.END)
        value = self.electron_dose_scale.get()
        self.electron_dose_entry.insert(0, str(value))
        self.electron_dose = value  
        self.config['Electrons']['dose'] = str(self.electron_dose)
        with open(self.config_file, 'w') as configfile:
            self.config.write(configfile)

    def update_electron_dose_entry(self, *args):
        value = int(self.electron_dose_entry.get())
        self.electron_dose_scale.set(value)
        self.electron_dose = value
        self.config['Electrons']['dose'] = str(self.electron_dose)
        with open(self.config_file, 'w') as configfile:
            self.config.write(configfile)  
        

    def update_energy_scale(self, *args):
        self.electron_energy_entry.delete(0, tk.END)
        value = self.electron_energy_scale.get()
        self.electron_energy_entry.insert(0, str(value))
        self.e_i = value  # Assuming you want to update self.e_i
        self.config['Energies']['e_i'] = str(self.e_i)
        with open(self.config_file, 'w') as configfile:
            self.config.write(configfile)

    def update_energy_entry(self, *args):
        value = int(self.electron_energy_entry.get())
        self.electron_energy_scale.set(value)
        self.e_i = value  # Assuming you want to update self.e_i
        self.config['Energies']['e_i'] = str(self.e_i)
        with open(self.config_file, 'w') as configfile:
            self.config.write(configfile)
            
            
    def update_tecount_scale(self, *args):
        self.tecount_entry.delete(0, tk.END)
        value = self.tecount_scale.get()
        self.tecount_entry.insert(0, str(value))
        self.tecount = value  # Assuming you want to update self.e_i
        self.config['Energies']['t_counting'] = str(self.tecount)
        with open(self.config_file, 'w') as configfile:
            self.config.write(configfile)

    def update_tecount_entry(self, *args):
        value = int(self.tecount_entry.get())
        self.tecount_scale.set(value)
        self.tecount = value  # Assuming you want to update self.e_i
        self.config['Energies']['t_counting'] = str(self.tecount)
        with open(self.config_file, 'w') as configfile:
            self.config.write(configfile)


    def update_simulation_type(self, event):
        self.simulation_type = event.widget.get()
        print(f"Updated simulation type to {self.simulation_type}")
    
    def update_proton_number(self, *args):
        # Get the value from the proton number entry box
        proton_number_value = self.proton_num_entry.get()
        # Check if the value is valid (you can add more validation as needed)
        if proton_number_value.isdigit():
            # Update the proton number variable
            self.proton_number = int(proton_number_value)
            # Update the config
            self.config['Material properties']['protonnum'] = str(self.proton_number)
            # Save the updated config to the file
            with open(self.config_file, 'w') as configfile:
                self.config.write(configfile)
        else:
            # Handle invalid input (e.g., show an error message)
            print("Invalid proton number input")

    
        
    def update_detector_settings(self, detector_threshold_entry, proton_num_entry, atomic_mass_entry, material_density_entry,
                             pixel_dimension_x_entry, pixel_dimension_y_entry, pixel_dimension_z_entry):
        # Read values from GUI widgets and convert them to the appropriate types
        detector_threshold = detector_threshold_entry.get()
        atomic_mass = atomic_mass_entry.get()
        proton_num = proton_num_entry.get()
        material_density = material_density_entry.get()
        pixel_dimension_x = pixel_dimension_x_entry.get()
        pixel_dimension_y = pixel_dimension_y_entry.get()
        pixel_dimension_z = pixel_dimension_z_entry.get()
        
         # Check if any of the required values are empty
        if any(value == '' for value in [detector_threshold, atomic_mass, proton_num, material_density,
                                     pixel_dimension_x, pixel_dimension_y, pixel_dimension_z]):
            # Show an error message if any of the values are empty
            messagebox.showerror("Error", "All values must be entered.")
            return

        # Update the config object with the new values
        self.config['Energies']['e_i'] = str(detector_threshold)
        self.config['Material properties']['atomicmass'] = str(atomic_mass)
        self.config['Material properties']['protonnum'] = str(proton_num)
        self.config['Material properties']['density'] = str(material_density)
        self.config['Material properties']['pixel_dimensions'] = f"{pixel_dimension_x}, {pixel_dimension_y}, {pixel_dimension_z}"

        # Save the updated config to the file
        with open(self.config_file, 'w') as configfile:
            self.config.write(configfile)

        # Optionally, print a message to confirm the update
        print("Detector settings updated and saved.")



    def check_system_specs(self):
        specs_window = tk.Toplevel(self.root)
        specs_window.title("System Specifications")
    
        os_info = platform.uname()
        cpu_count = psutil.cpu_count(logical=False)
        logical_cpu_count = psutil.cpu_count(logical=True)
    
        GPUs = GPUtil.getGPUs()
        num_of_gpus = len(GPUs)
    
        specs_text = tk.Text(specs_window, height=20, width=80)
        specs_text.pack()
    
        specs_text.insert(tk.END, f"System: {os_info.system}\n")
        specs_text.insert(tk.END, f"Node Name: {os_info.node}\n")
        specs_text.insert(tk.END, f"Release: {os_info.release}\n")
        specs_text.insert(tk.END, f"Version: {os_info.version}\n")
        specs_text.insert(tk.END, f"Machine: {os_info.machine}\n")
        specs_text.insert(tk.END, f"Processor: {os_info.processor}\n")
        specs_text.insert(tk.END, f"Physical Cores: {cpu_count}\n")
        specs_text.insert(tk.END, f"Logical Cores: {logical_cpu_count}\n")
        
    

    def user_guide_window():
        pass

    
    def detector_settings_widnow(self):
    
        detector_window = tk.Toplevel(self.root, bg='#E0EEEE', padx=20, pady=20)
        detector_window.title('InFluence: Detector Settings')

        # Load existing values from the config
        detector_threshold_default = self.config.get('Energies', 'e_i', fallback=0.0)
        material_density_default = self.config.get('Material properties', 'density', fallback=0.0)
        pixel_dimensions_default = self.config.get('Material properties', 'pixel_dimensions', fallback="0.0, 0.0, 0.0")
        proton_num_default = self.config.get('Material properties', 'protonnum', fallback=0)
        atomic_mass_default = self.config.get('Material properties', 'atomicmass', fallback=0.0)
    
        # Detector Threshold
        tk.Label(detector_window, text="Detector Threshold", bg='#E0EEEE', padx=4, pady=4).grid(row=0, column=0)
        detector_threshold_entry = tk.Entry(detector_window)
        detector_threshold_entry.insert(0, detector_threshold_default)  # Set the default value
        detector_threshold_entry.grid(row=0, column=1)
        tk.Label(detector_window, text="(keV)", bg='#E0EEEE').grid(row=0, column=2)

        # Material Density
        tk.Label(detector_window, text="Material Density", bg='#E0EEEE', padx=4, pady=4).grid(row=3, column=0)
        material_density_entry = tk.Entry(detector_window)
        material_density_entry.insert(0, material_density_default)  # Set the default value
        material_density_entry.grid(row=3, column=1)
        tk.Label(detector_window, text="(m^-3)", bg='#E0EEEE').grid(row=3, column=2)

        # Pixel Dimensions
        tk.Label(detector_window, text="Pixel Dimensions:", bg='#E0EEEE', padx=6, pady=6).grid(row=4, column=0)
        tk.Label(detector_window, text="x:", bg='#E0EEEE').grid(row=5, column=0)
        pixel_dimension_x_entry = tk.Entry(detector_window, width=10)
        pixel_dimension_x_entry.insert(0, pixel_dimensions_default.split(", ")[0])  # Set the default value
        pixel_dimension_x_entry.grid(row=5, column=1)

        tk.Label(detector_window, text="y:", bg='#E0EEEE', padx=6, pady=6).grid(row=6, column=0)
        pixel_dimension_y_entry = tk.Entry(detector_window, width=10)
        pixel_dimension_y_entry.insert(0, pixel_dimensions_default.split(", ")[1])  # Set the default value
        pixel_dimension_y_entry.grid(row=6, column=1)

        tk.Label(detector_window, text="z:", bg='#E0EEEE', padx=6, pady=6).grid(row=7, column=0)
        pixel_dimension_z_entry = tk.Entry(detector_window, width=10)
        pixel_dimension_z_entry.insert(0, pixel_dimensions_default.split(", ")[2])  # Set the default value
        pixel_dimension_z_entry.grid(row=7, column=1)
    
        # Proton Number Entry
        tk.Label(detector_window, text="Proton Number:", bg='#E0EEEE', padx=4, pady=4).grid(row=1, column=0)
        proton_num_entry = tk.Entry(detector_window)
        proton_num_entry.insert(0, proton_num_default)  # Set the default value
        proton_num_entry.grid(row=1, column=1)

        # Atomic Mass Entry
        tk.Label(detector_window, text="Atomic Mass:", bg='#E0EEEE', padx=4, pady=4).grid(row=2, column=0)
        atomic_mass_entry = tk.Entry(detector_window)
        atomic_mass_entry.insert(0, atomic_mass_default)  # Set the default value
        atomic_mass_entry.grid(row=2, column=1)

        # Elements, their proton numbers, and atomic masses
        element_data_dict = {
            'Silicon': {'Proton Number': 14, 'Atomic Mass': 28.085},
            'Germanium': {'Proton Number': 32, 'Atomic Mass': 72.63},
            'Beryllium': {'Proton Number': 4, 'Atomic Mass': 9.0122},
            'Gold': {'Proton Number': 79, 'Atomic Mass': 196.966569},
            'Platinum': {'Proton Number': 78, 'Atomic Mass': 195.084},
            'Aluminium': {'Proton Number': 13, 'Atomic Mass': 26.9815386}
        }



        # Create a Tkinter variable
        chosen_element = tk.StringVar(detector_window)
        chosen_element.set('Select Element')  # default value

        # Define the update_element_data function
        def update_element_data(*args):
            selected_element = chosen_element.get()
            if selected_element in element_data_dict:
                proton_num_entry.delete(0, tk.END)
                proton_num_entry.insert(0, element_data_dict[selected_element]['Proton Number'])

                atomic_mass_entry.delete(0, tk.END)
                atomic_mass_entry.insert(0, element_data_dict[selected_element]['Atomic Mass'])

        chosen_element.trace("w", update_element_data)

        # Create the dropdown menu
        tk.Label(detector_window, text="Choose an element:", bg='#E0EEEE', padx=4, pady=4).grid(row=0, column=2)
        element_menu = ttk.Combobox(detector_window, textvariable=chosen_element, values=list(element_data_dict.keys()), width=10)
        element_menu.grid(row=1, column=3)

        # Update Proton Number and Atomic Mass based on chosen element

        # Save detector settings
        save_button = tk.Button(detector_window, text="Save Settings",  bg='#C1CDCD', padx=6, pady=6, command=lambda: self.update_detector_settings(
        detector_threshold_entry, proton_num_entry, atomic_mass_entry, material_density_entry,
        pixel_dimension_x_entry, pixel_dimension_y_entry, pixel_dimension_z_entry))
        save_button.grid(row=8, column=4, columnspan=5)


    def menu_bar(self):
        # Create a frame for the menu bar with a different background color
        menu_frame = tk.Frame(self.root, bg='#0052cc')
        menu_frame.pack(side='top', fill='x')

        # Create a menu bar
        menu_bar = tk.Menu(menu_frame, bg='#0052cc', fg='#ffffff')
        self.root.config(menu=menu_bar)

        # Create the Options menu
        options_menu = tk.Menu(menu_bar, tearoff=0)
        options_menu.add_command(label='View System Spec', command=self.check_system_specs)
        options_menu.add_command(label='Plot Trajectories', command=self.plot)

        # Add Options to Menu Bar
        menu_bar.add_cascade(label='Options', menu=options_menu, background='#0052cc', foreground='#ffffff')

        # Detector Settings
        menu_bar.add_command(label='Detector Settings', command=self.detector_settings_widnow, background='#0052cc', foreground='#ffffff')

        # Create the Help menu
        help_menu = tk.Menu(menu_bar, tearoff=0)
        help_menu.add_command(label='User Tutorial', command=self.show_user_tutorial)
        help_menu.add_command(label='View Software Documentation')

        # Add Help to Menu Bar
        menu_bar.add_cascade(label='Help', menu=help_menu, background='#0052cc', foreground='#ffffff')


    def main_window(self):
    
        # Set background color
        self.root.configure(bg='#E0EEEE')
        
        # Create a frame
        self.frame = tk.Frame(self.root, bg='#E0EEEE', padx=20, pady=20)
        self.frame.pack(expand=True, fill='both')
        
        info_label = tk.Label(self.frame, text="Update parameters and press Enter on your keyboard to apply changes.", bg='#E0EEEE', padx=4, pady=4)
        info_label.grid(row=0, column=0, columnspan=3, sticky="W")

        # Simulation Type
        simulation_type_label = tk.Label(self.frame, text="Simulation Type", bg='#E0EEEE', padx=4, pady=4)
        simulation_type_label.grid(row=6, column=0)

        simulation_type_options = ["PySparse", "PyDense", "CPP"]
        simulation_type_combobox = ttk.Combobox(self.frame, values=simulation_type_options)
        simulation_type_combobox.grid(row=6, column=1)  # Removed columnspan
        simulation_type_combobox.bind("<<ComboboxSelected>>", self.update_simulation_type)
        
        # Counting Threshold Label
        tk.Label(self.frame, text="Counting Threshold", bg='#E0EEEE', padx=20, pady=3).grid(row=1, column=0)

        # Counting Threshold Entry
        self.tecount_entry = tk.Entry(self.frame)
        self.tecount_entry.grid(row=1, column=1)
        self.tecount_entry.bind('<Return>', self.update_tecount_entry)

        # Counting Threshold Scale
        self.tecount_scale = tk.Scale(self.frame, from_=5, to=80, orient='horizontal', length=300, bg='#E0EEEE', command=self.update_tecount_scale) 
        self.tecount_scale.grid(row=1, column=2, columnspan=2)


       
        # Electron Energy Label
        tk.Label(self.frame, text="Electron Energy", bg='#E0EEEE', padx=20, pady=3).grid(row=2, column=0)

        # Electron Energy Entry
        self.electron_energy_entry = tk.Entry(self.frame)
        self.electron_energy_entry.grid(row=2, column=1)
        self.electron_energy_entry.bind('<Return>', self.update_energy_entry)

        # Electron Energy Scale
        self.electron_energy_scale = tk.Scale(self.frame, from_=5, to=500, orient='horizontal', length=300, bg='#E0EEEE', command=self.update_energy_scale) 
        self.electron_energy_scale.grid(row=2, column=2, columnspan=2)

        # Electron Dose Label
        tk.Label(self.frame, text="Electron Dose", bg='#E0EEEE', padx=2, pady=3).grid(row=3, column=0)

        # Electron Dose Entry
        self.electron_dose_entry = tk.Entry(self.frame)
        self.electron_dose_entry.grid(row=3, column=1)
        self.electron_dose_entry.bind('<Return>', self.update_electron_dose_entry)

        # Electron Dose Scale
        self.electron_dose_scale = tk.Scale(self.frame, from_=0, to=3000000, orient='horizontal', length=300, bg='#E0EEEE', command=self.update_electron_dose_scale) 
        self.electron_dose_scale.grid(row=3, column=2, columnspan=2)
        tk.Label(self.frame, bg='#E0EEEE', width=10).grid(row=3, column=4) 

        # File path entry
        tk.Label(self.frame, text="Input Image File Path:", bg='#E0EEEE').grid(row=8, column=0)
        self.file_path_entry_in = tk.Entry(self.frame, width=50)
        self.file_path_entry_in.grid(row=8, column=1, columnspan=2)

        # File chooser button
        choose_file_button = tk.Button(self.frame, text="Choose File", command=self.choose_input_file, bg='#C1CDCD', padx=6, pady=6)
        choose_file_button.grid(row=8, column=3)

        # File path entry
        tk.Label(self.frame, text="Output Image File Path:", bg='#E0EEEE', padx=6, pady=6).grid(row=9, column=0)
        self.file_path_entry_out = tk.Entry(self.frame, width=50)
        self.file_path_entry_out.grid(row=9, column=1, columnspan=2)

        # File chooser button
        choose_file_button = tk.Button(self.frame, text="Choose File", command=self.choose_ouput_file, bg='#C1CDCD', padx=6, pady=6)
        choose_file_button.grid(row=9, column=3)
        tk.Label(self.frame, bg='#E0EEEE', width=10).grid(row=17, column=4) 
        
        # set paths
        set_paths_button = tk.Button(self.frame, text="Set File Paths", command=self.set_file_paths, bg='#C1CDCD', padx=6, pady=6)
        set_paths_button.grid(row=10, column=1)

        # Start Simulation Button
        start_button = tk.Button(self.frame, text="Start Simulation", command=self.start_simulation, bg='#32CD32', padx=6, pady=6)
        start_button.grid(row=15, column=2) # Changed column
        
        
    def build_gui(self):
        
        self.main_window()
        self.menu_bar()
        self.create_influence_display_image()
        self.root.mainloop()






    
