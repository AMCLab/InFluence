import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
from PIL import Image, ImageTk




from PIL import Image, ImageDraw, ImageFont
import numpy as np

# Create an image with a black background
width, height = 128, 128
image = Image.new('RGB', (width, height), 'black')
draw = ImageDraw.Draw(image)

# Define the color
color = '#E0EEEE'

# Load a font
font = ImageFont.load_default()

# Calculate text size and position to center the text
text = "InFluence"
text_width, text_height = draw.textsize(text, font)
text_width = text_width
text_height = text_height

text_x = (width - text_width) // 2
text_y = (height - text_height) // 2

# Draw the text
draw.text((text_x, text_y), text, fill=color, font=font)

influnce_image = image

def start_simulation():
    electron_energy = electron_energy_entry.get()
    atomic_number = atomic_number_entry.get()
    atomic_weight = atomic_weight_entry.get()
    material_density = material_density_entry.get()
    electron_dose = electron_dose_scale.get()
    print(f"Starting Simulation with Electron Energy: {electron_energy} keV, Atomic Number: {atomic_number}, Atomic Weight: {atomic_weight}, Material Density: {material_density} SI units, Electron Dose: {electron_dose}")


   
root = tk.Tk()
root.title("InFluence")

# Tooltip Label
def create_tooltip(widget, text):
    tooltip = tk.Toplevel(widget)
    tooltip.wm_overrideredirect(True)
    tooltip.wm_withdraw()
    tooltip_label = tk.Label(tooltip, text=text, relief='solid', borderwidth=1, background="#FFFFE0")
    tooltip_label.pack()

    def show_tooltip(event):
        tooltip.wm_geometry(f"+{event.x_root+20}+{event.y_root+20}")
        tooltip.wm_deiconify()

    def hide_tooltip(event):
        tooltip.wm_withdraw()

    widget.bind("<Enter>", show_tooltip)
    widget.bind("<Leave>", hide_tooltip)
    
    
def choose_file():
    file_path = filedialog.askopenfilename()
    file_path_entry.delete(0, tk.END)  # Remove current text in entry
    file_path_entry.insert(0, file_path)  # Insert the file path
    
    
def update_scale(*args):
    electron_dose_entry.delete(0, tk.END)
    electron_dose_entry.insert(0, str(electron_dose_scale.get()))

def update_entry(*args):
    value = int(electron_dose_entry.get())
    electron_dose_scale.set(value)
    

def update_energy_scale(*args):
    electron_energy_entry.delete(0, tk.END)
    electron_energy_entry.insert(0, str(electron_energy_scale.get()))


def update_energy_entry(*args):
    value = int(electron_energy_entry.get())
    electron_energy_scale.set(value)
    

    
# Set background color
root.configure(bg='#E0EEEE')

# Create a frame
frame = tk.Frame(root, bg='#E0EEEE', padx=20, pady=20)
frame.pack(expand=True, fill='both')



# Function to open Detector Settings window
def open_detector_settings():
    detector_window = tk.Toplevel(root, bg='#E0EEEE', padx=20, pady=20)
    detector_window.title('InFluence: Detector Settings')

    
    # Detector Threshold
    tk.Label(detector_window, text="Detector Threshold", bg='#E0EEEE', padx=4, pady=4).grid(row=0, column=0)
    detector_threshold_entry = tk.Entry(detector_window)
    detector_threshold_entry.grid(row=0, column=1)
    tk.Label(detector_window, text="(keV)", bg='#E0EEEE').grid(row=0, column=2)

    # Atomic Number
    tk.Label(detector_window, text="Proton Number", bg='#E0EEEE', padx=4, pady=4).grid(row=1, column=0)
    atomic_number_entry = tk.Entry(detector_window)
    atomic_number_entry.grid(row=1, column=1)

    # Atomic Weight
    tk.Label(detector_window, text="Atomic Mass", bg='#E0EEEE', padx=4, pady=4).grid(row=2, column=0)
    atomic_weight_entry = tk.Entry(detector_window)
    atomic_weight_entry.grid(row=2, column=1)

    # Material Density
    tk.Label(detector_window, text="Material Density", bg='#E0EEEE', padx=4, pady=4).grid(row=3, column=0)
    material_density_entry = tk.Entry(detector_window)
    material_density_entry.grid(row=3, column=1)
    tk.Label(detector_window, text="(m^-3)", bg='#E0EEEE').grid(row=3, column=2)


    # Pixel Dimensions
    tk.Label(detector_window, text="Pixel Dimensions:", bg='#E0EEEE', padx=6, pady=6).grid(row=4, column=0)
    tk.Label(detector_window, text="x:", bg='#E0EEEE').grid(row=5, column=0)
    pixel_dimension_x_entry = tk.Entry(detector_window, width=10)
    pixel_dimension_x_entry.grid(row=5, column=1)

    tk.Label(detector_window, text="y:", bg='#E0EEEE', padx=6, pady=6).grid(row=6, column=0)
    pixel_dimension_y_entry = tk.Entry(detector_window, width=10)
    pixel_dimension_y_entry.grid(row=6, column=1)

    tk.Label(detector_window, text="z:", bg='#E0EEEE', padx=6, pady=6).grid(row=7, column=0)
    pixel_dimension_z_entry = tk.Entry(detector_window, width=10)
    pixel_dimension_z_entry.grid(row=7, column=1)
    
    save_button = tk.Button(detector_window, text="Save Settings",  bg='#C1CDCD', padx=6, pady=6)
    save_button.grid(row=8, column = 4, columnspan=5)
    

def open_help_window():
    help_window = tk.Toplevel(root)
    help_window.title("Help")
    tk.Label(help_window, text="This is the help window. Add your content here.").pack()
    

# Create a frame for the menu bar with a different background color
menu_frame = tk.Frame(root, bg='#0052cc')
menu_frame.pack(side='top', fill='x')

# Create a menu bar
menu_bar = tk.Menu(menu_frame, bg='#0052cc', fg='#ffffff')
root.config(menu=menu_bar)

# Detector Settings
menu_bar.add_command(label='Detector Settings', command=open_detector_settings, background='#0052cc', foreground='#ffffff')

# Add Help to Menu
menu_bar.add_command(label='Help', command=open_help_window, background='#0052cc', foreground='#ffffff')






# Simulation Type
simulation_type_label = tk.Label(frame, text="Simulation Type", bg='#E0EEEE', padx=4, pady=4)
simulation_type_label.grid(row=9, column=0)

simulation_type_options = ["GPU", "CPU Python Sparse", "CPU Python Dense", "CPU Python Plot", "CPU C++"]
simulation_type_combobox = ttk.Combobox(frame, values=simulation_type_options)
simulation_type_combobox.grid(row=9, column=1)  # Removed columnspan

# Create a canvas for the question mark button
question_mark_canvas = tk.Canvas(frame, width=20, height=20, bg='#E0EEEE', highlightthickness=0)
question_mark_canvas.grid(row=9, column=2)  # Changed column

# Draw a circle in the canvas
question_mark_canvas.create_oval(2, 2, 18, 18, outline="black", fill='#E0EEEE')

# Place a label with the question mark inside the canvas
question_mark_label = tk.Label(question_mark_canvas, text="?", bg='#E0EEEE')
question_mark_canvas.create_window(10, 10, window=question_mark_label)

# Add the tooltip to the canvas
create_tooltip(question_mark_canvas, "Choose your simulation type.\n\nIf you have a GPU available, GPU execution is recommended.\nFill in other recommendations later.")


# Electron Energy Label
tk.Label(frame, text="Electron Energy", bg='#E0EEEE', padx=20, pady=3).grid(row=12, column=0)

# Electron Energy Entry
electron_energy_entry = tk.Entry(frame)
electron_energy_entry.grid(row=12, column=1)
electron_energy_entry.bind('<Return>', update_energy_entry)

# Electron Energy Scale
electron_energy_scale = tk.Scale(frame, from_=0, to=500, orient='horizontal', length=300, bg='#E0EEEE', command=update_energy_scale) 
electron_energy_scale.grid(row=12, column=2, columnspan=2)



# Electron Dose Label
tk.Label(frame, text="Electron Dose", bg='#E0EEEE', padx=20, pady=3).grid(row=14, column=0)

# Electron Dose Entry
electron_dose_entry = tk.Entry(frame)
electron_dose_entry.grid(row=14, column=1)
electron_dose_entry.bind('<Return>', update_entry)

# Electron Dose Scale
electron_dose_scale = tk.Scale(frame, from_=0, to=3000000, orient='horizontal', length=300, bg='#E0EEEE', command=update_scale) 
electron_dose_scale.grid(row=14, column=2, columnspan=2)
tk.Label(frame, bg='#E0EEEE', width=10).grid(row=14, column=4) 





# File path entry
tk.Label(frame, text="Input Image File Path:", bg='#E0EEEE').grid(row=18, column=0)
file_path_entry = tk.Entry(frame, width=50)
file_path_entry.grid(row=18, column=1, columnspan=2)

# File chooser button
choose_file_button = tk.Button(frame, text="Choose File", command=choose_file, bg='#C1CDCD', padx=6, pady=6)
choose_file_button.grid(row=18, column=3)

# File path entry
tk.Label(frame, text="Output Image File Path:", bg='#E0EEEE', padx=6, pady=6).grid(row=19, column=0)
file_path_entry = tk.Entry(frame, width=50)
file_path_entry.grid(row=19, column=1, columnspan=2)

# File chooser button
choose_file_button = tk.Button(frame, text="Choose File", command=choose_file, bg='#C1CDCD', padx=6, pady=6)
choose_file_button.grid(row=19, column=3)

tk.Label(frame, bg='#E0EEEE', width=10).grid(row=17, column=4) 

# Start Simulation Button
start_button = tk.Button(frame, text="Start Simulation", command=start_simulation, bg='#32CD32', padx=6, pady=6)
start_button.grid(row=20, column=6) # Changed column

# Stop Simulation Button
stop_button = tk.Button(frame, text="Stop Simulation", command=start_simulation, bg='#F08080', padx=6, pady=6) # Changed button name and command
stop_button.grid(row=21, column=6) # Changed column


tk.Label(frame, bg='#E0EEEE', width=10).grid(row=21, column=4) 

# Progress Bar
progress_label = tk.Label(frame, text="Simulation Progress", bg='#E0EEEE')
progress_label.grid(row=22, column=0) # Adjust the row and column as needed
progress_bar = ttk.Progressbar(frame, orient="horizontal", length=400, mode="determinate")
progress_bar.grid(row=22, column=1, columnspan=2)



# Convert the PIL Image to a PhotoImage


resized_image = influnce_image.resize((700,700), Image.ANTIALIAS)

influence_image_tk = ImageTk.PhotoImage(image=resized_image)


# Create a label to hold the image
image_label = tk.Label(frame, image=influence_image_tk, bg='#E0EEEE', padx=40, pady=40)
image_label.grid(row=0, column=6, rowspan=15)  

# Keep a reference to avoid garbage collection
image_label.image = influence_image_tk


root.mainloop()
