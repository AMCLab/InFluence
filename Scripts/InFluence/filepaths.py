import os
from pathlib import Path
import configparser

class Filepaths:
    # Get the directory where this Python script resides
    current_dir = Path(os.path.dirname(os.path.abspath(__file__)))
    
    # Define the location of the config.ini relative to the script's location
    config_path = current_dir / "config.ini"

    config = configparser.ConfigParser()
    config.read(config_path)
    
    InputImageLocation = config.get('Filepaths', 'img_in')
    SaveLocation = config.get('Filepaths', 'img_out')
