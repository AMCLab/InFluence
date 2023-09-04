import configparser


class Filepaths:
    config = configparser.ConfigParser()
    config.read('config.ini')

    InputImageLocation = config.get('Filepaths', 'img_in')
    SaveLocation = config.get('Filepaths', 'img_out')


