from setuptools import setup, find_packages

setup(
    name='InFluence',
    version='1.0',
    packages=find_packages(),  # Automatically discover and include all packages
    install_requires=[
        'argparse',
        'configparser',
        'inquirer',
        'matplotlib',
        'numpy',
    ],
    entry_points={
        'console_scripts': [
            'InFluence=interface_code.InFluence_CLI:main',  
        ],
    },
)
