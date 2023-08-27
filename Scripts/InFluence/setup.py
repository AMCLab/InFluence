from setuptools import setup

setup(
    name='InFluence',
    version='1.0',
    py_modules=['InFluence_CLI', 'main'],
    install_requires=[
        'argparse',
        'configparser',
        'inquirer',
        'matplotlib',
        'numpy',
    ],
    entry_points={
        'console_scripts': [
            'InFluence=InFluence_CLI:main',
        ],
    },
)
