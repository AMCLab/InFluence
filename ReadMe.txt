# InFluence




## Running from main.py

You can run influence from main.py by navigating into the InFluence directory and running python3 main.py.
You can folow the comments in the file to understand how to use the main.py file.

## Directly Running CLI 

You can navigate to the InFluence/interface_code and run python3 InFluenceCLI.py along with any commands listed in the CLI section below.



## Running GUI directly from source code:


You can navigate to the InFluence/interface_code and run python3 InFluenceCLI.py --gui and the gui will launch.




## Package  Installation

1. **Navigate to the InFluence Directory**: Open your terminal and navigate to the directory where you have InFluence installed. Use the `cd` command to change your working directory.

   ```bash
   cd /path/to/InFluence
Build the Wheel Distribution: To build the wheel distribution, run the following command:
bash
Copy code
python3 setup.py bdist_wheel
This command generates a wheel distribution file for InFluence.
Install the Package: Now, you can install the package using pip. Use the following command, replacing InFluence-1.0-py3-none-any.whl with the actual filename generated in Step 2:
bash
Copy code
pip install dist/InFluence-1.0-py3-none-any.whl
Usage

To use InFluence, you can simply launch the Command Line Interface (CLI) with the following command:

bash
Copy code
InFluence
This command will start the InFluence CLI, allowing you to interact with the package and its functionalities.


## Running tests 

To run tests run the command PyTest in the Tsts directory.






