This folder contains the scripts necessary to run a 4DSTEM simulation, add detector effects and finite beam fluences and then process the data with the single side band ptychographic reconstruction algorithm. 

Set up the conda environment with the my_env.yml file

conda-env create -n my_env -f= my_env.yml


The scripts should be run in the following order:

1. multislice_script.ipynb
2. InFluence_ver2.ipynb
3. ssbp_script.ipynb
4. ssb_ssim_script.ipynb

Read https://abtem.readthedocs.io/en/latest/ and  https://ptychography-4-0.github.io/ptychography/algorithms.html for more information on the multislice and single side band ptychography experiments. In general, the only things that need to be changed are the paths. The single side band requires manual changes to the parameters in cell 6. These must be worked out by inspection of the patterns. 
The user can tune the parameters defining the detector and electron probe in InFluence_ver2 script. 
