# Example: BALM SI Section Vina Docking

This repository guides you to recreate the data generated for the SI section of the BALM paper, which uses AutoDock Vina as a baseline for benchmarking binding affinity data. The guidance written here assumes that you have already cloned and installed this pipeline, as the Python scripts require the required packages to run correctly. 

The nomenclature of the Python scripts are self-explanatory as to which datasets they will recreate. To run them, simply open the command line and run the following command:

`` python python_script.py ``

For example, to recreate the results for Mpro with water, execute this command: 

`` python mpro_vina_water_v2_2.py ``

This operation will run for approximately 1.5 days on 20 cores of CPU, so to avoid losing progress, you can use the `` screen `` function so that it runs in the background:

```
# Opening a new screen session
screen -R docking session

# Activate the required environment
conda activate EZDock

# Make sure you're in the right folder, then run the docking
python mpro_vina_water_v2_2.py
```
Your output will be saved in the ``data`` folder and the free energy of binding results collected in a CSV file such as ``exp_data.csv``. 

Note: Do NOT run multiple docking scripts at the same time in the same folder. If you need to run multiple docking experiments in parallel, copy and paste this folder and run your subsequent docking experiments there. 