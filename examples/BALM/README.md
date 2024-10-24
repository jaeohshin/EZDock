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

## Input Files

The input files are stored in the following folders, which the individual scripts call from:

- Reference complexes (naming convention: `` {target}_protein.pdb ``) and ligands (`` {target}_ligands_{i}.sdf ``) for MCL1, SYK, and HIF2A: [mcl1_syk_hif2a](https://github.com/meyresearch/EZDock/tree/main/examples/BALM/mcl1_syk_hif2a/). These ligands come from [protein-ligand-benchmark](https://github.com/openforcefield/protein-ligand-benchmark/tree/main/data), and come in an SDF file for each target, each containing dozens of ligands. The docking script used for these 3 targets splits the SDF files into individual SDF files containing one ligand each, with the order they appear in the original file as their naming convention. (see the [split_sdf](https://github.com/meyresearch/EZDock/blob/676fabd72c8babe9e7f0f9f063a1e10ef9c98a51/examples/BALM/mcl1_syk_hif2a_docking_vina.py#L570) function on how this was done). 

    Alternatively, if you want to figure out which specific ligand was in the original SDF file, see the CSV files named ``{target}_name_mapping.csv`` to find out.
- Reference complex `` Mpro-protein.pdb `` for Mpro docking: [Mpro_complexes](https://github.com/meyresearch/EZDock/tree/main/examples/BALM/Mpro_complexes).
- Mpro ligands to dock to the protein: [Updated_Mpro_ligands](https://github.com/meyresearch/EZDock/tree/main/examples/BALM/Updated_Mpro_ligands).

For LP-PDBBind, the input data is stored in the CSV file [leakypdb.csv](https://github.com/meyresearch/EZDock/blob/main/examples/BALM/leakypdb_test.csv) as a list of PDB IDs, which the script calls from to download the necessary protein-ligand complex files before continuing to process the files accordingly in the same way that the other scripts do.

## Output Files

Your output will be saved in the ``data`` folder and the free energy of binding results collected in a CSV file such as ``exp_data.csv``. 

## Note

Do NOT run multiple docking scripts at the same time in the same folder. If you need to run multiple docking experiments in parallel, copy and paste this folder and run your subsequent docking experiments there. 