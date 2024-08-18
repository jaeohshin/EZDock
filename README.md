# EZDock: Automated Docking Pipeline

## Overview

This repository contains a comprehensive pipeline designed to evaluate protein-ligand interactions using molecular docking. The pipeline automates several critical steps in the docking process, making it easier and faster to conduct docking experiments. It supports multiple docking software tools, including **Vina** and **OpenEye (OE)**, and is designed for evaluating small molecules against biological targets. This is useful for a large compound library or molecules generated from a diffusion model, which this page will also use.

## Key Features

- **Automated Pipeline**: Automates the retrieval, preparation, and docking of protein-ligand complexes. Scripts used for docking are detailed in the ``pipeline_scripts`` folder.
- **Flexibility**: Supports docking with both Vina and OpenEye, providing flexibility depending on your needs.
- **Error Handling**: Includes robust error handling and logging to ensure smooth execution of the docking experiments.
- **Extensive Documentation**: Detailed descriptions of each part of the pipeline, including the methods for preparation, docking, and evaluation, detailed in the Jupyter Notebooks in the ``jupyter_notebooks`` folder.

## Pipeline Stages

### 1. Download and Split Protein-Ligand Complex (Optional)
- **Download**: Retrieves protein-ligand complexes from the PDB database.
- **Splitting**: Separates the complex into individual protein (receptor) and ligand files.
- **Fixing Receptor Structure**: Optional fixing of the receptor file for docking purposes, including adding missing residues and hydrogens.

### 2. Docking Preparation
- **Receptor Preparation**: Prepares the receptor file with the necessary modifications (hydrogen addition, metal ion preservation, etc.) for docking.
- **Ligand Preparation**: Ensures the ligand file is in the correct format and properly configured for interaction with the binding site.
- **Design Unit Creation (OE Only)**: Generates design units from the protein-ligand complex for OpenEye docking.

### 3. Docking and Scoring
- **Docking**: Executes the docking experiments using Vina or OE.
- **Scoring**: Extracts and records the binding scores, such as ΔG (Vina) or Chemgauss4 (OE), into a CSV file.
- **Data Storage**: Systematically stores all inputs, outputs, and intermediate files, ensuring reproducibility.

### 4. Evaluation
- **Error Rate Analysis**: Monitors the pipeline's performance by calculating the error rate.
- **Docking Results Evaluation**: Benchmarks the computed binding free energies against experimental values or cross-validates them against other docking platforms.

## Installation

1. **Clone the Repository**:
   ```
   git clone https://github.com/meyresearch/EZDock.git
   ```
   
2. **Install Dependencies**:
Make sure to install all required Python libraries as listed in requirements.txt.
    ```
    pip install -r requirements.txt
    ```

3. **If there is any missing packages**:
Go to [step_1b_docking_complex_to_dG.ipynb](https://github.com/meyresearch/EZDock/blob/main/jupyter_notebooks/step_1b_docking_complex_to_dG.ipynb) to see if anything is missing. Uncomment the notebook sections that are missing and execute to install them. 

## Using the pipeline 
The pipeline is done within the python script [EZDock_vina.py](https://github.com/meyresearch/EZDock/blob/main/EZDock_vina.py) and [EZDock_OE.py](https://github.com/meyresearch/EZDock/blob/main/EZDock_OE.py). 

To run this pipeline, you need to edit the python script and specify the paths for the ligand, receptor and reference complex you wish to use for docking, specifically these few lines:
```
all_ligands = glob.glob('Ligand_Folder/*.sdf') # change this to the path where the ligands are stored
receptor_file = 'receptor.pdb' # change this to the name of the receptor
complex_file = 'complex.pdb' # change this to the name of the reference complex
os.system(f'cp Complex_Folder/{complex_file} .') # change this to the path where the complexes are stored
os.system(f'cp Receptor_Folder/{receptor_file} .') # change this to the path where the receptors are stored
```

Then, you can specify whether you'd like the pipeline to preserve water, preserve metal, and specify the name of the output CSV file containing all the docking results, by changing the corresponding variables in this funcion：
```
oe_process_lig_prot(ligand_file, 
                    receptor_file, 
                    complex_file, 
                    preserve_water=True, # default is False
                    preserve_metal=True, # default is True
                    csv_out_file='oe_docking_data.csv') # output csv file
```
or in the Vina pipeline: 

```
vina_process_lig_prot(ligand_file, 
                        receptor_file, 
                        complex_file, 
                        preserve_water=True, # default is False
                        preserve_metal=True, # default is True
                        csv_out_file='oe_docking_data.csv') # output csv file
```

## Output

A selected number of results of the docking experiments are stored in the results/ directory, with each run creating a folder containing:

- Docking Scores: A CSV file with the docking scores and any errors encountered after running a high-throughput screening. There may be cleaned versions with only valid data kept, and updated with additional information like ligand SMILES.
- Starting Files: The files that used as input.
- Prepared Files: The prepared receptor and ligand files.
- Output Files: The docked ligand pose(s). (For Vina, this is {exp}_out.pdbqt; for OE, this is {exp}_docked.sdf and there might be 2 files, of which only one is the ligand in the correct place - only the score of the correct poses are kept in the resulting CSV file.)
- Logs: Detailed logs of the docking process for troubleshooting.

## Future Work

1. Use of ``argparse`` to allow users to define variables on the commandline and use the script as an executable. 

## Contributing
Contributions are welcome! Please fork this repository and submit a pull request with your improvements or bug fixes.

<span style="color:red">

## Data for MSc Project

[Click here](https://uoe-my.sharepoint.com/:u:/g/personal/s1732775_ed_ac_uk/EZn9Vb2VxoxOont6QMtDSg0BNgIGWPrm_rcKl6ZrWCAIGw?e=N2TKmb) to see the data files generated for the experiments conducted for robustness testing and for the diffusion model-generated molecules. You will need your UoE login. 

Note: this is 17 GB in a zipped file!

</span>

## Acknowledgments

[AutoDock Vina](https://vina.scripps.edu/)

[OpenEye Toolkits: Chemgauss4](https://docs.eyesopen.com/toolkits/cpp/dockingtk/scoring.html#section-scoring-chemgauss4)
