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

<font color='red'>
### Testing the pipeline

To test the robustness of this pipeline, we use the PDB IDs from LeakyPDB to obtain some free energy of binding results. The pipeline has error traps set up to capture and record any error encountered during the experiment. We will inspect the output files to investigate any causes of failure and seek to improve the pipeline this way.
</font>

### Electronic Lab Notebook (ELN)

I keep an detailed, active log of activities I've been up to with this repo in the folder [``ELN``](https://github.com/meyresearch/SILVR2/tree/evaluation/ELN) as supplementary information, though all relevant information will be updated to this README document. Feel free to check them out!

<font color='red'>
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
Set Up OpenEye and Vina:

Vina: Install and configure AutoDock Vina and ADFR suite.
OpenEye: Ensure that you have access to the OpenEye toolkit and have set up your license.
Usage
Running the Pipeline
To run the docking pipeline, execute the following command:

bash
复制代码
python mpro_oe_w_water_v2_error_handling.py --pdb_id <PDB_ID> --ligand_file <LIGAND_FILE>
Parameters
--pdb_id: The PDB ID of the protein-ligand complex to download.
--ligand_file: The path to the ligand file to be used in docking.
Example
bash
复制代码
python mpro_oe_w_water_v2_error_handling.py --pdb_id 1abc --ligand_file ligand.mol2
Output
The results of the docking experiments are stored in the results/ directory, with each run creating a timestamped folder containing:

Docking Scores: A CSV file with the docking scores and any errors.
Prepared Files: The prepared receptor and ligand files.
Logs: Detailed logs of the docking process for troubleshooting.
Contributing
Contributions are welcome! Please fork this repository and submit a pull request with your improvements or bug fixes.

License
This project is licensed under the MIT License - see the LICENSE file for details.

Acknowledgments
Vina: AutoDock Vina
OpenEye: OpenEye Scientific Software
References
Refer to the detailed documentation provided within the docs/ folder for a complete breakdown of each component in the pipeline.

</font>
