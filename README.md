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

## Output

A selected number of results of the docking experiments are stored in the results/ directory, with each run creating a folder containing:

- Docking Scores: A CSV file with the docking scores and any errors.
- Starting Files: The files that used as input.
- Prepared Files: The prepared receptor and ligand files.
- Output Files: The docked ligand pose(s). (For Vina, this is {exp}_out.pdbqt; for OE, this is {exp}_docked)
- Logs: Detailed logs of the docking process for troubleshooting.

## Contributing
Contributions are welcome! Please fork this repository and submit a pull request with your improvements or bug fixes.

<span style="color:red">

## Data for MSc Project

[Click here](https://uoe-my.sharepoint.com/:u:/g/personal/s1732775_ed_ac_uk/EZn9Vb2VxoxOont6QMtDSg0BNgIGWPrm_rcKl6ZrWCAIGw?e=N2TKmb) to see the data files generated for the experiments conducted for robustness testing and for the diffusion model-generated molecules.

</span>

## Acknowledgments

[AutoDock Vina](https://vina.scripps.edu/)

[OpenEye Toolkits: Chemgauss4](https://docs.eyesopen.com/toolkits/cpp/dockingtk/scoring.html#section-scoring-chemgauss4)
