import os
# giving permissions to run vina and ADFR scripts, change the directories of vina as required
os.chmod('./vina/vina_1.2.5_linux_x86_64', 0o755) # giving permissions to run vina scripts
os.chmod('./vina/vina_split_1.2.5_linux_x86_64', 0o755) # giving permissions to run vina scripts
os.chmod('./extract_ligand.sh', 0o755) # giving permissions to run extract_ligand.sh
os.chmod('./ADFRsuite/ADFRsuite_x86_64Linux_1.0', 0o755) # giving permissions to run ADFRsuite

import os
from Bio.PDB import PDBList, PDBParser, Select, PDBIO
from subprocess import Popen, PIPE
import re
import logging
from pathlib import Path
import contextlib
import subprocess
from rdkit import Chem
from rdkit.Chem import rdMolTransforms as rdmt
import MDAnalysis as mda
import numpy as np
import pandas as pd
from tqdm import tqdm 
import matplotlib.pyplot as plt
import seaborn as sns
from pdbfixer import PDBFixer
from openmm.app import PDBFile
from rdkit import Chem
from rdkit.Chem import SDWriter
from ase import Atoms
from ase.io.sdf import read_sdf
from ase.io import read
from iodata import load_one
from iodata.utils import angstrom

from openmm.app.element import zinc, iron, calcium  # Import other metals as needed
import openmm.unit as unit

# Set up the working directory
cwd = os.getcwd()

# Function to change directories
@contextlib.contextmanager
def set_directory(dirname: os.PathLike, mkdir: bool = False):
    pwd = os.getcwd()
    path = Path(dirname).resolve()
    if mkdir:
        path.mkdir(exist_ok=True, parents=True)
    os.chdir(path)
    yield path
    os.chdir(pwd)

# Function to run shell commands
def run_command(cmd, raise_error=True, input=None, timeout=None, **kwargs):
    """Run a shell command and handle possible errors."""
    # Popen is used to run the command
    sub = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE, **kwargs)
    if input is not None:
        sub.stdin.write(bytes(input, encoding='utf-8'))
    try:
        out, err = sub.communicate(timeout=timeout)
        return_code = sub.poll()
    # if the command times out, kill the process
    except subprocess.TimeoutExpired:
        sub.kill()
        print(f"Command {cmd} timeout after {timeout} seconds")
        return 999, "", ""  # 999 is a special return code for timeout
    out = out.decode('utf-8')
    err = err.decode('utf-8')
    if raise_error and return_code != 0:
        raise CommandExecuteError(f"Command {cmd} failed: \n{err}")
    return return_code, out, err

# Exception for command execution errors
class CommandExecuteError(Exception):
    def __init__(self, msg):
        self.msg = msg
    def __str__(self):
        return self.msg

# using Biopython to extract the receptor
class ReceptorSelect(Select):
    '''This class is used to select only the receptor residues from a PDB file, by excluding any HETATM and any water.'''

    def __init__(self, preserve_water=False, preserve_metal=False):
        self.preserve_water = preserve_water
        # List of common heavy metal ions in protein structures
        self.preserve_metal = preserve_metal
        # I deleted the heavy metal ions as they are not a valid AutoDock type! Be cautious here!
        self.heavy_metal_ions = ['ZN', 'FE', 'MG', 'CA', 'MN', 'CO', 'CU', 'NI', 'MO', 'W','YB', 'K', 'NA']

    def accept_residue(self, residue):
        # Exclude water molecules, assumed to have residue name 'HOH' or 'WAT'
        #TODO: What if there's water in the binding site?
        resname = residue.get_resname()
        # print(f"Processing residue: {resname}")

        if resname in ['HOH', 'WAT']:
            if self.preserve_water:
                # print("Preserving water molecule")
                return True
            else:
                # print("Excluding water molecule")
                return False

        # Include heavy metal ions
        if resname in self.heavy_metal_ions:
            if self.preserve_metal:
                # print("Preserving metal ion")
                return True
            else:
                # print("Excluding metal ion")
                return False
        
        # Check if any atom in the residue is from an ATOM record
        for atom in residue:
            '''in biopython, atom.get_full_id()[3] is accessing the fourth element of the full ID. This element represents the atom name, atom.get_full_id()[3][0] is accessing the first character of the atom name. This character represents the element symbol of the atom. The condition atom.get_full_id()[3][0] == ' ' checks whether the first character of the atom name is a space. If it is a space, then the atom is from an ATOM record, otherwise it is from a HETATM record.'''
            if atom.get_full_id()[3][0] == ' ':
                return True
            
        # print("Excluding residue")
        return False

def download_pdb_file(pdb_id):
    """Download PDB file using PDB ID."""
    pdbl = PDBList()
    filename = pdbl.retrieve_pdb_file(pdb_id, file_format='pdb', pdir='.', overwrite=True)
    return filename

def get_structure(filename):
    """Parse the structure from a PDB file."""
    parser = PDBParser()
    structure = parser.get_structure('structure', filename)
    return structure

def save_receptor(structure, receptor_file, preserve_water=False, preserve_metal=False):
    """Save the receptor part of the PDB file."""
    io = PDBIO()
    io.set_structure(structure)
    io.save(receptor_file, ReceptorSelect(preserve_water=preserve_water,
    preserve_metal = preserve_metal)) 

# Example usage:
# pdb_id = '1abc'
# filename = download_pdb_file(pdb_id)
# structure = get_structure(filename)
# save_receptor(structure, 'receptor.pdb', preserve_water=True)

def fix_receptor_file(input_receptor_file):
    '''This function must be used for complexes without crystal structurers on PDB, as there could be places where things should be fixed before docking. This must be done or else the docking software might report errors. '''

    #TODO: Fix this function - it unwantedly removes metal ions and water from the structure

    # Initialize the PDBFixer with the input file
    fixer = PDBFixer(filename=input_receptor_file)
    
    # Find and fix missing residues, atoms, and hydrogens
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    
    # Define the output file name
    output_receptor_file = input_receptor_file.replace('.pdb', '_fixed.pdb')
    
    # Save the fixed receptor file
    with open(output_receptor_file, 'w') as output_file:
        PDBFile.writeFile(fixer.topology, fixer.positions, output_file)
    
    print(f"Fixed receptor file saved as: {output_receptor_file}")
    
    return output_receptor_file


def rename_file(old_filename, new_filename):
    """Rename the downloaded file to a standard name."""
    if os.path.exists(old_filename) and not os.path.exists(new_filename):
        os.rename(old_filename, new_filename)
    return new_filename


class DockingPrepper:
    '''This class is used to prepare the ligand and protein files for docking using AutoDock Vina. It uses the ADFRsuite to convert the ligand and protein files to pdbqt format.
    # Example usage:
    prep = DockingPrepper(folder='path/to/folder', pdb_id = 'pdb_id.pdb', lig_file='ligand.pdb', prot_file='protein.pdb')
    
    # ---- AutoDock Vina Preparation ----
    prep.vina_process()

    '''
    def __init__(self, folder, lig_file, prot_file, lig_out_name=None, prot_out_name=None, add_h=True, preserve_water=False, preserve_metal=False):
        self.folder = folder
        self.lig_file = lig_file
        self.prot_file = prot_file
        self.lig_out_name = lig_out_name or f"{lig_file.split('.')[0]}.pdbqt"
        self.prot_out_name = prot_out_name or f"{prot_file.split('.')[0]}.pdbqt"
        self.add_h = add_h
        self.preserve_water = preserve_water
        self.preserve_metal = preserve_metal #TODO/WARMING: This is not implemented yet!

    # ----------------- AutoDock Vina Preparation -----------------
    def vina_prepare_ligand(self):
        """Prepare the ligand file for Vina."""
        if self.add_h:
            cmd = f'./ADFRsuite/ADFRsuite_x86_64Linux_1.0/bin/prepare_ligand -l {self.lig_file} -o {self.lig_out_name} -A hydrogens'
        else:
            cmd = f'./ADFRsuite/ADFRsuite_x86_64Linux_1.0/bin/prepare_ligand -l {self.lig_file} -o {self.lig_out_name}'
        self._run_command(cmd, f'The ligand file has been saved as {self.lig_out_name}')

    def vina_prepare_protein(self):
        """Prepare the protein file for Vina."""
        if self.preserve_water:
            cmd = f'./ADFRsuite/ADFRsuite_x86_64Linux_1.0/bin/prepare_receptor -r {self.prot_file} -o {self.prot_out_name} -A checkhydrogens -U nphs_lps_nonstdres'
            '''        
            [-U]  cleanup type:
             'nphs': merge charges and remove non-polar hydrogens
             'lps': merge charges and remove lone pairs
             'waters': remove water residues
             'nonstdres': remove chains composed entirely of residues of
                      types other than the standard 20 amino acids
             'deleteAltB': remove XX@B atoms and rename XX@A atoms->XX
             (default is 'nphs_lps_waters_nonstdres') 
             '''
        else:
            cmd = f'./ADFRsuite/ADFRsuite_x86_64Linux_1.0/bin/prepare_receptor -r {self.prot_file} -o {self.prot_out_name} -A checkhydrogens'
        self._run_command(cmd, f'The protein file {self.prot_file} has been converted to pdbqt format and saved in {self.prot_out_name}')

    def _run_command(self, cmd, success_message):
        """Run a shell command and handle possible errors."""
        try:
            subprocess.run(cmd, shell=True, check=True, timeout=999)
            print(success_message)
        except subprocess.CalledProcessError as e:
            print(f"Error occurred: {e.output}")
            raise e
        except subprocess.TimeoutExpired:
            print(f"Command '{cmd}' timed out.")
            raise e

    def vina_process(self):
        """Process the ligand and protein files."""
        with set_directory(self.folder):
        # if the lig_file and prot_file are not in pdb format, convert them to pdb format using obabel
            if self.lig_file.split('.')[-1] != 'pdb':
                lig_name = self.lig_file.split('.')[0]
                cmd = f'obabel {self.lig_file} -O {lig_name}.pdb'
                self._run_command(cmd, f'The ligand file has been converted to {lig_name}.pdb')
                self.lig_file = f'{lig_name}.pdb'

            if self.prot_file.split('.')[-1] != 'pdb':
                prot_name = self.prot_file.split('.')[0]
                cmd = f'obabel {self.prot_file} -O {prot_name}.pdb'
                self._run_command(cmd, f'The protein file has been converted to {prot_name}.pdb')
                self.prot_file = f'{prot_name}.pdb'
                
            self.vina_prepare_ligand()
            self.vina_prepare_protein()


def get_COM(file):
        if file.endswith('mol2') or file.endswith('xyz'):
            mol = load_one(file)
            ase_mol = Atoms(numbers=mol.atnums, positions=mol.atcoords / angstrom)
        elif file.endswith('sdf'):
            ase_mol = read_sdf(file)
        elif file.endswith('pdb'):
            ase_mol = read(file)
        else:
            raise NotImplementedError(f"file extension not supported for {file}")

        return ase_mol.get_center_of_mass()
    
def write_config_vina(lig_pdbqt,prot_pdbqt,center, config_fp = "config.txt", weights=None, boxsize=50, exhaustiveness=32, num_modes=1, energy_range=30, **kwargs):
    '''
    Write the config file for AutoDock Vina docking
    :param exhaustiveness: int, the exhaustiveness of the docking
    :param num_modes: int, the number of modes (conformations) to be generated
    :param energy_range: int, the energy range of the docking
    '''

    lines = ["receptor = {}".format(prot_pdbqt),
            "ligand = {}".format(lig_pdbqt),
            "scoring = vina",
            "",
            "center_x = {}".format(center[0]),
            "center_y = {}".format(center[1]),
            "center_z = {}".format(center[2]),
            "",
            "size_x = {}".format(boxsize),
            "size_y = {}".format(boxsize),
            "size_z = {}".format(boxsize),
            "",
            # "exhaustiveness = {}".format(exhaustiveness),
            # "num_modes = {}".format(num_modes),
            # "energy_range = {}".format(energy_range),
            ]
    if weights is not None:
        assert len(weights) == 6, "Autodock vina needs 6 weights"
        # --weight_gauss1 1 --weight_gauss2 0 --weight_repulsion 0  --weight_hydrophobic 0 --weight_hydrogen 0 --weight_rot 0"
        lines.extend([
            f"weight_gauss1 = {weights[0]}",
            f"weight_gauss2 = {weights[1]}",
            f"weight_repulsion = {weights[2]}",
            f"weight_hydrophobic = {weights[3]}",
            f"weight_hydrogen = {weights[4]}",
            f"weight_rot = {weights[5]}",
        ])
    with open(config_fp, "w") as f:
        f.write("\n".join(lines))

class DockingScorer:
    '''This class is used to score the ligand using AutoDock Vina. It uses the prepared ligand and protein files to run Vina and extract the docking score.
    Example usage:
    scorer = DockingScorer(folder='path/to/folder', lig_file='ligand.pdbqt', prot_file='protein.pdbqt', save_out_file=True)
    docking_score = scorer.vina_score_ligand()
    print(f"Docking Score: {docking_score}")
    '''
    def __init__(self, 
                 folder, 
                 lig_file, 
                 prot_file,
                 get_vina_poses=False,
                 weights=None, 
                 save_out_file=True, 
                 protein_sequence=None, 
                 smiles=None, 
                 exp_binding_affinity=None,
                 from_pdb=True,
                 csv_out_file='docking_data_playground_from_pdb.csv',
                 receptor_DU=None):
        self.folder = folder
        self.lig_file = lig_file
        self.prot_file = prot_file
        self.lig_name = os.path.splitext(lig_file)[0]
        self.prot_name = os.path.splitext(prot_file)[0]
        self.pdb_id = self.prot_name.split('_')[0]
        self.get_vina_poses = get_vina_poses
        self.weights = weights
        self.save_out_file = save_out_file
        #TODO: more work required here to extract the protein sequence, ligand SMILES string and exp dG value, from an online db (low priority)
        self.protein_sequence = protein_sequence
        self.smiles = smiles
        self.exp_binding_affinity = exp_binding_affinity
        self.from_pdb = from_pdb
        self.csv_out_file = csv_out_file
        self.receptor_DU = receptor_DU


    #TODO: add a functionality to use Vina GPU instead (low priority)
    
    # ----------------- AutoDock Vina Scoring -----------------
    def run_vina(self):
        '''Run Vina with the prepared files and return the output.'''
        
        with set_directory(self.folder):
            
            # reference ligand in this case is the ligand file as they are already in the binding site
            ref_ligand = self.lig_file
            
            # get the COM from the reference ligand so we can define a docking grid
            ligand_COM = get_COM(ref_ligand)
            write_config_vina(f'{self.lig_name}.pdbqt', f'{self.prot_name}.pdbqt',ligand_COM,config_fp = f"{self.lig_name}_{self.prot_name}_config.txt", weights=None)
            cmd = f"./vina/vina_1.2.5_linux_x86_64 --config {self.lig_name}_{self.prot_name}_config.txt --score_only" #TODO: make it so that this can also output poses

            if self.get_vina_poses == True:
                cmd = f"./vina/vina_1.2.5_linux_x86_64 --config {self.lig_name}_{self.prot_name}_config.txt"

            try:
                out, err = "", ""  # Initialize out and err
                code, out, err = run_command(cmd, timeout=100)
                if code != 0:
                    raise CommandExecuteError(f"Command failed with return code {code}")
                
                return out
            except CommandExecuteError as e:
                print(f"Error in {self.pdb_id}: {e}")
                print("out: ", out)
                print("err: ", err)
                raise e

    def extract_vina_score(self, out):
        '''Extract the docking score from Vina output.'''
        if not self.get_vina_poses:
            strings = re.split('Estimated Free Energy of Binding   :', out)
            line = strings[1].split('\n')[0]
            energy = float(line.strip().split()[0])
            return energy
        else:
            # TODO: if get_vina_poses is True, return the output and the energy
            energy = []
            for line in out.split('\n'):
                match = re.search(r'^\s*\d+\s+(-\d+\.\d+)', line)
                if match:
                    energy.append(float(match.group(1)))
            return min(energy) if energy else None
    
    def save_output(self, out, energy):
        '''Save the vina output to a file if required.'''
        if self.save_out_file:
            if self.from_pdb:
                with open(f"{self.pdb_id}.out", 'w') as f:
                    # Convert list to string if 'out' is a list
                    if isinstance(out, list):
                        out = "\n".join(out)
                    f.write(out)
                print(f"Output saved as {self.pdb_id}.out\n")

                # add data to the dataframe: Assign the energy value to a new column 'computed_dG' for the filtered rows
                df.loc[df['pdb_id'] == self.pdb_id, 'computed_dG'] = energy
                updated_rows = df[df['pdb_id'] == self.pdb_id]

                if not os.path.isfile(f'{self.csv_out_file}'):
                    updated_rows.to_csv(f'{self.csv_out_file}', index=False)
                else: # else it exists so append without writing the header
                    updated_rows.to_csv(f'{self.csv_out_file}', mode='a', header=False, index=False)
                print(f"Data saved to {self.csv_out_file}")
                return updated_rows
            else:
                if isinstance(energy, list): 
                    for i, e in enumerate(energy):
                        with open(f"{self.prot_name}_{self.lig_name}_{i}.out", 'w') as f:
                            # Convert list to string if 'out' is a list
                            if isinstance(out, list):
                                out = "\n".join(out)
                            f.write(out)
                        print(f"Output saved as {self.prot_name}_{self.lig_name}_{i}.out\n")  
                        
                        lig_name_with_index = f'{self.lig_name}_{i}'

                        df = pd.DataFrame({'ligand_name': [lig_name_with_index], 
                                           'protein_name': [self.prot_name], 
                                           'computed_dG': [e],
                                           'error_message': [None]})
                                           
                        if not os.path.isfile(f'{self.csv_out_file}'):
                            df.to_csv(f'{self.csv_out_file}', index=False)
                        else: # else it exists so append without writing the header
                            df.to_csv(f'{self.csv_out_file}', mode='a', header=False, index=False)
                        print(f"Data saved to {self.csv_out_file}")

                else: 
                    with open(f"{self.prot_name}_{self.lig_name}.out", 'w') as f:
                        # Convert list to string if 'out' is a list
                        if isinstance(out, list):
                            out = "\n".join(out)
                        f.write(out)
                    print(f"Output saved as {self.prot_name}_{self.lig_name}.out\n")
                    
                    df = pd.DataFrame({'ligand_name': [self.lig_name], 
                                    'protein_name': [self.prot_name], 
                                    'computed_dG': [energy],
                                    'error_message': [None]})

                    if not os.path.isfile(f'{self.csv_out_file}'):
                        df.to_csv(f'{self.csv_out_file}', index=False)
                    else: # else it exists so append without writing the header
                        df.to_csv(f'{self.csv_out_file}', mode='a', header=False, index=False)
                    print(f"Data saved to {self.csv_out_file}")
                return df


            
    ######################## FINAL FUNCTIONS ########################
    def vina_score_ligand(self):
        """Main method to score the ligand using Vina."""
        if self.from_pdb:
            self.check_files()
        try:
            out = self.run_vina()
            energy = self.extract_vina_score(out)
        except Exception as e:
            print(f"Error in {self.pdb_id}: {e}")
            out = f"Error in {self.pdb_id}: {e}"
            energy = None
            raise e
        
        print(f"{self.pdb_id}: Estimated Free Energy of Binding = {energy} kcal/mol")
        self.save_output(out, energy)
        return energy

   

# download the PDB file using the PDB ID
# filename = download_pdb_file(pdb_id)
# # split the protein and ligand from the PDB file
# receptor_file, ligand_file = pdb_to_prot_lig(pdb_id, filename)
# this is for when lig_files and prot_files are lists...
# for lig, prot in tqdm(zip(lig_files, prot_files)):
#         vina_process_lig_prot(lig, prot)

def vina_process_lig_prot(lig_file, prot_file, preserve_water=False, preserve_metal=False, csv_out_file='vina_docking_data.csv'):
    '''Final function to process each protein and ligand pair, prepare them for docking, score the ligand and save the output files.
    Both lig_file, prot_file and complex_file must be strings, not lists. If you want to process multiple ligands and proteins (i.e. LISTS), you should use this function in the following manner: 

    for lig, prot in tqdm(zip(lig_files, prot_files)):
        vina_process_lig_prot(lig, prot, complex_file)
        ...

    The reference complex file needs to be a protein-ligand complex, with a docked ligand in the binding site of the protein of interest.
    '''
    # TODO: handling lists..?
    if isinstance(lig_file, list) or isinstance(prot_file, list):
        raise ValueError("lig_file and prot_file must be strings, not lists. If you are trying to process multiple ligands and proteins, you should use this function in a loop.")
    
    # check if any of the input is none
    if lig_file is None or prot_file is None:
        raise ValueError("lig_file, prot_file must not be None. You risk deleting everything in your current directory if any of these is None!")

    lig_name = os.path.splitext(lig_file)[0]
    prot_name = os.path.splitext(prot_file)[0]

    try:     
        # prepare the ligand and receptor for Vina i.e. convert to pdbqt format
        prep = DockingPrepper('.',
                                lig_file=lig_file, 
                                prot_file=prot_file,
                                preserve_water=preserve_water,
                                preserve_metal=preserve_metal) # this is not really the pdb_id but more of a basename for the files
        # for vina preparation
        prep.vina_process()
                
        # using the DockingScorer class to get the docking score
        scorer = DockingScorer('.', 
                                lig_file, 
                                prot_file,
                                get_vina_poses=True,
                                save_out_file=True,
                                from_pdb=False,
                                csv_out_file=csv_out_file,
                                receptor_DU=None) # None because this is for Vina

        docking_score = scorer.vina_score_ligand()

    except Exception as e:
        # ------- START OF CODES: save the error message to the .out file, None as the energy score, and write both the energy and the error message to a csv file -------
        # Set up logging
        logging.basicConfig(filename=f'{lig_name}_{prot_name}error_log.txt', level=logging.ERROR, format='%(asctime)s - %(levelname)s - %(message)s')
        logging.error(f"Error processing {lig_file}: {str(e)}")
        print(f"Skipping {lig_file} due to an error: {str(e)}")
        
        # saving the error message to the .txt file and None as the energy score 
        energy = None
        df = pd.DataFrame({'ligand_name': [lig_file], 
                                'protein_name': [prot_file], 
                                'computed_dG': [energy],
                                'error_message' : [str(e)]}
                                )
    
        if not os.path.isfile(csv_out_file):
            df.to_csv(csv_out_file, index=True)
        else: # else it exists so append without writing the header
            df.to_csv(csv_out_file, mode='a', header=False, index=False)   
        print(f"Error message saved to {csv_out_file}") 
        # -------END OF CODES: save the error message to the .out file, None as the energy score, and write both the energy and the error message to a csv file -------
        pass

        
    new_dir = os.path.join(cwd, f'{lig_name}_{prot_name}')
    os.makedirs(new_dir, exist_ok=True)
    
    # copy the files to the new directory, this will overwrite the files if they already exist
    os.system(f'mv -f *{lig_name}* {new_dir}')
    os.system(f'mv -f *{prot_name}* {new_dir}')

    # if the new directory already exists, forcibly remove it so we can overwrite it
    if os.path.exists(f'data/{new_dir}'):
        os.system(f'rm -rf data/{new_dir}')
    # move the output dir to the data storage dir 'data'
    os.system(f'mv -f {new_dir} ./data/')
    print(f"Data saved in ./data/{new_dir}")

def split_sdf(input_sdf):
    # Load the multi-molecule SDF file
    supplier = Chem.SDMolSupplier(input_sdf)
    output_name = input_sdf.split('.')[0]
    
    # Count the number of molecules
    molecules = [mol for mol in supplier if mol is not None]
    num_molecules = len(molecules)
    
    output = []
    # Check if the SDF file contains more than one molecule
    if num_molecules > 1:
        # Loop over each molecule and write it to a separate file
        for i, mol in enumerate(molecules):
            output_sdf = f"{output_name}_{i+1}.sdf"
            output.append(output_sdf)

            writer = SDWriter(output_sdf)
            writer.write(mol)
            writer.close()
        print(f"{num_molecules} molecules have been successfully separated into individual SDF files.")
    else:
        print(f"The SDF file {input_sdf} contains one or no molecules. No splitting is required.")
        output = str(input_sdf)
    return output

################################ VINA PIPELINE ################################
hif2a_sdf = "mcl1_syk_hif2a/hif2a_ligands.sdf"
lig_files = split_sdf(hif2a_sdf)

mcl1_sdf = "mcl1_syk_hif2a/mcl1_ligands.sdf"
lig_files = split_sdf(mcl1_sdf)

syk_sdf = "mcl1_syk_hif2a/syk_ligands.sdf"
lig_files = split_sdf(syk_sdf)

import glob
# hif2a
all_ligands = glob.glob('mcl1_syk_hif2a/hif2a_ligands*.sdf') 

for ligand_file in all_ligands:
    
    os.system(f'cp {ligand_file} .')
    # ligand_file name without path name
    ligand_file = os.path.basename(ligand_file)

    receptor_file = 'hif2a_protein.pdb'
    os.system(f'cp mcl1_syk_hif2a/{receptor_file} .') 

    vina_process_lig_prot(ligand_file, 
                        receptor_file, 
                        preserve_water=True,  # this is optional, default is False
                        preserve_metal=True, 
                        csv_out_file='hif2a_vina_docking_data.csv')

# mcl1
all_ligands = glob.glob('mcl1_syk_hif2a/mcl1_ligands*.sdf')

for ligand_file in all_ligands:
    
    os.system(f'cp {ligand_file} .')
    # ligand_file name without path name
    ligand_file = os.path.basename(ligand_file)

    receptor_file = 'mcl1_protein.pdb'
    os.system(f'cp mcl1_syk_hif2a/{receptor_file} .') 

    vina_process_lig_prot(ligand_file, 
                        receptor_file, 
                        preserve_water=True,  # this is optional, default is False
                        preserve_metal=True, 
                        csv_out_file='mcl1_vina_docking_data.csv')

# syk
all_ligands = glob.glob('mcl1_syk_hif2a/syk_ligands*.sdf')

for ligand_file in all_ligands:
    
    os.system(f'cp {ligand_file} .')
    # ligand_file name without path name
    ligand_file = os.path.basename(ligand_file)

    receptor_file = 'syk_protein.pdb'
    os.system(f'cp mcl1_syk_hif2a/{receptor_file} .') 

    vina_process_lig_prot(ligand_file, 
                        receptor_file, 
                        preserve_water=True,  # this is optional, default is False
                        preserve_metal=True, 
                        csv_out_file='syk_vina_docking_data.csv')
