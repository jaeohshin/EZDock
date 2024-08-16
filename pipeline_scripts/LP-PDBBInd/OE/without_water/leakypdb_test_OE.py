import os
# giving permissions to run vina and ADFR scripts, change the directories of vina as required
os.chmod('./vina/vina_1.2.5_linux_x86_64', 0o755) # giving permissions to run vina scripts
os.chmod('./vina/vina_split_1.2.5_linux_x86_64', 0o755) # giving permissions to run vina scripts
os.chmod('./extract_ligand.sh', 0o755) # giving permissions to run extract_ligand.sh
os.chmod('./ADFRsuite/ADFRsuite_x86_64Linux_1.0', 0o755) # giving permissions to run ADFRsuite

# make sure to set the OE_LICENSE environment variable, the full path should be included, or else openeye will kill your kernel!
os.environ['OE_LICENSE'] = '/home/ian/oe_license.txt'
os.chmod('/home/ian/oe_license.txt', 0o755)

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
from asapdiscovery.data.backend.openeye import (
    oechem,
    oedocking,
    oegrid,
    oespruce,
    openeye_perceive_residues,
)
from asapdiscovery.modeling.schema import MoleculeComponent, MoleculeFilter
from asapdiscovery.modeling.modeling import split_openeye_mol, make_design_unit
import matplotlib.pyplot as plt
import seaborn as sns
import openeye.oechem as oechem

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
    
# from https://github.com/choderalab/perses/blob/main/examples/moonshot-mainseries/00-prep-receptor.py#L19-L35
def read_pdb_file(pdb_file):
    print(f'Reading receptor from {pdb_file}...')

    from openeye import oechem
    ifs = oechem.oemolistream()
    #ifs.SetFlavor(oechem.OEFormat_PDB, oechem.OEIFlavor_PDB_Default | oechem.OEIFlavor_PDB_DATA | oechem.OEIFlavor_PDB_ALTLOC)  # Causes extra protons on VAL73 for x1425
    ifs.SetFlavor(oechem.OEFormat_PDB, oechem.OEIFlavor_PDB_Default | oechem.OEIFlavor_PDB_DATA )

    if not ifs.open(pdb_file):
        oechem.OEThrow.Fatal("Unable to open %s for reading." % pdb_file)

    mol = oechem.OEGraphMol()
    if not oechem.OEReadMolecule(ifs, mol):
        oechem.OEThrow.Fatal("Unable to read molecule from %s." % pdb_file)
    ifs.close()

    return (mol)

# script from extract_ligand_oedu.py: https://docs.eyesopen.com/toolkits/python/oechemtk/oebio_examples/oebio_example_extract_ligand.html#section-example-oebio-extract-ligand
def ExtractLigandFromDU(du, ofs):
# @ <SNIPPET-ExtractLigandFromDesignUnit>
    lig = oechem.OEGraphMol()
    if not du.GetLigand(lig):
        oechem.OEThrow.Fatal("Error: Could not extract ligand from the OEDesignUnit.")
    oechem.OEWriteMolecule(ofs,lig) 
# @ </SNIPPET-ExtractLigandFromDesignUnit>

    ofs.close()

# Function to read a molecular complex/ligand/protein from a file
def oe_read_molecule(filename):
    ifs = oechem.oemolistream()
    if not ifs.open(filename):
        oechem.OEThrow.Fatal("Unable to open file %s" % filename)
    mol = oechem.OEGraphMol()
    oechem.OEReadMolecule(ifs, mol)
    ifs.close()
    return mol

def oe_split_complex(input_filename, ligand_output_filename, protein_output_filename):
    # Read the complex
    complex_molecule = oe_read_molecule(input_filename)

    # Split the complex
    protein = oechem.OEGraphMol()
    ligand = oechem.OEGraphMol()
    water = oechem.OEGraphMol()
    other = oechem.OEGraphMol()

    # Split the complex
    oechem.OESplitMolComplex(ligand, protein, water, other, complex_molecule)

    # Write the components to files
    oechem.OEWriteMolecule(oechem.oemolostream(protein_output_filename), protein)
    oechem.OEWriteMolecule(oechem.oemolostream(ligand_output_filename), ligand)
    
    # we dont need the water and other components
    # oechem.OEWriteMolecule(oechem.oemolostream(output_basename + "_water.pdb"), water)
    # oechem.OEWriteMolecule(oechem.oemolostream(output_basename + "_other.pdb"), other)

# Example usage
# input_filename = "./Mpro-x0072_0.pdb"  # Replace with your file path
# ligand_output_filename = "./Mpro-x0072_0_ligand.pdb"
# protein_output_filename = "./Mpro-x0072_0_protein.pdb"
# oe_split_complex(input_filename, ligand_output_filename, protein_output_filename)

# using Biopython to extract the receptor
class ReceptorSelect(Select):
    '''This class is used to select only the receptor residues from a PDB file, by excluding any HETATM and any water.'''
    # I deleted the heavy metal ions as they are not a valid AutoDock type!
    # def __init__(self):
    #     # List of common heavy metal ions in protein structures
    #     self.heavy_metal_ions = ['ZN', 'FE', 'MG', 'CA', 'MN', 'CO', 'CU', 'NI', 'MO', 'W','YB']

    def accept_residue(self, residue):
        # Exclude water molecules, assumed to have residue name 'HOH' or 'WAT'
        #TODO: What if there's water in the binding site?
        if residue.get_resname() in ['HOH', 'WAT']:
            return False
        
        # Include heavy metal ions
        # if residue.get_resname() in self.heavy_metal_ions:
        #     return True
        
        # Check if any atom in the residue is from an ATOM record
        for atom in residue:
            '''in biopython, atom.get_full_id()[3] is accessing the fourth element of the full ID. This element represents the atom name, atom.get_full_id()[3][0] is accessing the first character of the atom name. This character represents the element symbol of the atom. The condition atom.get_full_id()[3][0] == ' ' checks whether the first character of the atom name is a space. If it is a space, then the atom is from an ATOM record, otherwise it is from a HETATM record.'''
            if atom.get_full_id()[3][0] == ' ':
                return True
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

def save_receptor(structure, receptor_file):
    """Save the receptor part of the PDB file."""
    io = PDBIO()
    io.set_structure(structure)
    io.save(receptor_file, ReceptorSelect()) 

def rename_file(old_filename, new_filename):
    """Rename the downloaded file to a standard name."""
    if os.path.exists(old_filename) and not os.path.exists(new_filename):
        os.rename(old_filename, new_filename)
    return new_filename

def extract_ligand(filename, ligand_file):
    """Extract the ligand from the complex."""
    print(f'Extracting ligand from {filename}...')
    complex = read_pdb_file(filename)
    success, du = make_design_unit(complex)
    
    if success:
        oeoutfile = oechem.oemolostream(ligand_file)
        ExtractLigandFromDU(du, oeoutfile)
    else:
        print(f'Failed to extract ligand from {filename} using the OpenEye Toolkit,' 
              'trying a different script written using OESplitMolComplex...')
        # Extract ligand using shell script
        complex_id = filename.rsplit('.')[0]

        # script written using OESplitMolComplex:
        oe_split_complex(input_filename = filename, 
                         ligand_output_filename = f"{complex_id}_ligand.pdb", protein_output_filename = f"{complex_id}_receptor.pdb")
        
        # this is my shell script - not guaranteed to work for everything and hence not used!
        # run_command(f'./extract_ligand.sh {complex_id}') 

def pdb_to_prot_lig(pdb_id, filename):
    """Main function to handle the splitting of protein and ligand."""
    structure = get_structure(filename)

    # naming the output files
    receptor_file = f"{pdb_id}_receptor.pdb"
    ligand_file = f"{pdb_id}_ligand.pdb"
    filename = rename_file(filename, f"{pdb_id}.pdb")

    # ---old codes (uses Biopython for receptor extraction, and ASAP for ligand extraction - slower)---
    save_receptor(structure, receptor_file)
    extract_ligand(filename, ligand_file)

    # ---new codes (uses OE for ligand and receptor extraction)---
    # split the complex into protein and ligand
    # oe_split_complex(input_filename = filename, ligand_output_filename = ligand_file, protein_output_filename = receptor_file)
    
    print(f'The receptor file has been saved as {receptor_file}')
    print(f'The ligand file has been saved as {ligand_file}')
    return receptor_file, ligand_file

# Example usage:
# filename = download_pdb_file(pdb_id) # if starting from a PDB ID
# receptor_file, ligand_file = pdb_to_prot_lig(pdb_id, filename)

class DockingPrepper:
    '''This class is used to prepare the ligand and protein files for docking using AutoDock Vina. It uses the ADFRsuite to convert the ligand and protein files to pdbqt format.
    # Example usage:
    prep = DockingPrepper(folder='path/to/folder', pdb_id = 'pdb_id.pdb', lig_file='ligand.pdb', prot_file='protein.pdb')
    
    # ---- AutoDock Vina Preparation ----
    prep.vina_process()
    
    # Alternatively...
    # ---- OpenEye Preparation ----
    prep.oe_process()
    '''
    def __init__(self, folder, pdb_id, lig_file, prot_file, lig_out_name=None, prot_out_name=None, add_h=True):
        self.folder = folder
        self.pdb_id = pdb_id
        self.lig_file = lig_file
        self.prot_file = prot_file
        self.lig_out_name = lig_out_name or f"{lig_file.split('.')[0]}.pdbqt"
        self.prot_out_name = prot_out_name or f"{prot_file.split('.')[0]}.pdbqt"
        self.add_h = add_h

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
    
    # ----------------- OpenEye Preparation -----------------
    def oe_make_design_unit(self):
        ''' 
        This function makes design units from a protein-ligand complex.
        
        Input: input_file - the protein-ligand complex file
        Saved to file: {input_file_basename}_DU_{i}.oedu 
        Output: output_files - a list of the design unit files
        '''
        input_file = f"{self.pdb_id}.pdb"
    
        os.system(f"python ./OpenEye/make_design_units.py {input_file}")
        
        # this is the output pattern for the design unit files
        output_pattern = os.path.basename(input_file)[:-4] + "_DU_{}.oedu"
        # List the files in the current directory and filter for the expected output pattern
        output_files = []
        i = 0
        while True:
            output_file = output_pattern.format(i)
            if os.path.exists(output_file):
                output_files.append(output_file)
                i += 1
            else:
                break

        # output_files now contains the names of the output files
        print(f"Design unit was successfully made for {input_file}, output is saved to {output_files}.")
        return output_files

    def oe_make_receptor(self, input_file):

        oe_output_files = [] 
        for ifile in input_file:
            output_basename = os.path.basename(ifile)[:-5]
            os.system(f"python ./OpenEye/MakeReceptor.py -in {ifile} -out {output_basename}_receptor.oedu") # this will output {ifile}_receptor.oedu

            oe_output_files.append(f"{output_basename}_receptor.oedu")
            print(f'The receptor design unit file has been saved as {output_basename}_receptor.oedu')

        return oe_output_files 
    
    # Example usage:
    # input_file = "6EQ2.pdb"
    # complex_DU = oe_make_design_unit(input_file) # this will output ['6EQ2_DU_0.oedu']
    # receptor_DU = oe_make_receptor(complex_DU) # this will output ['6EQ2_DU_0_receptor.oedu']

    def oe_process(self):
        '''Process the ligand and protein files using OpenEye.'''
        print(f'Preparing the ligand and protein files using OpenEye Toolkits...')
        with set_directory(self.folder):
            complex_DU = self.oe_make_design_unit()
            oe_output_files = self.oe_make_receptor(input_file=complex_DU)
        
        # self.oe_output_files = [f"{self.pdb_id}_DU_0_receptor.oedu"]
        return oe_output_files
    
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
        self.lig_name = lig_file.split('.')[0]
        self.prot_name = prot_file.split('.')[0]
        self.pdb_id = self.prot_name.split('_')[0]
        self.weights = weights
        self.save_out_file = save_out_file
        #TODO: more work required here to extract the protein sequence, ligand SMILES string and exp dG value, from an online db (low priority)
        self.protein_sequence = protein_sequence
        self.smiles = smiles
        self.exp_binding_affinity = exp_binding_affinity
        self.from_pdb = from_pdb
        self.csv_out_file = csv_out_file
        self.receptor_DU = receptor_DU

    def check_files(self):
        '''Check and download the necessary files if needed'''
        if not os.path.isfile(f"{self.folder}/{self.prot_name}.pdbqt") or not os.path.isfile(f"{self.folder}/{self.lig_name}.pdbqt"):
            filename = download_pdb_file(self.pdb_id)
            self.prot_file, self.lig_file = pdb_to_prot_lig(self.pdb_id, filename)
    #TODO: add a functionality to use Vina GPU instead (low priority)
    
    # ----------------- AutoDock Vina Scoring -----------------
    def run_vina(self):
        '''Run Vina with the prepared files and return the output.'''
        cmd = f"./vina/vina_1.2.5_linux_x86_64 --receptor {self.prot_name}.pdbqt --ligand {self.lig_name}.pdbqt --autobox --score_only"
        with set_directory(self.folder):
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
        strings = re.split('Estimated Free Energy of Binding   :', out)
        line = strings[1].split('\n')[0]
        energy = float(line.strip().split()[0])
        return energy
    
    # ----------------- OpenEye Toolkits Scoring -----------------
    def oe_clean_then_dock(self):
        ''' This function docks a ligand to a receptor using OpenEye's CleanThenDockMolecules.py script.
        Input: lig_file - the ligand file to be docked
            receptor_DU - the receptor file to dock the ligand to
        Output: output_files - a list of the docked ligand files, which also contains the chemgauss4 score.'''
        
        output_files = [] 

        if self.receptor_DU is None:
            print("No receptor design unit file provided. Please provide a receptor file to dock the ligand.")
            raise ValueError("No receptor design unit file provided.")
        for receptor_du in self.receptor_DU:   
            output_basename = receptor_du.split('.')[0].replace('_receptor', '_ligand')
            os.system(f'python ./OpenEye/CleanThenDockMolecules.py -in {self.lig_file} -out {output_basename}_docked.sdf -receptor {receptor_du}') # output score is Chemgauss4, contained in the sdf file
            
            output_files.append(f"{output_basename}_docked.sdf")
        
        print(f'The docked ligand file has been saved as {output_files}')
        return output_files    
    
    # TODO: modify this function to calculate and extract the other scores from the docked ligands
    def extract_chemgauss4_scores(self, docked_ligand):
        ''' This function extracts the Chemgauss4 scores from the docked ligands.'''
        # docked_ligand should not be a list here.              
        if isinstance(docked_ligand, list):
            docked_ligand = docked_ligand[0]
        # Open the SDF file and create an SDMolSupplier object
        supplier = Chem.SDMolSupplier(docked_ligand)

        # energies = []
        # Iterate over each molecule in the SDF file
        for mol in supplier:
            if mol is None:
                continue
            # Extract the Chemgauss4 score if it exists
            if mol.HasProp('Chemgauss4'):
                energy = mol.GetProp('Chemgauss4')
                energy = float(energy)
                # energies.append(energy)
        return energy

    # # Example usage
    # docked_ligand = oe_clean_then_dock(ligand_file, receptor_DU) # this will output ['6EQ2_ligand_docked.sdf'], suppose ligand_file = '6EQ2_ligand.pdb'
    # Extracting Chemgauss4 scores from the docked ligands
    # scores = extract_chemgauss4_scores(docked_ligand)

    # ----------------- TODO: Putting things together -----------------

    
    # make the output into a pandas dataframe
    # add pdb_id, protein sequence, ligand SMILES string, dG, exp_dG as a column in the dataframe
    def extract_data_from_leakypdb(self, df):
        '''Extract the protein sequence, ligand SMILES string, and experimental binding affinity from the specified dataframe (leakypdb_test.csv).'''

        # # Read the DataFrame (must be done so that df is correctly recognised as a DataFrame, not a string.)    
        # df = pd.read_csv(df)

        # Filter the DataFrame for the specific PDB ID
        filtered_df = df[df['pdb_id'] == self.pdb_id] 

        # Check if the specific columns exist in the DataFrame
        required_columns = ['pdb_id','smiles', 'protein_sequence', 'binding_affinity']
        missing_columns = [col for col in required_columns if col not in df.columns]
        if missing_columns:
            return f"Missing columns in the data: {', '.join(missing_columns)}"
        # Extract the needed information
        result = filtered_df[required_columns]
        return result    

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

                df = self.extract_data_from_leakypdb(leakypdb) # define the dataframe to extract data from

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

    def oe_score_ligand(self):
        '''Main method to score the ligand using OpenEye Toolkits.'''
        # if self.from_pdb:
        #     self.check_files()
        try:
            docked_ligand = self.oe_clean_then_dock()
            out = docked_ligand
            energies = []
            if len(docked_ligand) > 1:
                for i, lig in enumerate(docked_ligand):
                    energy = self.extract_chemgauss4_scores(lig)
                    if energy is not None:
                        print(f"Docked ligand {lig}, molecule {i+1}: Chemgauss4 score = {energy:.2f}")
                        
                        self.save_output(out, energy)
                        energies.append(energy)
            else:
                energy = self.extract_chemgauss4_scores(docked_ligand)
                if energy is not None:
                    print(f"Docked ligand {docked_ligand}: Chemgauss4 score = {energy:.2f}")
                    self.save_output(out, energy)
                    energies = energy
            
            # debug line
            # print(energies)
        except Exception as e:
            print(f"Error in {self.pdb_id}: {e}")
            out = f"Error in {self.pdb_id}: {e}"
            energies = None
            print(f"{self.pdb_id}: No valid energy value available.")
            # self.save_output(out, energies)
            raise e
    
        return energies

def vina_process_pdb(pdb_ids, csv_out_file='vina_docking_data.csv'):
    '''Final function to process each PDB ID'''
    # capitalise the PDB IDs
    pdb_ids = [pdb_id.upper() for pdb_id in pdb_ids]
    if len(pdb_ids) == 1:
        destination_dir = f'data/{pdb_ids[0]}'
    else:
        destination_dir = f'data/{pdb_ids[0]}_to_{pdb_ids[-1]}'
    os.makedirs(destination_dir, exist_ok=True)
    
    for pdb_id in tqdm(pdb_ids):
        try:
            os.environ['OE_LICENSE'] = '/home/ian/oe_license.txt' # change this to your OE_LICENSE path
            filename = download_pdb_file(pdb_id)
            # split the protein and ligand from the PDB file
            receptor_file, ligand_file = pdb_to_prot_lig(pdb_id, filename)
            
            # prepare the ligand and receptor for Vina i.e. convert to pdbqt format
            prep = DockingPrepper('.',
                                    pdb_id,
                                   ligand_file, 
                                   receptor_file)
            # vina preparation
            prep.vina_process()
            
            # using the DockingScorer class to get the docking score
            scorer = DockingScorer('.', 
                                ligand_file, 
                                receptor_file, 
                                save_out_file=True,
                                from_pdb=True,
                                csv_out_file=csv_out_file,
                                receptor_DU=None) # None because this is for Vina, but for OE, this will be the receptor DU (required)
                                
            docking_score = scorer.vina_score_ligand()
            
            # make a new dir for each PDB ID
            new_dir = os.path.join(cwd, pdb_id)
            os.makedirs(new_dir, exist_ok=True)

            # move the files to the new directory, this will overwrite the files if they already exist
            subprocess.run(f'mv -f *{pdb_id}* {new_dir}', shell=True)
            # if the new directory already exists, forcibly remove it so we can overwrite it
            if os.path.exists(f'{destination_dir}/{new_dir}'):
                subprocess.run(f'rm -rf {destination_dir}/{new_dir}', shell=True)
            # move the output dir to the data storage dir 'data'
            subprocess.run(f'mv -f {new_dir} {destination_dir}', shell=True)
        except Exception as e:
            # ------- START OF CODES: save the error message to the .out filE, None as the energy score, and write both the energy and the error message to a csv file -------
            # Set up logging
            logging.basicConfig(filename=f'{pdb_id}_error_log.txt', level=logging.ERROR, format='%(asctime)s - %(levelname)s - %(message)s')
            logging.error(f"Error processing {pdb_id}: {str(e)}")
            print(f"Skipping {pdb_id} due to an error: {str(e)}")
        
            ## saving the error message to the .out file and None as the energy score 
            df = leakypdb
            # Filter the DataFrame for the specific PDB ID
            filtered_df = df[df['pdb_id'] == pdb_id] 

            # Check if the specific columns exist in the DataFrame
            required_columns = ['pdb_id','smiles', 'protein_sequence', 'binding_affinity']
            missing_columns = [col for col in required_columns if col not in df.columns]
            if missing_columns:
                return f"Missing columns in the data: {', '.join(missing_columns)}"
            # Extract the needed information
            df = filtered_df[required_columns]
            # set the energy to None
            df.loc[df['pdb_id'] == pdb_id, 'computed_dG'] = ''
            # add a new column for the error message
            df.loc[df['pdb_id'] == pdb_id, 'error_message'] = str(e)
            updated_rows = df[df['pdb_id'] == pdb_id]

            if not os.path.isfile(csv_out_file):
                updated_rows.to_csv(csv_out_file, index=True)
            else: # else it exists so append without writing the header
                updated_rows.to_csv(csv_out_file, mode='a', header=False, index=False)
            print(f"Error message saved to {csv_out_file}") 
            # -------END OF CODES: save the error message to the .out filE, None as the energy score, and write both the energy and the error message to a csv file -------

            # make a new dir for each PDB ID
            new_dir = os.path.join(cwd, pdb_id)
            os.makedirs(new_dir, exist_ok=True)
            # move the files to the new directory, this will overwrite the files if they already exist
            subprocess.run(f'mv -f *{pdb_id}* {new_dir}', shell=True)
            # if the new directory already exists, forcibly remove it so we can overwrite it
            if os.path.exists(f'{destination_dir}/{new_dir}'):
                subprocess.run(f'rm -rf {destination_dir}/{new_dir}', shell=True)
            # move the output dir to the data storage dir 'data'
            subprocess.run(f'mv -f {new_dir} {destination_dir}', shell=True)
            continue  # This will skip the current pdb_id and proceed with the next one

    print(f"Data saved in {destination_dir}")

def oe_process_pdb(pdb_ids, csv_out_file='docking_data.csv'):
    '''Final function to process each PDB ID'''
    # capitalise the PDB IDs
    pdb_ids = [pdb_id.upper() for pdb_id in pdb_ids]
    if len(pdb_ids) == 1:
        destination_dir = f'data/{pdb_ids[0]}'
    else:
        destination_dir = f'data/{pdb_ids[0]}_to_{pdb_ids[-1]}'
    os.makedirs(destination_dir, exist_ok=True)

    for pdb_id in tqdm(pdb_ids):
        try:
            os.environ['OE_LICENSE'] = '/home/ian/oe_license.txt' # change this to your OE_LICENSE path
            filename = download_pdb_file(pdb_id)
            # split the protein and ligand from the PDB file
            receptor_file, ligand_file = pdb_to_prot_lig(pdb_id, filename)
            
            # prepare the ligand and receptor for Vina i.e. convert to pdbqt format
            prep = DockingPrepper('.',
                                    pdb_id,
                                    ligand_file, 
                                    receptor_file)
            # vina preparation
            DU = prep.oe_process()
            
            # using the DockingScorer class to get the docking score
            scorer = DockingScorer('.', 
                                ligand_file, 
                                receptor_file, 
                                save_out_file=True,
                                from_pdb=True,
                                csv_out_file=csv_out_file,
                                receptor_DU=DU) # for OE, this will be the receptor DU produced from DockingPrepper.oe_process (required)
                                
            docking_score = scorer.oe_score_ligand()
            
            # make a new dir for each PDB ID
            new_dir = os.path.join(cwd, pdb_id)
            os.makedirs(new_dir, exist_ok=True)

            # move the files to the new directory, this will overwrite the files if they already exist
            subprocess.run(f'mv -f *{pdb_id}* {new_dir}', shell=True)

            # if the new directory already exists, forcibly remove it so we can overwrite it
            if os.path.exists(f'{destination_dir}/{new_dir}'):
                subprocess.run(f'rm -rf {destination_dir}/{new_dir}', shell=True)
            # move the output dir to the data storage dir 'data'
            subprocess.run(f'mv -f {new_dir} {destination_dir}', shell=True)
        except Exception as e:
            # ------- START OF CODES: save the error message to the .out filE, None as the energy score, and write both the energy and the error message to a csv file -------
            # Set up logging
            logging.basicConfig(filename=f'{pdb_id}_error_log.txt', level=logging.ERROR, format='%(asctime)s - %(levelname)s - %(message)s')
            logging.error(f"Error processing {pdb_id}: {str(e)}")
            print(f"Skipping {pdb_id} due to an error: {str(e)}")
        
            ## saving the error message to the .out file and None as the energy score 
            df = leakypdb
            # Filter the DataFrame for the specific PDB ID
            filtered_df = df[df['pdb_id'] == pdb_id] 

            # Check if the specific columns exist in the DataFrame
            required_columns = ['pdb_id','smiles', 'protein_sequence', 'binding_affinity']
            missing_columns = [col for col in required_columns if col not in df.columns]
            if missing_columns:
                return f"Missing columns in the data: {', '.join(missing_columns)}"
            # Extract the needed information
            df = filtered_df[required_columns]
            # set the energy to None
            df.loc[df['pdb_id'] == pdb_id, 'computed_dG'] = ''
            # add a new column for the error message
            df.loc[df['pdb_id'] == pdb_id, 'error_message'] = str(e)
            updated_rows = df[df['pdb_id'] == pdb_id]

            if not os.path.isfile(csv_out_file):
                updated_rows.to_csv(csv_out_file, index=True)
            else: # else it exists so append without writing the header
                updated_rows.to_csv(csv_out_file, mode='a', header=False, index=False)
            # -------END OF CODES: save the error message to the .out filE, None as the energy score, and write both the energy and the error message to a csv file -------

            # make a new dir for each PDB ID
            new_dir = os.path.join(cwd, pdb_id)
            os.makedirs(new_dir, exist_ok=True)
            # move the files to the new directory, this will overwrite the files if they already exist
            subprocess.run(f'mv -f *{pdb_id}* {new_dir}', shell=True)
            # if the new directory already exists, forcibly remove it so we can overwrite it
            if os.path.exists(f'{destination_dir}/{new_dir}'):
                subprocess.run(f'rm -rf {destination_dir}/{new_dir}', shell=True)
            # move the output dir to the data storage dir 'data'
            subprocess.run(f'mv -f {new_dir} {destination_dir}', shell=True)
            continue  # This will skip the current pdb_id and proceed with the next one

    print(f"Data saved in {destination_dir}")

    
# loading LeakyPDB test dataset
leakypdb = pd.read_csv('leakypdb_test.csv')
leakypdb.rename(columns={'Unnamed: 0': 'pdb_id',
                          'header': 'protein_family', 
                          'seq':'protein_sequence',
                          'kd/ki':'binding_affinity',
                          'value':'pkd/pki'}, inplace=True)
leakypdb['pdb_id'] = leakypdb['pdb_id'].str.upper()
leakypdb_ids=leakypdb['pdb_id'].tolist()
leakypdb_ids = [leakypdb_id.upper() for leakypdb_id in leakypdb_ids]

oe_process_pdb(pdb_ids = leakypdb_ids, csv_out_file='leakypdb_test_OE_data.csv')
