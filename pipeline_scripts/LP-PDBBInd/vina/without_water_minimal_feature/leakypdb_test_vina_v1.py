import os
# make sure to set the OE_LICENSE environment variable, the full path should be included, or else openeye will kill your kernel!
os.chmod('/home/ian/oe_license.txt', 0o755)
os.environ['OE_LICENSE'] = '/home/ian/oe_license.txt'
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
# this python script runs on the basis that you already have the following installed or in the current working directory or installed in your system:
# 1. AutoDock Vina
# 2. ADFR Suite
# 3. OpenEye Toolkits
# 4. ASAP Discovery
# 5. The leaky_pdb_test.csv file containing the PDB IDs

# Set up the working directory
cwd = os.getcwd()

# giving permissions to run vina and ADFR scripts, change the directories of vina as required
os.chmod('./vina/vina_1.2.5_linux_x86_64', 0o755) # giving permissions to run vina scripts
os.chmod('./vina/vina_split_1.2.5_linux_x86_64', 0o755) # giving permissions to run vina scripts
os.chmod('./extract_ligand.sh', 0o755) # giving permissions to run extract_ligand.sh
os.chmod('./ADFRsuite/ADFRsuite_x86_64Linux_1.0', 0o755) # giving permissions to run ADFRsuite

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

# script from extract_ligand_oedu.py
def ExtractLigandFromDU(du, ofs):
# @ <SNIPPET-ExtractLigandFromDesignUnit>
    lig = oechem.OEGraphMol()
    if not du.GetLigand(lig):
        oechem.OEThrow.Fatal("Error: Could not extract ligand from the OEDesignUnit.")
    oechem.OEWriteMolecule(ofs,lig) 
# @ </SNIPPET-ExtractLigandFromDesignUnit>

    ofs.close()

# using Biopython to extract the receptor
class ReceptorSelect(Select):
    '''This class is used to select only the r`eceptor residues from a PDB file, by excluding any HETATM and any water.'''
    # I deleted the heavy metal ions as they are not a valid AutoDock type!
    # def __init__(self):
    #     # List of common heavy metal ions in protein structures
    #     self.heavy_metal_ions = ['ZN', 'FE', 'MG', 'CA', 'MN', 'CO', 'CU', 'NI', 'MO', 'W','YB']

    def accept_residue(self, residue):
        # Exclude water molecules, assumed to have residue name 'HOH' or 'WAT'
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
    
class ProteinLigandSplitter:
    def __init__(self, pdb_id):
        self.pdb_id = pdb_id
        self.pdbl = PDBList()
        self.parser = PDBParser()
        self.filename = None
        self.structure = None
        self.receptor_file = f"{self.pdb_id}_receptor.pdb"
        self.ligand_file = f"{self.pdb_id}_ligand.pdb"

    def download_pdb_file(self, pdb_id):
        self.filename = self.pdbl.retrieve_pdb_file(pdb_id, file_format='pdb', pdir='.', overwrite=True)
        self.structure = self.parser.get_structure('structure', self.filename)
    
    def extract_receptor(self):
        io = PDBIO()
        io.set_structure(self.structure)
        io.save(self.receptor_file, ReceptorSelect())

    def rename_downloaded_file(self):
        new_filename = f"{self.pdb_id}.pdb"
        if os.path.exists(self.filename) and not os.path.exists(new_filename):
            os.rename(self.filename, new_filename)
        self.filename = new_filename

    def extract_ligand(self):
        # debug line
        print(f'Extracting ligand from {self.filename}...')
        complex = read_pdb_file(self.filename)
        # du: OEDesignUnit
        # TODO: error trap for 'success'
        success, du = make_design_unit(complex)
        # use script from extract_ligand_oedu.py to extract ligand from complex
        oeoutfile = oechem.oemolostream(self.ligand_file)
        ExtractLigandFromDU(du, oeoutfile)

    def split_protein_ligand(self):
        self.download_pdb_file(self.pdb_id)
        self.extract_receptor()
        self.rename_downloaded_file()
        self.extract_ligand()
        print(f'The receptor file has been saved as {self.receptor_file}')
        print(f'The ligand file has been saved as {self.ligand_file}')
        return self.receptor_file, self.ligand_file

# # Example usage:
# splitter = ProteinLigandSplitter('1a2b')
# splitter.split_protein_ligand()

class VinaPrepper:
    '''This class is used to prepare the ligand and protein files for docking using AutoDock Vina. It uses the ADFRsuite to convert the ligand and protein files to pdbqt format.
    # Example usage:
    prep = VinaPrepper(folder='path/to/folder', lig_file='ligand.pdb', prot_file='protein.pdb')
    prep.process()
    '''
    def __init__(self, folder, lig_file, prot_file, lig_out_name=None, prot_out_name=None, add_h=True):
        self.folder = folder
        self.lig_file = lig_file
        self.prot_file = prot_file
        self.lig_out_name = lig_out_name or f"{lig_file.split('.')[0]}.pdbqt"
        self.prot_out_name = prot_out_name or f"{prot_file.split('.')[0]}.pdbqt"
        self.add_h = add_h

    def prepare_ligand(self):
        """Prepare the ligand file for Vina."""
        if self.add_h:
            cmd = f'./ADFRsuite/ADFRsuite_x86_64Linux_1.0/bin/prepare_ligand -l {self.lig_file} -o {self.lig_out_name} -A hydrogens'
        else:
            cmd = f'./ADFRsuite/ADFRsuite_x86_64Linux_1.0/bin/prepare_ligand -l {self.lig_file} -o {self.lig_out_name}'
        self._run_command(cmd, f'The ligand file has been saved as {self.lig_out_name}')

    def prepare_protein(self):
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
            raise
        except subprocess.TimeoutExpired:
            print(f"Command '{cmd}' timed out.")
            raise

    def process(self):
        """Process the ligand and protein files."""
        with set_directory(self.folder):
            self.prepare_ligand()
            self.prepare_protein()

class VinaScorer:
    '''This class is used to score the ligand using AutoDock Vina. It uses the prepared ligand and protein files to run Vina and extract the docking score.
    # Example usage:
    scorer = VinaScorer(folder='path/to/folder', lig_file='ligand.pdbqt', prot_file='protein.pdbqt', save_out_file=True)
    docking_score = scorer.score_ligand()
    print(f"Docking Score: {docking_score}")
    '''
    def __init__(self, folder, lig_file, prot_file, weights=None, save_out_file=True, protein_sequence=None, smiles=None, exp_binding_affinity=None):
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

    def prepare_files(self):
        '''Check and prepare the necessary files if needed'''
        if not os.path.isfile(f"{self.folder}/{self.prot_name}.pdbqt") or not os.path.isfile(f"{self.folder}/{self.lig_name}.pdbqt"):
            splitter = ProteinLigandSplitter(self.pdb_id)
            self.prot_file, self.lig_file = splitter.split_protein_ligand()
    #TODO: add a functionality to use Vina GPU instead (low priority)
    #TODO: add a functionality to use OpenEye instead (high priority)
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
                raise

    def get_docking_score(self, out):
        '''Extract the docking score from Vina output.'''
        strings = re.split('Estimated Free Energy of Binding   :', out)
        line = strings[1].split('\n')[0]
        energy = float(line.strip().split()[0])
        return energy
    
    # make the output into a pandas dataframe
    # add pdb_id, protein sequence, ligand SMILES string, dG, exp_dG as a column in the dataframe
    def extract_data(self, df):
        '''Extract the protein sequence, ligand SMILES string, and experimental binding affinity from the specified dataframe.'''

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
            with open(f"{self.pdb_id}.out", 'w') as f:
                f.write(out)
            print(f"Output saved as {self.pdb_id}.out\n")

            df = self.extract_data(leakypdb) # define the dataframe to extract data from
            # add data to the dataframe
            # Assign the energy value to a new column 'computed_dG' for the filtered rows
            df.loc[df['pdb_id'] == self.pdb_id, 'computed_dG'] = energy
            updated_rows = df[df['pdb_id'] == self.pdb_id]

            if not os.path.isfile('leakypdb_test_AutoDock_data.csv'):
                updated_rows.to_csv('leakypdb_test_AutoDock_data.csv', index=False)
            else: # else it exists so append without writing the header
                updated_rows.to_csv('leakypdb_test_AutoDock_data.csv', mode='a', header=False, index=False)
            print(f"Data saved to leakypdb_test_AutoDock_data.csv")
            return updated_rows

    def score_ligand(self):
        """Main method to score the ligand using Vina."""
        self.prepare_files()
        try:
            out = self.run_vina()
            energy = self.get_docking_score(out)
        except Exception as e:
            print(f"Error in {self.pdb_id}: {e}")
            energy = None
        print(f"{self.pdb_id}: Estimated Free Energy of Binding = {energy} kcal/mol")
        self.save_output(out, energy)
        return energy

def process_pdb(pdb_ids):
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
            # split the protein and ligand from the PDB file
            splitter = ProteinLigandSplitter(pdb_id)
            receptor_file, ligand_file = splitter.split_protein_ligand()
            
            # prepare the ligand and receptor for Vina i.e. convert to pdbqt format
            prep = VinaPrepper('.', ligand_file, receptor_file)
            prep.process()
            
            # using the VinaScorer class to get the docking score
            scorer = VinaScorer('.', ligand_file, receptor_file, save_out_file=True)
            docking_score = scorer.score_ligand()
            
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
            if not os.path.isfile('leakypdb_test_AutoDock_data.csv'):
                updated_rows.to_csv('leakypdb_test_AutoDock_data.csv', index=False)
            else: # else it exists so append without writing the header
                updated_rows.to_csv('leakypdb_test_AutoDock_data.csv', mode='a', header=False, index=False)
            
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

# loading LeakyPDB dataset
leakypdb = pd.read_csv('leaky_pdb_test.csv')
leakypdb.rename(columns={'Unnamed: 0': 'pdb_id',
                          'header': 'protein_family', 
                          'seq':'protein_sequence',
                          'kd/ki':'binding_affinity'}, inplace=True)
leakypdb['pdb_id'] = leakypdb['pdb_id'].str.upper()
#TODO: kd/ki column may not be all experimentally determined, some may be predicted, need to check this
leakypdb_ids=leakypdb['pdb_id'].tolist()
leakypdb_ids = [leakypdb_id.upper() for leakypdb_id in leakypdb_ids]

# make sure to set the OE_LICENSE environment variable, the full path should be included, or else openeye will kill your kernel!
# os.environ['OE_LICENSE'] = '/home/ian/oe_license.txt' # change this to your OE_LICENSE path
# os.chmod('/home/ian/oe_license.txt', 0o755) # change this to your OE_LICENSE path
# test the pipeline for the LeakyPDB dataset
process_pdb(leakypdb_ids)
#TODO: fix issues with not getting peptide ligands (low priority)