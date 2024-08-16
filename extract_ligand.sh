#!/bin/bash

# error trap: exit if the script was not provided with the required argument, exit with a status code of 1
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <pdb_id>"
    exit 1
fi

# the shell script takes one argument: the pdb_id of the PDB file
pdb_id=$1
pdb_file="${pdb_id}.pdb"
output_file="${pdb_id}_ligand.pdb"

# error trap: check if the PDB file exists
if [ ! -f "$pdb_file" ]; then
    echo "PDB file ${pdb_file} not found!"
    exit 1
fi

# # Define standard residues, water and solvents
standard_residues="ALA|CYS|ASP|GLU|PHE|GLY|HIS|ILE|LYS|LEU|MET|ASN|PRO|GLN|ARG|SER|THR|VAL|TRP|TYR|DA|DT|DG|DC|A|T|G|C|U|HOH|WAT|DOD|K|NA|CL|MG|CA|SO4|PO4|EDO|ACT|PEG|MPD|GOL|PG4|PE4|MES|BME|TRS|CIT|HEP|EPE|ECS|TAR|IMD|DMS|ZN|FE|MG|CA|MN|CO|CU|NI|MO|W|YB"

# Define standard residues, water and solvents with word boundaries
# standard_residues="\\bALA\\b|\\bCYS\\b|\\bASP\\b|\\bGLU\\b|\\bPHE\\b|\\bGLY\\b|\\bHIS\\b|\\bILE\\b|\\bLYS\\b|\\bLEU\\b|\\bMET\\b|\\bASN\\b|\\bPRO\\b|\\bGLN\\b|\\bARG\\b|\\bSER\\b|\\bTHR\\b|\\bVAL\\b|\\bTRP\\b|\\bTYR\\b|\\bDA\\b|\\bDT\\b|\\bDG\\b|\\bDC\\b|\\bA\\b|\\bT\\b|\\bG\\b|\\bC\\b|\\bU\\b|\\bHOH\\b|\\bWAT\\b|\\bDOD\\b|\\bK\\b|\\bNA\\b|\\bMG\\b|\\bCA\\b|\\bEDO\\b|\\bACT\\b|\\bPEG\\b|\\bMPD\\b|\\bGOL\\b|\\bPG4\\b|\\bPE4\\b|\\bMES\\b|\\bBME\\b|\\bTRS\\b|\\bCIT\\b|\\bHEP\\b|\\bEPE\\b|\\bECS\\b|\\bTAR\\b|\\bIMD\\b|\\bZN\\b|\\bFE\\b|\\bMG\\b|\\bCA\\b|\\bMN\\b|\\bCO\\b|\\bCU\\b|\\bNI\\b|\\bMO\\b|\\bW\\b|\\bYB\\b"

# Extract HETATM lines that do not correspond to standard residues or solvents
grep "^HETATM" "$pdb_file" | awk -v std_residues="$standard_residues" '
{
    resname = substr($0, 18, 3);
    if (resname == "LIG" || resname !~ ("\\b" std_residues "\\b")) {
        print;
    }
}' > "$output_file"

echo "Ligand extraction completed. Output saved to ${output_file}"
