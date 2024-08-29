#!/usr/bin/env python
# (C) 2022 Cadence Design Systems, Inc. (Cadence) 
# All rights reserved.
# TERMS FOR USE OF SAMPLE CODE The software below ("Sample Code") is
# provided to current licensees or subscribers of Cadence products or
# SaaS offerings (each a "Customer").
# Customer is hereby permitted to use, copy, and modify the Sample Code,
# subject to these terms. Cadence claims no rights to Customer's
# modifications. Modification of Sample Code is at Customer's sole and
# exclusive risk. Sample Code may require Customer to have a then
# current license or subscription to the applicable Cadence offering.
# THE SAMPLE CODE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED.  OPENEYE DISCLAIMS ALL WARRANTIES, INCLUDING, BUT
# NOT LIMITED TO, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
# PARTICULAR PURPOSE AND NONINFRINGEMENT. In no event shall Cadence be
# liable for any damages or liability in connection with the Sample Code
# or its use.

#############################################################################
# This program demonstrates how to extract the ligand from a DesignUnit.
#############################################################################
import sys
from openeye import oechem


def ExtractLigandFromDU(du, ofs):
# @ <SNIPPET-ExtractLigandFromDesignUnit>
    lig = oechem.OEGraphMol()
    if not du.GetLigand(lig):
        oechem.OEThrow.Fatal("Error: Could not extract ligand from the OEDesignUnit.")
    oechem.OEWriteMolecule(ofs, lig)
# @ </SNIPPET-ExtractLigandFromDesignUnit>

    ofs.close()


def main(argv=[__name__]):
    if len(sys.argv) != 3:
        oechem.OEThrow.Usage("%s [<oedu infile>] [<ligand outfile>]" % argv[0])

    du = oechem.OEDesignUnit()
    if not oechem.OEReadDesignUnit(argv[1], du):
        oechem.OEThrow.Fatal("Unable to open %s for reading OEDesignUnit" % argv[1])

    ofs = oechem.oemolostream()
    if not ofs.open(argv[2]):
        oechem.OEThrow.Fatal("Unable to open %s for writing" % argv[2])

    ExtractLigandFromDU(du, ofs)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
