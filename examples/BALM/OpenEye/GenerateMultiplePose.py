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

import sys
from openeye import oechem
from openeye import oedocking


class MyOptions(oedocking.OEPositOptions):
    def __init__(self):
        oedocking.OEPositOptions.__init__(self)
        param1 = oechem.OEUIntParameter("-numPoses", 1)
        param1.AddLegalRange("1", "20")
        param1.SetBrief("Number of poses to generate")
        self._param1 = self.AddParameter(param1)
        pass

    def GetNumPoses(self):
        if self._param1.GetHasValue():
            return int(self._param1.GetStringValue())
        return int(self._param1.GetStringDefault())

def main(argv=[__name__]):
    positOpts = MyOptions()
    opts = oechem.OERefInputAppOptions(positOpts, "PoseMolecules", oechem.OEFileStringType_Mol3D,
                                       oechem.OEFileStringType_DU, oechem.OEFileStringType_DU, "-receptor")
    if oechem.OEConfigureOpts(opts, argv, False) == oechem.OEOptsConfigureStatus_Help:
        return 0
    positOpts.UpdateValues(opts)

    ifs = oechem.oemolistream()
    if not ifs.open(opts.GetInFile()):
        oechem.OEThrow.Fatal("Unable to open %s for reading" % opts.GetInFile())

    rfs = oechem.oeifstream()
    if not rfs.open(opts.GetRefFile()):
        oechem.OEThrow.Fatal("Unable to open %s for reading" % opts.GetRefFile())

    ofs = oechem.oeofstream()
    if not ofs.open(opts.GetOutFile()):
        oechem.OEThrow.Fatal("Unable to open %s for writing" % opts.GetOutFile())

    poser = oedocking.OEPosit(positOpts)
    du = oechem.OEDesignUnit()
    count = 0
    while oechem.OEReadDesignUnit(rfs, du):
        if not du.HasReceptor():
            oechem.OEThrow.Fatal("Design unit %s does not contain a receptor" % du.GetTitle())
        poser.AddReceptor(du)
        count += 1
    if count == 0:
        oechem.OEThrow.Fatal("Receptor input does not contain any design unit")

    oechem.OEThrow.Info("Number of conformers: %d" % positOpts.GetNumPoses())
    oechem.OEThrow.Info("Best Receptor pose flag: %s" % positOpts.GetBestReceptorPoseOnly())
    
    for mcmol in ifs.GetOEMols():
        oechem.OEThrow.Info("posing %s" % mcmol.GetTitle())
        results = oedocking.OEPositResults()
        ret_code = poser.Dock(results, mcmol, positOpts.GetNumPoses())

        if ret_code == oedocking.OEDockingReturnCode_Success:
            for result in results.GetSinglePoseResults():
                posedDU = result.GetDesignUnit()
                posedDU.SetDoubleData(poser.GetName(), result.GetProbability())
                oechem.OEThrow.Info("Receptor used: %s pose probability: %f" % (posedDU.GetTitle(), result.GetProbability()))
                oechem.OEWriteDesignUnit(ofs, posedDU)
        else:
            errMsg = oedocking.OEDockingReturnCodeGetName(ret_code)
            oechem.OEThrow.Warning("%s: %s" % (mcmol.GetTitle(), errMsg))
    return 0

    

if __name__ == "__main__":
    sys.exit(main(sys.argv))
