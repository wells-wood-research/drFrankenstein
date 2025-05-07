import os
from os import path as p
import pandas as pd
from subprocess import call, PIPE
from shutil import copy


from OperatingTools import file_parsers
################################################################################
def copy_final_frcmod(config):
    finalCreationDir = config["runtimeInfo"]["madeByCreator"]["finalCreationDir"]
    proposedFrcmod = config["runtimeInfo"]["madeByStitching"]["proposedFrcmod"]
    moleculeName = config["moleculeInfo"]["moleculeName"]

    finalFrcmod = p.join(finalCreationDir, f"{moleculeName}.frcmod")

    copy(proposedFrcmod, finalFrcmod)
    
################################################################################
def get_capping_atom_ids(config):
    chargeGroups = config["runtimeInfo"]["madeByCharges"]["chargeGroups"]

    cappingAtomIds = []
    for chargeGroupName, chargeGroupData in chargeGroups.items():
        if chargeGroupName.startswith("NME_cap") or chargeGroupName.startswith("ACE_cap"):
            cappingAtomIds.extend(chargeGroupData["indexes"])

    return cappingAtomIds

################################################################################

def create_final_lib_and_mol2(cappingAtomIds, config):
    cappedMol2 = config["runtimeInfo"]["madeByStitching"]["moleculeMol2"]

    finalCreationDir = config["runtimeInfo"]["madeByCreator"]["finalCreationDir"]
    moleculeName = config["moleculeInfo"]["moleculeName"]
    finalMol2 = p.join(finalCreationDir, f"{moleculeName}.mol2")
    finalLib= p.join(finalCreationDir, f"{moleculeName}.lib")

    finalTleapInput = p.join(finalCreationDir, "tleap.in")
    
    with open(finalTleapInput, "w") as f:
        f.write("source leaprc.gaff2 \n")
        f.write(f"{moleculeName}  = loadmol2  {cappedMol2}\n")
        f.write(f"set {moleculeName} head {moleculeName}.1.N \n")
        f.write(f"set {moleculeName} tail {moleculeName}.1.C \n")
        for atomId in cappingAtomIds:
            f.write(f"remove {moleculeName} {moleculeName}.1.{atomId}\n")
        f.write(f"check {moleculeName} \n")
        f.write(f"saveoff  {moleculeName} {finalLib} \n")
        f.write(f"saveMol2 {moleculeName} {finalMol2} 1\n")
        f.write("quit")

    tleapOutput = p.join(finalCreationDir, f"tleap.out")

    tleapCommand: list = ["tleap", "-f", finalTleapInput, ">", tleapOutput]

    os.chdir(finalCreationDir)
    call(tleapCommand, stdout=PIPE)
################################################################################
def get_capping_proton_ids(cappingHeteroAtomNames, atomDf, bondDf):
    cappingProtonIds = []
    for cappingHeteroAtomName in cappingHeteroAtomNames:
        heteroAtomId = atomDf[atomDf["ATOM_NAME"]==cappingHeteroAtomName]["ATOM_ID"].to_list()[0]

        bondedToHeteroAtom_A = bondDf[(bondDf["ATOM_A_ID"] == heteroAtomId)]["ATOM_B_ID"].to_list()
        bondedToHeteroAtom_B = bondDf[(bondDf["ATOM_B_ID"] == heteroAtomId)]["ATOM_A_ID"].to_list()

        bondedToHeteroAtomIds = bondedToHeteroAtom_A + bondedToHeteroAtom_B
        protonIds = atomDf[atomDf["ATOM_ID"].isin(bondedToHeteroAtomIds) & atomDf["ATOM_TYPE"].str.startswith("h")]["ATOM_ID"].values
        cappingProtonIds.extend(protonIds)

    return cappingProtonIds