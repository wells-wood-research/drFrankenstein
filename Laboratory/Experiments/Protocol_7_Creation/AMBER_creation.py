import os
from os import path as p
import pandas as pd
from subprocess import call, PIPE
from shutil import copy
from copy import copy as object_copy
import sys

import parmed as pmd
from OperatingTools import file_parsers


def duplicate_capping_parameters(config: dict) -> dict:
    """
    Loads a frcmod file, identifies parameters containing the 'CX' atom type,
    duplicates them by creating new, distinct parameter objects for 'CT' and 'XC',
    and writes a new frcmod file.

    Args:
        config (dict): A dictionary containing the configuration, including the
                       path to the input frcmod file.
    """
    finalFrcMod = config["runtimeInfo"]["madeByCreator"]["finalFrcmod"]

    params = pmd.load_file(finalFrcMod)

    newParamSections = {
        'bond': params.bond_types,
        'angle': params.angle_types,
        'dihedral': params.dihedral_types,
        'improper': params.improper_types,
    }

    # --- 1. Identify and duplicate parameters containing 'CX' ---
    for section_name, param_dict in newParamSections.items():
        newParams = {}

        for atomTypes, paramObject in param_dict.items():
            if "CX" in atomTypes:
                for newAtomType in ["CT", "XC"]:
                    newAtomNames = tuple(
                        [newAtomType if atomType == "CX" else atomType for atomType in atomTypes]
                    )
                    
                    newParams[newAtomNames] = object_copy(paramObject)

        # --- 2. Update the main parameter set with the new ones ---
        if newParams:
            param_dict.update(newParams)

    # --- 3. Write the new frcmod file ---

    params.write(finalFrcMod, style="frcmod")

    
    return config


################################################################################
def copy_final_frcmod(config: dict) -> dict:
    """
    Copies the final frcmod file from the Stitching directory to the Creator directory.

    Args:
        config (dict): The configuration dictionary.

    Returns:
        config (dict): The updated configuration dictionary.
    """
    finalCreationDir = config["runtimeInfo"]["madeByCreator"]["finalCreationDir"]
    moleculeFrcmod = config["runtimeInfo"]["madeByStitching"]["moleculeFrcmod"]
    moleculeName = config["moleculeInfo"]["moleculeName"]

    finalFrcmod = p.join(finalCreationDir, f"{moleculeName}.frcmod")
    config["runtimeInfo"]["madeByCreator"]["finalFrcmod"] = finalFrcmod

    copy(moleculeFrcmod, finalFrcmod)

    return config
    
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
    cappedMol2 = config["runtimeInfo"]["madeByAssembly"]["cappedMol2"]
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

        bondedToHeteroAtomA = bondDf[(bondDf["ATOM_A_ID"] == heteroAtomId)]["ATOM_B_ID"].to_list()
        bondedToHeteroAtomB = bondDf[(bondDf["ATOM_B_ID"] == heteroAtomId)]["ATOM_A_ID"].to_list()

        bondedToHeteroAtomIds = bondedToHeteroAtomA + bondedToHeteroAtomB
        protonIds = atomDf[atomDf["ATOM_ID"].isin(bondedToHeteroAtomIds) & atomDf["atomType"].str.startswith("h")]["ATOM_ID"].values
        cappingProtonIds.extend(protonIds)

    return cappingProtonIds