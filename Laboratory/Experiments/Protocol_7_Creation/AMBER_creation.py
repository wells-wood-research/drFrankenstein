import os
from os import path as p
import pandas as pd
from subprocess import run
from shutil import copy
from copy import copy as object_copy
import sys
from pdbUtils import pdbUtils
import parmed as pmd
from OperatingTools import file_parsers


## CLEAN CODE CLASSES ##
class FilePath:
    pass
class DirectoryPath:
    pass
from typing import Dict, List, Tuple

OPTIONAL_BACKBONE_ALIASES = {"H", "HA"}

def extract_backbone_types_from_prmtop(moleculePrmtop: FilePath, backboneAliases: Dict[str, List[str]]) -> dict:
    """Extract AMBER backbone atom types from a PRMTOP file."""

    ## extract backbone atom types from prmtop
    parmedPrmtop = pmd.load_file(moleculePrmtop)
    prmtopDf = parmedPrmtop.to_dataframe()
    prmtopDf["name"] = prmtopDf["name"].astype(str).str.strip()

    ## gather backbone atom types from prmtop
    backboneAtomTypes = {}
    missingAliases = {}
    for canonicalName, aliases in backboneAliases.items():
        if len(aliases) == 0:
            if canonicalName in OPTIONAL_BACKBONE_ALIASES:
                continue
            missingAliases[canonicalName] = aliases
            continue

        matchedAtomType = None
        for alias in aliases:
            matchedTypes = prmtopDf.loc[prmtopDf["name"] == alias, "type"].values
            if len(matchedTypes) > 0:
                matchedAtomType = matchedTypes[0]
                break

        if matchedAtomType is None:
            if canonicalName in OPTIONAL_BACKBONE_ALIASES:
                continue
            missingAliases[canonicalName] = aliases
            continue

        ## keep all configured aliases usable downstream
        for alias in aliases:
            backboneAtomTypes[alias] = matchedAtomType

    if missingAliases:
        missingText = ", ".join([f"{name}:{aliases}" for name, aliases in missingAliases.items()])
        raise ValueError(
            f"Could not find configured backboneAliases in PRMTOP {moleculePrmtop}: {missingText}"
        )
      
    return backboneAtomTypes

def add_ncaa_ncaa_dihedrals(config: dict) -> None:
    """Add ncAA-ncAA dihedral parameters to the final frcmod."""
    ## unpack config
    finalFrcMod = config["runtimeInfo"]["madeByCreator"]["finalFrcmod"]
    moleculePrmtop = config["runtimeInfo"]["madeByStitching"]["moleculePrmtop"]
    backboneAliases = config["moleculeInfo"]["backboneAliases"]

    ## extract backbone atom types from prmtop
    backboneAtomTypes = extract_backbone_types_from_prmtop(moleculePrmtop, backboneAliases)
    ## init mapping between canonical atom Names to amber types
    amberBackboneTypes = {"O": "O",
                          "C": "C",
                          "N": "N",
                          "CA": "CT",
                          "H": "H",
                          "HA": "HC"}
    ## create a mapping between amber types to gaff2 (made by us) atom types
    amber2gaffTypes = {}
    for canonicalName, customNames in backboneAliases.items():
        if len(customNames) == 0:
            continue
        amberType = amberBackboneTypes[canonicalName]
        mappedAlias = next((alias for alias in customNames if alias in backboneAtomTypes), None)
        if mappedAlias is None:
            continue
        amber2gaffTypes[amberType] = backboneAtomTypes[mappedAlias]
    ## make ncaa-ncaa dihedrals
    parmedFrcmod = pmd.load_file(finalFrcMod)
    newParams = {}
    for atomTypes, paramObject in parmedFrcmod.dihedral_types.items():
        if any(atomType in amber2gaffTypes for atomType in atomTypes):
            newAtomNames = tuple(amber2gaffTypes.get(atomType, atomType) for atomType in atomTypes)
            print(newAtomNames)
            newParams[newAtomNames] = object_copy(paramObject)

    ## update frcmod
    if newParams:
        parmedFrcmod.dihedral_types.update(newParams)

    parmedFrcmod.write(finalFrcMod, style="frcmod")

    return None
def add_wildcard_dihedrals(config: dict) -> None:
    """Add wildcard dihedral parameters to the final frcmod."""
    ## unpack config
    finalFrcMod = config["runtimeInfo"]["madeByCreator"]["finalFrcmod"]
    moleculePrmtop = config["runtimeInfo"]["madeByStitching"]["moleculePrmtop"]
    backboneAliases = config["moleculeInfo"]["backboneAliases"]

    backboneAtomTypes = extract_backbone_types_from_prmtop(moleculePrmtop, backboneAliases)


    parmedFrcmod = pmd.load_file(finalFrcMod)
    newParams = {}
    for atomTypes, paramObject in parmedFrcmod.dihedral_types.items():
        if any(atomType in backboneAtomTypes for atomType in atomTypes):
            newAtomNames = tuple(
                ["X" if not atomType in backboneAtomTypes.keys() else atomType for atomType in atomTypes])

            newParams[newAtomNames] = object_copy(paramObject)


    ## load frcmod to be updated
    parmedFrcmod = pmd.load_file(finalFrcMod)

    newParams = {}
    # --- 1. Create a new parameter set with wildcard dihedrals ---
    for atomTypes, paramObject in parmedFrcmod.dihedral_types.items():
        if any(atomType in backboneAtomTypes for atomType in atomTypes):
            newAtomNames = tuple(
                ["X" if not atomType in backboneAtomTypes else atomType for atomType in atomTypes]
                )
           
            newParams[newAtomNames] = object_copy(paramObject)

    # --- 2. Update the main parameter set with the new ones ---
    if newParams:
        parmedFrcmod.dihedral_types.update(newParams)

    # --- 3. Write the new frcmod file ---
    parmedFrcmod.write(finalFrcMod, style="frcmod")
    return None


def duplicate_capping_parameters(config: dict) -> dict:
    """Duplicate capping parameters for terminal and capped ncAA cases."""
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
            if "CT" in atomTypes:
                for newAtomType in ["CX", "XC"]:
                    newAtomNames = tuple(
                        [newAtomType if atomType == "CT" else atomType for atomType in atomTypes]
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
    """Copy the stitched frcmod file into the final creation directory."""
    finalCreationDir = config["runtimeInfo"]["madeByCreator"]["finalCreationDir"]
    moleculeFrcmod = config["runtimeInfo"]["madeByStitching"]["moleculeFrcmod"]
    moleculeName = config["moleculeInfo"]["moleculeName"]

    finalFrcmod = p.join(finalCreationDir, f"{moleculeName}.frcmod")
    config["runtimeInfo"]["madeByCreator"]["finalFrcmod"] = finalFrcmod

    copy(moleculeFrcmod, finalFrcmod)

    return config
    
################################################################################
def get_capping_atom_ids(config: dict) -> list[int]:
    """Return the atom IDs that belong to the capping groups."""

    cappedMol2 = config["runtimeInfo"]["madeByStitching"]["moleculeMol2"]
    creationDir = config["runtimeInfo"]["madeByCreator"]["finalCreationDir"]
    tmpPdb = p.join(creationDir, "tmp.pdb")
    obabelCommand = ["obabel", "-i", "mol2", cappedMol2, "-O", tmpPdb]
    run(obabelCommand, capture_output=True)
    cappingAtomNames = ["N_N", "H_N", "C_N", "H1_N", "H2_N", "H3_N", ## NME cap
                        "C_C", "O_C", "C2_C", "H1_C", "H2_C", "H3_C"] ## ACE cap

    tmpdf = pdbUtils.pdb2df(tmpPdb)
    cappingAtomIds = tmpdf[tmpdf["ATOM_NAME"].isin(cappingAtomNames)]["ATOM_ID"].to_list()

    os.remove(tmpPdb)

    return cappingAtomIds

################################################################################

def create_final_lib_and_mol2(cappingAtomIds: list[int], config: dict) -> None:
    """Create the final AMBER lib and mol2 files."""
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
    run(tleapCommand, capture_output=True)
################################################################################
def get_capping_proton_ids(cappingHeteroAtomNames: list[str], atomDf: pd.DataFrame, bondDf: pd.DataFrame) -> list[int]:
    """Return proton IDs bound to the capping hetero atoms."""
    cappingProtonIds = []
    for cappingHeteroAtomName in cappingHeteroAtomNames:
        heteroAtomId = atomDf[atomDf["ATOM_NAME"]==cappingHeteroAtomName]["ATOM_ID"].to_list()[0]

        bondedToHeteroAtomA = bondDf[(bondDf["ATOM_A_ID"] == heteroAtomId)]["ATOM_B_ID"].to_list()
        bondedToHeteroAtomB = bondDf[(bondDf["ATOM_B_ID"] == heteroAtomId)]["ATOM_A_ID"].to_list()

        bondedToHeteroAtomIds = bondedToHeteroAtomA + bondedToHeteroAtomB
        protonIds = atomDf[atomDf["ATOM_ID"].isin(bondedToHeteroAtomIds) & atomDf["atomType"].str.startswith("h")]["ATOM_ID"].values
        cappingProtonIds.extend(protonIds)

    return cappingProtonIds
