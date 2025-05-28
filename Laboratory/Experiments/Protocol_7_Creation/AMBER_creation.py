import os
from os import path as p
import pandas as pd
from subprocess import call, PIPE
from shutil import copy


from OperatingTools import file_parsers
from typing import List, Dict # Added for type hinting

# Placeholder classes (extend if needed)
class FilePath:
    pass
class DirectoryPath:
    pass

################################################################################
def copy_final_frcmod(config: dict) -> None:
    """
    Copies the proposed FRCMOD file to the final creation directory.

    Args:
        config: Configuration dictionary containing paths and molecule information.
    """
    finalCreationDir = config["runtimeInfo"]["madeByCreator"]["finalCreationDir"]
    proposedFrcmod = config["runtimeInfo"]["madeByStitching"]["proposedFrcmod"]
    moleculeName = config["moleculeInfo"]["moleculeName"]

    finalFrcmod = p.join(finalCreationDir, f"{moleculeName}.frcmod")

    copy(proposedFrcmod, finalFrcmod)
    
################################################################################
def get_capping_atom_ids(config: dict) -> List[int]:
    """
    Retrieves a list of atom IDs belonging to capping groups (NME or ACE).

    Args:
        config: Configuration dictionary containing charge group information.

    Returns:
        A list of integer atom IDs for capping group atoms.
    """
    chargeGroups: Dict[str, Dict] = config["runtimeInfo"]["madeByCharges"]["chargeGroups"]

    cappingAtomIds: List[int] = []
    for chargeGroupName, chargeGroupData in chargeGroups.items():
        if chargeGroupName.startswith("NME_cap") or chargeGroupName.startswith("ACE_cap"):
            cappingAtomIds.extend(chargeGroupData["indexes"])

    return cappingAtomIds

################################################################################

def create_final_lib_and_mol2(capping_atom_ids: List[int], config: dict) -> None:
    """
    Creates final .lib and .mol2 files for the uncapped molecule using TLEAP.

    Args:
        capping_atom_ids: A list of atom IDs to be removed (capping groups).
        config: Configuration dictionary containing paths and molecule information.
    """
    cappedMol2 = config["runtimeInfo"]["madeByStitching"]["moleculeMol2"]

    finalCreationDir: DirectoryPath = config["runtimeInfo"]["madeByCreator"]["finalCreationDir"] # type: ignore
    moleculeName: str = config["moleculeInfo"]["moleculeName"]
    finalMol2: FilePath = p.join(finalCreationDir, f"{moleculeName}.mol2") # type: ignore
    finalLib: FilePath = p.join(finalCreationDir, f"{moleculeName}.lib") # type: ignore

    finalTleapInput: FilePath = p.join(finalCreationDir, "tleap.in") # type: ignore
    
    with open(finalTleapInput, "w") as f:
        f.write("source leaprc.gaff2 \n")
        f.write(f"{moleculeName}  = loadmol2  {cappedMol2}\n")
        f.write(f"set {moleculeName} head {moleculeName}.1.N \n")
        f.write(f"set {moleculeName} tail {moleculeName}.1.C \n")
        for atomId in capping_atom_ids:
            f.write(f"remove {moleculeName} {moleculeName}.1.{atomId}\n")
        f.write(f"check {moleculeName} \n")
        f.write(f"saveoff  {moleculeName} {finalLib} \n")
        f.write(f"saveMol2 {moleculeName} {finalMol2} 1\n")
        f.write("quit")

    tleapOutput: FilePath = p.join(finalCreationDir, f"tleap.out") # type: ignore

    tleapCommand: list = ["tleap", "-f", finalTleapInput, ">", tleapOutput]

    os.chdir(finalCreationDir) # type: ignore
    call(tleapCommand, stdout=PIPE)
################################################################################
def get_capping_proton_ids(capping_hetero_atom_names: List[str], 
                            atom_df: pd.DataFrame, 
                            bond_df: pd.DataFrame) -> List[int]:
    """
    Identifies proton IDs that are bonded to specified capping heteroatoms.

    Args:
        capping_hetero_atom_names: A list of names for heteroatoms in capping groups.
        atom_df: DataFrame containing atom information (ATOM_ID, ATOM_NAME, ATOM_TYPE).
        bond_df: DataFrame containing bond information (ATOM_A_ID, ATOM_B_ID).

    Returns:
        A list of integer atom IDs for protons bonded to the capping heteroatoms.
    """
    cappingProtonIds: List[int] = []
    for cappingHeteroAtomName in capping_hetero_atom_names:
        heteroAtomId = atom_df[atom_df["ATOM_NAME"]==cappingHeteroAtomName]["ATOM_ID"].to_list()[0]

        bondedToHeteroAtom_A = bond_df[bond_df["ATOM_A_ID"] == heteroAtomId]["ATOM_B_ID"].to_list()
        bondedToHeteroAtom_B = bond_df[bond_df["ATOM_B_ID"] == heteroAtomId]["ATOM_A_ID"].to_list()

        bondedToHeteroAtomIds = bondedToHeteroAtom_A + bondedToHeteroAtom_B
        protonIds = atom_df[atom_df["ATOM_ID"].isin(bondedToHeteroAtomIds) & atom_df["ATOM_TYPE"].str.startswith("h")]["ATOM_ID"].values
        cappingProtonIds.extend(protonIds)

    return cappingProtonIds