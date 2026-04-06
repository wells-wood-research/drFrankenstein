import os
from os import path as p
import warnings
import parmed
from parmed.exceptions import ParameterWarning

# Suppress ParameterWarning
warnings.filterwarnings('ignore', category=ParameterWarning)

from .. import Assembly_Assistant

from typing import Tuple
class FilePath:
    pass
class DirectoryPath:
    pass

def replace_parameters(moleculeFrcmod: FilePath,
                        parameterDataset: FilePath) -> None:
    """
    Reads though a FRCMOD file and replaces parameters with those
    from another FRCMOD / DAT file
    Overwrites original file

    Args:
        moleculeFrcmod (FilePath): path to original FRCMOD file
        paraneterDataset (FilePath): path to FRCMOD / DAT file
    Returns:
        None
    """

    ## load parameter files into parmed
    molParams = parmed.load_file(moleculeFrcmod)
    ffParams = parmed.amber.AmberParameterSet(parameterDataset)

    # Replace BOND parameters
    for bondKey in molParams.bond_types:
        if bondKey in ffParams.bond_types:
            bondKey = ffParams.bond_types[bondKey]

    # Replace ANGLE parameters
    for angleKey in molParams.angle_types:
        if angleKey in ffParams.angle_types:
            angleKey = ffParams.angle_types[angleKey]

    # Replace DIHEDRAL parameters
    for dihedralKey in molParams.dihedral_types:
        wildCardDihedralKey = ("X", dihedralKey[1], dihedralKey[2], "X")
        if dihedralKey in ffParams.dihedral_types:
            dihedralKey = ffParams.dihedral_types[dihedralKey]
        elif wildCardDihedralKey in ffParams.dihedral_types:
            dihedralKey = ffParams.dihedral_types[wildCardDihedralKey]
            
    ## hard-coded changes to PHI and PSI using "CX" wildcard ##
    ## manual change of PHI
    molParams.dihedral_types[("C", "N", "CA", "C")] = ffParams.dihedral_types[("C", "N", "CX", "C")] 
    ## manual change of PSI
    molParams.dihedral_types[("N", "CA", "C", "N")] = ffParams.dihedral_types[("N", "CX", "C", "N")] 


    for improperKey in molParams.improper_types:
        if improperKey in ffParams.improper_types:
            improperKey = ffParams.improper_types[improperKey]
    
    molParams.write(moleculeFrcmod, style="frcmod")
    return None

def change_capping_types_amber(mol2File: FilePath, config: dict) -> None:
    """
    Uses a map to change capping atom types from auto-assigned gaff2
    to protein-style amber types

    Args:
        mol2File (FilePath): path to mol2 file
        config (dict): config dict

    Returns:
        None
    """
    # Define AMBER19ff atom type mapping
    amberDefaultMap = {
                    "N_N": "N",
                    "H_N": "H",
                    "C_N": "CX", 
                    "H1_N": "H3",
                    "H2_N": "H3",
                    "H2_N": "H3",
                    "C_C": "C",
                    "O_C": "O",
                    "CM_C": "CX", 
                    "H1_C": "H3",
                    "H2_C": "H3",
                    "H3_C": "H3"
    }
    
    parmedStructure = parmed.load_file(mol2File)
    for atom in parmedStructure.atoms:
        if atom.name in amberDefaultMap:
            atom.type = amberDefaultMap[atom.name]

    parmedStructure.save(mol2File, overwrite=True)

    return None

def change_backbone_types_amber(mol2File: FilePath, config: dict) -> None:
    """
    Uses a map to change backbone atom types from auto-assigned gaff2
    to protein-style amber types

    Args:
        mol2File (FilePath): path to mol2 file
        config (dict): config dict

    Returns:
        None

    """

    # Define AMBER19ff atom type mapping
    amberDefaultMap = {
                    "N" : "N",
                    "HN": "H",
                    "CA": "CT",
                    "HA": "H1",
                    "C" : "C",
                    "O" : "O",
                    "N_N": "N",
                    "H_N": "H",
                    "C_N": "CT", 
                    "H1_N": "H3",
                    "H2_N": "H3",
                    "H3_N": "H3",
                    "C_C": "C",
                    "O_C": "O",
                    "C2_C": "CT", 
                    "H1_C": "H3",
                    "H2_C": "H3",
                    "H3_C": "H3"
    }

    parmedStructure = parmed.load_file(mol2File)
    for atom in parmedStructure.atoms:
        if atom.name in amberDefaultMap:
            atom.type = amberDefaultMap[atom.name]

    parmedStructure.save(mol2File, overwrite=True)



