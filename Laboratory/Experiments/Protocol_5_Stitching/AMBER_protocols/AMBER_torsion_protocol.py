## BASIC LIBRARIES ##
from os import path as p
import pandas as pd

## drFRANKENSTEIN LIBRARIES ##
from . import AMBER_helper_functions

## CLEAN CODE CLASSES ##
from typing import Tuple
class FilePath:
    pass
class DirectoryPath:
    pass

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def get_MM_torsion_energies(config: dict, torsionTag: str) -> Tuple[dict, dict]:    
    """
    Gets MM[torsion] energy for each torsion we have scanned
    This is done by extracting torsion parameters from FRCMOD file
    Then reconstructing energy as a sum of the cosine functions described

    Args:
        config (dict): config containing all run information
        torsionTag (str): tag for torsion we are interested in

    Returns:
        mmTorsionEnergies (dict): energies 
        mmCosineComponents (dict): cosine components
    """
    ## get torsion parameters from FRCMOD file
    mmTorsionParameters = extract_torsion_parameters_from_frcmod(config, torsionTag)

    ## reconstruct torsion energies from parameters
    mmTorsionEnergies, mmCosineComponents = AMBER_helper_functions.construct_MM_torsion_energies(mmTorsionParameters)

    return  mmTorsionEnergies, mmCosineComponents

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def extract_torsion_parameters_from_frcmod(config: dict, torsionTag: str) -> dict:
    """
    Reads through a frcmod file 
    Finds torsion parameters for each torsion that we have scanned
    Returns a dict with the torsion tag as the key and the torsion parameters as the value

    Args:
        molFrcmod (FilePath): frcmod file
        atomTypeMap (dict): dict mapping atom names to atom types
        config (dict): config dict

    Returns:
        mmTorsionParameters (dict): dict with the torsion tag as the key and the torsion parameters as the value
    """

    molFrcmod = config["runtimeInfo"]["madeByStitching"]["moleculeFrcmod"]
    parsedFrcmod: dict = AMBER_helper_functions.parse_frcmod(molFrcmod)

    atomTypeMap = config["runtimeInfo"]["madeByStitching"]["atomTypeMap"]

    ## get torsion atom names
    torsionAtoms: list = torsionTag.split("-")
    ## get torsion atom types
    torsionAtomTypes: list = [atomTypeMap[atom] for atom in torsionAtoms]
    ## find torsion parameters for these atom types (account for reverse torsion order)
    torsionAtomTypes_reversed = torsionAtomTypes[::-1]
    lookForAtoms = [torsionAtomTypes, torsionAtomTypes_reversed]
    try:
        torsionParameters = [entry for entry in parsedFrcmod["DIHEDRALS"] if entry["atoms"] in lookForAtoms]

    except:
        raise Exception(f"Torsion {torsionTag} not found in frcmod file")

    return torsionParameters
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
