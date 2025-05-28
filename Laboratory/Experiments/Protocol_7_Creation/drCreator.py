import os
from os import path as p
import pandas as pd
from subprocess import call, PIPE
from shutil import copy

from . import CHARMM_creation
from . import AMBER_creation


from typing import Dict, List # Added for type hinting

# Placeholder classes (extend if needed)
class FilePath:
    pass
class DirectoryPath:
    pass
################################################################################

def create_the_monster(config: dict) -> dict:
    """
    Main protocol to create the final molecular representation (AMBER or CHARMM).
    This involves calling specific creation routines based on the chosen force field,
    handling capping groups, and generating final topology/parameter files.

    Args:
        config: Configuration dictionary containing all necessary paths, 
                molecule information, and previous protocol outputs.

    Returns:
        The updated configuration dictionary with paths to the final created files
        and a checkpoint flag indicating completion.
    """
    ## create a runtimeInfo 
    config["runtimeInfo"]["madeByCreator"] = {}

    ## unpack config
    forcefield: str = config["parameterFittingInfo"]["forceField"]

    ## make a dir
    outputDir: DirectoryPath = config["pathInfo"]["outputDir"] # type: ignore
    finalCreationDir: DirectoryPath = p.join(outputDir, "07_final_creation") # type: ignore
    os.makedirs(finalCreationDir, exist_ok=True)
    config["runtimeInfo"]["madeByCreator"]["finalCreationDir"] = finalCreationDir

    if forcefield == "AMBER":
        cappingAtomIds: List[int] = AMBER_creation.get_capping_atom_ids(config)
        AMBER_creation.create_final_lib_and_mol2(capping_atom_ids=cappingAtomIds, config=config)
        AMBER_creation.copy_final_frcmod(config=config)
    elif forcefield == "CHARMM":
        config = CHARMM_creation.get_donor_acceptors(config=config)
        CHARMM_creation.create_final_rtf(config=config)
        CHARMM_creation.copy_final_prm(config=config)
    config["checkpointInfo"]["finalCreationComplete"] = True
    return config



################################################################################

if __name__ == "__main__":
    raise NotImplementedError