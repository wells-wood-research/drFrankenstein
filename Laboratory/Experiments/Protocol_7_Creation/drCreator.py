import os
from os import path as p
import pandas as pd
from subprocess import call, PIPE
from shutil import copy

from . import CHARMM_creation
from . import AMBER_creation


################################################################################

def create_the_monster(config):
    ## create a runtimeInfo 
    config["runtimeInfo"]["madeByCreator"] = {}


    ## unpack config
    forcefield = config["parameterFittingInfo"]["forceField"]

    ## make a dir
    outputDir = config["pathInfo"]["outputDir"]
    finalCreationDir = p.join(outputDir, "07_final_creation")
    os.makedirs(finalCreationDir, exist_ok=True)
    config["runtimeInfo"]["madeByCreator"]["finalCreationDir"] = finalCreationDir



    ##TODO: this could be done during the capping step instead

    if forcefield == "AMBER":
        cappingAtomIds = AMBER_creation.get_capping_atom_ids(config)
        AMBER_creation.create_final_lib_and_mol2(cappingAtomIds, config)
        AMBER_creation.copy_final_frcmod(config)
    elif forcefield == "CHARMM":
        config = CHARMM_creation.get_donor_acceptors(config)
        CHARMM_creation.create_final_rtf(config)
        CHARMM_creation.copy_final_prm(config)
    exit()
    config["checkpointInfo"]["finalCreationComplete"] = True
    return config



################################################################################

if __name__ == "__main__":
    raise NotImplementedError