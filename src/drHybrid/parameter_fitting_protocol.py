
## BASIC LIBRARIES ##
import os
from os import path as p
import sys
import random

## CLEAN CODE CLASSES ##
class FilePath:
    pass
class DirectoryPath:
    pass
## ADD SRC TO PATH ##
currentFilePath: FilePath = os.path.abspath(__file__)
currentDir: DirectoryPath = os.path.dirname(currentFilePath)
srcDir: DirectoryPath = os.path.dirname(currentDir)
sys.path.append(srcDir)

## drFRANKENSTEIN LIBRARIES ##
import drInputs
import drFourier
from drHelper import print_dict
from drHybrid import MM_torsion_protocol, MM_total_protocol, QMMM_fitting_protocol, shared_utils

####################################################################
def dummy_inputs():

    ## get config yaml file from argpass command-line-argument
    configYaml = drInputs.get_config_input_arg()
    ## load into dict, check for bad formatting
    ## TODO: write config_checker for bad args
    config = drInputs.read_input_yaml(configYaml)

    moleculeName = "_".join(p.basename(config['moleculeInfo']['moleculePdb']).split(".")[0].split("_"))

    config["pathInfo"]["cappedPdb"] = p.join(config["pathInfo"]["outputDir"], "01_capped_amino_acids", f"{moleculeName}_capped.pdb")

    config["pathInfo"]["chargeFittingDir"] = p.join(config["pathInfo"]["outputDir"], "03_charge_calculations", "charge_fitting")
    config["pathInfo"]["torsionScanningDir"] = p.join(config["pathInfo"]["outputDir"], "02_torsion_scanning")

    ## TODO: get torsion scanning protocol to add this
    torsionTags = []
    for dirName in os.listdir(config["pathInfo"]["torsionScanningDir"]):
        if not dirName.startswith("torsion"):
            continue
        torsionTags.append(dirName.split("_")[1])

    config["torsionScanInfo"]["torsionTags"] = torsionTags

    config["pathInfo"]["gaff2Dat"] = "/home/esp/anaconda3/envs/Igor/dat/leap/parm/gaff2.dat"

    return config

####################################################################
def torsion_fitting_protocol(config):
    
    config = shared_utils.sort_out_directories(config)
    config = MM_torsion_protocol.get_generate_initial_frcmod(config)


    ##TODO: get torsion scanning to collate final_scan_energies to one csv, put that in config
    torsionTags = config["torsionScanInfo"]["torsionTags"]
    
    nShuffles = 10

    for shuffleIndex in range(nShuffles):
        print(f"Shuffle round {shuffleIndex+1}")
        shuffledTorsionTags = shared_utils.shuffle_torsion_tags(torsionTags)
        for torsionTag in shuffledTorsionTags:
            mmTotalEnergy = MM_total_protocol.get_MM_total_energies(config, torsionTag)
            mmTorsionEnergy = MM_torsion_protocol.get_MM_torsion_energies(config, torsionTag)
            torsionParameterDf = QMMM_fitting_protocol.fit_torsion_parameters(config, torsionTag, mmTotalEnergy, mmTorsionEnergy)
            config = shared_utils.update_frcmod(config, torsionTag, torsionParameterDf)
    
    config["checkpointInfo"]["torsionFittingComplete"] = True
    return config

####################################################################
if __name__ == "__main__":
    config = dummy_inputs()
    torsion_fitting_protocol(config)