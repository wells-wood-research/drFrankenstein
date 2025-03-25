
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
from drHelper import print_dict
from drHybrid import MM_torsion_protocol, MM_total_protocol, QMMM_fitting_protocol, shared_utils, Plotter



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
            mmTorsionEnergy, mmCosineComponents = MM_torsion_protocol.get_MM_torsion_energies(config, torsionTag)
            torsionParameterDf = QMMM_fitting_protocol.fit_torsion_parameters(config, torsionTag, mmTotalEnergy, mmTorsionEnergy, shuffleIndex, mmCosineComponents)
            config = shared_utils.update_frcmod(config, torsionTag, torsionParameterDf)

    for torsionTag in torsionTags:
        fittingGif = p.join(config["pathInfo"]["qmmmParameterFittingDir"], torsionTag, f"torsion_fitting.gif")
        torsionFittingDir = p.join(config["pathInfo"]["qmmmParameterFittingDir"], torsionTag)
        Plotter.make_gif(torsionFittingDir, fittingGif)
    
    config["checkpointInfo"]["torsionFittingComplete"] = True
    return config

####################################################################
if __name__ == "__main__":
    raise NotImplementedError