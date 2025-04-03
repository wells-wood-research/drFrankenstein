
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

## drFRANKENSTEIN LIBRARIES ##
from . import MM_torsion_protocol
from . import MM_total_protocol
from . import QMMM_fitting_protocol
from . import Stitching_Assistant
from . import Stitching_Plotter
from OperatingTools import Timer


# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
@Timer.time_function()
def torsion_fitting_protocol(config: dict) -> dict:
    """
    Main protocol for torsion fitting
    For each torsion that has had QM scans performed (creates QM[total]):
        1. Get MM total energies for all conformers using OpenMM
        2. Get MM torsion energies for all conformers straight from the parameters
        3. Calculate QM[Torsion] by {QM[torsion] = QM[total] - MM[total] + MM[torsion]}
        4. Use the fourier transform  method to fit cosine parameters to QM[torsion]
        5. Update parameters (this changes the MM[total] and MM[torsion] for subsequent torsions)

    The above protocol is repeated several times, each time the order of the torsions is shuffled

    Args:
        config (dict) : the drFrankenstein config containing all run information
    Returns:
        config (dict): updated config
    
    """
    ## create output directories
    config = Stitching_Assistant.sort_out_directories(config)
    ## create a basic frcmod using antechamber
    config = MM_torsion_protocol.get_generate_initial_frcmod(config)
    ## get torsion tags from config
    torsionTags = config["torsionScanInfo"]["torsionTags"]
    
    nShuffles = 10   ## TODO: put in config
    ## run the torsion fitting protocol, each time, shuffle the torsion order
    for shuffleIndex in range(nShuffles):
        print(f"Shuffle round {shuffleIndex+1}") ## TODO: SPLASH
        ## shuffle the order of our torsion tags
        shuffledTorsionTags = Stitching_Assistant.shuffle_torsion_tags(torsionTags)
        ## run the torsion fitting protocol for each torsion
        for torsionTag in shuffledTorsionTags:
            mmTotalEnergy = MM_total_protocol.get_MM_total_energies(config, torsionTag)
            mmTorsionEnergy, mmCosineComponents = MM_torsion_protocol.get_MM_torsion_energies(config, torsionTag)
            torsionParameterDf = QMMM_fitting_protocol.fit_torsion_parameters(config, torsionTag, mmTotalEnergy, mmTorsionEnergy, shuffleIndex, mmCosineComponents)
            config = Stitching_Assistant.update_frcmod(config, torsionTag, torsionParameterDf)
    ## make a gif for each torsion being fitted - so satisfying!
    for torsionTag in torsionTags:
        fittingGif = p.join(config["pathInfo"]["qmmmParameterFittingDir"], torsionTag, f"torsion_fitting.gif")
        torsionFittingDir = p.join(config["pathInfo"]["qmmmParameterFittingDir"], torsionTag)
        Stitching_Plotter.make_gif(torsionFittingDir, fittingGif)
    ## update config checkpoint flag
    config["checkpointInfo"]["torsionFittingComplete"] = True
    return config

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲####
if __name__ == "__main__":
    raise NotImplementedError