
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
from .AMBER_protocols import AMBER_torsion_protocol
from .AMBER_protocols import AMBER_total_protocol
from .AMBER_protocols import AMBER_helper_functions
from .CHARMM_protocols import CHARMM_helper_functions
from .CHARMM_protocols import CHARMM_total_protocol
from .CHARMM_protocols import CHARMM_torsion_protocol
from . import QMMM_fitting_protocol
from . import Stitching_Assistant
from . import Stitching_Plotter
from OperatingTools import Timer


# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
@Timer.time_function()
def torsion_fitting_protocol_AMBER(config: dict) -> dict:
    """
    Main protocol for torsion fitting for AMBER parameters
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
    ## create runtimeInfo entry for stitching
    config["runtimeInfo"]["madeByStitching"] = {}

    ## create output directories
    config = Stitching_Assistant.sort_out_directories(config)
    ## create a basic frcmod using antechamber
    config = AMBER_helper_functions.get_generate_initial_frcmod(config)
    ## get torsion tags from config
    torsionTags = config["runtimeInfo"]["madeByTwisting"]["torsionTags"]
    
    nShuffles = 10   ## TODO: put in config
    ## run the torsion fitting protocol, each time, shuffle the torsion order
    for shuffleIndex in range(nShuffles):
        print(f"Shuffle round {shuffleIndex+1}") ## TODO: SPLASH
        ## shuffle the order of our torsion tags
        shuffledTorsionTags = Stitching_Assistant.shuffle_torsion_tags(torsionTags)
        ## run the torsion fitting protocol for each torsion
        for torsionTag in shuffledTorsionTags:
            mmTotalEnergy = AMBER_total_protocol.get_MM_total_energies(config, torsionTag)
            mmTorsionEnergy, mmCosineComponents = AMBER_torsion_protocol.get_MM_torsion_energies(config, torsionTag)
            torsionParameterDf = QMMM_fitting_protocol.fit_torsion_parameters(config, torsionTag, mmTotalEnergy, mmTorsionEnergy, shuffleIndex, mmCosineComponents)
            config = AMBER_helper_functions.update_frcmod(config, torsionTag, torsionParameterDf)
    ## make a gif for each torsion being fitted - so satisfying!
    for torsionTag in torsionTags:
        fittingGif = p.join(config["runtimeInfo"]["madeByStitching"]["qmmmParameterFittingDir"], torsionTag, f"torsion_fitting.gif")
        torsionFittingDir = p.join(config["runtimeInfo"]["madeByStitching"]["qmmmParameterFittingDir"], torsionTag)
        Stitching_Plotter.make_gif(torsionFittingDir, fittingGif)
    ## update config checkpoint flag
    config["checkpointInfo"]["torsionFittingComplete"] = True
    return config




@Timer.time_function()
def torsion_fitting_protocol_CHARMM(config: dict) -> dict:
    """
    Main protocol for torsion fitting for CHARMM parameters
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

    ## create runtimeInfo entry for stitching
    config["runtimeInfo"]["madeByStitching"] = {}
    ## create output directories
    config = Stitching_Assistant.sort_out_directories(config)
    ## create RTF, PRM and PSF files 
    config = CHARMM_helper_functions.create_initial_CHARMM_parameters(config)
    ## get torsion tags from config
    torsionTags = config["runtimeInfo"]["madeByTwisting"]["torsionTags"]
    nShuffles = 10   ## TODO: put in config
    ## run the torsion fitting protocol, each time, shuffle the torsion order
    for shuffleIndex in range(nShuffles):
        print(f"Shuffle round {shuffleIndex+1}") ## TODO: SPLASH
        ## shuffle the order of our torsion tags
        shuffledTorsionTags = Stitching_Assistant.shuffle_torsion_tags(torsionTags)
        ## run the torsion fitting protocol for each torsion
        for torsionTag in shuffledTorsionTags:
            mmTotalEnergy = CHARMM_total_protocol.get_MM_total_energies(config, torsionTag)
            mmTorsionEnergy, mmCosineComponents = CHARMM_torsion_protocol.get_MM_torsion_energies(config, torsionTag)
            
            torsionParameterDf = QMMM_fitting_protocol.fit_torsion_parameters(config, torsionTag, mmTotalEnergy, mmTorsionEnergy, shuffleIndex, mmCosineComponents)
            config = CHARMM_helper_functions.update_prm(config, torsionTag, torsionParameterDf)
    ## make a gif for each torsion being fitted - so satisfying!
    for torsionTag in torsionTags:
        fittingGif = p.join(config["runtimeInfo"]["madeByStitching"]["qmmmParameterFittingDir"], torsionTag, f"torsion_fitting.gif")
        torsionFittingDir = p.join(config["runtimeInfo"]["madeByStitching"]["qmmmParameterFittingDir"], torsionTag)
        Stitching_Plotter.make_gif(torsionFittingDir, fittingGif)
    ## update config checkpoint flag
    config["checkpointInfo"]["torsionFittingComplete"] = True
    return config

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲####
if __name__ == "__main__":
    raise NotImplementedError