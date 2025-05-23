
## BASIC LIBRARIES ##
import os
from os import path as p
import sys
import random
from tqdm import tqdm

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
from OperatingTools import Timer, cleaner


# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
@Timer.time_function("Parameter Fitting", "PARAMETER FITTING")
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
    config = AMBER_helper_functions.copy_assembled_parameters(config)

    AMBER_helper_functions.edit_mol2_partial_charges(config)
    AMBER_helper_functions.run_tleap_to_make_params(config)
    # Minimal test

    tqdmBarOptions = {
        "desc": f"\033[32mRunning Parameter Fitting\033[0m",
        "ascii": "-🗲→",    
        "colour": "yellow",
        "unit":  "scan",
        "dynamic_ncols": False, 
        "ncols": 102,
        "leave": True,
        "position": 1
    }

    ## get torsion tags from config
    torsionTags = config["runtimeInfo"]["madeByTwisting"]["torsionTags"]
    
    nShuffles = config["parameterFittingInfo"]["nShuffles"]
    shuffledTorsionTags = []
    for i in range(nShuffles):
        random.shuffle(torsionTags)
        shuffledTorsionTags.extend(torsionTags)

    currentParameters = {}
    shuffleIndex = 0
    ## run the torsion fitting protocol, each time, shuffle the torsion order
    for  torsionTag in tqdm(shuffledTorsionTags, **tqdmBarOptions):
        if not shuffleIndex == 0:
            config["runtimeInfo"]["madeByStitching"]["moleculeFrcmod"] = config["runtimeInfo"]["madeByStitching"]["proposedFrcmod"]
            AMBER_helper_functions.run_tleap_to_make_params(config)
        mmTotalEnergy = AMBER_total_protocol.get_MM_total_energies(config, torsionTag)
        mmTorsionEnergy, mmCosineComponents = AMBER_torsion_protocol.get_MM_torsion_energies(config, torsionTag)
        torsionParameterDf = QMMM_fitting_protocol.fit_torsion_parameters(config, torsionTag, mmTotalEnergy, mmTorsionEnergy, shuffleIndex, mmCosineComponents)
        currentParameters[torsionTag] = torsionParameterDf.to_dict()
        config = AMBER_helper_functions.update_frcmod(config, torsionTag, torsionParameterDf, shuffleIndex)
        shuffleIndex += 1

    config["runtimeInfo"]["madeByStitching"]["moleculeFrcmod"] = config["runtimeInfo"]["madeByStitching"]["proposedFrcmod"]
    config["runtimeInfo"]["madeByStitching"]["finalParameters"] = currentParameters
    ## make gif out of the torsion fitting
    for torsionTag in torsionTags:
        fittingGif = p.join(config["runtimeInfo"]["madeByStitching"]["qmmmParameterFittingDir"], torsionTag, f"torsion_fitting.gif")
        torsionFittingDir = p.join(config["runtimeInfo"]["madeByStitching"]["qmmmParameterFittingDir"], torsionTag)
        Stitching_Plotter.make_gif(torsionFittingDir, fittingGif)
    ## clean up 
    cleaner.clean_up_stitching(config)

    ## update config checkpoint flag
    config["checkpointInfo"]["torsionFittingComplete"] = True
    return config
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
@Timer.time_function("Parameter Fitting", "PARAMETER_FITTING")
def torsion_fitting_protocol_CHARMM(config: dict, debug = False) -> dict:
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
    config = CHARMM_helper_functions.copy_assembled_parameters(config)

    ## get torsion tags from config
    torsionTags = config["runtimeInfo"]["madeByTwisting"]["torsionTags"]
    
    nShuffles = config["parameterFittingInfo"]["nShuffles"]
    shuffledTorsionTags = []
    for i in range(nShuffles):
        random.shuffle(torsionTags)
        shuffledTorsionTags.extend(torsionTags)

    tqdmBarOptions = {
        "desc": f"\033[32mRunning Parameter Fitting\033[0m",
        "ascii": "-🗲→",    
        "colour": "yellow",
        "unit":  "scan",
        "dynamic_ncols": False, 
        "ncols": 102,
        "leave": True,
        "position": 1
    }
    currentParameters = {}
    counter = 1
    shuffleIndex = 0
    ## run the torsion fitting protocol, each time, shuffle the torsion order
    for  torsionTag in tqdm(shuffledTorsionTags, **tqdmBarOptions):
        if not counter == 1:
            config["runtimeInfo"]["madeByStitching"]["moleculePrm"] = config["runtimeInfo"]["madeByStitching"]["proposedPrm"]
        ## run the torsion fitting protocol for each torsion
        mmTotalEnergy = CHARMM_total_protocol.get_MM_total_energies(config, torsionTag, debug)
        mmTorsionEnergy, mmCosineComponents = CHARMM_torsion_protocol.get_MM_torsion_energies(config, torsionTag, debug)
        torsionParameterDf = QMMM_fitting_protocol.fit_torsion_parameters(config, torsionTag, mmTotalEnergy, mmTorsionEnergy, shuffleIndex, mmCosineComponents, debug)
        config = CHARMM_helper_functions.update_prm(config, torsionTag, torsionParameterDf, shuffleIndex)
        currentParameters[torsionTag] = torsionParameterDf.to_dict(orient = "records")
        if counter % len(torsionTags) == 0:
            shuffleIndex += 1
        counter += 1
    ## update config moleculePrm
    config["runtimeInfo"]["madeByStitching"]["moleculePrm"] = config["runtimeInfo"]["madeByStitching"]["proposedPrm"]
    config["runtimeInfo"]["madeByStitching"]["finalParameters"] = currentParameters
    ## make a gif for each torsion being fitted - so satisfying!
    for torsionTag in torsionTags:
        fittingGif = p.join(config["runtimeInfo"]["madeByStitching"]["qmmmParameterFittingDir"], torsionTag, f"torsion_fitting.gif")
        torsionFittingDir = p.join(config["runtimeInfo"]["madeByStitching"]["qmmmParameterFittingDir"], torsionTag)
        Stitching_Plotter.make_gif(torsionFittingDir, fittingGif)
        
    ## clean up 
    cleaner.clean_up_stitching(config)

    ## update config checkpoint flag
    config["checkpointInfo"]["torsionFittingComplete"] = True
    return config

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲####
if __name__ == "__main__":
    raise NotImplementedError