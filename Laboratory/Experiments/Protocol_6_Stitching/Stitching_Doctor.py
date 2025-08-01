
## BASIC LIBRARIES ##
import os
from os import path as p
import sys
import random
from tqdm import tqdm
import pandas as pd

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
def torsion_fitting_protocol_AMBER(config: dict, debug=False) -> dict:
    """
    Main protocol for torsion fitting for AMBER parameters.
    For each torsion that has had QM scans performed (creates QM[total]):
        1. Get MM total energies for all conformers using OpenMM.
        2. Get MM torsion energies for all conformers straight from the parameters.
        3. Calculate QM[Torsion] by {QM[torsion] = QM[total] - MM[total] + MM[torsion]}.
        4. Use the Fourier transform method to fit cosine parameters to QM[torsion].
        5. Update parameters (this changes the MM[total] and MM[torsion] for subsequent torsions).

    The above protocol is repeated several times, each time the order of the torsions is shuffled.

    Args:
        config (dict): The drFrankenstein config containing all run information.
        debug (bool): If True, enables debug mode for additional logging or checks. Defaults to False.
    Returns:
        config (dict): Updated config with fitted parameters and runtime information.
    """
    ## Initialize runtimeInfo entry for stitching
    config["runtimeInfo"]["madeByStitching"] = {}

    ## Create output directories
    config = Stitching_Assistant.sort_out_directories(config)
    ## Create a basic frcmod using antechamber
    config = AMBER_helper_functions.copy_assembled_parameters(config)
    ## Update mol2 file with partial charges and run tleap to generate parameter files
    AMBER_helper_functions.edit_mol2_partial_charges(config)
    AMBER_helper_functions.run_tleap_to_make_params(config)

    ## Get options for tqdm loading bar
    tqdmBarOptions = Stitching_Assistant.init_tqdm_bar_options()

    ## Unpack config
    torsionTags = config["runtimeInfo"]["madeByTwisting"]["torsionTags"]
    maxShuffles = config["parameterFittingInfo"]["maxShuffles"]
    minShuffles = config["parameterFittingInfo"]["minShuffles"]
    maeTolTotal = config["parameterFittingInfo"]["maeTolTotal"]
    maeTolTorsion = config["parameterFittingInfo"]["maeTolTorsion"]

    ## Remove torsions that failed QM scans
    torsionTags = Stitching_Assistant.remove_exploded_torsions(config)

    ## Create a repeating list of the torsion tags, with shuffled order
    shuffledTorsionTags = Stitching_Assistant.shuffle_torsion_tags(torsionTags, maxShuffles)

    ## Initialize dictionaries to store current parameters and mean average errors
    currentParameters = {}
    meanAverageErrorTorsion = {torsionTag: [] for torsionTag in torsionTags}
    meanAverageErrorTotal = {torsionTag: [] for torsionTag in torsionTags}
    rmsMaeTorsion = []
    rmsMaeTotal = []
    counter = 1
    shuffleIndex = 1
    converged = False

    ### START OF FITTING LOOP ###
    ## Run the torsion fitting protocol, shuffling the torsion order each iteration
    for torsionTag in tqdm(shuffledTorsionTags, **tqdmBarOptions):
        ## Update the config with the current parameters, unless first iteration
        if counter != 1:
            config["runtimeInfo"]["madeByStitching"]["moleculeFrcmod"] = config["runtimeInfo"]["madeByStitching"]["proposedFrcmod"]
            AMBER_helper_functions.run_tleap_to_make_params(config)
        
        ## Use OpenMM to get MM total energies using scan geometries
        mmTotalEnergy = AMBER_total_protocol.get_MM_total_energies(config, torsionTag, debug)
        ## Extract torsion energies from parameters
        mmTorsionEnergy, mmCosineComponents = AMBER_torsion_protocol.get_MM_torsion_energies(config, torsionTag)
        ## Fit torsion parameters using Fourier Transform
        torsionParameterDf, maeTorsion, maeTotal = QMMM_fitting_protocol.fit_torsion_parameters(
            config, torsionTag, mmTotalEnergy, mmTorsionEnergy, shuffleIndex, mmCosineComponents, debug
        )

        ## Store mean average errors
        meanAverageErrorTorsion[torsionTag].append(maeTorsion)
        meanAverageErrorTotal[torsionTag].append(maeTotal)

        ## Update parameter file
        config = AMBER_helper_functions.update_frcmod(config, torsionTag, torsionParameterDf, shuffleIndex)
        ## Store current parameters
        currentParameters[torsionTag] = torsionParameterDf.to_dict(orient="records")

        ## At the end of one shuffle, check for convergence
        if counter % len(torsionTags) == 0:
            ## Calculate RMS of MAE for torsion and total fits
            rmsMaeTorsion.append(Stitching_Assistant.rms_of_mae_dict(meanAverageErrorTorsion))
            rmsMaeTotal.append(Stitching_Assistant.rms_of_mae_dict(meanAverageErrorTotal))
            if Stitching_Assistant.check_mae_convergence(
                rmsMaeTorsion[-1], rmsMaeTotal[-1], maeTolTorsion, maeTolTotal
            ) and shuffleIndex >= minShuffles:
                config["runtimeInfo"]["madeByStitching"]["shufflesCompleted"] = shuffleIndex
                converged = True
                break
            shuffleIndex += 1
        counter += 1
    ### END OF FITTING LOOP ###

    ## If the protocol does not converge, update the config with max shuffles
    if not converged:
        config["runtimeInfo"]["madeByStitching"]["shufflesCompleted"] = maxShuffles

    ## Update config with final parameters
    config["runtimeInfo"]["madeByStitching"]["moleculeFrcmod"] = config["runtimeInfo"]["madeByStitching"]["proposedFrcmod"]
    config["runtimeInfo"]["madeByStitching"]["finalParameters"] = currentParameters

    ## Run plotting protocols
    for torsionTag in torsionTags:
        fittingGif = p.join(config["runtimeInfo"]["madeByStitching"]["qmmmParameterFittingDir"], torsionTag, "torsion_fitting.gif")
        torsionFittingDir = p.join(config["runtimeInfo"]["madeByStitching"]["qmmmParameterFittingDir"], torsionTag)
        ## Create a GIF for each torsion being fitted
        Stitching_Plotter.make_gif(torsionFittingDir, fittingGif)
        ## Plot mean average errors
        Stitching_Plotter.plot_mean_average_error(torsionFittingDir, meanAverageErrorTorsion[torsionTag], meanAverageErrorTotal[torsionTag])

    ## Save mean average error data to DataFrames
    maeTorsionDf = pd.DataFrame.from_dict(meanAverageErrorTorsion)
    maeTorsionDf["All_Torsions"] = rmsMaeTorsion
    maeTotalDf = pd.DataFrame.from_dict(meanAverageErrorTotal)
    maeTotalDf["All_Torsions"] = rmsMaeTotal

    ## Save MAE data and plot run-wide errors
    maeTorsionDf.to_csv(p.join(config["runtimeInfo"]["madeByStitching"]["qmmmParameterFittingDir"], "mean_average_error.csv"), index=False)
    Stitching_Plotter.plot_run_mean_average_error(config["runtimeInfo"]["madeByStitching"]["qmmmParameterFittingDir"], maeTorsionDf, maeTotalDf)

    ## Clean up temporary files
    cleaner.clean_up_stitching(config)

    ## Update config checkpoint flag
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

    ## unpack config
    torsionTags = config["runtimeInfo"]["madeByTwisting"]["torsionTags"]
    maxShuffles = config["parameterFittingInfo"]["maxShuffles"]
    minShuffles = config["parameterFittingInfo"]["minShuffles"]
    converganceTolTotal = config["parameterFittingInfo"]["converganceTolTotal"]
    converganceTolTorsion = config["parameterFittingInfo"]["converganceTolTorsion"]

    ## create a repeating list of the torsion tags, with shuffled order
    shuffledTorsionTags = Stitching_Assistant.shuffle_torsion_tags(torsionTags, maxShuffles)
    ## get options for tqdm loading bar
    tqdmBarOptions = Stitching_Assistant.init_tqdm_bar_options()

    ## init dict to store current parameters
    currentParameters = {}
    ## init counters for shuffling
    counter = 1
    shuffleIndex = 1
    ## set up dict to store the mean average error for TOTAL andTORSION
    meanAverageErrorTorsion = {torsionTag: [] for torsionTag in torsionTags}
    meanAverageErrorTotal = {torsionTag: [] for torsionTag in torsionTags}
    rmsMaeTorsion = []
    rmsMaeTotal = []
    converged = False
        ### START OF FITTING LOOP ###
    ## run the torsion fitting protocol, each time, shuffle the torsion order
    for torsionTag in tqdm(shuffledTorsionTags, **tqdmBarOptions):
        ## update the config with the current parameters, unless first iteration (no params to update!)
        if not counter == 1:
            config["runtimeInfo"]["madeByStitching"]["moleculePrm"] = config["runtimeInfo"]["madeByStitching"]["proposedPrm"]
        ## Use OpemMM to get MM level single-points using scan geometries
        mmTotalEnergy = CHARMM_total_protocol.get_MM_total_energies(config, torsionTag, debug)
        ## Extract torsion parameters from PRM
        mmTorsionEnergy, mmCosineComponents = CHARMM_torsion_protocol.get_MM_torsion_energies(config, torsionTag, debug)
        ## fit torsion parameters using Fourier Transform
        torsionParameterDf, maeTorsion, maeTotal = QMMM_fitting_protocol.fit_torsion_parameters(config,
                                                                                                 torsionTag,
                                                                                                   mmTotalEnergy,
                                                                                                     mmTorsionEnergy,
                                                                                                       shuffleIndex,
                                                                                                         mmCosineComponents,
                                                                                                           debug)
        ## store mean average errors
        meanAverageErrorTorsion[torsionTag].append(maeTorsion)
        meanAverageErrorTotal[torsionTag].append(maeTotal)

        ## update parameter file
        config = CHARMM_helper_functions.update_prm(config, torsionTag, torsionParameterDf, shuffleIndex)
        ## store current parameters
        currentParameters[torsionTag] = torsionParameterDf.to_dict(orient = "records")
        ## At the end of one shuffle, check for convergence, if not converged update shuffleIndex
        if counter % len(torsionTags) == 0:
            ## caclulate RMS of MAE for torsion and total fits
            rmsMaeTorsion.append(Stitching_Assistant.rms_of_mae_dict(meanAverageErrorTorsion))
            rmsMaeTotal.append(Stitching_Assistant.rms_of_mae_dict(meanAverageErrorTotal))
            if Stitching_Assistant.check_mae_convergence(rmsMaeTorsion[-1], rmsMaeTotal[-1], converganceTolTorsion, converganceTolTotal) and shuffleIndex > minShuffles:
                config["runtimeInfo"]["madeByStitching"]["shufflesCompleted"] = shuffleIndex
                converged = True
                break
            shuffleIndex += 1
        counter += 1
        ### END OF FITTING LOOP ###
    ## in the case where the protocol does not converge, update the config
    if not converged:
        config["runtimeInfo"]["madeByStitching"]["shufflesCompleted"] = maxShuffles

    ## update config moleculePrm
    config["runtimeInfo"]["madeByStitching"]["moleculePrm"] = config["runtimeInfo"]["madeByStitching"]["proposedPrm"]
    config["runtimeInfo"]["madeByStitching"]["finalParameters"] = currentParameters
    ## run plotting protocols
    for torsionTag in torsionTags:
        fittingGif = p.join(config["runtimeInfo"]["madeByStitching"]["qmmmParameterFittingDir"], torsionTag, f"torsion_fitting.gif")
        torsionFittingDir = p.join(config["runtimeInfo"]["madeByStitching"]["qmmmParameterFittingDir"], torsionTag)
        ## make a gif for each torsion being fitted - so satisfying!
        Stitching_Plotter.make_gif(torsionFittingDir, fittingGif)
        ## plot mean average errors
        Stitching_Plotter.plot_mean_average_error(torsionFittingDir, meanAverageErrorTorsion[torsionTag], meanAverageErrorTotal[torsionTag])

    maeTotalDf = pd.DataFrame.from_dict(meanAverageErrorTotal)
    maeTotalDf["All_Torsions"] = rmsMaeTotal
    maeTorsionDf = pd.DataFrame.from_dict(meanAverageErrorTorsion)
    maeTorsionDf["All_Torsions"] = rmsMaeTorsion

    Stitching_Plotter.plot_run_mean_average_error(config["runtimeInfo"]["madeByStitching"]["qmmmParameterFittingDir"], maeTorsionDf, maeTotalDf)
    ## clean up 
    cleaner.clean_up_stitching(config)

    ## update config checkpoint flag
    config["checkpointInfo"]["torsionFittingComplete"] = True
    return config

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲####
if __name__ == "__main__":
    raise NotImplementedError