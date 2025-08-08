## BASIC LIBRARIES ##
import os
from os import path as p
import sys
import random
from tqdm import tqdm
import pandas as pd
from collections import defaultdict

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
def torsion_fitting_protocol(config: dict, debug=False) -> dict:
    """
    Main protocol for torsion fitting for AMBER or CHARMM parameters.
    For each torsion that has had QM scans performed (creates QM[total]):
        1. Get MM total energies for all conformers using OpenMM.
        2. Get MM torsion energies for all conformers straight from the parameters.
        3. Calculate QM[Torsion] by {QM[torsion] = QM[total] - MM[total] + MM[torsion]}.
        4. Use the Fourier transform method to fit cosine parameters to QM[torsion].
        5. Update parameters (this changes the MM[total] and MM[torsion] for subsequent torsions).

    The above protocol is repeated several times, each time the order of the torsions is shuffled.
    Memory usage is optimized by writing MAE results to a CSV file after each shuffle and clearing them from memory.

    Args:
        config (dict): The drFrankenstein config containing all run information.
        forcefield (str): The force field being used, either "AMBER" or "CHARMM".
        debug (bool): If True, enables debug mode for additional logging or checks. Defaults to False.
    Returns:
        config (dict): Updated config with fitted parameters and runtime information.
    """

    forcefield = config["parameterFittingInfo"]["forceField"]
    # Set up forcefield-specific helpers and keys
    if forcefield == "AMBER":
        helper_functions = AMBER_helper_functions
        total_protocol = AMBER_total_protocol
        torsion_protocol = AMBER_torsion_protocol
        param_key = "moleculeFrcmod"
        proposed_param_key = "proposedFrcmod"
        update_param_func = helper_functions.update_frcmod
    else: # CHARMM
        helper_functions = CHARMM_helper_functions
        total_protocol = CHARMM_total_protocol
        torsion_protocol = CHARMM_torsion_protocol
        param_key = "moleculePrm"
        proposed_param_key = "proposedPrm"
        update_param_func = helper_functions.update_prm

    ## Initialize runtimeInfo entry for stitching
    config["runtimeInfo"]["madeByStitching"] = {}

    ## Create output directories
    config = Stitching_Assistant.sort_out_directories(config)
    
    ## Create initial parameter files
    config = helper_functions.copy_assembled_parameters(config)
    if forcefield == "AMBER":
        helper_functions.edit_mol2_partial_charges(config)
        helper_functions.run_tleap_to_make_params(config)

    ## Unpack config
    torsionTags = config["runtimeInfo"]["madeByTwisting"]["torsionTags"]
    maxShuffles = config["parameterFittingInfo"]["maxShuffles"]
    minShuffles = config["parameterFittingInfo"]["minShuffles"]
    # Standardize on maeTol keys, assuming this is the consistent naming in the config
    maeTolTotal = config["parameterFittingInfo"].get("maeTolTotal", config["parameterFittingInfo"].get("converganceTolTotal"))
    maeTolTorsion = config["parameterFittingInfo"].get("maeTolTorsion", config["parameterFittingInfo"].get("converganceTolTorsion"))


    ## Remove torsions that failed QM scans and shuffle the rest
    torsionTags = Stitching_Assistant.remove_exploded_torsions(config)
    shuffledTorsionTags = Stitching_Assistant.shuffle_torsion_tags(torsionTags, maxShuffles)

    ## Get options for tqdm loading bar and initialize containers
    tqdmBarOptions = Stitching_Assistant.init_tqdm_bar_options()
    currentParameters = {}
    meanAverageErrorTorsion = defaultdict(list)
    meanAverageErrorTotal = defaultdict(list)
    rmsMaeTorsion = []
    rmsMaeTotal = []
    counter = 1
    shuffleIndex = 1
    converged = False

    ## Set up Mean Average Error CSV for memory-efficient logging
    maeCsv = p.join(config["runtimeInfo"]["madeByStitching"]["qmmmParameterFittingDir"], "mean_average_errors.csv")
    config["runtimeInfo"]["madeByStitching"]["maeCsv"] = maeCsv
    with open(maeCsv, "w") as f:
        f.write("shuffle,torsion_tag,mae_torsion,mae_total\n")

    ### START OF FITTING LOOP ###
    ## Run the torsion fitting protocol, shuffling the torsion order each iteration
    for torsionTag in tqdm(shuffledTorsionTags, **tqdmBarOptions):
        ## Update the config with the current parameters, unless first iteration
        if counter > 1:
            config["runtimeInfo"]["madeByStitching"][param_key] = config["runtimeInfo"]["madeByStitching"][proposed_param_key]
            if forcefield == "AMBER":
                helper_functions.run_tleap_to_make_params(config)
        
        ## Use OpenMM to get MM total energies using scan geometries
        mmTotalEnergy = total_protocol.get_MM_total_energies(config, torsionTag, debug)
        ## Extract torsion energies from parameters
        mmTorsionEnergy, mmCosineComponents = torsion_protocol.get_MM_torsion_energies(config, torsionTag, debug)
        ## Fit torsion parameters using Fourier Transform
        torsionParameterDf, maeTorsion, maeTotal = QMMM_fitting_protocol.fit_torsion_parameters(
            config, torsionTag, mmTotalEnergy, mmTorsionEnergy, shuffleIndex, mmCosineComponents, debug
        )

        ## Store mean average errors for the current shuffle
        meanAverageErrorTorsion[torsionTag].append(maeTorsion)
        meanAverageErrorTotal[torsionTag].append(maeTotal)

        ## Update parameter file and store current parameters
        config = update_param_func(config, torsionTag, torsionParameterDf, shuffleIndex)
        currentParameters[torsionTag] = torsionParameterDf.to_dict(orient="records")

        ## At the end of one shuffle, check for convergence
        if counter % len(torsionTags) == 0:
            ## Calculate RMS of MAE for torsion and total fits
            rmsMaeTorsion.append(Stitching_Assistant.rms_of_mae_dict(meanAverageErrorTorsion))
            rmsMaeTotal.append(Stitching_Assistant.rms_of_mae_dict(meanAverageErrorTotal))

            ## Write buffer to CSV file and clear it to save memory
            with open(maeCsv, "a") as f:
                for tag in meanAverageErrorTorsion:
                    mae_t = meanAverageErrorTorsion[tag][-1]
                    mae_tot = meanAverageErrorTotal[tag][-1]
                    f.write(f"{shuffleIndex},{tag},{mae_t},{mae_tot}\n")
            
            meanAverageErrorTorsion.clear()
            meanAverageErrorTotal.clear()
            
            ## Check for convergence
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
    config["runtimeInfo"]["madeByStitching"][param_key] = config["runtimeInfo"]["madeByStitching"][proposed_param_key]
    config["runtimeInfo"]["madeByStitching"]["finalParameters"] = currentParameters

    ## Run plotting protocols
    for torsionTag in torsionTags:
        fittingGif = p.join(config["runtimeInfo"]["madeByStitching"]["qmmmParameterFittingDir"], torsionTag, "torsion_fitting.gif")
        torsionFittingDir = p.join(config["runtimeInfo"]["madeByStitching"]["qmmmParameterFittingDir"], torsionTag)
        ## Create a GIF for each torsion being fitted
        Stitching_Plotter.make_gif(torsionFittingDir, fittingGif)
        ## Plot mean average errors from the CSV data
        Stitching_Plotter.plot_mean_average_error(torsionFittingDir, maeCsv, torsionTag)

    ## Load MAE data from CSV for final run-wide analysis and plotting
    if os.path.exists(maeCsv):
        full_mae_df = pd.read_csv(maeCsv)
        maeTorsionDf = full_mae_df.pivot(index='shuffle', columns='torsion_tag', values='mae_torsion').reset_index(drop=True)
        if len(rmsMaeTorsion) == len(maeTorsionDf):
             maeTorsionDf["All_Torsions"] = rmsMaeTorsion
       
        maeTotalDf = full_mae_df.pivot(index='shuffle', columns='torsion_tag', values='mae_total').reset_index(drop=True)
        if len(rmsMaeTotal) == len(maeTotalDf):
            maeTotalDf["All_Torsions"] = rmsMaeTotal

        Stitching_Plotter.plot_run_mean_average_error(config["runtimeInfo"]["madeByStitching"]["qmmmParameterFittingDir"], maeTorsionDf, maeTotalDf)

    ## Clean up temporary files
    cleaner.clean_up_stitching(config)

    ## Update config checkpoint flag
    config["checkpointInfo"]["torsionFittingComplete"] = True
    return config

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲####
if __name__ == "__main__":
    raise NotImplementedError