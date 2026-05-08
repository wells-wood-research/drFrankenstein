## BASIC LIBRARIES ##
import os
from os import path as p
import sys
from tqdm import tqdm
import pandas as pd
from collections import defaultdict
import gc
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
from OperatingTools import Timer, cleaner, drLogger

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
@drLogger.experiment_logger("Torsion Parameter Fitting")
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
        param_extraction_function = AMBER_torsion_protocol.extract_torsion_parameters_from_frcmod
        param_key = "moleculeFrcmod"
        update_param_func = helper_functions.update_frcmod
    else: # CHARMM
        helper_functions = CHARMM_helper_functions
        total_protocol = CHARMM_total_protocol
        torsion_protocol = CHARMM_torsion_protocol
        param_extraction_function = CHARMM_torsion_protocol.extract_torsion_parameters_from_prm
        param_key = "moleculePrm"
        update_param_func = helper_functions.update_prm
        paramFile = config["runtimeInfo"]["madeByAssembly"]["assembledPrm"]

    ## Initialize runtimeInfo entry for stitching
    config["runtimeInfo"]["madeByStitching"] = {}

    ## Create output directories
    config = Stitching_Assistant.sort_out_directories(config)
    
    ## Create initial parameter files ## TODO: move some of this to Assembly
    config = helper_functions.copy_assembled_parameters(config)
    if forcefield == "AMBER":
        helper_functions.edit_mol2_partial_charges(config)
        paramFile = config["runtimeInfo"]["madeByAssembly"]["assembledFrcmod"]
        helper_functions.run_tleap_to_make_params(paramFile, config)
    elif forcefield == "CHARMM":
        helper_functions.edit_rtf_charges(config)

        
    ## Unpack config
    torsionTags = config["runtimeInfo"]["madeByTwisting"]["torsionTags"]
    maxShuffles = config["parameterFittingInfo"]["maxShuffles"]
    minShuffles = config["parameterFittingInfo"]["minShuffles"]
    maxCosineFunctions = config["parameterFittingInfo"]["maxCosineFunctions"]
    fittingProtocol = config["parameterFittingInfo"].get("fittingProtocol", "ROBUST_SPARSE_HARMONIC")
    progressiveCosineFitting = config["parameterFittingInfo"].get("progressiveCosineFitting", False)
    freezeConvergedTorsions = config["parameterFittingInfo"].get("freezeConvergedTorsions", False)
    seed = config["miscInfo"]["seed"]
    converganceTolerance = config["parameterFittingInfo"].get("converganceTolerance", None)
    torsionConvergenceTolerance = config["parameterFittingInfo"].get("torsionConvergenceTolerance", converganceTolerance)


    ## Remove torsions that failed QM scans and shuffle the rest
    torsionTags = Stitching_Assistant.remove_exploded_torsions(config)
    ## Get options for tqdm loading bar and initialize containers
    tqdmBarOptions = Stitching_Assistant.init_tqdm_bar_options()

    ## Set up Mean Average Error CSV for memory-efficient logging
    maeCsv = p.join(config["runtimeInfo"]["madeByStitching"]["qmmmParameterFittingDir"], "mean_average_errors.csv")
    config["runtimeInfo"]["madeByStitching"]["maeCsv"] = maeCsv
    with open(maeCsv, "w") as f:
        f.write("shuffle,torsion_tag,mae_torsion,mae_total\n")


    activeTorsionTags = list(torsionTags)
    frozenTorsionTags = []
    totalShufflesCompleted = 0
    firstFittingStep = True

    globalJointFitting = fittingProtocol == "GLOBAL_JOINT_HARMONIC"

    # OVERFIT_PRUNE_HARMONIC and GLOBAL_JOINT_HARMONIC already avoid
    # order-dependent staged growth in their own design.
    if fittingProtocol in ["OVERFIT_PRUNE_HARMONIC", "GLOBAL_JOINT_HARMONIC"]:
        progressiveCosineFitting = False
    if globalJointFitting:
        # Global joint mode is a single global round (all torsions fit from
        # the same parameter state, then all updates applied together).
        maxShuffles = 1
        minShuffles = 1

    if progressiveCosineFitting:
        cosineSchedule = list(range(1, maxCosineFunctions + 1))
    else:
        cosineSchedule = [maxCosineFunctions]

    for cosineLimit in cosineSchedule:
        if not activeTorsionTags:
            break

        stageConverged = False
        for stageShuffleIndex in tqdm(range(1, maxShuffles + 1), **tqdmBarOptions):
            if not activeTorsionTags:
                stageConverged = True
                break

            meanAverageErrorTorsion = defaultdict(list)
            meanAverageErrorTotal = defaultdict(list)

            globalShuffleIndex = totalShufflesCompleted + stageShuffleIndex
            workingTorsionTags = list(activeTorsionTags)
            if not globalJointFitting:
                derivedSeed = seed * globalShuffleIndex if seed != 0 else globalShuffleIndex
                random.Random(derivedSeed).shuffle(workingTorsionTags)
            resetAppliedByTorsion = {}
            pendingTorsionUpdates = {}

            if forcefield == "AMBER" and not firstFittingStep:
                helper_functions.run_tleap_to_make_params(paramFile, config)

            for torsionTag in workingTorsionTags:
                if forcefield == "AMBER" and not firstFittingStep and not globalJointFitting:
                    helper_functions.run_tleap_to_make_params(paramFile, config)

                mmTotalEnergy = total_protocol.get_MM_total_energies(config, torsionTag, paramFile, debug)
                mmTorsionEnergy, mmCosineComponents = torsion_protocol.get_MM_torsion_energies(config, torsionTag, paramFile)
                torsionParameterDf, maeTorsion, maeTotal, resetApplied = QMMM_fitting_protocol.fit_torsion_parameters(
                    config,
                    torsionTag,
                    mmTotalEnergy,
                    mmTorsionEnergy,
                    globalShuffleIndex,
                    mmCosineComponents,
                    maxCosineFunctions=cosineLimit,
                    debug=debug,
                )

                meanAverageErrorTorsion[torsionTag].append(maeTorsion)
                meanAverageErrorTotal[torsionTag].append(maeTotal)
                resetAppliedByTorsion[torsionTag] = resetApplied
                if globalJointFitting:
                    pendingTorsionUpdates[torsionTag] = torsionParameterDf
                else:
                    paramFile = update_param_func(paramFile, config, torsionTag, torsionParameterDf, globalShuffleIndex)
                    firstFittingStep = False

            if globalJointFitting:
                for torsionTag in workingTorsionTags:
                    torsionParameterDf = pendingTorsionUpdates[torsionTag]
                    paramFile = update_param_func(paramFile, config, torsionTag, torsionParameterDf, globalShuffleIndex)
                if len(workingTorsionTags) > 0:
                    firstFittingStep = False

            rmsMaeTorsion = Stitching_Assistant.rms_of_mae_dict(meanAverageErrorTorsion)
            rmsMaeTotal = Stitching_Assistant.rms_of_mae_dict(meanAverageErrorTotal)

            with open(maeCsv, "a") as f:
                for tag in meanAverageErrorTorsion:
                    maeTorsionTag = meanAverageErrorTorsion[tag][-1]
                    maeTotalTag = meanAverageErrorTotal[tag][-1]
                    f.write(f"{globalShuffleIndex},{tag},{maeTorsionTag},{maeTotalTag}\n")
                f.write(f"{globalShuffleIndex},All_Torsions,{rmsMaeTorsion},{rmsMaeTotal}\n")

            if (
                freezeConvergedTorsions
                and torsionConvergenceTolerance is not None
                and stageShuffleIndex >= minShuffles
            ):
                for torsionTag in list(activeTorsionTags):
                    latestMae = meanAverageErrorTorsion.get(torsionTag, [None])[-1]
                    if resetAppliedByTorsion.get(torsionTag, False):
                        continue
                    if latestMae is not None and latestMae < torsionConvergenceTolerance:
                        activeTorsionTags.remove(torsionTag)
                        if torsionTag not in frozenTorsionTags:
                            frozenTorsionTags.append(torsionTag)

            if not activeTorsionTags:
                stageConverged = True
                totalShufflesCompleted = globalShuffleIndex
                break

            if (
                Stitching_Assistant.check_mae_convergence(
                    rmsMaeTorsion, rmsMaeTotal, converganceTolerance
                )
                and stageShuffleIndex >= minShuffles
            ):
                stageConverged = True
                totalShufflesCompleted = globalShuffleIndex
                break

            if stageShuffleIndex % 10 == 0 and stageShuffleIndex >= minShuffles:
                if Stitching_Assistant.check_mae_flatline(maeCsv):
                    stageConverged = True
                    totalShufflesCompleted = globalShuffleIndex
                    break

            totalShufflesCompleted = globalShuffleIndex

            meanAverageErrorTorsion.clear()
            meanAverageErrorTotal.clear()
            del rmsMaeTorsion
            del rmsMaeTotal
            gc.collect()

        if not stageConverged:
            break
        if not progressiveCosineFitting:
            break

    config["runtimeInfo"]["madeByStitching"]["shufflesCompleted"] = totalShufflesCompleted
    config["runtimeInfo"]["madeByStitching"]["frozenTorsionTags"] = frozenTorsionTags

    ## Update config with final parameters
    config["runtimeInfo"]["madeByStitching"][param_key] = paramFile

    ##TODO: construct final parameters
    config["runtimeInfo"]["madeByStitching"]["finalParameters"] = Stitching_Assistant.construct_final_params(config, paramFile, param_extraction_function, torsionTags)
    

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
        Stitching_Plotter.plot_run_mean_average_error(config["runtimeInfo"]["madeByStitching"]["qmmmParameterFittingDir"], maeCsv )
    ## Clean up temporary files
    cleaner.clean_up_stitching(config)
    ## Update config checkpoint flag
    config["checkpointInfo"]["torsionFittingComplete"] = True
    return config





# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲####
if __name__ == "__main__":
    raise NotImplementedError
