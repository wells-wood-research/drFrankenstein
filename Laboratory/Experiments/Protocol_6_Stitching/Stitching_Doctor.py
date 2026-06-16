## BASIC LIBRARIES ##
import os
from os import path as p
import sys
from tqdm import tqdm
import pandas as pd
from collections import defaultdict
import gc
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
from OperatingTools import Timer, cleaner, drLogger, drSplash


def _run_fitting_loop(
    config: dict,
    forcefield: str,
    helper_functions,
    total_protocol,
    torsion_protocol,
    update_param_func,
    paramFile: FilePath,
    torsionTags: list,
    shuffledTorsionTags: list,
    tqdmBarOptions: dict,
    minShuffles: int,
    converganceTolerance,
    fittingScoresCsv: FilePath,
    debug: bool,
    statusTableRows: list | None = None,
    statusTableTags: list | None = None,
    statusTableColumns: int = 4,
    statusTableNcols: int | None = None,
    displayConvergedTags: set | None = None,
    displayScores: dict | None = None,
) -> tuple[FilePath, bool]:
    """Run iterative torsion fitting until convergence or shuffle limit."""


    ## init empties for storing data, counters, and flags
    meanAverageErrorTorsion = defaultdict(list)
    meanAverageErrorTotal = defaultdict(list)
    fitScoreTorsion = defaultdict(list)
    fitScoreTotal = defaultdict(list)
    torsionMetricsByTag = {}
    totalMetricsByTag = {}
    counter = 1
    shuffleIndex = 1
    converged = False

    convergedTags = set()
    flatlinedTags = set()
    if displayConvergedTags is None:
        displayConvergedTags = set()
    if displayScores is None:
        displayScores = {}
    if not torsionTags:
        config["runtimeInfo"]["madeByStitching"]["convergedTags"] = []
        return paramFile, True

    ## Run the torsion fitting protocol, shuffling the torsion order each iteration
    for torsionTag in tqdm(shuffledTorsionTags, **tqdmBarOptions):

        ## Skip fitting for this torsion if it has already converged or flatlined
        if torsionTag in convergedTags:
            continue
        elif torsionTag in flatlinedTags:
            continue


        ## Update the config with the current parameters, unless first iteration
        if counter > 1:
            if forcefield == "AMBER":
                helper_functions.run_tleap_to_make_params(paramFile, config)

        ## RESET TORSION IN QUESTION TO ZERO BEFORE FITTING TO AVOID PARAMETER INTERFERENCE
        flattenedTorsionDf =  pd.DataFrame({"Amplitude": [0.0], "Phase": [180], "Period": [1.0]})
        paramFile = update_param_func(paramFile, config, torsionTag, flattenedTorsionDf, shuffleIndex)
        if forcefield == "AMBER":
            helper_functions.run_tleap_to_make_params(paramFile, config)

        ## Use OpenMM to get MM total energies using scan geometries
        mmTotalEnergy = total_protocol.get_MM_total_energies(config, torsionTag, paramFile, debug)
        ## Extract torsion energies from parameters
        mmTorsionEnergy, mmCosineComponents = torsion_protocol.get_MM_torsion_energies(config, torsionTag, paramFile)
        ## Fit torsion parameters using Fourier Transform
        torsionParameterDf, torsionMetrics, totalMetrics, maeTorsion, maeTotal, torsionConverged = QMMM_fitting_protocol.fit_torsion_parameters(
            config, torsionTag, mmTotalEnergy, mmTorsionEnergy, shuffleIndex, mmCosineComponents, debug
        )                

        ## Store mean average errors for the current shuffle
        meanAverageErrorTorsion[torsionTag].append(maeTorsion)
        meanAverageErrorTotal[torsionTag].append(maeTotal)
        fitScoreTorsion[torsionTag].append(torsionMetrics["composite_score"])
        fitScoreTotal[torsionTag].append(totalMetrics["composite_score"])
        torsionMetricsByTag[torsionTag] = torsionMetrics
        totalMetricsByTag[torsionTag] = totalMetrics
        displayScores[torsionTag] = torsionMetrics["composite_score"]

        # Record the composite score and all subcomponents immediately after this fit
        try:
            with open(fittingScoresCsv, "a") as f:
                f.write(
                    f"{shuffleIndex},{config['runtimeInfo']['madeByStitching']['maxTorsions']},{torsionTag},"
                    f"{maeTorsion},{maeTotal},{torsionMetrics['composite_score']},{totalMetrics['composite_score']},"
                    f"{torsionMetrics['location_score']},{torsionMetrics['amplitude_score']},"
                    f"{torsionMetrics['stationary_count_score']},{torsionMetrics['normalized_mae_score']},"
                    f"{totalMetrics['location_score']},{totalMetrics['amplitude_score']},"
                    f"{totalMetrics['stationary_count_score']},{totalMetrics['normalized_mae_score']}\n"
                )
        except Exception:
            # Don't let logging failures break the fitting loop
            pass

        ## Update parameter file
        paramFile = update_param_func(paramFile, config, torsionTag, torsionParameterDf, shuffleIndex)

        ## Freeze torsions as soon as the torsion score is within tolerance.
        if torsionConverged:
            convergedTags.add(torsionTag)
            displayConvergedTags.add(torsionTag)
        if statusTableRows and statusTableTags:
            drSplash.update_torsion_status_table(
                statusTableRows,
                statusTableTags,
                displayConvergedTags,
                scores=displayScores,
                columns=statusTableColumns,
                ncols=statusTableNcols,
            )

        ##################### END OF SHUFFLE ####################
        ## At the end of one shuffle, check for convergence
        if counter % len(torsionTags) == 0:
            ## Calculate RMS of MAE for torsion and total fits
            rmsMaeTorsion = Stitching_Assistant.rms_of_mae_dict(meanAverageErrorTorsion)
            rmsMaeTotal = Stitching_Assistant.rms_of_mae_dict(meanAverageErrorTotal)
            rmsScoreTorsion = Stitching_Assistant.rms_of_mae_dict(fitScoreTorsion)
            rmsScoreTotal = Stitching_Assistant.rms_of_mae_dict(fitScoreTotal)

            # Write only the All_Torsions summary line to CSV (per-fit rows are written immediately after each fit)
            tagCount = max(len(torsionMetricsByTag), 1)
            summaryTorsionLocation = sum(m["location_score"] for m in torsionMetricsByTag.values()) / tagCount
            summaryTorsionAmplitude = sum(m["amplitude_score"] for m in torsionMetricsByTag.values()) / tagCount
            summaryTorsionCount = sum(m["stationary_count_score"] for m in torsionMetricsByTag.values()) / tagCount
            summaryTorsionNormMae = sum(m["normalized_mae_score"] for m in torsionMetricsByTag.values()) / tagCount
            summaryTotalLocation = sum(m["location_score"] for m in totalMetricsByTag.values()) / tagCount
            summaryTotalAmplitude = sum(m["amplitude_score"] for m in totalMetricsByTag.values()) / tagCount
            summaryTotalCount = sum(m["stationary_count_score"] for m in totalMetricsByTag.values()) / tagCount
            summaryTotalNormMae = sum(m["normalized_mae_score"] for m in totalMetricsByTag.values()) / tagCount
            try:
                with open(fittingScoresCsv, "a") as f:
                    f.write(
                        f"{shuffleIndex},{config['runtimeInfo']['madeByStitching']['maxTorsions']},All_Torsions,"
                        f"{rmsMaeTorsion},{rmsMaeTotal},{rmsScoreTorsion},{rmsScoreTotal},"
                        f"{summaryTorsionLocation},{summaryTorsionAmplitude},{summaryTorsionCount},{summaryTorsionNormMae},"
                        f"{summaryTotalLocation},{summaryTotalAmplitude},{summaryTotalCount},{summaryTotalNormMae}\n"
                    )
            except Exception:
                # Don't let logging failures break the fitting loop
                pass

            ## every 3 shuffles, check to see if composite scores are flatlined
            if shuffleIndex % 3 == 0 and shuffleIndex >= minShuffles:
                thisRoundFlatlinedTorsions = Stitching_Assistant.check_scores_for_flatline(torsionTags, convergedTags, fittingScoresCsv)
                if len(thisRoundFlatlinedTorsions) > 0:
                   flatlinedTags.update(thisRoundFlatlinedTorsions)


            ## at the end of each shuffle, clear memory
            meanAverageErrorTorsion.clear()
            meanAverageErrorTotal.clear()
            fitScoreTorsion.clear()
            fitScoreTotal.clear()
            torsionMetricsByTag.clear()
            totalMetricsByTag.clear()
            del rmsMaeTorsion
            del rmsMaeTotal
            del rmsScoreTorsion
            del rmsScoreTotal
            gc.collect()

            shuffleIndex += 1
        counter += 1

    config["runtimeInfo"]["madeByStitching"]["convergedTags"] = sorted(convergedTags)
    return paramFile, converged


# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
@drLogger.experiment_logger("Torsion Parameter Fitting")
@Timer.time_function("Parameter Fitting", "PARAMETER FITTING")
def torsion_fitting_protocol(config: dict, debug: bool = False) -> dict:
    """Run the full torsion fitting protocol for the chosen force field."""

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
    seed = config["miscInfo"]["seed"]
    # Standardize on maeTol keys, assuming this is the consistent naming in the config
    converganceTolerance = config["parameterFittingInfo"].get("converganceTolerance", None)

    ## duplicate maxTorsions to runtimeInfo for easy access
    config["runtimeInfo"]["madeByStitching"]["maxTorsions"] = config["parameterFittingInfo"]["maxCosineFunctions"]


    ## Remove torsions that failed QM scans and shuffle the rest
    allTorsionTags = Stitching_Assistant.remove_exploded_torsions(config)
    torsionTags = list(allTorsionTags)
    shuffledTorsionTags = Stitching_Assistant.shuffle_torsion_tags(torsionTags, maxShuffles, seed)
    ## Get options for tqdm loading bar and initialize containers
    tqdmBarOptions = Stitching_Assistant.init_tqdm_bar_options()
    statusLine = None
    bannerLine = None
    statusTableRows = None
    statusTableTags = list(allTorsionTags)
    statusTableColumns = 4
    statusTableNcols = tqdmBarOptions.get("ncols", 102)
    displayConvergedTags: set[str] = set()
    displayScores: dict[str, float] = {}
    
    if not tqdmBarOptions.get("disable", False):
        statusLine = tqdm(
            total=1,
            initial=1,
            desc="",
            bar_format="{desc}",
            position=1,
            leave=True,
            ncols=statusTableNcols,
            mininterval=0,
        )
        statusLine.refresh()
        drSplash.set_status_line(statusLine)
    tableStartPosition = 2 if statusLine else 1
    if statusTableTags and not tqdmBarOptions.get("disable", False):
        bannerLine = tqdm(
            total=1,
            initial=1,
            desc=drSplash.build_torsion_status_banner(),
            bar_format="{desc}",
            position=tableStartPosition,
            leave=True,
            ncols=statusTableNcols,
            mininterval=0,
        )
        bannerLine.refresh()
        tableStartPosition += 1
        statusTableRows = drSplash.init_torsion_status_table(
            statusTableTags,
            convergedTags=displayConvergedTags,
            scores=displayScores,
            columns=statusTableColumns,
            ncols=statusTableNcols,
            position=tableStartPosition,
            tol = config["parameterFittingInfo"].get("converganceTolerance", 0.1)
        )

    ## Set up Mean Average Error CSV for memory-efficient logging
    fittingScoresCsv = p.join(config["runtimeInfo"]["madeByStitching"]["qmmmParameterFittingDir"], "mean_average_errors.csv")
    config["runtimeInfo"]["madeByStitching"]["fittingScoresCsv"] = fittingScoresCsv
    with open(fittingScoresCsv, "w") as f:
        f.write(
            "shuffle,n_cosines,torsion_tag,mae_torsion,mae_total,fit_score_torsion,fit_score_total,"
            "torsion_location_score,torsion_amplitude_score,torsion_stationary_count_score,torsion_normalized_mae_score,"
            "total_location_score,total_amplitude_score,total_stationary_count_score,total_normalized_mae_score\n"
        )


    cosineScheme = range(2, config["parameterFittingInfo"]["maxCosineFunctions"]+1)
    cosineScheme = [config["parameterFittingInfo"]["maxCosineFunctions"]]
    for nCosines in cosineScheme:
        config["runtimeInfo"]["madeByStitching"]["maxTorsions"] = nCosines
        paramFile, converged = _run_fitting_loop(
            config=config,
            forcefield=forcefield,
            helper_functions=helper_functions,
            total_protocol=total_protocol,
            torsion_protocol=torsion_protocol,
            update_param_func=update_param_func,
            paramFile=paramFile,
            torsionTags=torsionTags,
            shuffledTorsionTags=shuffledTorsionTags,
            tqdmBarOptions=tqdmBarOptions,
            minShuffles=minShuffles,
            converganceTolerance=converganceTolerance,
            fittingScoresCsv=fittingScoresCsv,
            debug=debug,
            statusTableRows=statusTableRows,
            statusTableTags=statusTableTags,
            statusTableColumns=statusTableColumns,
            statusTableNcols=statusTableNcols,
            displayConvergedTags=displayConvergedTags,
            displayScores=displayScores,
        )

        convergedTags = set(config["runtimeInfo"]["madeByStitching"].get("convergedTags", []))
        torsionTags = [tag for tag in torsionTags if tag not in convergedTags]
        # torsionTags = Stitching_Assistant._remove_converged_torsion_tags(fittingScoresCsv, converganceTolerance, torsionTags, topN = 2)
        shuffledTorsionTags = Stitching_Assistant.shuffle_torsion_tags(torsionTags, maxShuffles, seed) if torsionTags else []

        print("Remaining torsions after convergence check:", torsionTags)
        ## run again with the remaining torsions to try to get more to converge, but only if some have not converged and we have not exceeded max shuffles

        if not torsionTags:
            break

    ## If the protocol does not converge, update the config with max shuffles
    if not converged:
        config["runtimeInfo"]["madeByStitching"]["shufflesCompleted"] = maxShuffles

    if statusTableRows:
        drSplash.close_torsion_status_table(statusTableRows)
    if bannerLine:
        bannerLine.close()
    if statusLine:
        drSplash.clear_status_line()
        statusLine.close()

    ## Update config with final parameters
    config["runtimeInfo"]["madeByStitching"][param_key] = paramFile

    ##TODO: construct final parameters
    config["runtimeInfo"]["madeByStitching"]["finalParameters"] = Stitching_Assistant.construct_final_params(config, paramFile, param_extraction_function, allTorsionTags)
    

    ## Run plotting protocols
    for torsionTag in allTorsionTags:
        fittingGif = p.join(config["runtimeInfo"]["madeByStitching"]["qmmmParameterFittingDir"], torsionTag, "torsion_fitting.gif")
        torsionFittingDir = p.join(config["runtimeInfo"]["madeByStitching"]["qmmmParameterFittingDir"], torsionTag)
        ## Create a GIF for each torsion being fitted
        Stitching_Plotter.make_gif(torsionFittingDir, fittingGif)
        ## Plot mean average errors from the CSV data
        # Stitching_Plotter.plot_mean_average_error(torsionFittingDir, fittingScoresCsv, torsionTag)

    ## Load MAE data from CSV for final run-wide analysis and plotting
    if os.path.exists(fittingScoresCsv):
        Stitching_Plotter.plot_run_mean_average_error(config["runtimeInfo"]["madeByStitching"]["qmmmParameterFittingDir"], fittingScoresCsv )
    ## Clean up temporary files
    cleaner.clean_up_stitching(config)
    ## Update config checkpoint flag
    config["checkpointInfo"]["torsionFittingComplete"] = True
    return config





# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲####
if __name__ == "__main__":
    raise NotImplementedError
