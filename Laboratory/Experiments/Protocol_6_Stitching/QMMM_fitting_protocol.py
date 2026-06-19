import os
from os import path as p
import numpy as np
import pandas as pd
import sys
from scipy.signal import savgol_filter

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
from . import drFourier
from . import Stitching_Plotter
from . import Stitching_Assistant

def _resolve_savgol_settings(smoothingConfig: object, signalLength: int) -> dict | None:
    """Resolve and validate Savitzky-Golay settings from config."""
    if smoothingConfig in (None, False):
        return None

    if smoothingConfig is True:
        settings = {"window_length": 9, "polyorder": 2, "mode": "interp"}
    elif isinstance(smoothingConfig, dict):
        settings = {
            "window_length": int(smoothingConfig.get("windowLength", 9)),
            "polyorder": int(smoothingConfig.get("polyorder", 2)),
            "mode": str(smoothingConfig.get("mode", "interp")),
        }
    else:
        raise TypeError("parameterFittingInfo.sagvolSmoothing must be bool, dict, or null.")

    windowLength = settings["window_length"]
    if windowLength % 2 == 0:
        windowLength += 1
    if signalLength % 2 == 0 and windowLength >= signalLength:
        windowLength = signalLength - 1
    elif windowLength > signalLength:
        windowLength = signalLength
    if windowLength < 3:
        return None
    if settings["polyorder"] >= windowLength:
        settings["polyorder"] = windowLength - 1
    if settings["polyorder"] < 1:
        return None

    settings["window_length"] = windowLength
    return settings

##############################################################################
def fit_torsion_parameters(config: dict,
                            torsionTag: str,
                              mmTotalEnergy: np.ndarray,
                                mmTorsionEnergy: np.ndarray,
                                  shuffleIndex: int,
                                    mmCosineComponents: dict,
                                      debug: bool = False):
    """Fit torsion parameters to QM data and return fit metrics."""
    ## unpack config
    qmmmFittingDir = config["runtimeInfo"]["madeByStitching"]["qmmmParameterFittingDir"] 
    sagvolSmoothing = config["parameterFittingInfo"]["sagvolSmoothing"]
    converganceTolerance = config["parameterFittingInfo"].get("converganceTolerance", None)

    ## make the dir if it doesn't exist
    qmmmTorsionFittingDir = p.join(qmmmFittingDir, torsionTag)
    os.makedirs(qmmmTorsionFittingDir,exist_ok=True)

    ## retrieve QM total energies from config (from torsion scanning)
    qmTotalEnergy = get_qm_scan_energies(config, torsionTag)
    ## normalise QM total energy
    qmTotalEnergy = qmTotalEnergy - qmTotalEnergy.min()
    qmTotalEnergyRaw = qmTotalEnergy.copy()
    mmTotalEnergy = mmTotalEnergy - mmTotalEnergy.min()
    mmTotalEnergyRaw = mmTotalEnergy.copy()

    ## smooth QM and MM total profiles with the same settings before deriving torsion profile
    savgolSettings = _resolve_savgol_settings(sagvolSmoothing, signalLength=len(qmTotalEnergy))
    if savgolSettings is not None:
        qmTotalEnergy = savgol_filter(qmTotalEnergy, **savgolSettings)
        qmTotalEnergy = qmTotalEnergy - qmTotalEnergy.min()
        mmTotalEnergy = savgol_filter(mmTotalEnergy, **savgolSettings)
        mmTotalEnergy = mmTotalEnergy - mmTotalEnergy.min()

    ## calculate QM torsion energy to fit parameters to 
    qmTorsionEnergy = qmTotalEnergy - mmTotalEnergy + mmTorsionEnergy
    ## normalise QM Torsion energy
    qmTorsionEnergy = qmTorsionEnergy - qmTorsionEnergy.min()

    torsionParametersDf, cosineComponents, reconstructedSignal, meanAverageErrorTorsion, fitScoreTorsion = drFourier.fourier_transform_protocol(
        qmTorsionEnergy,
        torsionTag,
        qmmmTorsionFittingDir,
        config=config,
        qmTotalEnergy=qmTotalEnergy,
        mmTotalEnergy=mmTotalEnergy,
        mmTorsionEnergy=mmTorsionEnergy,
    )
    mmFittedTotalEnergy = mmTotalEnergy - mmTorsionEnergy + reconstructedSignal
    mmFittedTotalEnergy = mmFittedTotalEnergy - mmFittedTotalEnergy.min()
    meanAverageErrorTotal = np.mean(np.abs(qmTotalEnergy - mmFittedTotalEnergy))
    totalMetrics = Stitching_Assistant.calculate_profile_fit_metrics(qmTotalEnergy, mmFittedTotalEnergy)
    torsionMetrics = Stitching_Assistant.calculate_profile_fit_metrics(qmTorsionEnergy, reconstructedSignal)
    torsionConverged = Stitching_Assistant.check_torsion_convergence(torsionMetrics["composite_score"], totalMetrics["composite_score"], converganceTolerance)

    Stitching_Plotter.plot_qmmm_energies(qmTotalEnergy,
                                          qmTorsionEnergy,
                                            mmFittedTotalEnergy,
                                              reconstructedSignal,
                                                cosineComponents,
                                                  torsionMetrics,
                                                    totalMetrics,
                                                    meanAverageErrorTorsion,
                                                    meanAverageErrorTotal,
                                                    config["runtimeInfo"]["madeByStitching"]["maxTorsions"],
                                                    torsionConverged,
                                                    qmmmTorsionFittingDir,
                                                      shuffleIndex,
                                                      tol=converganceTolerance,
                                                      qmTotalEnergyRaw=qmTotalEnergyRaw,
                                                      mmTotalEnergyRaw=mmTotalEnergyRaw)

    return torsionParametersDf, torsionMetrics, totalMetrics, meanAverageErrorTorsion, meanAverageErrorTotal, torsionConverged

##############################################################################
def get_qm_scan_energies(config: dict, torsionTag: str) -> np.array:
    """Load the QM scan energies for a torsion tag."""
    qmScanEnergyCsv = config["runtimeInfo"]["madeByTwisting"]["finalScanEnergies"][torsionTag]
    if not p.exists(qmScanEnergyCsv):
        raise FileNotFoundError(f"Couldn't find QM scan energies for torsion {torsionTag}.")

    qmScanEnergyDf = pd.read_csv(qmScanEnergyCsv)

    
    return qmScanEnergyDf[torsionTag].to_numpy()
