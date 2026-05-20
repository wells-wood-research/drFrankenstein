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

##############################################################################
def fit_torsion_parameters(config: dict,
                            torsionTag: str,
                              mmTotalEnergy: np.ndarray,
                                mmTorsionEnergy: np.ndarray,
                                  shuffleIndex: int,
                                    mmCosineComponents: dict,
                                      debug: bool = False):
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
    ## calculate QM torsion energy to fit parameters to 
    qmTorsionEnergy = qmTotalEnergy - mmTotalEnergy + mmTorsionEnergy
    ## normalise QM Torsion energy
    qmTorsionEnergy = qmTorsionEnergy - qmTorsionEnergy.min()

    if not sagvolSmoothing is None:
        qmTorsionEnergy = savgol_filter(qmTorsionEnergy, window_length=9, polyorder=2)

    meanAverageErrorTotal = np.mean(np.abs(qmTotalEnergy - mmTotalEnergy))
    totalMetrics = Stitching_Assistant.calculate_profile_fit_metrics(qmTotalEnergy, mmTotalEnergy)
    
    torsionParametersDf, cosineComponents, reconstructedSignal, meanAverageErrorTorsion, fitScoreTorsion = drFourier.fourier_transform_protocol(qmTorsionEnergy,
                                                                                  torsionTag,
                                                                                    qmmmTorsionFittingDir,
                                                                                    config=config)
    torsionMetrics = Stitching_Assistant.calculate_profile_fit_metrics(qmTorsionEnergy, reconstructedSignal)
    torsionConverged = Stitching_Assistant.check_torsion_convergence(torsionMetrics["composite_score"], converganceTolerance)

    Stitching_Plotter.plot_qmmm_energies(qmTotalEnergy,
                                          qmTorsionEnergy,
                                            mmTotalEnergy,
                                              reconstructedSignal,
                                                cosineComponents,
                                                  torsionMetrics,
                                                    totalMetrics,
                                                    meanAverageErrorTorsion,
                                                    meanAverageErrorTotal,
                                                    config["runtimeInfo"]["madeByStitching"]["maxTorsions"],
                                                    torsionConverged,
                                                    qmmmTorsionFittingDir,
                                                      shuffleIndex)

    return torsionParametersDf, torsionMetrics, totalMetrics, meanAverageErrorTorsion, meanAverageErrorTotal, torsionConverged

##############################################################################
def get_qm_scan_energies(config: dict, torsionTag: str) -> np.array:
    qmScanEnergyCsv = config["runtimeInfo"]["madeByTwisting"]["finalScanEnergies"][torsionTag]
    if not p.exists(qmScanEnergyCsv):
        raise FileNotFoundError(f"Couldn't find QM scan energies for torsion {torsionTag}.")

    qmScanEnergyDf = pd.read_csv(qmScanEnergyCsv)

    
    return qmScanEnergyDf[torsionTag].to_numpy()
