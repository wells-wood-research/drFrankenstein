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

def _get_energy_amplitude(energy: np.ndarray) -> float:
    return float(np.max(energy) - np.min(energy))


def _reconstruct_torsion_energy_from_params(torsionParametersDf: pd.DataFrame, nPoints: int = 36) -> np.ndarray:
    angle = np.radians(np.linspace(0, 360, nPoints, endpoint=False))
    reconstructedEnergy = np.zeros_like(angle)
    for _, row in torsionParametersDf.iterrows():
        amplitude = float(row["Amplitude"])
        periodicity = abs(float(row["Period"]))
        phase = np.radians(float(row["Phase"]))
        reconstructedEnergy += amplitude * (1 + np.cos(periodicity * angle - phase))
    reconstructedEnergy = reconstructedEnergy - reconstructedEnergy.min()
    return reconstructedEnergy


def _should_reset_to_flatline(qmTorsionEnergy: np.ndarray,
                              fittedTorsionParametersDf: pd.DataFrame,
                              explosionRatioThreshold: float | None) -> tuple[bool, float]:
    if explosionRatioThreshold is None:
        return False, 0.0
    qmAmplitude = _get_energy_amplitude(qmTorsionEnergy)
    if qmAmplitude <= 1e-8:
        return False, 0.0
    fittedTorsionEnergy = _reconstruct_torsion_energy_from_params(fittedTorsionParametersDf, len(qmTorsionEnergy))
    fittedAmplitude = _get_energy_amplitude(fittedTorsionEnergy)
    amplitudeRatio = fittedAmplitude / qmAmplitude
    return amplitudeRatio > explosionRatioThreshold, amplitudeRatio


def _build_flatline_parameter_df() -> pd.DataFrame:
    return pd.DataFrame([{
        "Amplitude": 0.0,
        "Period": 1,
        "Phase": 0.0,
    }])

##############################################################################
def fit_torsion_parameters(config: dict,
                            torsionTag: str,
                              mmTotalEnergy: np.ndarray,
                                mmTorsionEnergy: np.ndarray,
                                  shuffleIndex: int,
                                    mmCosineComponents: dict,
                                       maxCosineFunctions: int | None = None,
                                        debug: bool = False):
    ## unpack config
    qmmmFittingDir = config["runtimeInfo"]["madeByStitching"]["qmmmParameterFittingDir"] 
    sagvolSmoothing = config["parameterFittingInfo"]["sagvolSmoothing"]
    amplitudeExplosionRatioThreshold = config["parameterFittingInfo"].get("amplitudeExplosionRatioThreshold", None)

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
        qmTorsionEnergy = savgol_filter(qmTorsionEnergy, window_length=5, polyorder=2)

    
    Stitching_Plotter.plot_qmmm_energies(qmTotalEnergy,
                                          qmTorsionEnergy,
                                            mmTotalEnergy,
                                              mmTorsionEnergy,
                                                mmCosineComponents,
                                                  qmmmTorsionFittingDir,
                                                    shuffleIndex)

    torsionParametersDf, cosineComponents, meanAverageErrorTorsion = drFourier.fourier_transform_protocol(qmTorsionEnergy,
                                                                                   torsionTag,
                                                                                     qmmmTorsionFittingDir,
                                                                                     config=config,
                                                                                     maxCosineFunctions=maxCosineFunctions)
    shouldResetToFlatline, amplitudeRatio = _should_reset_to_flatline(
        qmTorsionEnergy, torsionParametersDf, amplitudeExplosionRatioThreshold
    )
    if shouldResetToFlatline:
        torsionParametersDf = _build_flatline_parameter_df()
        meanAverageErrorTorsion = float(np.mean(np.abs(qmTorsionEnergy)))
        meanAverageErrorTotal = float(np.mean(np.abs(qmTotalEnergy - mmTotalEnergy)))
        resetTorsions = config["runtimeInfo"]["madeByStitching"].setdefault("resetToFlatlineTorsions", [])
        resetTorsions.append({
            "torsionTag": torsionTag,
            "shuffleIndex": shuffleIndex,
            "amplitudeRatio": amplitudeRatio,
        })
        return torsionParametersDf, meanAverageErrorTorsion, meanAverageErrorTotal, True
    
    meanAverageErrorTotal = np.mean(np.abs(qmTotalEnergy - mmTotalEnergy))

    return torsionParametersDf, meanAverageErrorTorsion, meanAverageErrorTotal, False

##############################################################################
def get_qm_scan_energies(config: dict, torsionTag: str) -> np.array:
    qmScanEnergyCsv = config["runtimeInfo"]["madeByTwisting"]["finalScanEnergies"][torsionTag]
    if not p.exists(qmScanEnergyCsv):
        raise FileNotFoundError(f"Couldn't find QM scan energies for torsion {torsionTag}.")

    qmScanEnergyDf = pd.read_csv(qmScanEnergyCsv)

    
    return qmScanEnergyDf[torsionTag].to_numpy()
