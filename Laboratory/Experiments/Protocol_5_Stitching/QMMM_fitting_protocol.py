import os
from os import path as p
import numpy as np
import pandas as pd

import sys
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

##############################################################################
def fit_torsion_parameters(config: dict,
                            torsionTag: str,
                              mmTotalEnergy: np.ndarray,
                                mmTorsionEnergy: np.ndarray,
                                  shuffleIndex: int,
                                    mmCosineComponents: dict,
                                      debug: bool = False):
    if debug:
        print(f"Fitting torsion parameters for {torsionTag}")
    qmmmFittingDir = config["runtimeInfo"]["madeByStitching"]["qmmmParameterFittingDir"] 
    qmmmTorsionFittingDir = p.join(qmmmFittingDir, torsionTag)
    os.makedirs(qmmmTorsionFittingDir,exist_ok=True)
    if debug:
        print(f"Getting QM scan energies for {torsionTag}")
    qmTotalEnergy = get_qm_scan_energies(config, torsionTag)

    qmTotalEnergy = qmTotalEnergy - qmTotalEnergy.min()

    qmTorsionEnergy = qmTotalEnergy - mmTotalEnergy + mmTorsionEnergy

    qmTorsionEnergy = qmTorsionEnergy - qmTorsionEnergy.min()


    Stitching_Plotter.plot_qmmm_energies(qmTotalEnergy, qmTorsionEnergy, mmTotalEnergy, mmTorsionEnergy, mmCosineComponents, qmmmTorsionFittingDir, shuffleIndex)
    if debug:
        print(f"Running fourier transform for {torsionTag}")
    torsionParametersDf, cosineComponents = drFourier.fourier_transform_protocol(qmTorsionEnergy, torsionTag, qmmmTorsionFittingDir)

    return torsionParametersDf

##############################################################################
def get_qm_scan_energies(config: dict, torsionTag: str) -> np.array:
    qmScanEnergyCsv = config["runtimeInfo"]["madeByTwisting"]["finalScanEnergies"][torsionTag]
    if not p.exists(qmScanEnergyCsv):
        raise FileNotFoundError(f"Couldn't find QM scan energies for torsion {torsionTag}.")

    qmScanEnergyDf = pd.read_csv(qmScanEnergyCsv)

    
    return qmScanEnergyDf[torsionTag].to_numpy()