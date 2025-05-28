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
                            torsion_tag: str,
                              mm_total_energy: np.ndarray,
                                mm_torsion_energy: np.ndarray,
                                  shuffle_index: int,
                                    mm_cosine_components: dict,
                                      debug: bool = False) -> pd.DataFrame:
    """
    Fits torsion parameters by comparing QM and MM energies.

    Args:
        config: Configuration dictionary.
        torsion_tag: Identifier for the torsion being fitted.
        mm_total_energy: NumPy array of MM total energies.
        mm_torsion_energy: NumPy array of MM torsion energies for the specific dihedral.
        shuffle_index: Index used for shuffling or selecting a specific conformer/scan.
        mm_cosine_components: Dictionary of MM cosine components for the dihedral.
        debug: Boolean flag for debug mode.

    Returns:
        DataFrame containing the fitted torsion parameters.
    """
    qmmmFittingDir = config["runtimeInfo"]["madeByStitching"]["qmmmParameterFittingDir"] 
    qmmmTorsionFittingDir = p.join(qmmmFittingDir, torsion_tag)
    os.makedirs(qmmmTorsionFittingDir,exist_ok=True)

    qmTotalEnergy = get_qm_scan_energies(config, torsion_tag)

    qmTotalEnergy = qmTotalEnergy - qmTotalEnergy.min()

    qmTorsionEnergy = qmTotalEnergy - mm_total_energy + mm_torsion_energy

    qmTorsionEnergy = qmTorsionEnergy - qmTorsionEnergy.min()

    Stitching_Plotter.plot_qmmm_energies(qmTotalEnergy, qmTorsionEnergy, mm_total_energy, mm_torsion_energy, mm_cosine_components, qmmmTorsionFittingDir, shuffle_index)

    torsionParametersDf, cosineComponents = drFourier.fourier_transform_protocol(qmTorsionEnergy,
                                                                                  torsion_tag,
                                                                                    qmmmTorsionFittingDir,
                                                                                      forcefeild=config["parameterFittingInfo"]["forceField"])

    return torsionParametersDf

##############################################################################
def get_qm_scan_energies(config: dict, torsion_tag: str) -> np.ndarray:
    """
    Retrieves QM scan energies for a specific torsion from a CSV file.

    Args:
        config: Configuration dictionary.
        torsion_tag: Identifier for the torsion.

    Returns:
        NumPy array of QM scan energies.
    """
    qmScanEnergyCsv = config["runtimeInfo"]["madeByTwisting"]["finalScanEnergies"][torsion_tag]
    if not p.exists(qmScanEnergyCsv):
        raise FileNotFoundError(f"Couldn't find QM scan energies for torsion {torsion_tag}.")

    qmScanEnergyDf = pd.read_csv(qmScanEnergyCsv)

    
    return qmScanEnergyDf[torsion_tag].to_numpy()