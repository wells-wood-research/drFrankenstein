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
from drHelper import print_dict
import drHybrid.drFourier as drFourier
from drHybrid import Plotter
##############################################################################
def dummy_inputs(config):
    torsionScanDir = "/home/esp/scriptDevelopment/drFrankenstein/ALA_outputs/02_torsion_scanning"

    torsionDirs = [p.join(torsionScanDir, dir) for dir in os.listdir(torsionScanDir)]

    torsionTags = [dir.split("_")[1] for dir in os.listdir(torsionScanDir)]

    config["torsionScanInfo"]["finalScanEnergies"] = {}

    for torsionDir, torsionTag in zip(torsionDirs, torsionTags):
        scanCsv = p.join(torsionDir,"scan_data","final_scan_energies.csv")
        config["torsionScanInfo"]["finalScanEnergies"][torsionTag] = scanCsv

    return config

##############################################################################
def fit_torsion_parameters(config, torsionTag, mmTotalEnergy, mmTorsionEnergy, shuffleIndex, mmCosineComponents):
    qmmmFittingDir = config["pathInfo"]["qmmmParameterFittingDir"] 
    qmmmTorsionFittingDir = p.join(qmmmFittingDir, torsionTag)
    os.makedirs(qmmmTorsionFittingDir,exist_ok=True)

    qmTotalEnergy = get_qm_scan_energies(config, torsionTag)

    qmTotalEnergy = qmTotalEnergy - qmTotalEnergy.min()

    qmTorsionEnergy = qmTotalEnergy - mmTotalEnergy + mmTorsionEnergy

    qmTorsionEnergy = qmTorsionEnergy - qmTorsionEnergy.min()


    Plotter.plot_qmmm_energies(qmTotalEnergy, qmTorsionEnergy, mmTotalEnergy, mmTorsionEnergy, mmCosineComponents, qmmmTorsionFittingDir, shuffleIndex)

    torsionParametersDf, cosineComponents = drFourier.fourier_transform_protocol(qmTorsionEnergy, torsionTag, qmmmTorsionFittingDir)

    return torsionParametersDf

##############################################################################
def get_qm_scan_energies(config: dict, torsionTag: str) -> np.array:
    qmScanEnergyCsv = config["torsionScanInfo"]["finalScanEnergies"][torsionTag]

    qmScanEnergyDf = pd.read_csv(qmScanEnergyCsv)
    
    return qmScanEnergyDf[torsionTag].to_numpy()