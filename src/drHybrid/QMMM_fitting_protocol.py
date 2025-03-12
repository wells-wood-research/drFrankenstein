import os
from os import path as p
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Set the Agg backend before importing pyplot
import matplotlib.pyplot as plt
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
import drFourier
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
def fit_torsion_parameters(config, torsionTag, mmTotalEnergy, mmTorsionEnergy):
    qmmmFittingDir = p.join(config["pathInfo"]["parameterFittingTopDir"], "qm-mm_parameter_fitting")
    config["pathInfo"]["qmmmParameterFittingDir"] = qmmmFittingDir
    os.makedirs(qmmmFittingDir,exist_ok=True)


    qmmmTorsionFittingDir = p.join(qmmmFittingDir, torsionTag)
    os.makedirs(qmmmTorsionFittingDir,exist_ok=True)

    qmTotalEnergy = get_qm_scan_energies(config, torsionTag)

    qmTorsionEnergy = qmTotalEnergy - mmTotalEnergy + mmTorsionEnergy

    qmTorsionEnergy = qmTorsionEnergy - qmTorsionEnergy.min()

    plot_qmmm_energies(qmTotalEnergy, qmTorsionEnergy, mmTotalEnergy, mmTorsionEnergy, qmmmTorsionFittingDir)

    torsionParametersDf = drFourier.fourier_transform_protocol(qmTorsionEnergy, torsionTag, qmmmTorsionFittingDir)

    return torsionParametersDf



##############################################################################
def plot_qmmm_energies(qmTotalEnergy, qmTorsionEnergy, mmTotalEnergy, 
                       mmTorsionEnergy, outDir):
    angles = np.linspace(0, 360, len(qmTotalEnergy))

    plt.figure()
    plt.plot(angles, qmTotalEnergy, label='QM Total Energy', 
             linestyle='--', color='orange')
    plt.plot(angles, qmTorsionEnergy, label='QM Torsion Energy', 
             linewidth=2, color='red')
    plt.plot(angles, mmTotalEnergy, label='MM Total Energy', 
             linestyle='--', color='green')
    plt.plot(angles, mmTorsionEnergy, label='MM Torsion Energy', 
             linestyle='--', color='blue')
    plt.legend()
    plt.xlabel('Torison Angle')
    plt.ylabel('Energy (Kcal / mol)')
    plt.title('QM/MM Energies')
    plt.savefig(p.join(outDir, "QM-MM_torsion_energies.png"))
    plt.close()
##############################################################################
def get_qm_scan_energies(config: dict, torsionTag: str) -> np.array:
    qmScanEnergyCsv = config["torsionScanInfo"]["finalScanEnergies"][torsionTag]

    qmScanEnergyDf = pd.read_csv(qmScanEnergyCsv)
    
    return qmScanEnergyDf[torsionTag].to_numpy()