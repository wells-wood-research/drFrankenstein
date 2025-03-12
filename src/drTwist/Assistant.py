"""
This script contains:

- set_up_directories
- identify_rotatable_bonds [TODO: remove rdkit dependency]
- process_scan_data 

"""
import os
from os import path as p
import sys
import pandas as pd
import numpy as np
import glob
import re
from scipy.signal import argrelextrema

## CLEAN CODE CLASSES ##
class FilePath:
    pass
class DirectoryPath:
    pass
from typing import Dict, List, Tuple
## ADD SRC TO PATH ##
currentFilePath: FilePath = os.path.abspath(__file__)
currentDir: DirectoryPath = os.path.dirname(currentFilePath)
srcDir: DirectoryPath = os.path.dirname(currentDir)
sys.path.append(srcDir)

## RDKIT IMPORTS ##
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem.rdmolfiles import MolToPDBFile
RDLogger.DisableLog('rdApp.warning')

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def read_singlepoint_energy(spOrcaOutput):
    with open(spOrcaOutput, "r") as f:
        for line in f:
            if line.startswith("FINAL SINGLE POINT ENERGY"):
                singlePointEnergy = float(line.split()[-1])
                return singlePointEnergy
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def add_mid_points(indexes: np.array):
    newIndexes = []
    for i, indexA in enumerate(indexes[:-1]):
        newIndexes.append(indexA)
        indexB = indexes[i+1]
        midPointIndex = (indexA + indexB) // 2
        newIndexes.append(midPointIndex)
    newIndexes.append(indexes[-1])
    
    return newIndexes
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def find_local_extrema(energies):
    # Convert the energies to a numpy array
    energiesArray = energies.to_numpy()

    # Extend the array to handle periodicity
    extendedEnergies = np.concatenate(
        (energiesArray[-1:], energiesArray, energiesArray[:1])
    )

    # Find local minima
    localMinimaIndexes = argrelextrema(extendedEnergies, np.less)[0] - 1

    # Find local maxima
    localMaximaIndexes = argrelextrema(extendedEnergies, np.greater)[0] - 1

    # Filter out-of-bound indices
    localMinimaIndexes = localMinimaIndexes[
        (localMinimaIndexes >= 0) & (localMinimaIndexes < len(energiesArray))
    ]
    localMaximaIndexes = localMaximaIndexes[
        (localMaximaIndexes >= 0) & (localMaximaIndexes < len(energiesArray))
    ]

    combinedExtremaIndexes = sorted(
        np.concatenate((localMinimaIndexes, localMaximaIndexes))
    )
    return [int(index) for index in combinedExtremaIndexes]

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def find_scan_xyz_files(scanDir: DirectoryPath, expectedNumberOfFiles: int):
    scanXyzs = sorted([xyzFile for xyzFile in glob.glob(p.join(scanDir, 'orca_scan.[0-9][0-9][0-9].xyz'))])
    if not len(scanXyzs) == expectedNumberOfFiles:
        raise(ValueError(f"Number of scan xyz files ({len(scanXyzs)}) does not match expected ({expectedNumberOfFiles})"))

    return scanXyzs
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def process_energy_outputs(energyDf):
    energyDf = rescale_and_sort_energy_angles(energyDf)
    energyDf = take_min_duplicate_angles(energyDf)

    energyDf["Energy"] = energyDf["Energy"].apply(hartree_to_kcal_per_mol)
    return energyDf
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def hartree_to_kcal_per_mol(energy):
    return energy * 627.5095
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def rescale_torsion_angles(angle):
    angle = angle % 360  # reduce the angle to the 0-360 range
    return angle
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def rescale_and_sort_energy_angles(energyDf):
    energyDf["Energy"] = energyDf["Energy"] - energyDf["Energy"].min()
    energyDf["Angle"] = energyDf["Angle"].apply(rescale_torsion_angles)
    energyDf = energyDf.sort_values(by="Angle", ascending=True)
    return energyDf
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def take_min_duplicate_angles(energyDf):
    energyDf['Angle'] = energyDf['Angle'].round(0)
    # Find the index of the minimum energy for each unique Angle
    minEnergyIndexes = energyDf.groupby('Angle')["Energy"].idxmin()

    minEnergyIndexes = minEnergyIndexes.dropna()

    # Select the rows with the minimum energy for each Angle
    return energyDf.loc[minEnergyIndexes].reset_index(drop=True)
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def read_scan_energy_data(conformertorsionScanDir):
    scanDat: FilePath = p.join(conformertorsionScanDir, "orca_scan.relaxscanact.dat")
    if not os.path.exists(scanDat):
        raise FileNotFoundError(f"Scan data not found in {conformertorsionScanDir}")
    scanDf: pd.DataFrame = pd.read_csv(scanDat, sep='\s+', header=None)
    return scanDf
    
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def find_final_xyz(conformertorsionScanDir):
    allXyzFiles = glob.glob(p.join(conformertorsionScanDir, "*.xyz"))
    scanXYZFiles = sorted([f for f in allXyzFiles if re.match(r'orca_scan\.\d+\.xyz$', os.path.basename(f))])

    finalXyzFile = scanXYZFiles[-1]
    return finalXyzFile
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def measure_current_torsion_angle(conformerXyz, torsionIndexes):
    atomCoords = xyz_to_np_array(conformerXyz)
    torsionAngle = calculate_torsion_angle(atomCoords, torsionIndexes)
    roundedToTenAngle = round(torsionAngle, -1)
    return roundedToTenAngle
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def xyz_to_np_array(filePath):
    with open(filePath, 'r') as file:
        lines = file.readlines()

    # Skip the first two lines (atom count and comment)
    dataLines = lines[2:]

    # Parse the coordinates
    coords = []
    for line in dataLines:
        parts = line.split()
        if len(parts) >= 4:
            x, y, z = map(float, parts[1:4])
            coords.append([x, y, z])

    return np.array(coords)
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def calculate_torsion_angle(coords, torsionIndexes):
    # Vectors between points
    b1 = coords[torsionIndexes[1]] - coords[torsionIndexes[0]]
    b2 = coords[torsionIndexes[2]] - coords[torsionIndexes[1]]
    b3 = coords[torsionIndexes[3]] - coords[torsionIndexes[2]]

    # Normal vectors to the planes
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)

    # Normalize the normal vectors
    n1 /= np.linalg.norm(n1)
    n2 /= np.linalg.norm(n2)

    # Calculate the angle
    angleInRadians = np.arctan2(np.dot(np.cross(n1, n2), b2 / np.linalg.norm(b2)), np.dot(n1, n2))

    # Convert from radians to degrees
    angleIdDegrees = np.degrees(angleInRadians)

    return angleIdDegrees

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

def set_up_directories(config: dict) -> dict:
    outputDir = config["pathInfo"]["outputDir"]
    torsionTopDir = p.join(outputDir, "03_torsion_scanning")
    os.makedirs(torsionTopDir, exist_ok=True)
    config["pathInfo"]["torsionTopDir"] = torsionTopDir

    return config

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

def identify_rotatable_bonds(inputPdb) -> List[Tuple[int,int,int,int]]:
    # Load the molecule from a PDB file
    mol = Chem.MolFromPDBFile(inputPdb, removeHs=False)
    # Identify torsion angles for rotatable bonds
    torsionAngles = []
    for bond in mol.GetBonds():
        if bond.IsInRing():
            continue
        if bond.GetBondType() == Chem.BondType.SINGLE:
            atom2 = bond.GetBeginAtom()
            atom3 = bond.GetEndAtom()

            ## get atom names for begin and end atoms
            atom2Name = atom2.GetPDBResidueInfo().GetName().strip()
            atom3Name = atom3.GetPDBResidueInfo().GetName().strip()

            ## dont scan amide bonds
            if atom2Name in ["NN", "C"] and atom3Name in ["NN", "C"]:
                continue
            if atom2Name in ["N", "CC1"] and atom3Name in ["N", "CC1"]:
                continue

            
            if not (atom2.IsInRing() or atom3.IsInRing()):
                # Find neighboring atoms for torsion angle
                neighborsBegin = [a for a in atom2.GetNeighbors() if a.GetIdx() != atom3.GetIdx()]
                neighborsEnd = [a for a in atom3.GetNeighbors() if a.GetIdx() != atom2.GetIdx()]
                
                for atom1 in neighborsBegin:
                    for atom4 in neighborsEnd:

                        atom1Name = atom1.GetPDBResidueInfo().GetName().strip()
                        atom4Name = atom4.GetPDBResidueInfo().GetName().strip()

                        ## dont scan bonds with non-polar hydrogens at either as atoms 1 or 4
                        if atom1Name.startswith("H"):
                            if atom2Name.startswith("C"):
                                continue
                        if atom4Name.startswith("H"):
                            if atom3Name.startswith("C"):
                                continue
                        ## add torsion data to list
                        torsionAngles.append({
                            'atoms': (atom1Name, atom2Name, atom3Name, atom4Name),
                            'indices': (atom1.GetIdx(), atom2.GetIdx(), atom3.GetIdx(), atom4.GetIdx())
                        })


    return torsionAngles
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

def process_scan_data(scanDfs: List[pd.DataFrame],
                       torsionTopDir: DirectoryPath,
                       torsionTag: str):

    ## make a dir to store output csv files
    dataDir = p.join(torsionTopDir, "scan_data")
    os.makedirs(dataDir, exist_ok=True)

    ## merge scan dataframes
    # scanDfs = [process_energy_outputs(scanDf) for scanDf in scanDfs]
    mergedDf = merge_scan_dfs(scanDfs)
 
    ## remove cols with large jumps in energy
    mergedDf = detect_jumps_in_data(mergedDf)

    ## write to csv
    mergedScanCsv = p.join(dataDir, f"scan_energies.csv")
    mergedDf.to_csv(mergedScanCsv, index=False)

    scanAverageDf = pd.DataFrame()
    scanAverageDf["Angle"] = mergedDf["Angle"]

    ## calculate averages
    finalScanEnergiesCsv = p.join(dataDir, "final_scan_energies.csv")
    scanAverageDf[torsionTag] = mergedDf.drop(columns="Angle").mean(axis=1)
    scanAverageDf.to_csv(finalScanEnergiesCsv, index=False)

    return finalScanEnergiesCsv, scanAverageDf
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

def merge_scan_dfs(scanDfs):
    mergedDf = pd.DataFrame()
    mergedDf["Angle"] = np.arange(0,360,10)
    for colIndex, scanDf in enumerate(scanDfs):
        mergedDf = mergedDf.merge(scanDf[["Angle", "Energy"]], on = "Angle", how= "left")
        mergedDf.rename(columns={"Energy":f"Energy_{colIndex + 1}"}, inplace=True)
    return mergedDf
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

def detect_jumps_in_data(df):
    # Calculate differences between consecutive values
    diffDf = df.drop(columns='Angle').diff().abs()
    
    # Identify columns with any difference greater than the threshold
    jumpyCols = diffDf.columns[((diffDf > 10).any())]
    
    # Drop these columns from the original DataFrame
    cleanDf = df.drop(columns=jumpyCols)
    
    return cleanDf
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
