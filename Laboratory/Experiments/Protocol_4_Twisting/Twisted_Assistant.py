"""
This script contains:

- set_up_directories
- identify_rotatable_bonds [TODO: remove rdkit dependency]
- process_scan_data 

"""
import os
from os import path as p
import re
from collections import defaultdict
import pandas as pd
import numpy as np
import glob
import re
from shutil import rmtree
from scipy.signal import argrelextrema
import random
import parmed
from parmed.charmm import CharmmParameterSet, CharmmPsfFile

## CLEAN CODE CLASSES ##
class FilePath:
    pass
class DirectoryPath:
    pass
from typing import List, Tuple



from OperatingTools import drOrca
from OperatingTools import file_parsers



def choose_torsions_to_scan(config: dict) -> dict:
    """
    Uses the "runScansOn" entry in the "torsionScanInfo" entry in config to choose what torsions to scan

    Args:
        config (dict): the config dict
    
    Returns:    
        config (dict): the config dict updated
    """
    ## unpack config
    runScansOn = config["torsionScanInfo"]["runScansOn"]
    allDihedrals = config["runtimeInfo"]["madeByTwisting"]["allDihedrals"]

    torsionsToScan = allDihedrals["remainingDihedrals"]

    if runScansOn["phiPsi"]:
        torsionsToScan = {**torsionsToScan, **allDihedrals["phiPsiDihedrals"]}
    if runScansOn["nonPolarProtons"]:
        torsionsToScan = {**torsionsToScan, **allDihedrals["nonPolarProtonsDihedrals"]}
    if runScansOn["polarProtons"]:
        torsionsToScan = {**torsionsToScan, **allDihedrals["polarProtonDihedrals"]}
    if runScansOn["nonAromaticRings"]:
        torsionsToScan = {**torsionsToScan, **allDihedrals["nonAromaticRingDihedrals"]}


    config["runtimeInfo"]["madeByTwisting"]["torsionsToScan"] = torsionsToScan
    return config

def get_unique_bonds(adjacencyMatrix: defaultdict):
    """
    Processes the adjacency matrix to get a stable list of unique bond pairs
    (atomIdx1, atomIdx2) assuming 0-based indexing consistent with xyz2df.
    """
    uniqueBondPairs = [] 
    processedBondsSet = set()

    sortedAtom1Indices = sorted(adjacencyMatrix.keys())

    for atom1 in sortedAtom1Indices:
        sortedNeighbors = sorted(adjacencyMatrix[atom1])
        for atom2 in sortedNeighbors:
            canonicalPair = tuple(sorted((atom1, atom2))) # e.g. (0,1) not (1,0)
            if canonicalPair not in processedBondsSet:
                processedBondsSet.add(canonicalPair)
                uniqueBondPairs.append(canonicalPair) # Store the (smaller_idx, larger_idx) pair
    return uniqueBondPairs

def calculate_bond_lengths_for_frame_vectorized(
    uniqueBondPairs: list[tuple[int, int]], 
    coordsDf: pd.DataFrame
) -> dict[str, float]:
    """
    Calculates bond lengths for a single geometry frame using vectorized operations.
    """
    if not uniqueBondPairs:
        return {}

    atom1Indices = [pair[0] for pair in uniqueBondPairs]
    atom2Indices = [pair[1] for pair in uniqueBondPairs]

    coords1 = coordsDf.loc[atom1Indices, ['x', 'y', 'z']].values
    coords2 = coordsDf.loc[atom2Indices, ['x', 'y', 'z']].values
    
    diffSq = (coords1 - coords2)**2
    distances = np.sqrt(np.sum(diffSq, axis=1))
    
    bondLengthsDict = {}
    for i, pair in enumerate(uniqueBondPairs):
        # Create bond label from the canonical pair used for indexing
        bondLabel = f"{pair[0]}-{pair[1]}" 
        bondLengthsDict[bondLabel] = distances[i]
        
    return bondLengthsDict

def have_scans_exploded(scanDir: DirectoryPath, configDict: dict, tolerance: float = 0.5):
    adjacencyMatrix = configDict["runtimeInfo"]["madeByTwisting"]["adjacencyMatrix"]
    orcaScanXyzs = sorted(glob.glob(p.join(scanDir, "orca_scan.[0-9][0-9][0-9].xyz")))

    if not orcaScanXyzs:
        return False 

    # Assumes adjacencyMatrix uses 0-based indexing for atoms
    uniqueBondPairs = get_unique_bonds(adjacencyMatrix)

    if not uniqueBondPairs:
        return False 

    minLengthsSoFar = {}
    maxLengthsSoFar = {}

    # Process first file to initialize min/max
    firstXyzDf = file_parsers.xyz2df(orcaScanXyzs[0])
    currentBondLengths = calculate_bond_lengths_for_frame_vectorized(uniqueBondPairs, firstXyzDf)
        
    for pair in uniqueBondPairs: # Iterate using the definitive list of bonds
        bondLabel = f"{pair[0]}-{pair[1]}"
        length = currentBondLengths.get(bondLabel)
        if length is not None:
            minLengthsSoFar[bondLabel] = length
            maxLengthsSoFar[bondLabel] = length
    
    if len(orcaScanXyzs) < 2: # Need at least two frames to check for variation
         return False

    for xyzFile in orcaScanXyzs[1:]:
        xyzDf = file_parsers.xyz2df(xyzFile)
        currentBondLengths = calculate_bond_lengths_for_frame_vectorized(uniqueBondPairs, xyzDf)

        for bondLabelInFrame, length in currentBondLengths.items():
            # Only process bonds that were successfully defined and calculated from uniqueBondPairs
            if bondLabelInFrame not in minLengthsSoFar: 
                minLengthsSoFar[bondLabelInFrame] = length
                maxLengthsSoFar[bondLabelInFrame] = length
            else:
                minLengthsSoFar[bondLabelInFrame] = min(minLengthsSoFar[bondLabelInFrame], length)
                maxLengthsSoFar[bondLabelInFrame] = max(maxLengthsSoFar[bondLabelInFrame], length)

            minL = minLengthsSoFar[bondLabelInFrame]
            maxL = maxLengthsSoFar[bondLabelInFrame]
            
            # Check for explosion condition for this bond immediately
            if maxL > minL * (1.0 + tolerance):
                return True # Explosion detected
    
    return False # No explosion detected after checking all files and bonds




def gather_scan_data(averagesDf: pd.DataFrame, torsionTag: str, config: dict) -> dict:
    """
    Gets the global minimum angle and barrier height for the torsion in question

    Args:
        averagesDf (pd.DataFrame): dataframe containing single point energies
        torsionTag (str): torsion tag
        config (dict): config containing all run information
    
    Returns:
        config (dict): updated config
    """
    ## If no data was found for this torsion, skip
    if averagesDf.empty or any(pd.isna(averagesDf[torsionTag])):

        print(f"Skipping torsion {torsionTag} because no data was found.")
        config["runtimeInfo"]["madeByTwisting"]["torsionsToScan"][torsionTag]["globalMinimaAngle"] = None
        config["runtimeInfo"]["madeByTwisting"]["torsionsToScan"][torsionTag]["barrierHeight"] = None
        return config


    globalMinimaAngle = averagesDf["Angle"].loc[averagesDf[torsionTag].idxmin()]
    globalMinimaEnergy = averagesDf[torsionTag].loc[averagesDf[torsionTag].idxmin()]
    globalMaximaEnergy = averagesDf[torsionTag].loc[averagesDf[torsionTag].idxmax()]



    barrierHeight = round(globalMaximaEnergy - globalMinimaEnergy, 3)

    config["runtimeInfo"]["madeByTwisting"]["torsionsToScan"][torsionTag]["globalMinimaAngle"] = int(globalMinimaAngle)
    config["runtimeInfo"]["madeByTwisting"]["torsionsToScan"][torsionTag]["barrierHeight"] = float(barrierHeight)


    return config

def load_amber_params(config):
    ## unpack config
    assembledPrmtop = config["runtimeInfo"]["madeByAssembly"]["assembledPrmtop"]

    parmedPrmtop = parmed.load_file(assembledPrmtop)

    return parmedPrmtop

def load_charmm_params(config):
    ## unpack config
    assembledPsf = config["runtimeInfo"]["madeByAssembly"]["assembledPsf"]
    assembledPrm = config["runtimeInfo"]["madeByAssembly"]["assembledPrm"]
    assembledRtf = config["runtimeInfo"]["madeByAssembly"]["assembledRtf"]

    parmedPsf = CharmmPsfFile(assembledPsf)
    parmedPsf.load_parameters(CharmmParameterSet(assembledPrm, assembledRtf))

    return parmedPsf

def identify_rotatable_bonds(config: dict, mode: str = "AMBER") -> List[Tuple[int,int,int,int]]:

    if mode == "AMBER":
        moleculeParams = load_amber_params(config)
    elif mode == "CHARMM":
        moleculeParams = load_charmm_params(config)
    else:
        raise ValueError(f"Unknown mode: {mode}")
    ## init empty lists
    aromaticDihedrals = []
    nonAromaticRingDihedrals = []       
    nonPolarProtonsDihedrals = []         ## eg N - C - C - H
    polarProtonDihedrals = []           ## eg N - C - N - H
    remainingDihedrals = []             ## store remaining rotatable bonds
    phiPsiDihedrals = []

    rings = find_ring_atoms(moleculeParams)
    adjacencyMatrix = _construct_adjacency_matrix(moleculeParams)

    aromaticRings, nonAromaticRings = classify_rings_aromatic(rings, moleculeParams, adjacencyMatrix)

    for dihedral in moleculeParams.dihedrals:
        atomNames = extract_dihedral_atom_names(dihedral)
        atomTypes = extract_dihedral_atom_types(dihedral)
        atomIndexes = extract_dihedral_atom_indexes(dihedral)
        ## skip amide dihedrals
        if _is_amide_dihedral(dihedral):
            continue
        ## store dihedrals that end with a non-polar proton
        if _is_non_polar_proton_dihedral(dihedral):
            nonPolarProtonsDihedrals.append((atomTypes, atomNames, atomIndexes))
            continue
        ## store dihedrals that end with a polar proton
        if _is_polar_proton_dihedral(dihedral):
            polarProtonDihedrals.append((atomTypes, atomNames, atomIndexes))
        ## deal with rings
        if _is_a_ring_dihedral(dihedral, rings):
            ## store aromatic dihedrals
            if _is_a_ring_dihedral(dihedral, aromaticRings):
                aromaticDihedrals.append((atomTypes, atomNames, atomIndexes))
                continue
            else:
                ## store non-aromatic ring dihedrals
                nonAromaticRingDihedrals.append((atomTypes, atomNames, atomIndexes))
                continue

        ## deal with phi / psi dihedrals
        backboneAliases = config["moleculeInfo"].get("backboneAliases", None)
        if backboneAliases is not None:
            if _is_a_phi_dihedral(dihedral, backboneAliases):
                ## store phi dihedrals
                phiPsiDihedrals.append((atomTypes, atomNames, atomIndexes))
                continue
            elif _is_a_psi_dihedral(dihedral, backboneAliases):
                ## store psi dihedrals
                phiPsiDihedrals.append((atomTypes, atomNames, atomIndexes))
                continue

        ## store remaining rotatable dihedrals
        remainingDihedrals.append((atomTypes, atomNames, atomIndexes))
    ## assign torsion tags to the dihedrals
    taggedAromaticDihedrals = assign_torsion_tags(aromaticDihedrals)
    taggedNonAromaticRingDihedrals = assign_torsion_tags(nonAromaticRingDihedrals)
    taggedPolarProtonDihedrals = assign_torsion_tags(polarProtonDihedrals)
    taggedNonPolarProtonDihedrals = assign_torsion_tags(nonPolarProtonsDihedrals)
    taggedPhiPsiDihedrals = assign_torsion_tags(phiPsiDihedrals)
    taggedRemainingDihedrals = assign_torsion_tags(remainingDihedrals)

    ## update config
    config["runtimeInfo"]["madeByTwisting"]["adjacencyMatrix"]          = dict(adjacencyMatrix)

    allDihedrals = {}
    allDihedrals["aromaticDihedrals"] = taggedAromaticDihedrals
    allDihedrals["nonAromaticRingDihedrals"] = taggedNonAromaticRingDihedrals
    allDihedrals["polarProtonDihedrals"] = taggedPolarProtonDihedrals
    allDihedrals["nonPolarProtonsDihedrals"] = taggedNonPolarProtonDihedrals
    allDihedrals["phiPsiDihedrals"] = taggedPhiPsiDihedrals
    allDihedrals["remainingDihedrals"] = taggedRemainingDihedrals

    config["runtimeInfo"]["madeByTwisting"]["allDihedrals"]             = allDihedrals

    return config

def _is_a_phi_dihedral(dihedral, backboneAliases):
    if not ("N" in backboneAliases.keys() and "CA" in backboneAliases.keys()):
        return False
    dihedralAtomNames = extract_dihedral_atom_names(dihedral)
    if dihedralAtomNames[1] in backboneAliases["N"] and dihedralAtomNames[2] in backboneAliases["CA"]:
        return True
    elif dihedralAtomNames[2] in backboneAliases["N"] and dihedralAtomNames[1] in backboneAliases["CA"]:
        return True
    else:
        return False
    
def _is_a_psi_dihedral(dihedral, backboneAliases):
    if not ("CA" in backboneAliases.keys() and "C" in backboneAliases.keys()):
        return False
    dihedralAtomNames = extract_dihedral_atom_names(dihedral)
    if dihedralAtomNames[1] in backboneAliases["CA"] and dihedralAtomNames[2] in backboneAliases["C"]:
        return True
    elif dihedralAtomNames[2] in backboneAliases["CA"] and dihedralAtomNames[1] in backboneAliases["C"]:
        return True
    else:
        return False

def exclude_backbone_torsions(config: dict) -> dict:
    """
    Removes Phi and Psi Angles from the rotatable bonds list

    Args:
        config (dict): the config dict
    
    Returns:
        config (dict): the config dict updated
    """
    ## unpack config
    uniqueRotatableDihedrals = config["runtimeInfo"]["madeByTwisting"]["torsionsToScan"]
    forceField = config["parameterFittingInfo"]["forceField"]
    if forceField == "CHARMM":
        phiCenterTypes = ("NH1", "CT1") ## C N CA C
        psiCenterTypes = ("CT1", "C") ## N CA C N
    elif forceField == "AMBER":
        phiCenterTypes = ("N", "CT") ## C N CA C
        psiCenterTypes = ("CT", "C") ## N CA C N    

    tagsToRemove = []
    for torsionTag, dihedralData in uniqueRotatableDihedrals.items():
        dihedralCenterTypes = dihedralData["ATOM_TYPES"][1:3]
        if dihedralCenterTypes == phiCenterTypes or dihedralCenterTypes == phiCenterTypes[::-1]:
            tagsToRemove.append(torsionTag)
        elif dihedralCenterTypes == psiCenterTypes or dihedralCenterTypes == psiCenterTypes[::-1]:
            tagsToRemove.append(torsionTag)

    for torsionTag in tagsToRemove:
        uniqueRotatableDihedrals.pop(torsionTag)
    config["runtimeInfo"]["madeByTwisting"]["uniqueRotatableDihedrals"] = uniqueRotatableDihedrals
    return config

        
        



# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def get_conformer_xyzs(config, seed=1818, temperature=300):
    ## get conformer XYZ files
    conformerXyzs = config["runtimeInfo"]["madeByConformers"]["conformerXyzs"]
    conformerEnergies = config["runtimeInfo"]["madeByConformers"]["conformerEnergies"]
    nConformers = config["torsionScanInfo"]["nConformers"]

    if nConformers == -1 or nConformers >= len(conformerXyzs):
        return conformerXyzs  # Return all conformers if nConformers is -1 or too large

    # Boltzmann sampling
    random.seed(seed)
    kBoltzmann = 0.0019872041  # Boltzmann constant in kcal/mol·K
    kT = kBoltzmann * temperature

    # Get energies and compute Boltzmann weights
    energies = [conformerEnergies[key] for key in sorted(conformerEnergies.keys())]
    boltzmannFactors = [np.exp(-energy / kT) for energy in energies]
    totalWeight = sum(boltzmannFactors)
    probabilities = [factor / totalWeight for factor in boltzmannFactors]

    # Ensure conformerXyzs is aligned with sorted energies
    sortedKeys = sorted(conformerEnergies.keys())
    sortedConformerXyzs = [conformerXyzs[list(conformerEnergies.keys()).index(key)] for key in sortedKeys]

    # Sample conformers based on Boltzmann probabilities
    selectedIndices = np.random.choice(
        len(sortedConformerXyzs),
        size=nConformers,
        replace=False,  # No replacement to avoid duplicates
        p=probabilities
    )
    selectedConformerXyzs = [sortedConformerXyzs[i] for i in selectedIndices]
    print(selectedConformerXyzs)

    return selectedConformerXyzs

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

def create_orca_terminated_flag(orcaDir, orcaOut):

    if drOrca.did_orca_finish_normallly(orcaOut):
        with open(p.join(orcaDir, "ORCA_FINISHED_NORMALLY"), "w") as f:
            f.write("ORCA FINISHED NORMALLY")
    else:
        with open(p.join(orcaDir, "ORCA_CRASHED"), "w") as f:
            f.write("ORCA CRASHED")


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
def read_scan_energy_data(conformerTorsionScanDir):
    scanDat: FilePath = p.join(conformerTorsionScanDir, "orca_scan.relaxscanact.dat")
    if not os.path.exists(scanDat):
        raise FileNotFoundError(f"Scan data not found in {conformerTorsionScanDir}")
    scanDf: pd.DataFrame = pd.read_csv(scanDat, sep='\s+', header=None)
    return scanDf
    
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def find_final_xyz(conformerTorsionScanDir):
    allXyzFiles = glob.glob(p.join(conformerTorsionScanDir, "*.xyz"))
    scanXyzFiles = sorted([f for f in allXyzFiles if re.match(r'orca_scan\.\d+\.xyz$', os.path.basename(f))])

    finalXyzFile = scanXyzFiles[-1]
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
    angleInDegrees = np.degrees(angleInRadians)

    return angleInDegrees

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

def set_up_directories(config: dict) -> dict:
    outputDir = config["pathInfo"]["outputDir"]
    torsionTopDir = p.join(outputDir, "04_torsion_scanning")
    os.makedirs(torsionTopDir, exist_ok=True)
    config["runtimeInfo"]["madeByTwisting"]["torsionDir"] = torsionTopDir

    return config

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def assign_torsion_tags(dihedrals: list[tuple[tuple[str], tuple[str]]]) -> dict[tuple[str], tuple[str], tuple[str]]:
    taggedDihedrals = {}
    seenTags = set()
    for dihedralData in dihedrals:
        torsionTag = "-".join(dihedralData[0])
        if torsionTag not in seenTags:
            seenTags.add(torsionTag)
            taggedDihedrals[torsionTag] = {"ATOM_TYPES": dihedralData[0], "ATOM_NAMES": dihedralData[1], "ATOM_INDEXES": dihedralData[2]}
        else:
           continue
    return taggedDihedrals



def get_unique_dihedrals(dihedralGroup: list[tuple[tuple[str], tuple[str]]]):
    uniqueDihedrals = {}
    for dihedralData in dihedralGroup:
        if dihedralData[0] not in uniqueDihedrals.keys():
            uniqueDihedrals[dihedralData[0]] = [dihedralData[1]]
        else:
            uniqueDihedrals[dihedralData[0]].append(dihedralData[1])

    return uniqueDihedrals


def _is_amide_dihedral(dihedral: parmed.topologyobjects.Dihedral) -> bool:
    atomTypes = extract_dihedral_atom_types(dihedral)
    if atomTypes[1] in ["N","NH1"] and atomTypes[2] == "C":
        return True
    elif atomTypes[2] in ["N","NH1"] and atomTypes[1] == "C":
        return True
    return False
def _is_a_ring_dihedral(dihedral: parmed.topologyobjects.Dihedral, rings: List[set[int]]) -> bool:
    atomIndexes = extract_dihedral_atom_indexes(dihedral)
    for ring in rings:
        if atomIndexes[1] in ring and atomIndexes[2] in ring:
            return True
    return False


def classify_rings_aromatic(rings: List[set[int]], parmedPsf: CharmmPsfFile, adjacencyMatrix: defaultdict) -> Tuple[List[set[int]], List[set[int]]]:
    aromaticValenceCheck = {
                        6 : [3],        ## Carbons MUST have 3 bonded atoms
                        7 : [2,3],      ## Nitrogens MUST have 2 or 3 bonded atoms
                        8 : [2]         ## Oxygens MUST have 2 bonded atoms

    }
    
    aromaticRings = []
    nonAromaticRings = []
    for ringIndex in rings:
        ringElements = [parmedPsf.atoms[idx].element for idx in ringIndex]
        ringValences = [len(adjacencyMatrix[idx]) for idx in ringIndex]  

        # possibleAromaticValences = [aromaticValenceCheck.get(element, (0)) for element in ringElements]

        passedCheckAtoms = [valence in aromaticValenceCheck.get(element, (0))
                            for valence, element in zip(ringValences, ringElements)]
        if all(passedCheckAtoms):
            aromaticRings.append(ringIndex)
        else:
            nonAromaticRings.append(ringIndex)
    return aromaticRings, nonAromaticRings

def find_ring_atoms(structure):
    """
    Identify atoms in rings within a ParmEd Structure object and return as a list of sets.
    Each set contains the indices of atoms in a distinct ring.
    Args:
        structure: ParmEd Structure object (e.g., loaded from a PSF file).
    Returns:
        list: List of sets, where each set contains atom indices (atom.idx) for a ring.
    """
    def dfs(current_atom, parent_atom, visited, path, rings, start_atom):
        """
        DFS to detect cycles and collect atoms in rings.
        """
        visited.add(current_atom)
        path.append(current_atom)

        for bond in current_atom.bonds:
            neighbor = bond.atom1 if bond.atom2 == current_atom else bond.atom2
            if neighbor == parent_atom:
                continue
            # If neighbor is in path and not the immediate parent, we found a cycle
            if neighbor in path:
                cycle_start_idx = path.index(neighbor)
                # Extract atoms from path back to neighbor
                cycle_atoms = path[cycle_start_idx:]
                # Store as a set of atom indices
                ring = set(atom.idx for atom in cycle_atoms)
                # Only add if not already found (avoid duplicates)
                if ring not in rings:
                    rings.append(ring)
            elif neighbor not in visited:
                # Continue DFS if neighbor hasn't been visited
                dfs(neighbor, current_atom, visited, path, rings, start_atom)

        path.pop()
        # Allow revisiting atoms for detecting multiple rings
        visited.remove(current_atom) # Uncomment if needed for complex graphs

    rings = []
    visited = set()

    # Run DFS from each atom to find all cycles
    for atom in structure.atoms:
        if atom not in visited:
            dfs(atom, None, visited.copy(), [], rings, atom)

    return rings


def _construct_adjacency_matrix(parmedPsf: CharmmPsfFile) -> np.ndarray:
    # Build adjacency list for graph representation
    adjacency = defaultdict(list)
    for bond in parmedPsf.bonds:
        idx1, idx2 = bond.atom1.idx, bond.atom2.idx
        adjacency[idx1].append(idx2)
        adjacency[idx2].append(idx1)

    return adjacency


def _is_polar_proton_dihedral(dihedral: parmed.topologyobjects.Dihedral) -> bool:

    atomElements = extract_dihedral_atom_elements(dihedral)
    if atomElements[0] == 1 and atomElements[1] != 6:
        return True
    elif atomElements[2] != 6 and atomElements[3] == 1:
        return True
    else:
        return False

def _is_non_polar_proton_dihedral(dihedral: parmed.topologyobjects.Dihedral) -> bool:

    atomElements = extract_dihedral_atom_elements(dihedral)
    if atomElements[0] == 1 and atomElements[1] == 6:
        return True
    elif atomElements[2] == 6 and atomElements[3] == 1:
        return True
    else:
        return False

def extract_dihedral_atom_types(dihedral: parmed.topologyobjects.Dihedral) -> Tuple[str,str,str,str]:
    return dihedral.atom1.type, dihedral.atom2.type, dihedral.atom3.type, dihedral.atom4.type

def extract_dihedral_atom_elements(dihedral: parmed.topologyobjects.Dihedral) -> Tuple[str,str,str,str]:
    return dihedral.atom1.element, dihedral.atom2.element, dihedral.atom3.element, dihedral.atom4.element

def extract_dihedral_atom_indexes(dihedral: parmed.topologyobjects.Dihedral) -> Tuple[str,str,str,str]:
    return dihedral.atom1.idx, dihedral.atom2.idx, dihedral.atom3.idx, dihedral.atom4.idx

def extract_dihedral_atom_names(dihedral: parmed.topologyobjects.Dihedral) -> Tuple[str,str,str,str]:
    return dihedral.atom1.name, dihedral.atom2.name, dihedral.atom3.name, dihedral.atom4.name

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

def boltzmann_weighted_average(row: pd.Series):

    ## set constants
    T = 300
    kB = 0.0019872041  # Boltzmann constant in kcal/mol·K
    beta = 1 / (kB * T)

    energies = row.drop(columns="Angle").values
    boltzmannFactors = np.exp(-beta * energies)
    weightedAverages = energies * boltzmannFactors  # Compute E_i * exp(-βE_i)
    boltzmannAverage = np.sum(weightedAverages) / np.sum(boltzmannFactors)
    return boltzmannAverage
    

def process_scan_data(scanDfs: List[pd.DataFrame],
                       torsionTopDir: DirectoryPath,
                       torsionTag: str):

    ## make a dir to store output csv files
    dataDir = p.join(torsionTopDir, "scan_data")
    os.makedirs(dataDir, exist_ok=True)

    ## merge scan dataframes
    # scanDfs = [process_energy_outputs(scanDf) for scanDf in scanDfs]
    mergedDf = merge_scan_dfs(scanDfs)

    ## write to csv
    mergedScanCsv = p.join(dataDir, f"scan_energies.csv")
    mergedDf.to_csv(mergedScanCsv, index=False)

    scanAverageDf = pd.DataFrame()
    scanAverageDf["Angle"] = mergedDf["Angle"]

    ## calculate averages
    finalScanEnergiesCsv = p.join(dataDir, "final_scan_energies.csv")
    ## swap out arithmetic mean for boltzmann-weighted average
    # scanAverageDf[torsionTag] = mergedDf.drop(columns="Angle").mean(axis=1)
    scanAverageDf[torsionTag] = mergedDf.apply(boltzmann_weighted_average, axis=1)
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




