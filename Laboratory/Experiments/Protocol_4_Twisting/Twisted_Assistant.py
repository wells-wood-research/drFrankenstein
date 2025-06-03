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

class FilePath:
    pass

class DirectoryPath:
    pass
from typing import List, Tuple
from OperatingTools import drOrca

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
    globalMinimaAngle = averagesDf['Angle'].loc[averagesDf[torsionTag].idxmin()]
    globalMinimaEnergy = averagesDf[torsionTag].loc[averagesDf[torsionTag].idxmin()]
    globalMaximaEnergy = averagesDf[torsionTag].loc[averagesDf[torsionTag].idxmax()]
    barrierHeight = round(globalMaximaEnergy - globalMinimaEnergy, 3)
    config['runtimeInfo']['madeByTwisting']['rotatableDihedrals'][torsionTag]['globalMinimaAngle'] = int(globalMinimaAngle)
    config['runtimeInfo']['madeByTwisting']['rotatableDihedrals'][torsionTag]['barrierHeight'] = float(barrierHeight)
    return config

def load_amber_params(config):
    assembledPrmtop = config['runtimeInfo']['madeByAssembly']['assembledPrmtop']
    parmedPrmtop = parmed.load_file(assembledPrmtop)
    return parmedPrmtop

def load_charmm_params(config):
    assembledPsf = config['runtimeInfo']['madeByAssembly']['assembledPsf']
    assembledPrm = config['runtimeInfo']['madeByAssembly']['assembledPrm']
    assembledRtf = config['runtimeInfo']['madeByAssembly']['assembledRtf']
    parmedPsf = CharmmPsfFile(assembledPsf)
    parmedPsf.load_parameters(CharmmParameterSet(assembledPrm, assembledRtf))
    return parmedPsf

def identify_rotatable_bonds(config: dict, mode: str='AMBER') -> List[Tuple[int, int, int, int]]:
    if mode == 'AMBER':
        moleculeParams = load_amber_params(config)
    else:
        moleculeParams = load_charmm_params(config)
    rotatableDihedrals = []
    aromaticDihedrals = []
    nonAromaticRingDihedrals = []
    terminalDihedrals = []
    rings = find_ring_atoms(moleculeParams)
    adjacencyMatrix = _construct_adjacency_matrix(moleculeParams)
    (aromaticRings, nonAromaticRings) = classify_rings_aromatic(rings, moleculeParams, adjacencyMatrix)
    for dihedral in moleculeParams.dihedrals:
        atomNames = extract_dihedral_atom_names(dihedral)
        atomTypes = extract_dihedral_atom_types(dihedral)
        atomIndexes = extract_dihedral_atom_indexes(dihedral)
        if _is_amide_dihedral(dihedral):
            continue
        if _is_terminal_non_polar_dihedral(dihedral):
            terminalDihedrals.append((atomTypes, atomNames, atomIndexes))
            continue
        if _is_a_ring_dihedral(dihedral, rings):
            if _is_a_ring_dihedral(dihedral, aromaticRings):
                aromaticDihedrals.append((atomTypes, atomNames, atomIndexes))
                continue
            else:
                nonAromaticRingDihedrals.append((atomTypes, atomNames, atomIndexes))
                continue
        rotatableDihedrals.append((atomTypes, atomNames, atomIndexes))
    taggedRotatableDihedrals = assign_torsion_tags(rotatableDihedrals)
    taggedAromaticDihedrals = assign_torsion_tags(aromaticDihedrals)
    taggedNonAromaticRingDihedrals = assign_torsion_tags(nonAromaticRingDihedrals)
    taggedTerminalDihedrals = assign_torsion_tags(terminalDihedrals)
    config['runtimeInfo']['madeByTwisting']['rotatableDihedrals'] = taggedRotatableDihedrals
    config['runtimeInfo']['madeByTwisting']['aromaticDihedrals'] = taggedAromaticDihedrals
    config['runtimeInfo']['madeByTwisting']['nonAromaticRingDihedrals'] = taggedNonAromaticRingDihedrals
    config['runtimeInfo']['madeByTwisting']['terminalDihedrals'] = taggedTerminalDihedrals
    return config

def exclude_backbone_torsions(config: dict) -> dict:
    """
    Removes Phi and Psi Angles from the rotatable bonds list

    Args:
        config (dict): the config dict
    
    Returns:
        config (dict): the config dict updated
    """
    uniqueRotatableDihedrals = config['runtimeInfo']['madeByTwisting']['rotatableDihedrals']
    forceFeild = config['parameterFittingInfo']['forceField']
    if forceFeild == 'CHARMM':
        phiCenterTypes = ('NH1', 'CT1')
        psiCenterTypes = ('CT1', 'C')
    elif forceFeild == 'AMBER':
        phiCenterTypes = ('N', 'CT')
        psiCenterTypes = ('CT', 'C')
    tagsToRemove = []
    for (torsionTag, dihedralData) in uniqueRotatableDihedrals.items():
        dihedralCenterTypes = dihedralData['ATOM_TYPES'][1:3]
        if dihedralCenterTypes == phiCenterTypes or dihedralCenterTypes == phiCenterTypes[::-1]:
            tagsToRemove.append(torsionTag)
        elif dihedralCenterTypes == psiCenterTypes or dihedralCenterTypes == psiCenterTypes[::-1]:
            tagsToRemove.append(torsionTag)
    for torsionTag in tagsToRemove:
        uniqueRotatableDihedrals.pop(torsionTag)
    config['runtimeInfo']['madeByTwisting']['uniqueRotatableDihedrals'] = uniqueRotatableDihedrals
    return config

def get_conformer_xyzs(config, seed=1818):
    conformerXyzs = config['runtimeInfo']['madeByConformers']['conformerXyzs']
    nConformers = config['torsionScanInfo']['nConformers']
    if nConformers == -1 or nConformers > len(conformerXyzs):
        pass
    elif nConformers < len(conformerXyzs):
        random.seed(seed)
        conformerXyzs = random.sample(conformerXyzs, nConformers)
    return conformerXyzs

def create_orca_terminated_flag(orcaDir, orcaOut):
    if drOrca.did_orca_finish_normallly(orcaOut):
        with open(p.join(orcaDir, 'ORCA_FINISHED_NORMALLY'), 'w') as f:
            f.write('ORCA FINISHED NORMALLY')
    else:
        with open(p.join(orcaDir, 'ORCA_CRASHED'), 'w') as f:
            f.write('ORCA CRASHED')

def read_singlepoint_energy(spOrcaOutput):
    with open(spOrcaOutput, 'r') as f:
        for line in f:
            if line.startswith('FINAL SINGLE POINT ENERGY'):
                singlePointEnergy = float(line.split()[-1])
                return singlePointEnergy

def add_mid_points(indexes: np.array):
    newIndexes = []
    for (i, indexA) in enumerate(indexes[:-1]):
        newIndexes.append(indexA)
        indexB = indexes[i + 1]
        midPointIndex = (indexA + indexB) // 2
        newIndexes.append(midPointIndex)
    newIndexes.append(indexes[-1])
    return newIndexes

def find_local_extrema(energies):
    energiesArray = energies.to_numpy()
    extendedEnergies = np.concatenate((energiesArray[-1:], energiesArray, energiesArray[:1]))
    localMinimaIndexes = argrelextrema(extendedEnergies, np.less)[0] - 1
    localMaximaIndexes = argrelextrema(extendedEnergies, np.greater)[0] - 1
    localMinimaIndexes = localMinimaIndexes[(localMinimaIndexes >= 0) & (localMinimaIndexes < len(energiesArray))]
    localMaximaIndexes = localMaximaIndexes[(localMaximaIndexes >= 0) & (localMaximaIndexes < len(energiesArray))]
    combinedExtremaIndexes = sorted(np.concatenate((localMinimaIndexes, localMaximaIndexes)))
    return [int(index) for index in combinedExtremaIndexes]

def find_scan_xyz_files(scanDir: DirectoryPath, expectedNumberOfFiles: int):
    scanXyzs = sorted([xyzFile for xyzFile in glob.glob(p.join(scanDir, 'orca_scan.[0-9][0-9][0-9].xyz'))])
    if not len(scanXyzs) == expectedNumberOfFiles:
        raise ValueError(f'Number of scan xyz files ({len(scanXyzs)}) does not match expected ({expectedNumberOfFiles})')
    return scanXyzs

def process_energy_outputs(energyDf):
    energyDf = rescale_and_sort_energy_angles(energyDf)
    energyDf = take_min_duplicate_angles(energyDf)
    energyDf['Energy'] = energyDf['Energy'].apply(hartree_to_kcal_per_mol)
    return energyDf

def hartree_to_kcal_per_mol(energy):
    return energy * 627.5095

def rescale_torsion_angles(angle):
    angle = angle % 360
    return angle

def rescale_and_sort_energy_angles(energyDf):
    energyDf['Energy'] = energyDf['Energy'] - energyDf['Energy'].min()
    energyDf['Angle'] = energyDf['Angle'].apply(rescale_torsion_angles)
    energyDf = energyDf.sort_values(by='Angle', ascending=True)
    return energyDf

def take_min_duplicate_angles(energyDf):
    energyDf['Angle'] = energyDf['Angle'].round(0)
    minEnergyIndexes = energyDf.groupby('Angle')['Energy'].idxmin()
    minEnergyIndexes = minEnergyIndexes.dropna()
    return energyDf.loc[minEnergyIndexes].reset_index(drop=True)

def read_scan_energy_data(conformertorsionScanDir):
    scanDat: FilePath = p.join(conformertorsionScanDir, 'orca_scan.relaxscanact.dat')
    if not os.path.exists(scanDat):
        raise FileNotFoundError(f'Scan data not found in {conformertorsionScanDir}')
    scanDf: pd.DataFrame = pd.read_csv(scanDat, sep='\\s+', header=None)
    return scanDf

def find_final_xyz(conformertorsionScanDir):
    allXyzFiles = glob.glob(p.join(conformertorsionScanDir, '*.xyz'))
    scanXYZFiles = sorted([f for f in allXyzFiles if re.match('orca_scan\\.\\d+\\.xyz$', os.path.basename(f))])
    finalXyzFile = scanXYZFiles[-1]
    return finalXyzFile

def measure_current_torsion_angle(conformerXyz, torsionIndexes):
    atomCoords = xyz_to_np_array(conformerXyz)
    torsionAngle = calculate_torsion_angle(atomCoords, torsionIndexes)
    roundedToTenAngle = round(torsionAngle, -1)
    return roundedToTenAngle

def xyz_to_np_array(filePath):
    with open(filePath, 'r') as file:
        lines = file.readlines()
    dataLines = lines[2:]
    coords = []
    for line in dataLines:
        parts = line.split()
        if len(parts) >= 4:
            (x, y, z) = map(float, parts[1:4])
            coords.append([x, y, z])
    return np.array(coords)

def calculate_torsion_angle(coords, torsionIndexes):
    b1 = coords[torsionIndexes[1]] - coords[torsionIndexes[0]]
    b2 = coords[torsionIndexes[2]] - coords[torsionIndexes[1]]
    b3 = coords[torsionIndexes[3]] - coords[torsionIndexes[2]]
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    n1 /= np.linalg.norm(n1)
    n2 /= np.linalg.norm(n2)
    angleInRadians = np.arctan2(np.dot(np.cross(n1, n2), b2 / np.linalg.norm(b2)), np.dot(n1, n2))
    angleIdDegrees = np.degrees(angleInRadians)
    return angleIdDegrees

def set_up_directories(config: dict) -> dict:
    outputDir = config['pathInfo']['outputDir']
    torsionTopDir = p.join(outputDir, '04_torsion_scanning')
    os.makedirs(torsionTopDir, exist_ok=True)
    config['runtimeInfo']['madeByTwisting']['torsionDir'] = torsionTopDir
    return config

def assign_torsion_tags(dihedrals: list[tuple[tuple[str], tuple[str]]]) -> dict[tuple[str], tuple[str], tuple[str]]:
    taggedDihedrals = {}
    seenTags = set()
    for dihedralData in dihedrals:
        torsionTag = '-'.join(dihedralData[0])
        if torsionTag not in seenTags:
            seenTags.add(torsionTag)
            taggedDihedrals[torsionTag] = {'ATOM_TYPES': dihedralData[0], 'ATOM_NAMES': dihedralData[1], 'ATOM_INDEXES': dihedralData[2]}
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
    if atomTypes[1] in ['N', 'NH1'] and atomTypes[2] == 'C':
        return True
    elif atomTypes[2] in ['N', 'NH1'] and atomTypes[1] == 'C':
        return True
    return False

def _is_a_ring_dihedral(dihedral: parmed.topologyobjects.Dihedral, rings: List[set[int]]) -> bool:
    atomIndexes = extract_dihedral_atom_indexes(dihedral)
    for ring in rings:
        if atomIndexes[1] in ring and atomIndexes[2] in ring:
            return True
    return False

def classify_rings_aromatic(rings: List[set[int]], parmedPsf: CharmmPsfFile, adjacencyMatrix: defaultdict) -> Tuple[List[set[int]], List[set[int]]]:
    aromaticValenceCheck = {6: [3], 7: [2, 3], 8: [2]}
    aromaticRings = []
    nonAromaticRings = []
    for ringIndex in rings:
        ringElements = [parmedPsf.atoms[idx].element for idx in ringIndex]
        ringValences = [len(adjacencyMatrix[idx]) for idx in ringIndex]
        passedCheckAtoms = [valence in aromaticValenceCheck.get(element, 0) for (valence, element) in zip(ringValences, ringElements)]
        if all(passedCheckAtoms):
            aromaticRings.append(ringIndex)
        else:
            nonAromaticRings.append(ringIndex)
    return (aromaticRings, nonAromaticRings)

def find_ring_atoms(structure):
    """
    Identify atoms in rings within a ParmEd Structure object and return as a list of sets.
    Each set contains the indices of atoms in a distinct ring.
    Args:
        structure: ParmEd Structure object (e.g., loaded from a PSF file).
    Returns:
        list: List of sets, where each set contains atom indices (atom.idx) for a ring.
    """

    def dfs(currentAtom, parentAtom, visited, path, rings, startAtom):
        """
        DFS to detect cycles and collect atoms in rings.
        """
        visited.add(currentAtom)
        path.append(currentAtom)
        for bond in currentAtom.bonds:
            neighbor = bond.atom1 if bond.atom2 == currentAtom else bond.atom2
            if neighbor == parentAtom:
                continue
            if neighbor in path:
                cycleStartIdx = path.index(neighbor)
                cycleAtoms = path[cycleStartIdx:]
                ring = set((atom.idx for atom in cycleAtoms))
                if ring not in rings:
                    rings.append(ring)
            elif neighbor not in visited:
                dfs(neighbor, currentAtom, visited, path, rings, startAtom)
        path.pop()
        visited.remove(currentAtom)
    rings = []
    visited = set()
    for atom in structure.atoms:
        if atom not in visited:
            dfs(atom, None, visited.copy(), [], rings, atom)
    return rings

def _construct_adjacency_matrix(parmedPsf: CharmmPsfFile) -> np.ndarray:
    adjacency = defaultdict(list)
    for bond in parmedPsf.bonds:
        (idx1, idx2) = (bond.atom1.idx, bond.atom2.idx)
        adjacency[idx1].append(idx2)
        adjacency[idx2].append(idx1)
    return adjacency

def _is_terminal_non_polar_dihedral(dihedral: parmed.topologyobjects.Dihedral) -> bool:
    atomElements = extract_dihedral_atom_elements(dihedral)
    if atomElements[0] == 1 and atomElements[1] == 6:
        return True
    elif atomElements[2] == 6 and atomElements[3] == 1:
        return True
    else:
        return False

def extract_dihedral_atom_types(dihedral: parmed.topologyobjects.Dihedral) -> Tuple[str, str, str, str]:
    return (dihedral.atom1.type, dihedral.atom2.type, dihedral.atom3.type, dihedral.atom4.type)

def extract_dihedral_atom_elements(dihedral: parmed.topologyobjects.Dihedral) -> Tuple[str, str, str, str]:
    return (dihedral.atom1.element, dihedral.atom2.element, dihedral.atom3.element, dihedral.atom4.element)

def extract_dihedral_atom_indexes(dihedral: parmed.topologyobjects.Dihedral) -> Tuple[str, str, str, str]:
    return (dihedral.atom1.idx, dihedral.atom2.idx, dihedral.atom3.idx, dihedral.atom4.idx)

def extract_dihedral_atom_names(dihedral: parmed.topologyobjects.Dihedral) -> Tuple[str, str, str, str]:
    return (dihedral.atom1.name, dihedral.atom2.name, dihedral.atom3.name, dihedral.atom4.name)

def process_scan_data(scanDfs: List[pd.DataFrame], torsionTopDir: DirectoryPath, torsionTag: str):
    dataDir = p.join(torsionTopDir, 'scan_data')
    os.makedirs(dataDir, exist_ok=True)
    mergedDf = merge_scan_dfs(scanDfs)
    mergedDf = detect_jumps_in_data(mergedDf)
    mergedScanCsv = p.join(dataDir, f'scan_energies.csv')
    mergedDf.to_csv(mergedScanCsv, index=False)
    scanAverageDf = pd.DataFrame()
    scanAverageDf['Angle'] = mergedDf['Angle']
    finalScanEnergiesCsv = p.join(dataDir, 'final_scan_energies.csv')
    scanAverageDf[torsionTag] = mergedDf.drop(columns='Angle').mean(axis=1)
    scanAverageDf.to_csv(finalScanEnergiesCsv, index=False)
    return (finalScanEnergiesCsv, scanAverageDf)

def merge_scan_dfs(scanDfs):
    mergedDf = pd.DataFrame()
    mergedDf['Angle'] = np.arange(0, 360, 10)
    for (colIndex, scanDf) in enumerate(scanDfs):
        mergedDf = mergedDf.merge(scanDf[['Angle', 'Energy']], on='Angle', how='left')
        mergedDf.rename(columns={'Energy': f'Energy_{colIndex + 1}'}, inplace=True)
    return mergedDf

def detect_jumps_in_data(df):
    diffDf = df.drop(columns='Angle').diff().abs()
    jumpyCols = diffDf.columns[(diffDf > 10).any()]
    cleanDf = df.drop(columns=jumpyCols)
    return cleanDf