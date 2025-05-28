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

    globalMinimaAngle = averagesDf["Angle"].loc[averagesDf[torsionTag].idxmin()]
    globalMinimaEnergy = averagesDf[torsionTag].loc[averagesDf[torsionTag].idxmin()]
    globalMaximaEnergy = averagesDf[torsionTag].loc[averagesDf[torsionTag].idxmax()]



    barrierHeight = round(globalMaximaEnergy - globalMinimaEnergy, 3)

    config["runtimeInfo"]["madeByTwisting"]["rotatableDihedrals"][torsionTag]["globalMinimaAngle"] = int(globalMinimaAngle)
    config["runtimeInfo"]["madeByTwisting"]["rotatableDihedrals"][torsionTag]["barrierHeight"] = float(barrierHeight)


    return config

def load_amber_params(config: dict) -> parmed.amber.AmberParameterSet:
    ## unpack config
    assembledPrmtop = config["runtimeInfo"]["madeByAssembly"]["assembledPrmtop"]

    parmedPrmtop = parmed.load_file(assembledPrmtop)

    return parmedPrmtop

def load_charmm_params(config: dict) -> CharmmPsfFile:
    ## unpack config
    assembledPsf = config["runtimeInfo"]["madeByAssembly"]["assembledPsf"]
    assembledPrm = config["runtimeInfo"]["madeByAssembly"]["assembledPrm"]
    assembledRtf = config["runtimeInfo"]["madeByAssembly"]["assembledRtf"]

    parmedPsf = CharmmPsfFile(assembledPsf)
    parmedPsf.load_parameters(CharmmParameterSet(assembledPrm, assembledRtf))

    return parmedPsf

def identify_rotatable_bonds(config: dict, mode: str = "AMBER") -> dict:

    if mode == "AMBER":
        moleculeParams = load_amber_params(config)
    else:
        moleculeParams = load_charmm_params(config)

    rotatableDihedrals = []
    aromaticDihedrals = []
    nonAromaticRingDihedrals = []
    terminalDihedrals = []

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
        ## store terminal non-polar dihedrals
        if _is_terminal_non_polar_dihedral(dihedral):
            terminalDihedrals.append((atomTypes, atomNames, atomIndexes))
            continue
        ## store aromatic dihedrals
        if _is_a_ring_dihedral(dihedral, rings):
            if _is_a_ring_dihedral(dihedral, aromaticRings):
                aromaticDihedrals.append((atomTypes, atomNames, atomIndexes))
                continue
            else:
                ## store non-aromatic ring dihedrals
                nonAromaticRingDihedrals.append((atomTypes, atomNames, atomIndexes))
                continue
        ## store rotatable dihedrals
        rotatableDihedrals.append((atomTypes, atomNames, atomIndexes))
        
    taggedRotatableDihedrals = assign_torsion_tags(rotatableDihedrals)
    taggedAromaticDihedrals = assign_torsion_tags(aromaticDihedrals)
    taggedNonAromaticRingDihedrals = assign_torsion_tags(nonAromaticRingDihedrals)
    taggedTerminalDihedrals = assign_torsion_tags(terminalDihedrals)

    config["runtimeInfo"]["madeByTwisting"]["rotatableDihedrals"] = taggedRotatableDihedrals
    config["runtimeInfo"]["madeByTwisting"]["aromaticDihedrals"] = taggedAromaticDihedrals
    config["runtimeInfo"]["madeByTwisting"]["nonAromaticRingDihedrals"] = taggedNonAromaticRingDihedrals
    config["runtimeInfo"]["madeByTwisting"]["terminalDihedrals"] = taggedTerminalDihedrals

    return config


def exclude_backbone_torsions(config: dict) -> dict:
    """
    Removes Phi and Psi Angles from the rotatable bonds list

    Args:
        config (dict): the config dict
    
    Returns:
        config (dict): the config dict updated
    """
    ## unpack config
    uniqueRotatableDihedrals = config["runtimeInfo"]["madeByTwisting"]["rotatableDihedrals"]
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
def get_conformer_xyzs(config: dict, seed: int = 1818) -> list[FilePath]:
    ## get conformer XYZ files
    conformerXyzs = config["runtimeInfo"]["madeByConformers"]["conformerXyzs"]
    nConformers = config["torsionScanInfo"]["nConformers"]

    if nConformers == -1 or nConformers > len(conformerXyzs):
        pass
    elif nConformers < len(conformerXyzs):
        random.seed(seed)
        conformerXyzs = random.sample(conformerXyzs, nConformers)
        
    return conformerXyzs

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

def create_orca_terminated_flag(orcaDir: DirectoryPath, orcaOut: FilePath) -> None:

    if drOrca.did_orca_finish_normallly(orcaOut):
        with open(p.join(orcaDir, "ORCA_FINISHED_NORMALLY"), "w") as f:
            f.write("ORCA FINISHED NORMALLY")
    else:
        with open(p.join(orcaDir, "ORCA_CRASHED"), "w") as f:
            f.write("ORCA CRASHED")


# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def read_singlepoint_energy(spOrcaOutput: FilePath) -> float:
    with open(spOrcaOutput, "r") as f:
        for line in f:
            if line.startswith("FINAL SINGLE POINT ENERGY"):
                singlePointEnergy = float(line.split()[-1])
                return singlePointEnergy
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def add_mid_points(indexes: np.ndarray) -> list[int]:
    newIndexes = []
    for i, indexA in enumerate(indexes[:-1]):
        newIndexes.append(indexA)
        indexB = indexes[i+1]
        midPointIndex = (indexA + indexB) // 2
        newIndexes.append(midPointIndex)
    newIndexes.append(indexes[-1])
    
    return newIndexes
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def find_local_extrema(energies: pd.Series) -> list[int]:
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
def find_scan_xyz_files(scanDir: DirectoryPath, expectedNumberOfFiles: int) -> list[FilePath]:
    scanXyzs = sorted([xyzFile for xyzFile in glob.glob(p.join(scanDir, 'orca_scan.[0-9][0-9][0-9].xyz'))])
    if not len(scanXyzs) == expectedNumberOfFiles:
        raise(ValueError(f"Number of scan xyz files ({len(scanXyzs)}) does not match expected ({expectedNumberOfFiles})"))

    return scanXyzs
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def process_energy_outputs(energyDf: pd.DataFrame) -> pd.DataFrame:
    energyDf = rescale_and_sort_energy_angles(energyDf)
    energyDf = take_min_duplicate_angles(energyDf)

    energyDf["Energy"] = energyDf["Energy"].apply(hartree_to_kcal_per_mol)
    return energyDf
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def hartree_to_kcal_per_mol(energy: float) -> float:
    return energy * 627.5095
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def rescale_torsion_angles(angle: float) -> float:
    angle = angle % 360  # reduce the angle to the 0-360 range
    return angle
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def rescale_and_sort_energy_angles(energyDf: pd.DataFrame) -> pd.DataFrame:
    energyDf["Energy"] = energyDf["Energy"] - energyDf["Energy"].min()
    energyDf["Angle"] = energyDf["Angle"].apply(rescale_torsion_angles)
    energyDf = energyDf.sort_values(by="Angle", ascending=True)
    return energyDf
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def take_min_duplicate_angles(energyDf: pd.DataFrame) -> pd.DataFrame:
    energyDf['Angle'] = energyDf['Angle'].round(0)
    # Find the index of the minimum energy for each unique Angle
    minEnergyIndexes = energyDf.groupby('Angle')["Energy"].idxmin()

    minEnergyIndexes = minEnergyIndexes.dropna()

    # Select the rows with the minimum energy for each Angle
    return energyDf.loc[minEnergyIndexes].reset_index(drop=True)
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def read_scan_energy_data(conformerTorsionScanDir: DirectoryPath) -> pd.DataFrame:
    scanDat: FilePath = p.join(conformerTorsionScanDir, "orca_scan.relaxscanact.dat")
    if not os.path.exists(scanDat):
        raise FileNotFoundError(f"Scan data not found in {conformerTorsionScanDir}")
    scanDf: pd.DataFrame = pd.read_csv(scanDat, sep='\s+', header=None)
    return scanDf
    
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def find_final_xyz(conformerTorsionScanDir: DirectoryPath) -> FilePath:
    allXyzFiles = glob.glob(p.join(conformerTorsionScanDir, "*.xyz"))
    scanXYZFiles = sorted([f for f in allXyzFiles if re.match(r'orca_scan\.\d+\.xyz$', os.path.basename(f))])

    finalXyzFile = scanXYZFiles[-1]
    return finalXyzFile
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def measure_current_torsion_angle(conformerXyz: FilePath, torsionIndexes: list[int]) -> float:
    atomCoords = xyz_to_np_array(conformerXyz)
    torsionAngle = calculate_torsion_angle(atomCoords, torsionIndexes)
    roundedToTenAngle = round(torsionAngle, -1)
    return roundedToTenAngle
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def xyz_to_np_array(filePath: FilePath) -> np.ndarray:
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
def calculate_torsion_angle(coords: np.ndarray, torsionIndexes: list[int]) -> float:
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
    torsionTopDir = p.join(outputDir, "04_torsion_scanning")
    os.makedirs(torsionTopDir, exist_ok=True)
    config["runtimeInfo"]["madeByTwisting"]["torsionDir"] = torsionTopDir

    return config

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

def get_non_symmetric_rotatable_bonds(rotatable_bonds: list[dict], config: dict) -> list[dict]:
    ## unpack config
    cappedPdb = config["runtimeInfo"]["madeByCapping"]["cappedPdb"]
    symmetryData = config["runtimeInfo"]["madeByCapping"]["symmetryData"]

    ## get groups of atoms that are symmetrically equivalent
    symmetricAtoms = [symmetryGroup for symmetryGroup in symmetryData if len(symmetryGroup) > 1]
    ## create a dict that assigns an int to each group of symmetrical atoms
    enumeratedSymmetryGroups = {i: atoms for i, atoms in enumerate(symmetricAtoms)}
    ## create a dict that maps each atom to its assigned int
    atomToSymmetryGroupIndex = {}
    for key, atoms in enumeratedSymmetryGroups.items():
        for atom in atoms:
            atomToSymmetryGroupIndex[atom] = key

    ## Replace atom names with assigned int
    for bond in rotatable_bonds:
        bond['atoms'] = tuple(atomToSymmetryGroupIndex.get(atom, atom) for atom in bond['atoms'])

    ## get unique rotatable  bonds
    uniqueBonds = []
    seenAtoms = set()
    for bond in rotatable_bonds:
        atoms = bond['atoms']
        if atoms not in seenAtoms:
            seenAtoms.add(atoms)
            uniqueBonds.append(bond)

    ## made a dict that maps assigned int back to the first atom name in each group
    symmetryGroupIndexToFirstAtom = {}
    for key, atoms in enumeratedSymmetryGroups.items():
        firstAtom = sorted(atoms)[0]
        for atom in atoms:
            symmetryGroupIndexToFirstAtom[key] = firstAtom
    ## replace assigned int with first atom in each group
    for bond in uniqueBonds:
        bond['atoms'] = tuple(symmetryGroupIndexToFirstAtom.get(atom,atom) for atom in bond['atoms'])

    return uniqueBonds



 

def assign_torsion_tags(dihedrals: list[tuple[tuple[str], tuple[str], tuple[int]]]) -> dict[str, dict[str, tuple]]:
    taggedDihedrals = {}
    seenTags = set()
    for dihedralData in dihedrals:
        torsionTag = "-".join(dihedralData[0]) # Atom types for tag
        if torsionTag not in seenTags:
            seenTags.add(torsionTag)
            taggedDihedrals[torsionTag] = {"ATOM_TYPES": dihedralData[0], 
                                           "ATOM_NAMES": dihedralData[1], 
                                           "ATOM_INDEXES": dihedralData[2]}
        else:
           continue
    return taggedDihedrals



def get_unique_dihedrals(dihedral_group: list[tuple[tuple[str], tuple[str], tuple[int]]]) -> dict[tuple[str], list[tuple[str]]]:
    uniqueDihedrals = {}
    # dihedralData is expected to be (atom_types_tuple, atom_names_tuple, atom_indexes_tuple)
    for atom_types_tuple, atom_names_tuple, _ in dihedral_group: # atom_indexes_tuple is ignored for now
        if atom_types_tuple not in uniqueDihedrals.keys():
            uniqueDihedrals[atom_types_tuple] = [atom_names_tuple]
        else:
            uniqueDihedrals[atom_types_tuple].append(atom_names_tuple)

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


def classify_rings_aromatic(rings: List[set[int]], parmedPsf: CharmmPsfFile, adjacency_matrix: defaultdict) -> Tuple[List[set[int]], List[set[int]]]:
    aromaticValenceCheck = {
                        6 : [3],        ## Carbons MUST have 3 bonded atoms
                        7 : [2,3],      ## Nitrogens MUST have 2 or 3 bonded atoms
                        8 : [2]         ## Oxygens MUST have 2 bonded atoms

    }
    
    aromaticRings = []
    nonAromaticRings = []
    for ringIndex in rings:
        ringElements = [parmedPsf.atoms[idx].element for idx in ringIndex]
        ringValences = [len(adjacency_matrix[idx]) for idx in ringIndex]  

        # possibleAromaticValences = [aromaticValenceCheck.get(element, (0)) for element in ringElements]

        passedCheckAtoms = [valence in aromaticValenceCheck.get(element, (0))
                            for valence, element in zip(ringValences, ringElements)]
        if all(passedCheckAtoms):
            aromaticRings.append(ringIndex)
        else:
            nonAromaticRings.append(ringIndex)
    return aromaticRings, nonAromaticRings

def find_ring_atoms(structure: parmed.Structure) -> List[set[int]]:
    """
    Identify atoms in rings within a ParmEd Structure object and return as a list of sets.
    Each set contains the indices of atoms in a distinct ring.
    Args:
        structure: ParmEd Structure object (e.g., loaded from a PSF file).
    Returns:
        list: List of sets, where each set contains atom indices (atom.idx) for a ring.
    """
    def dfs(current_atom: parmed.topologyobjects.Atom, 
            parent_atom: parmed.topologyobjects.Atom | None, 
            visited: set[parmed.topologyobjects.Atom], 
            path: list[parmed.topologyobjects.Atom], 
            rings: List[set[int]], 
            start_atom: parmed.topologyobjects.Atom) -> None:
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
        # visited.remove(current_atom) # This line was commented out, keeping it as is.

    rings_collector: List[set[int]] = []
    visited_atoms: set[parmed.topologyobjects.Atom] = set()

    # Run DFS from each atom to find all cycles
    for atom_obj in structure.atoms:
        if atom_obj not in visited_atoms:
            # For each new DFS traversal, we need a fresh visited set for that traversal's context,
            # but the overall visited_atoms set ensures we don't restart DFS from already processed components.
            dfs(atom_obj, None, set(), [], rings_collector, atom_obj)
            # Add all atoms visited in this DFS traversal to the overall visited_atoms set
            # This part is tricky because the 'visited' set in dfs is modified.
            # A simpler approach for ensuring atoms are not re-processed by initiating DFS:
            # After a DFS completes from 'atom_obj', all atoms reachable from 'atom_obj' could be considered 'processed'
            # for the purpose of initiating new DFS traversals.
            # However, the original code allows restarting DFS for multi-ring detection,
            # so we might not need to add to visited_atoms here if rings_collector is the main output.
            # The original `visited.remove(current_atom)` being commented out also affects this.
            # For now, let's assume `rings_collector` correctly accumulates rings without needing complex visited management here.
            pass


    return rings_collector


def _construct_adjacency_matrix(parmedPsf: CharmmPsfFile) -> defaultdict[int, list[int]]:
    # Build adjacency list for graph representation
    adjacency = defaultdict(list)
    for bond in parmedPsf.bonds:
        idx1, idx2 = bond.atom1.idx, bond.atom2.idx
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

def extract_dihedral_atom_types(dihedral: parmed.topologyobjects.Dihedral) -> Tuple[str,str,str,str]:
    return dihedral.atom1.type, dihedral.atom2.type, dihedral.atom3.type, dihedral.atom4.type

def extract_dihedral_atom_elements(dihedral: parmed.topologyobjects.Dihedral) -> Tuple[str,str,str,str]:
    return dihedral.atom1.element, dihedral.atom2.element, dihedral.atom3.element, dihedral.atom4.element

def extract_dihedral_atom_indexes(dihedral: parmed.topologyobjects.Dihedral) -> Tuple[str,str,str,str]:
    return dihedral.atom1.idx, dihedral.atom2.idx, dihedral.atom3.idx, dihedral.atom4.idx

def extract_dihedral_atom_names(dihedral: parmed.topologyobjects.Dihedral) -> Tuple[str,str,str,str]:
    return dihedral.atom1.name, dihedral.atom2.name, dihedral.atom3.name, dihedral.atom4.name

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

def process_scan_data(scan_dfs: List[pd.DataFrame],
                       torsion_top_dir: DirectoryPath,
                       torsion_tag: str) -> Tuple[FilePath, pd.DataFrame]:

    ## make a dir to store output csv files
    dataDir = p.join(torsion_top_dir, "scan_data")
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
    scanAverageDf[torsion_tag] = mergedDf.drop(columns="Angle").mean(axis=1)
    scanAverageDf.to_csv(finalScanEnergiesCsv, index=False)

    return finalScanEnergiesCsv, scanAverageDf
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

def merge_scan_dfs(scan_dfs: List[pd.DataFrame]) -> pd.DataFrame:
    mergedDf = pd.DataFrame()
    mergedDf["Angle"] = np.arange(0,360,10)
    for colIndex, scanDf in enumerate(scan_dfs):
        mergedDf = mergedDf.merge(scanDf[["Angle", "Energy"]], on = "Angle", how= "left")
        mergedDf.rename(columns={"Energy":f"Energy_{colIndex + 1}"}, inplace=True)
    return mergedDf
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

def detect_jumps_in_data(df: pd.DataFrame) -> pd.DataFrame:
    # Calculate differences between consecutive values
    diffDf = df.drop(columns='Angle').diff().abs()
    
    # Identify columns with any difference greater than the threshold
    jumpyCols = diffDf.columns[((diffDf > 10).any())]
    
    # Drop these columns from the original DataFrame
    cleanDf = df.drop(columns=jumpyCols)
    
    return cleanDf
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲



