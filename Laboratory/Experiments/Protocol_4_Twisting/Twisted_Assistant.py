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
## RDKIT IMPORTS ##
from rdkit import Chem, RDLogger


## CLEAN CODE CLASSES ##
class FilePath:
    pass
class DirectoryPath:
    pass
from typing import List, Tuple


RDLogger.DisableLog('rdApp.warning')

from OperatingTools import drOrca
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def get_conformer_xyzs(config, seed = 1818):
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
    config["runtimeInfo"]["madeByTwisting"]["torsionDir"] = torsionTopDir

    return config

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

def get_non_symmetric_rotatable_bonds(rotatableBonds, config):
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
    for bond in rotatableBonds:
        bond['atoms'] = tuple(atomToSymmetryGroupIndex.get(atom, atom) for atom in bond['atoms'])

    ## get unique rotatable  bonds
    uniqueBonds = []
    seenAtoms = set()
    for bond in rotatableBonds:
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


def identify_rotatable_bonds_CHARMM(config) -> List[Tuple[int,int,int,int]]:
    ## unpack config
    assembledPsf = config["runtimeInfo"]["madeByAssembly"]["assembledPsf"]
    assembledPrm = config["runtimeInfo"]["madeByAssembly"]["assembledPrm"]
    assembledRtf = config["runtimeInfo"]["madeByAssembly"]["assembledRtf"]

    parmedPsf = CharmmPsfFile(assembledPsf)
    parmedPsf.load_parameters(CharmmParameterSet(assembledPrm, assembledRtf))
    # Create bond lookup dictionary
    bond_lookup = {}
    for bond in parmedPsf.bonds:
        # Use sorted indices as key to handle undirected bonds
        key = tuple(sorted([bond.atom1.idx, bond.atom2.idx]))
        bond_lookup[key] = bond

    rotatableDihedrals = []
    aromaticDihedrals = []
    nonAromaticRingDihedrals = []
    terminalDihedrals = []

    rings = find_ring_atoms(parmedPsf)
    adjacencyMatrix = _construct_adjacency_matrix(parmedPsf)

    aromaticRings, nonAromaticRings = classify_rings_aromatic(rings, parmedPsf, adjacencyMatrix)

    for dihedral in parmedPsf.dihedrals:
        atomNames = extract_dihedral_atom_names(dihedral)
        atomTypes = extract_dihedral_atom_types(dihedral)

        if _is_terminal_non_polar_dihedral(dihedral):
            terminalDihedrals.append((atomTypes, atomNames))
            continue

        if _is_a_ring_dihedral(dihedral, rings):
            if _is_a_ring_dihedral(dihedral, aromaticRings):
                aromaticDihedrals.append((atomTypes, atomNames))
                continue
            else:
                nonAromaticRingDihedrals.append((atomTypes, atomNames))
                continue
        rotatableDihedrals.append((atomTypes, atomNames))
    uniqueRotatableDihedrals = get_unique_dihedrals(rotatableDihedrals)
    for key, value in uniqueRotatableDihedrals.items():
        print(key, value, len(value))

    exit()
 
def get_unique_dihedrals(dihedralGroup: list[tuple[tuple[str], tuple[str]]]):
    uniqueDihedrals = {}

    for dihedral in dihedralGroup:
        if not dihedral[0] in uniqueDihedrals.keys():
            uniqueDihedrals[dihedral[0]] = [dihedral[1]]
        else:
            uniqueDihedrals[dihedral[0]].append(dihedral[1])

    return uniqueDihedrals
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

        print(ringElements, ringValences)
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


def _is_terminal_non_polar_dihedral(dihedral: parmed.topologyobjects.Dihedral) -> bool:

    atomElements = extract_dihedral_atom_elements(dihedral)
    print(atomElements)
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

def identify_rotatable_bonds(config) -> List[Tuple[int,int,int,int]]:
    ## unpack config
    cappedPdb = config["runtimeInfo"]["madeByCapping"]["cappedPdb"]
    # Load the molecule from a PDB file
    mol = Chem.MolFromPDBFile(cappedPdb, removeHs=False)
    # Identify torsion angles for rotatable bonds
    rotatableBonds = []
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
            nTerminalAtomNames = config["moleculeInfo"]["nTermini"]
            cTerminalAtomNames = config["moleculeInfo"]["cTermini"]
            nTerminalAmideAtoms = nTerminalAtomNames + ["CC1"]
            cTerminalAmideAtoms = cTerminalAtomNames + ["NN"]
            if atom2Name in nTerminalAmideAtoms and atom3Name in nTerminalAmideAtoms:
                continue
            if atom2Name in cTerminalAmideAtoms and atom3Name in cTerminalAmideAtoms:
                continue
            
            if not (atom2.IsInRing() or atom3.IsInRing()):
                # Find neighboring atoms for torsion angle
                neighborsBegin = [a for a in atom2.GetNeighbors() if a.GetIdx() != atom3.GetIdx()]
                neighborsEnd = [a for a in atom3.GetNeighbors() if a.GetIdx() != atom2.GetIdx()]
                
                for atom1 in neighborsBegin:
                    for atom4 in neighborsEnd:

                        atom1Name = atom1.GetPDBResidueInfo().GetName().strip()
                        atom4Name = atom4.GetPDBResidueInfo().GetName().strip()

                        # dont scan bonds with non-polar hydrogens at either as atoms 1 or 4
                        if atom1Name.startswith("H"):
                            if atom2Name.startswith("C"):
                                continue
                        if atom4Name.startswith("H"):
                            if atom3Name.startswith("C"):
                                continue
                        ## add torsion data to list
                        rotatableBonds.append({
                            'atoms': (atom1Name, atom2Name, atom3Name, atom4Name),
                            'indices': (atom1.GetIdx(), atom2.GetIdx(), atom3.GetIdx(), atom4.GetIdx())
                        })

    return rotatableBonds
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
def clean_up_opt_dir(optDir, cleanUpLevel):
    """
    For an opt dir:
    vital = orca_opt.xyz (for use in forwards)
    nice2have = orca_opt.out (for debugging)
    
    OPTIONS: off, basic, full
 
    """
    if cleanUpLevel == "off":
        return
    
    keepFiles = ["orca_opt.xyz", "ORCA_FINISHED_NORMALLY", "ORCA_CRASHED"]

    for file in os.listdir(optDir):
        if file in keepFiles:
            continue
        elif cleanUpLevel in ["basic", "full"]:
            os.remove(p.join(optDir, file))
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def clean_up_scan_dir(scanDir, cleanUpLevel):
    """
    For a scan dir
    vital = orca_scan_XYZ.xyz, orca_scan.relaxscanact.dat, orca_scan.out
    """
    if cleanUpLevel == "off":
        return

    scanXyzs = [file for file in os.listdir(scanDir) if re.match(r'orca_scan\.\d\d\d\.xyz$', file)]

    keepFiles = ["orca_scan.relaxscanact.dat",  ## need for getting energies
                 "orca_scan.inp",               ## need for getting angles
                 "ORCA_FINISHED_NORMALLY",      ## need as a positive flag for charge fitting later
                 "ORCA_CRASHED"]                ## need as a negative flag for charge fitting later

    for file in os.listdir(scanDir):
        if file in scanXyzs:
            continue
        elif file in keepFiles:
            continue
        elif cleanUpLevel in ["basic", "full"]:
            os.remove(p.join(scanDir, file))
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

def clean_up_singlepoint_dir(spDir, cleanUpLevel):
    """
    For a singlepoint dir
    vital = None
    nice2have = orca_singlepoint.out
    """

    keepFiles = ["orca_opt.xyz", "ORCA_FINISHED_NORMALLY", "ORCA)CRASHED"]

    if cleanUpLevel == "off":
        return
    if cleanUpLevel == "basic":
        for file in os.listdir(spDir):
            if file in keepFiles:
                continue
            else:
                os.remove(p.join(spDir, file))

    elif cleanUpLevel == "full":
        rmtree(spDir)


