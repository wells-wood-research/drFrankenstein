import os
from os import path as p
import pandas as pd
import numpy as np

## CLEAN CODE CLASSES ##
class FilePath:
    pass
class DirectoryPath:
    pass
from typing import List, Tuple
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def reorder_atom_ids(pdbDf: pd.DataFrame) -> pd.DataFrame:
    """
    Reorders the ATOM_ID column in a pdb DataFrame, starting with 1

    Args:
        pdbDf (pd.DataFrame): a PDB file loaded into a DataFrame using pdbUtils
    Returns:
        pdbDf (pd.DataFrame): the same DataFrame with reordered ATOM_ID column

    """
    pdbDf["ATOM_ID"] = np.arange(len(pdbDf["ATOM_ID"])) + 1
    return pdbDf

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def find_bonded_atoms(molDf: pd.DataFrame,
                       atomName: str,
                         distanceCutoff: float = 1.6) -> list:
    """
    Finds all atoms within a certain distance of a given atom

    Args:
        molDf (pd.DataFrame): dataframe of the molecule
        atomName (str): name of the atom to find bonded atoms for
        distanceCutoff (float): distance cutoff
    
    Returns:
        bondedAtoms (list): list of bonded atoms
    
    """
    ## make a copy of input dataframe
    tmpDf = molDf.copy()
    ## make a dataframe containing only our atom of interest
    atomDf = molDf[molDf["ATOM_NAME"] == atomName]
    ## get atom coords
    atomX, atomY, atomZ = atomDf.loc[:,['X', 'Y', 'Z']].astype(float).iloc[0]
    ## calculate distance to all other atoms
    tmpDf['distance_to_atom'] = tmpDf.apply(calculate_distance,
                                             axis=1, atomX=atomX, atomY=atomY, atomZ=atomZ)
    ## find all atoms within distance cutoff
    bondedDf = tmpDf[tmpDf['distance_to_atom'] < distanceCutoff]
    ## convert to list
    bondedAtoms = bondedDf["ATOM_NAME"].to_list()
    ## remove atom of interest
    bondedAtoms = [atom for atom in bondedAtoms if atom != atomName]
    return bondedAtoms
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def decide_atom_to_delete_C_termini(atomNames: list) -> str:
    """
    From a list of atom names, decide which one to delete in relation to the C-Termini

    Args:
        atomNames (list): list of atom names
    Returns:
        atomToDelete (str): name of the atom to delete (can be "None")
    
    """
    atomNames = sorted([atomName for atomName in atomNames if atomName.startswith("O")])
    if len(atomNames) == 2:
        atomToDelete = atomNames[-1]
    else: 
        atomToDelete = "None"

    return atomToDelete
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def decide_atom_to_delete_N_termini(atomNames: list) -> str:
    """
    From a list of atom names, decide which one to delete in relation to the N-Termini

    Args:
        atomNames (list): list of atom names
    Returns:
        atomToDelete (str): name of the atom to delete (can be "None")
    
    """
    atomNames = sorted([atomName for atomName in atomNames if atomName.startswith("H")])
    if len(atomNames) > 1:
        atomToDelete = atomNames[-1]
    else:
        atomToDelete = "None"
    return atomToDelete

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

def  decide_atom_to_delete_C_terminal_proton(atomNames: list) -> str:
    """
    From a list of atom names, decide which one to delete in relation to the C-Terminal proton

    Args:
        atomNames (list): list of atom names
    Returns:
        atomToDelete (str): name of the atom to delete (can be "None")
    
    """   
    atomNames = sorted([atomName for atomName in atomNames if atomName.startswith("H")])
    if len(atomNames) == 1:
        atomToDelete = atomNames[0]
    elif len(atomNames) == 2:
        atomToDelete = atomNames[-1]
    else: 
        atomToDelete = "None"

    return atomToDelete


#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def calculate_distance(row: pd.Series,
                        atomX: float,
                          atomY: float,
                            atomZ: float) -> float:
    """
    Calculates distance between a row of a pdb dataframe and a set of coords
    Used in an apply function
    
    Args:
        row (pd.Series): row of a pdb dataframe
        atomX (float): x coord
        atomY (float): y coord
        atomZ (float): z coord
    Returns:
        distance (float): distance
    """
    return np.sqrt((row['X'] - atomX) ** 2 +
                   (row['Y'] - atomY) ** 2 +
                   (row['Z'] - atomZ) ** 2)

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def find_capping_pdbs() -> Tuple[FilePath, FilePath]:
    """
    Looks through drFrankenstein's src directory to find the NME and ACE PDBs

    Args:
        None
    Returns:
        nmePdb (FilePath): path to NME.pdb
        acePdb (FilePath): path to ACE.pdb
    """
    ## get location of this file
    drCapperDir = p.dirname(p.abspath(__file__))
    nmePdb = p.join(drCapperDir, "NME.pdb")
    acePdb = p.join(drCapperDir, "ACE.pdb")
    if not p.isfile(nmePdb):
        raise FileNotFoundError(f"Could not find NME.pdb in {drCapperDir}")
    if not p.isfile(acePdb):
        raise FileNotFoundError(f"Could not find ACE.pdb in {drCapperDir}")
    return nmePdb, acePdb

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def apply_transformation(pdbDf: pd.DataFrame,
                          rotationMatrix: np.ndarray,
                            translationVector: np.ndarray) -> pd.DataFrame:
    """
    Applies rotation and translation to a PDB DataFrame

    Args:
        pdbDf (pd.DataFrame): PDB DataFrame
        rotationMatrix (np.ndarray): rotation matrix
        translationVector (np.ndarray): translation vector

    Returns:
        pdbDf (pd.DataFrame): PDB DataFrame
    """
    coords = pdbDf[['X', 'Y', 'Z']].values
    transformedCoords = np.dot(coords, rotationMatrix) + translationVector
    pdbDf[['X', 'Y', 'Z']] = transformedCoords
    return pdbDf

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

def calculate_transformation(originalDf: pd.DataFrame,
                              targetDf: pd.DataFrame,
                                alignmentAtoms: list) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculates the transformation and rotations required to transform the original molecule to the target molecule,
    Used for moving the methyl protons of the N-Me cap

    Args:
        originalDf (pd.DataFrame): DataFrame of the original molecule
        targetDf (pd.DataFrame): DataFrame of the target molecule
        alignmentAtoms (list): list of atom names to be used to transform the rest of the molecule

    Returns:
        rotationMatrix (np.ndarray): rotation matrix
        translationVector (np.ndarray): translation vector
    
    """

    # Extract coordinates for specified atoms
    originalCoords = originalDf[originalDf['ATOM_NAME'].isin(alignmentAtoms)][['X', 'Y', 'Z']].values
    targetCoords = targetDf[targetDf['ATOM_NAME'].isin(alignmentAtoms)][['X', 'Y', 'Z']].values

    # Calculate centroids
    originalCentroid = np.mean(originalCoords, axis=0)
    targetCentroid = np.mean(targetCoords, axis=0)

    # Center the coordinates
    originalCoordsCentered = originalCoords - originalCentroid
    targetCoordsCentered = targetCoords - targetCentroid

    # Compute covariance matrix
    covarianceMatrix = np.dot(originalCoordsCentered.T, targetCoordsCentered)

    # Singular Value Decomposition
    V, S, Wt = np.linalg.svd(covarianceMatrix)
    rotationMatrix = np.dot(V, Wt)

    # Ensure a proper rotation (determinant should be 1)
    if np.linalg.det(rotationMatrix) < 0:
        V[:, -1] = -V[:, -1]
        rotationMatrix = np.dot(V, Wt)

    # Calculate translation
    translationVector = targetCentroid - np.dot(originalCentroid, rotationMatrix)

    return rotationMatrix, translationVector
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

def transform_whole(originalDf: pd.DataFrame,
                     targetDf: pd.DataFrame,
                       atomNames: list) -> pd.DataFrame:
    """
    Transforms the original molecule to the target molecule, given at least three atom names
    Used to move the methyl protons of the N-Me and Ace capping groups

    Args:
        originalDf (pd.DataFrame): DataFrame of the original molecule
        targetDf (pd.DataFrame): DataFrame of the target molecule
        atomNames (list): list of atom names to be used to transform the rest of the molecule (must be length == 3)
    Returns:
        transformedDf (pd.DataFrame): DataFrame of the transformed molecule
    """
    if not len(atomNames) == 3:
        raise ValueError("atomNames must be length == 3")

    # Calculate transformation
    rotationMatrix, translationVector = calculate_transformation(originalDf, targetDf, atomNames)


    transformedDf = apply_transformation(originalDf.copy(), rotationMatrix, translationVector)

    return transformedDf

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def place_capping_atom_by_vector(capDf: pd.DataFrame,
                                  vectorAtomCoordsA: np.ndarray,
                                    vectorAtomCoordsB: np.ndarray,
                                      bondedAtomCoords: np.ndarray,
                                        atomNameToPlace: str,
                                          bondLength: float) -> pd.DataFrame:
    """
    Places a capping atom in its DataFrame by:
    1. Calculate vector from atom A to B
    2. Use vector and bond length to place capping atom, relative to bonded atom

    Args:
        capDf (pd.DataFrame): DataFrame of capping atoms
        vectorAtomCoordsA (np.ndarray): coordinates of atom A
        vectorAtomCoordsB (np.ndarray): coordinates of atom B
        bondedAtomCoords (np.ndarray): coordinates of bonded atom
        atomNameToPlace (str): name of atom to place
        bondLength (float): length of bond
    Returns:
        capDf (pd.DataFrame): DataFrame with capping atom placed
    """
    ## calculate a vector from atom A to B
    vectorAtoB =  vectorAtomCoordsB - vectorAtomCoordsA
    ## normalise the vector
    vectorAtoB /= np.linalg.norm(vectorAtoB)
    ## claculate new coordinates for capping atom
    newCoords = bondedAtomCoords + vectorAtoB * bondLength
    ## place capping atom in capDf
    capDf.loc[capDf["ATOM_NAME"] == atomNameToPlace, ["X", "Y", "Z"]] = newCoords

    return capDf
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def place_capping_atom_by_average_of_vectors(capDf: pd.DataFrame,
                                              vectorAtomCoordsA: np.ndarray,
                                                vectorAtomCoordsB: np.ndarray,
                                                  vectorAtomCoordsC: np.ndarray,
                                                    vectorAtomCoordsD: np.ndarray,
                                                      bondedAtomCoords: np.ndarray,
                                                        atomNameToPlace: str,
                                                          bondLength: float) -> pd.DataFrame:
    """
    Places a capping atom in its DataFrame by:
    1. Calculate vector from atom A to B
    2. Calculate vector from atom C to D
    3. Calculate average of vectors A->B and C->D
    4. Use average vector and bond length to place capping atom, relative to bonded atom
    
    Args:
        capDf (pd.DataFrame): DataFrame of capping atoms
        vectorAtomCoordsA (np.ndarray): coordinates of atom A
        vectorAtomCoordsB (np.ndarray): coordinates of atom B
        vectorAtomCoordsC (np.ndarray): coordinates of atom C
        vectorAtomCoordsD (np.ndarray): coordinates of atom D
        bondedAtomCoords (np.ndarray): coordinates of bonded atom
        atomNameToPlace (str): name of atom to place
        bondLength (float): length of bond
    Returns:
        capDf (pd.DataFrame): DataFrame with capping atom placed
    """
    ## calculate a vector from atom A to B
    vectorAtoB =  vectorAtomCoordsB - vectorAtomCoordsA
    ## calculate a vector from atom C to D
    vectorCtoD =   vectorAtomCoordsD - vectorAtomCoordsC
    ## calculate average direction
    averageDirection = (vectorAtoB + vectorCtoD) / 2
    ## normalise the vector
    averageDirection /= np.linalg.norm(averageDirection)  # Normalize
    ## create new vector with dimensions from bondLength
    newVector = - bondLength * averageDirection
    ## claculate new coordinates for capping atom
    newCoords = bondedAtomCoords - newVector
    ## update capDf with new position for capping atom
    capDf.loc[capDf["ATOM_NAME"] == atomNameToPlace, ["X", "Y", "Z"]] = newCoords
    return capDf