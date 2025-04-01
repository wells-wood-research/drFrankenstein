## BASIC IMPORTS ##
import os
from os import path as p
import numpy as np
import pandas as pd

## CUSTOM LIBRARIES ##
from pdbUtils import pdbUtils

## CLEAN CODE CLASSES ##
class FilePath:
    pass
class DirectoryPath:
    pass
from typing import List, Tuple
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def capping_protocol(config: dict) -> dict:
    """
    Main protocol for adding capping groups to N and C termini of our non-canonical amino acid
    This is important as it allows for better calculation of charges and torsion angles 
    and gives us a nice closed-shell system with no unpaired electrons

    Args:
        config (dict) : the drFrankenstein config containing all run information
    Returns:
        config (dict): updated config
    
    """
    ## unpack config ##
    outputDir = config["pathInfo"]["outputDir"]
    molPdb = config["moleculeInfo"]["moleculePdb"]
    moleculeName = config["moleculeInfo"]["moleculeName"]

    ## make a dir for capping, update config ##
    cappingDir = p.join(outputDir, "01_termini_capping")
    os.makedirs(cappingDir, exist_ok=True)
    config["pathInfo"]["cappingDir"] = outputDir

    ## load pdb file into a dataframe
    molDf = pdbUtils.pdb2df(molPdb)
    ## remove any extra atoms bound to the termini to make room for capping groups
    molDf = trim_termini(molDf, config)
    ## get capping PDBs from the drCapper directory
    nmePdb, acePdb = find_capping_pdbs()
    ## add N-Methyl (NCH3) caps to each C-Terminus
    nmeCappedDf = add_nMethyl_caps(molDf, nmePdb, config)
    ## add Acetyl (COCH3) caps to each N-Terminus
    aceCappedDf = add_acetyl_caps(nmeCappedDf, acePdb, config)
    ## edit the ATOM_ID column to make it 1 -> N
    cappedDf = reorder_atom_ids(aceCappedDf)
    ## write to PDB file, update config
    cappedPdb = p.join(cappingDir, f"{moleculeName}_capped.pdb")
    pdbUtils.df2pdb(cappedDf, cappedPdb)
    config["moleculeInfo"]["cappedPdb"] = cappedPdb
    ## update config, this lets drMD know that capping is complete and lets it skip this step next time
    config["checkpointInfo"]["cappingComplete"] = True


    return config
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
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
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def add_nMethyl_caps(molDf: pd.DataFrame,
                      nmePdb: FilePath,
                        config: dict) -> pd.DataFrame:
    """
    Protocol for adding N-Methyl (NCH3) caps to each C-Terminus

    Args:
        molDf (pd.DataFrame): dataframe of the molecule
        nmePdb (FilePath): path to the N-Methyl (NCH3) PDB file
        config (dict): drFrankenstein config

    Returns:
        nmeCappedDf (pd.DataFrame): dataframe of NMethyl capped molecule with updated coordinates
    
    """
    cappedDf = molDf.copy()
    cTerminalAtoms = config["moleculeInfo"]["cTermini"]
    nmeDf = pdbUtils.pdb2df(nmePdb)
    orginalNmeDf = nmeDf.copy()

    dfsToConcat = [molDf]

    for i, cTerminalAtom in enumerate(cTerminalAtoms):
        tmpNmeDf = place_NN(cappedDf, cTerminalAtom, nmeDf)
        tmpNmeDf = place_HNN1(cappedDf, cTerminalAtom, tmpNmeDf)
        tmpNmeDf = place_CN(cappedDf, cTerminalAtom, tmpNmeDf)
        tmpNmeDf = transform_whole(orginalNmeDf, tmpNmeDf, atomNames=["NN", "HNN1", "CN"])
        dfsToConcat.append(tmpNmeDf)

    nmeCappedDf = pd.concat(dfsToConcat, axis=0)
    return nmeCappedDf
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def add_acetyl_caps(molDf: pd.DataFrame,
                     acePdb: FilePath,
                       config: dict) -> pd.DataFrame:
    """
    Protocol for adding Acetyl (CCH3) caps to each N-Terminus

    Args:
        molDf (pd.DataFrame): dataframe of the molecule
        acePdb (FilePath): path to the Acetyl (CCH3) PDB file
        config (dict): drFrankenstein config

    Returns:
        aceCappedDf (pd.DataFrame): dataframe of Acetyl capped molecule with updated coordinates
    
    """
    cappedDf = molDf.copy()
    nTerminalAtoms = config["moleculeInfo"]["nTermini"]
    aceDf = pdbUtils.pdb2df(acePdb)

    dfsToConcat = [molDf]

    for i, nTerminalAtom in enumerate(nTerminalAtoms):
        target_aceDf = place_CC1(cappedDf, nTerminalAtom, aceDf)
        target_aceDf = place_OC(cappedDf, nTerminalAtom, aceDf)
        target_aceDf = place_CC2(cappedDf, nTerminalAtom, aceDf)
        aceDf = pdbUtils.pdb2df(acePdb)
        transformed_aceDf = transform_whole(aceDf, target_aceDf, atomNames=["CC1", "OC", "CC2"])
        dfsToConcat.append(transformed_aceDf)

    aceCappedDf = pd.concat(dfsToConcat, axis=0)
    return aceCappedDf
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
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
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
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

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def place_CC2(cappedDf: pd.DataFrame,
               nTerminalAtom: str,
                 aceDf: pd.DataFrame) -> pd.DataFrame:
    """
    Places the CC2 capping atom in the Acetyl capping group
    
    Args:
        cappedDf (pd.DataFrame): DataFrame of capped molecule
        nTerminalAtom (str): name of N-Terminus atom
        aceDf (pd.DataFrame): DataFrame of Acetyl capping group
    Returns:
        aceDf (pd.DataFrame): DataFrame with CC2 capping atom placed
    """
    ## get coordinates
    nTerminalCoords = cappedDf[cappedDf["ATOM_NAME"] == nTerminalAtom].loc[:, ["X", "Y", "Z"]].values[0]
    cc1Coords = aceDf[aceDf["ATOM_NAME"] == "CC1"].loc[:, ["X", "Y", "Z"]].values[0]
    ocCoords = aceDf[aceDf["ATOM_NAME"] == "OC"].loc[:, ["X", "Y", "Z"]].values[0]
    ## place the CC2 capping atom
    aceDf = place_capping_atom_by_average_of_vectors(capDf = aceDf,
                                                      vectorAtomCoordsA = nTerminalCoords,
                                                       vectorAtomCoordsB = cc1Coords,
                                                        vectorAtomCoordsC = ocCoords,
                                                         vectorAtomCoordsD = cc1Coords, 
                                                          bondedAtomCoords = cc1Coords,
                                                           atomNameToPlace = "CC2",
                                                            bondLength= 1.45)
    return aceDf
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def place_CC1(cappedDf: pd.DataFrame,
               nTerminalAtom: str,
                 aceDf: pd.DataFrame) -> pd.DataFrame:
    """
    Places the CC1 capping atom in the Acetyl capping group
    
    Args:
        cappedDf (pd.DataFrame): DataFrame of capped molecule
        nTerminalAtom (str): name of N-Terminus atom
        aceDf (pd.DataFrame): DataFrame of Acetyl capping group
    Returns:
        aceDf (pd.DataFrame): DataFrame with CC1 capping atom placed
    """
    
    ## get coordinates for N-Terminus
    nTerminalCoords = cappedDf[cappedDf["ATOM_NAME"] == nTerminalAtom].loc[:, ["X", "Y", "Z"]].values[0]
    ##  find atoms bound to N-Terminus
    bondedAtoms = find_bonded_atoms(cappedDf, nTerminalAtom)
    ## get coordinates for proton bound to N-Terminus
    hAtomNames = [atom for atom in bondedAtoms if atom.startswith("H")]
    if len(hAtomNames)  == 0:
        hAtomName = [atom for atom in bondedAtoms if atom.startswith("C") and not atom == "CA"][0]
    else:
        hAtomName = hAtomNames[0]
    hCoords = cappedDf[cappedDf["ATOM_NAME"] == hAtomName].loc[:, ["X", "Y", "Z"]].values[0]
    ## get coordinates for alpha-carbon or equivalent
    caAtomNames = [atom for atom in bondedAtoms if atom.startswith("C")]
    if len(caAtomNames) == 0:
        caName = [atom for atom in bondedAtoms if not atom.startswith("H")][0]
    elif "CA" in caAtomNames:
        caName = "CA"
    else:
        caName = caAtomNames[0]
    caCoords = cappedDf[cappedDf["ATOM_NAME"] == caName].loc[:, ["X", "Y", "Z"]].values[0]
    ## place the CC1 capping atom
    aceDf = place_capping_atom_by_average_of_vectors(capDf = aceDf,
                                                      vectorAtomCoordsA = hCoords,
                                                       vectorAtomCoordsB = nTerminalCoords,
                                                        vectorAtomCoordsC = caCoords,
                                                          vectorAtomCoordsD =  nTerminalCoords,
                                                             bondedAtomCoords = nTerminalCoords,
                                                              atomNameToPlace = "CC1",
                                                               bondLength = 1.4)

    return aceDf
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def place_OC(cappedDf: pd.DataFrame,
              nTerminalAtom: str,
                aceDf: pd.DataFrame) -> pd.DataFrame:
    """
    Places the OC capping atom in the Acetyl capping group
    
    Args:
        cappedDf (pd.DataFrame): DataFrame of capped molecule
        nTerminalAtom (str): name of N-Terminus atom
        aceDf (pd.DataFrame): DataFrame of Acetyl capping group
    Returns:
        aceDf (pd.DataFrame): DataFrame with OC capping atom placed
    """

    ## get coordinates for N-Terminus
    nTerminalCoords = cappedDf[cappedDf["ATOM_NAME"] == nTerminalAtom].loc[:, ["X", "Y", "Z"]].values[0]
    ## find atoms bound to N-Terminus
    bondedAtoms = find_bonded_atoms(cappedDf, nTerminalAtom)
    ## find coordinates for proton bound to N-Terminus
    hAtomNames = [atom for atom in bondedAtoms if atom.startswith("H")]
    if len(hAtomNames)  == 0:
        hAtomName = [atom for atom in bondedAtoms if atom.startswith("C") and not atom == "CA"][0]
    else:
        hAtomName = hAtomNames[0]
    hCoords = cappedDf[cappedDf["ATOM_NAME"] == hAtomName].loc[:, ["X", "Y", "Z"]].values[0]
    ## get coordinates for CC1 in the Acetyl capping group
    cc1Coords = aceDf[aceDf["ATOM_NAME"] == "CC1"].loc[:, ["X", "Y", "Z"]].values[0]
    ## place the OC capping atom
    aceDf = place_capping_atom_by_vector(capDf=aceDf,
                                        vectorAtomCoordsA = hCoords,
                                        vectorAtomCoordsB = nTerminalCoords,
                                        bondedAtomCoords = cc1Coords,
                                        atomNameToPlace = "OC",
                                        bondLength = 1.2)

    return aceDf
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

def transform_whole(originalDf: pd.DataFrame,
                     targetDf: pd.DataFrame,
                       atomNames: list) -> pd.DataFrame:
    """
    Transforms the original molecule to the target molecule

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
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

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
    covariance_matrix = np.dot(originalCoordsCentered.T, targetCoordsCentered)

    # Singular Value Decomposition
    V, S, Wt = np.linalg.svd(covariance_matrix)
    rotationMatrix = np.dot(V, Wt)

    # Ensure a proper rotation (determinant should be 1)
    if np.linalg.det(rotationMatrix) < 0:
        V[:, -1] = -V[:, -1]
        rotationMatrix = np.dot(V, Wt)

    # Calculate translation
    translationVector = targetCentroid - np.dot(originalCentroid, rotationMatrix)

    return rotationMatrix, translationVector
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

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
    transformed_coords = np.dot(coords, rotationMatrix) + translationVector
    pdbDf[['X', 'Y', 'Z']] = transformed_coords
    return pdbDf

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def place_CN(cappedDf: pd.DataFrame,
              cTerminalAtom: str,
                nmeDf: pd.DataFrame) -> pd.DataFrame:
    """
    Places the CN capping atom in the N-Methyl capping group
    
    Args:
        cappedDf (pd.DataFrame): DataFrame of capped molecule
        cTerminalAtom (str): name of C-Terminus atom
        nmeDf (pd.DataFrame): DataFrame of N-Methyl capping group
    Returns:
        nmeDf (pd.DataFrame): DataFrame with CN capping atom placed
    """
    ## get coordinates for C-Terminus, NN and HNN atoms
    cTerminalCoords = cappedDf[cappedDf["ATOM_NAME"] == cTerminalAtom].loc[:, ["X", "Y", "Z"]].values[0]
    nnCoords = nmeDf[nmeDf["ATOM_NAME"] == "NN"].loc[:, ["X", "Y", "Z"]].values[0]
    hnnCoords = nmeDf[nmeDf["ATOM_NAME"] == "HNN1"].loc[:, ["X", "Y", "Z"]].values[0]
    ## place CN
    nmeDf = place_capping_atom_by_average_of_vectors(capDf = nmeDf,
                                             vectorAtomCoordsA = cTerminalCoords,
                                               vectorAtomCoordsB = nnCoords,
                                                 vectorAtomCoordsC = hnnCoords,
                                                   vectorAtomCoordsD = nnCoords,
                                                     bondedAtomCoords = nnCoords,
                                                       atomNameToPlace = "CN",
                                                         bondLength = 1.4)
    
    return nmeDf
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def place_HNN1(cappedDf: pd.DataFrame,
                cTerminalAtom: str,
                  nmeDf: pd.DataFrame) -> pd.DataFrame:
    """
    Places the HNN1 capping atom in the N-Methyl capping group
    
    Args:
        cappedDf (pd.DataFrame): DataFrame of capped molecule
        cTerminalAtom (str): name of C-Terminus atom
        nmeDf (pd.DataFrame): DataFrame of N-Methyl capping group
    Returns:
        nmeDf (pd.DataFrame): DataFrame with HNN1 capping atom placed
    """
    # Get coordinates from protein
    bondedAtoms = find_bonded_atoms(cappedDf, cTerminalAtom)

    oEquivalentName = [atomName for atomName in bondedAtoms if atomName.startswith("O")][0]
    
    cTerminalCoords = cappedDf[cappedDf["ATOM_NAME"] == cTerminalAtom].loc[:, ["X", "Y", "Z"]].values[0]
    oCoords = cappedDf[cappedDf["ATOM_NAME"] == oEquivalentName].loc[:, ["X", "Y", "Z"]].values[0]
    
    # Get original NN coordinates
    nnCoords = nmeDf[nmeDf["ATOM_NAME"] == "NN"].loc[:, ["X", "Y", "Z"]].values[0]

    nmeDf = place_capping_atom_by_vector(capDf=nmeDf,
                                        vectorAtomCoordsA = oCoords,
                                        vectorAtomCoordsB = cTerminalCoords,
                                        bondedAtomCoords = nnCoords,
                                        atomNameToPlace = "HNN1",
                                        bondLength = 1.0)

    return nmeDf
###################################################################
def place_NN(cappedDf: pd.DataFrame,
              cTerminalAtom: str,
                nmeDf: pd.DataFrame) -> pd.DataFrame:
    """
    Places the NN capping atom in the N-Methyl capping group
    
    Args:
        cappedDf (pd.DataFrame): DataFrame of capped molecule
        cTerminalAtom (str): name of C-Terminus atom
        nmeDf (pd.DataFrame): DataFrame of N-Methyl capping group
    Returns:
        nmeDf (pd.DataFrame): DataFrame with NN capping atom placed
    """    
    # find atoms bonded to C-Terminus
    bondedAtoms = find_bonded_atoms(cappedDf, cTerminalAtom)
    ## get coordinates for CA and O atoms, or best equivalent
    caEquivalentName = [atomName for atomName in bondedAtoms if atomName.startswith("C")][0]
    oEquivalentName = [atomName for atomName in bondedAtoms if atomName.startswith("O")][0]

    cTerminalCoords = cappedDf[cappedDf["ATOM_NAME"] == cTerminalAtom].loc[:, ["X", "Y", "Z"]].values[0]
    caCoords = cappedDf[cappedDf["ATOM_NAME"] == caEquivalentName].loc[:, ["X", "Y", "Z"]].values[0]
    oCoords = cappedDf[cappedDf["ATOM_NAME"] == oEquivalentName].loc[:, ["X", "Y", "Z"]].values[0]
    ## place NN
    nmeDf = place_capping_atom_by_average_of_vectors(capDf=nmeDf,
                                                     vectorAtomCoordsA=oCoords,
                                                      vectorAtomCoordsB=cTerminalCoords,
                                                       vectorAtomCoordsC=caCoords,
                                                        vectorAtomCoordsD=cTerminalCoords,
                                                         bondedAtomCoords=cTerminalCoords,
                                                          atomNameToPlace="NN",
                                                           bondLength=1.4)

    
    return nmeDf
#####################################################################

def trim_termini(molDf: pd.DataFrame,
                  config: dict) -> pd.DataFrame:
    """
    Removes any extra atoms bound to the termini to make room for capping groups
    1. Remove extra protons from N-Termini
    2. Remove extra oxygens and protons from C-Termini

    Args:
        molDf (pd.DataFrame): dataframe of the molecule
        config (dict): drFrankenstein config
    
    """
    ## unpack config ##
    nTerminalAtoms = config["moleculeInfo"]["nTermini"]
    cTerminalAtoms = config["moleculeInfo"]["cTermini"]
    ## get the proton (or equivalent) atoms bound to N-Terminus
    allAtomsToDelete = []
    for nTerminalAtom in nTerminalAtoms:
        bondedAtoms = find_bonded_atoms(molDf, nTerminalAtom)
        nTerminalProton = decide_atom_to_delete_N_termini(bondedAtoms)
        allAtomsToDelete.append(nTerminalProton)

    ## get the oxygen (or equivalent) atoms bound to C-Terminus, and anything bound to that oxygen
    for cTerminalAtom in cTerminalAtoms:
        bondedAtoms = find_bonded_atoms(molDf, cTerminalAtom)
        cTerminalOxygen = decide_atom_to_delete_C_termini(bondedAtoms)
        allAtomsToDelete.append(cTerminalOxygen)
        if not cTerminalOxygen == "None":
            bondedAtoms = find_bonded_atoms(molDf, cTerminalOxygen)
            cTerminalProton = decide_atom_to_delete_C_terminal_proton(bondedAtoms)
            allAtomsToDelete.append(cTerminalProton)

    molDf = molDf[~molDf["ATOM_NAME"].isin(allAtomsToDelete)]

    return molDf

#####################################################################
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
#####################################################################
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
#####################################################################
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

#####################################################################
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

#####################################################################
def calculate_distance(row: pd.Series,
                        atomX: float,
                          atomY: float,
                            atomZ: float):
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

#####################################################################


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


#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
if __name__ == "__main__":
    raise NotImplementedError