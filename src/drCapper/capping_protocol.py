## BASIC IMPORTS ##
import os
from os import path as p
import yaml
import numpy as np
import pandas as pd

## CUSTOM LIBRARIES ##
from pdbUtils import pdbUtils
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def capping_protocol(config: dict):

    outputDir = config["pathInfo"]["outputDir"]
    molPdb = config["moleculeInfo"]["moleculePdb"]
    moleculeName = config["moleculeInfo"]["moleculeName"]

    cappingDir = p.join(outputDir, "01_termini_capping")
    os.makedirs(cappingDir, exist_ok=True)
    config["pathInfo"]["cappingDir"] = outputDir

    
    molDf = pdbUtils.pdb2df(molPdb)

    molDf = trim_termini(molDf, config)

    nmePdb, acePdb = find_capping_pdbs()

    nmeCappedDf = add_nMethyl_caps(molDf, nmePdb, config)
    aceCappedDf = add_acetyl_caps(nmeCappedDf, acePdb, config)

    aceCappedDf["ATOM_ID"] = np.arange(len(aceCappedDf["ATOM_ID"])) + 1

    cappedPdb = p.join(cappingDir, f"{moleculeName}_capped.pdb")

    pdbUtils.df2pdb(aceCappedDf, cappedPdb)

    config["moleculeInfo"]["cappedPdb"] = cappedPdb

    config["checkpointInfo"]["cappingComplete"] = True

    return config

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def add_acetyl_caps(molDf, acePdb, config):
    tmpDf = molDf.copy()
    nTerminalAtoms = config["moleculeInfo"]["nTermini"]
    aceDf = pdbUtils.pdb2df(acePdb)

    dfsToConcat = [molDf]

    maxResId = molDf["RES_ID"].max()

    for i, nTerminalAtom in enumerate(nTerminalAtoms):
        target_aceDf = place_CC1(tmpDf, nTerminalAtom, aceDf)
        pdbUtils.df2pdb(target_aceDf,"cc1.pdb")
        target_aceDf = place_OC(tmpDf, nTerminalAtom, aceDf)
        pdbUtils.df2pdb(target_aceDf,"oc.pdb")
        target_aceDf = place_CC2(tmpDf, nTerminalAtom, aceDf)
        pdbUtils.df2pdb(target_aceDf, "cc2.pdb")
        pdbUtils.df2pdb(aceDf, "original.pd")
        aceDf = pdbUtils.pdb2df(acePdb)
        transformed_aceDf = transform_whole(aceDf, target_aceDf, atomNames=["CC1", "OC", "CC2"])
        ## set resID (cope with 0 indexing)
        transformed_aceDf["RES_ID"] = maxResId + i + 1

        pdbUtils.df2pdb(transformed_aceDf, "final_ace.pdb")
        dfsToConcat.append(transformed_aceDf)

    aceCappedDf = pd.concat(dfsToConcat, axis=0)
    return aceCappedDf
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

def place_OC(tmpDf, nTerminalAtom, aceDf):
    ## get coords
    nTerminalCoords = tmpDf[tmpDf["ATOM_NAME"] == nTerminalAtom].loc[:, ["X", "Y", "Z"]].values[0]

    bondedAtoms = find_bonded_atoms(tmpDf, nTerminalAtom)
    hAtomNames = [atom for atom in bondedAtoms if atom.startswith("H")]
    if len(hAtomNames)  == 0:
        hAtomName = [atom for atom in bondedAtoms if atom.startswith("C") and not atom == "CA"][0]
    else:
        hAtomName = hAtomNames[0]

    hCoords = tmpDf[tmpDf["ATOM_NAME"] == hAtomName].loc[:, ["X", "Y", "Z"]].values[0]
    cc1Coords = aceDf[aceDf["ATOM_NAME"] == "CC1"].loc[:, ["X", "Y", "Z"]].values[0]

    vectorHtoN = get_vector(hCoords, nTerminalCoords)

    bondLength = 1.2

    newCoordsOC = cc1Coords + vectorHtoN * bondLength
    
    ocIndex = aceDf[aceDf["ATOM_NAME"]=="OC"].index

    aceDf.loc[ocIndex,["X", "Y", "Z"]] = newCoordsOC

    return aceDf

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def place_CC2(tmpDf, nTerminalAtom, aceDf):
    nTerminalCoords = tmpDf[tmpDf["ATOM_NAME"] == nTerminalAtom].loc[:, ["X", "Y", "Z"]].values[0]
    cc1Coords = aceDf[aceDf["ATOM_NAME"] == "CC1"].loc[:, ["X", "Y", "Z"]].values[0]
    ocCoords = aceDf[aceDf["ATOM_NAME"] == "OC"].loc[:, ["X", "Y", "Z"]].values[0]

    # Define vectors
    vectorNtoCC1 = nTerminalCoords - cc1Coords
    vectorOCtoCC1 = ocCoords - cc1Coords
    
    # Calculate average direction
    averageDirection = (vectorNtoCC1 + vectorOCtoCC1) / 2
    averageDirection /= np.linalg.norm(averageDirection)  # Normalize
    
    bondLength = 1.45  # C-C distance
    
    # Calculate new NN position
    vectorCC1toCC2 = - bondLength * averageDirection
    newCoordsCC2 = cc1Coords + vectorCC1toCC2

    cc2Index = aceDf[aceDf["ATOM_NAME"]=="CC2"].index
    aceDf.loc[cc2Index,["X", "Y", "Z"]] = newCoordsCC2

    return aceDf
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

def place_CC1(tmpDf, nTerminalAtom, aceDf):
    ## get coords
    nTerminalCoords = tmpDf[tmpDf["ATOM_NAME"] == nTerminalAtom].loc[:, ["X", "Y", "Z"]].values[0]

    bondedAtoms = find_bonded_atoms(tmpDf, nTerminalAtom)
    hAtomNames = [atom for atom in bondedAtoms if atom.startswith("H")]
    if len(hAtomNames)  == 0:
        hAtomName = [atom for atom in bondedAtoms if atom.startswith("C") and not atom == "CA"][0]
    else:
        hAtomName = hAtomNames[0]

    
    hCoords = tmpDf[tmpDf["ATOM_NAME"] == hAtomName].loc[:, ["X", "Y", "Z"]].values[0]
    caName = [atom for atom in bondedAtoms if atom.startswith("CA")][0]
    caCoords = tmpDf[tmpDf["ATOM_NAME"] == caName].loc[:, ["X", "Y", "Z"]].values[0]

    # Define vectors
    vectorCAtoN = caCoords - nTerminalCoords
    vectorHtoN = hCoords - nTerminalCoords
    
    # Calculate average direction
    averageDirection = (vectorCAtoN + vectorHtoN) / 2
    averageDirection /= np.linalg.norm(averageDirection)  # Normalize
    
    bondLength = 1.4  # CN-NN distance
    
    # Calculate new NN position
    vectorNNtoCN = - bondLength * averageDirection
    newCoordsCN = nTerminalCoords + vectorNNtoCN

    cnIndex = aceDf[aceDf["ATOM_NAME"]=="CC1"].index
    aceDf.loc[cnIndex,["X", "Y", "Z"]] = newCoordsCN

    return aceDf


#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def add_nMethyl_caps(molDf, nmePdb, config):
    tmpDf = molDf.copy()
    cTerminalAtoms = config["moleculeInfo"]["cTermini"]
    nmeDf = pdbUtils.pdb2df(nmePdb)

    maxResId = molDf["RES_ID"].max()

    dfsToConcat = [molDf]

    for i, cTerminalAtom in enumerate(cTerminalAtoms):
        target_nmeDf = place_NN(tmpDf, cTerminalAtom, nmeDf)
        target_nmeDf = place_HNN1(tmpDf, cTerminalAtom, target_nmeDf)
        target_nmeDf= place_CN(tmpDf, cTerminalAtom, target_nmeDf)
        transformed_nmeDf = transform_whole(nmeDf, target_nmeDf, atomNames=["NN", "HNN1", "CN"])
        ## fix resID (cope with 0 indexing)
        transformed_nmeDf["RES_ID"] = maxResId + i + 1
    dfsToConcat.append(transformed_nmeDf)

    nmeCappedDf = pd.concat(dfsToConcat, axis=0)
    return nmeCappedDf
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

def transform_whole(originalDf, targetDf, atomNames):

    # Calculate transformation
    rotation_matrix, translation_vector = calculate_transformation(originalDf, targetDf, atomNames)

    transformedDf = apply_transformation(originalDf.copy(), rotation_matrix, translation_vector)

    return transformedDf
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

def calculate_transformation(df1, df2, atom_names):
    # Extract coordinates for specified atoms
    coords1 = df1[df1['ATOM_NAME'].isin(atom_names)][['X', 'Y', 'Z']].values

    coords2 = df2[df2['ATOM_NAME'].isin(atom_names)][['X', 'Y', 'Z']].values

    # Calculate centroids
    centroid1 = np.mean(coords1, axis=0)
    centroid2 = np.mean(coords2, axis=0)

    # Center the coordinates
    centered_coords1 = coords1 - centroid1
    centered_coords2 = coords2 - centroid2

    # Compute covariance matrix
    covariance_matrix = np.dot(centered_coords1.T, centered_coords2)

    # Singular Value Decomposition
    V, S, Wt = np.linalg.svd(covariance_matrix)
    rotation_matrix = np.dot(V, Wt)

    # Ensure a proper rotation (determinant should be 1)
    if np.linalg.det(rotation_matrix) < 0:
        V[:, -1] = -V[:, -1]
        rotation_matrix = np.dot(V, Wt)

    # Calculate translation
    translation_vector = centroid2 - np.dot(centroid1, rotation_matrix)

    return rotation_matrix, translation_vector
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

def apply_transformation(df, rotation_matrix, translation_vector):
    coords = df[['X', 'Y', 'Z']].values
    transformed_coords = np.dot(coords, rotation_matrix) + translation_vector
    df[['X', 'Y', 'Z']] = transformed_coords
    return df

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

def place_CN(tmpDf, cTerminalAtom, nmeDf):
    ## get coords
    cTerminalCoords = tmpDf[tmpDf["ATOM_NAME"] == cTerminalAtom].loc[:, ["X", "Y", "Z"]].values[0]
    nnCoords = nmeDf[nmeDf["ATOM_NAME"] == "NN"].loc[:, ["X", "Y", "Z"]].values[0]
    hnnCoords = nmeDf[nmeDf["ATOM_NAME"] == "HNN1"].loc[:, ["X", "Y", "Z"]].values[0]

    # Define vectors
    vectorCtoNN = cTerminalCoords - nnCoords
    vectorHNNtoNN = hnnCoords - nnCoords
    
    # Calculate average direction
    averageDirection = (vectorCtoNN + vectorHNNtoNN) / 2
    averageDirection /= np.linalg.norm(averageDirection)  # Normalize
    
    bondLength = 1.4  # CN-NN distance
    
    # Calculate new NN position
    vectorNNtoCN = - bondLength * averageDirection
    newCoordsCN = nnCoords + vectorNNtoCN

    cnIndex = nmeDf[nmeDf["ATOM_NAME"]=="CN"].index
    nmeDf.loc[cnIndex,["X", "Y", "Z"]] = newCoordsCN

    return nmeDf
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲


def place_HNN1(tmpDf, cTerminalAtom, nmeDf):
    # Get coordinates from protein
    bondedAtoms = find_bonded_atoms(tmpDf, cTerminalAtom)
    if "CA" not in bondedAtoms or "O" not in tmpDf["ATOM_NAME"].values:
        raise ValueError("Missing CA or O bonded to C-Terminus!")
    
    cTerminalCoords = tmpDf[tmpDf["ATOM_NAME"] == cTerminalAtom].loc[:, ["X", "Y", "Z"]].values[0]
    oCoords = tmpDf[tmpDf["ATOM_NAME"] == "O"].loc[:, ["X", "Y", "Z"]].values[0]
    
    # Get original NN coordinates
    nnCoords = nmeDf[nmeDf["ATOM_NAME"] == "NN"].loc[:, ["X", "Y", "Z"]].values[0]

    vectorOtoC = get_vector(oCoords, cTerminalCoords)

    bondLength = 1.0

    vectorFromNN = bondLength * vectorOtoC

    newHNN1Coords = nnCoords + vectorFromNN

    hnn1Index = nmeDf[nmeDf["ATOM_NAME"]=="HNN1"].index
    nmeDf.loc[hnn1Index,["X", "Y", "Z"]] = newHNN1Coords

    return nmeDf
###################################################################
def place_NN(tmpDf, cTerminalAtom, nmeDf):
    # Get coordinates from protein
    bondedAtoms = find_bonded_atoms(tmpDf, cTerminalAtom)
    if not any([item.startswith("CA") for item in bondedAtoms]):
        raise ValueError("Missing CA C-Terminus!")
    
    cTerminalCoords = tmpDf[tmpDf["ATOM_NAME"] == cTerminalAtom].loc[:, ["X", "Y", "Z"]].values[0]
    caCoords = tmpDf[tmpDf["ATOM_NAME"] == "CA"].loc[:, ["X", "Y", "Z"]].values[0]
    oCoords = tmpDf[tmpDf["ATOM_NAME"] == "O"].loc[:, ["X", "Y", "Z"]].values[0]
    
    # Get original NN coordinates
    nnCoords_orig = nmeDf[nmeDf["ATOM_NAME"] == "NN"].loc[:, ["X", "Y", "Z"]].values[0]
    
    # Define vectors
    vectorCtoO = oCoords - cTerminalCoords
    vectorCtoCA = caCoords - cTerminalCoords
    
    # Calculate average direction
    averageDirection = (vectorCtoO + vectorCtoCA) / 2
    averageDirection /= np.linalg.norm(averageDirection)  # Normalize
    
    bondLength = 1.4  # C-NN distance
    
    # Calculate new NN position
    vectorCtoNN = - bondLength * averageDirection
    newCoordsNN = cTerminalCoords + vectorCtoNN
    translation = newCoordsNN - nnCoords_orig
    
    # Apply translation
    nmeTransformedDf = nmeDf.copy()
    nmeTransformedDf.loc[:, ["X", "Y", "Z"]] = nmeDf.apply(
        lambda row: apply_translation(row, translation),
        axis=1
    )
    
    return nmeTransformedDf

###############################################################################
def apply_translation(row, translation):
    position = row[['X', 'Y', 'Z']].values
    new_position = position + translation
    return pd.Series(new_position, index=['X', 'Y', 'Z'])

def apply_rotation(row, rotationMatrix, pivot):
    position = row[['X', 'Y', 'Z']].values
    relative_position = position - pivot
    rotated_position = np.dot(rotationMatrix, relative_position) + pivot
    return pd.Series(rotated_position, index=['X', 'Y', 'Z'])

def calculate_rotation_matrix(axis, angle):
    axis = axis / np.linalg.norm(axis)
    a = np.cos(angle / 2.0)
    b, c, d = -axis * np.sin(angle / 2.0)
    return np.array([[a*a + b*b - c*c - d*d, 2*(b*c - a*d), 2*(b*d + a*c)],
                     [2*(b*c + a*d), a*a + c*c - b*b - d*d, 2*(c*d - a*b)],
                     [2*(b*d - a*c), 2*(c*d + a*b), a*a + d*d - b*b - c*c]])

def get_perpendicular_vector(vector):
    perpVector = np.cross(vector, [1, 0, 0])
    if np.linalg.norm(perpVector) == 0:
        perpVector = np.cross(vector, [0, 1, 0])
    perpVector /= np.linalg.norm(perpVector)
    return perpVector

def get_vector(coordsA, coordsB):
    vectorAtoB = coordsB - coordsA
    vectorAtoB /= np.linalg.norm(vectorAtoB)
    return vectorAtoB
#####################################################################

def trim_termini(molDf, config):
    nTerminalAtoms = config["moleculeInfo"]["nTermini"]
    cTerminalAtoms = config["moleculeInfo"]["cTermini"]

    allAtomsToDelete = []
    for nTerminalAtom in nTerminalAtoms:
        bondedAtoms = find_bonded_atoms(molDf, nTerminalAtom)
        atomToDelete = decide_atom_to_delete_N_termini(bondedAtoms)
        allAtomsToDelete.append(atomToDelete)
        bondedAtoms = find_bonded_atoms(molDf, atomToDelete)
        bondedAtoms.remove(nTerminalAtom)
        allAtomsToDelete.extend(bondedAtoms)

    for cTerminalAtom in cTerminalAtoms:
        bondedAtoms = find_bonded_atoms(molDf, cTerminalAtom)
        atomToDelete = decide_atom_to_delete_C_termini(bondedAtoms)
        bondedAtoms = find_bonded_atoms(molDf, atomToDelete)
        bondedAtoms.remove(cTerminalAtom)
        allAtomsToDelete.extend(bondedAtoms)

    molDf = molDf[~molDf["ATOM_NAME"].isin(allAtomsToDelete)]

    return molDf
#####################################################################
def decide_atom_to_delete_C_termini(atomNames):
    atomNames = sorted([atomName for atomName in atomNames if atomName.startswith("O")])
    atomToDelete = atomNames[-1]

    return atomToDelete
#####################################################################
def decide_atom_to_delete_N_termini(atomNames):
    atomNames = sorted([atomName for atomName in atomNames if atomName.startswith("H")])
    atomToDelete = atomNames[-1]
    return atomToDelete

#####################################################################
def find_bonded_atoms(molDf, atomName, distanceCutoff=1.6):
    tmpDf = molDf.copy()

    atomDf = molDf[molDf["ATOM_NAME"] == atomName]
    atomX, atomY, atomZ = atomDf.loc[:,['X', 'Y', 'Z']].astype(float).iloc[0]

    tmpDf['distance_to_atom'] = tmpDf.apply(calculate_distance,
                                             axis=1, atom_x=atomX, atom_y=atomY, atom_z=atomZ)

    bondedDf = tmpDf[tmpDf['distance_to_atom'] < distanceCutoff]
    bondedAtoms = bondedDf["ATOM_NAME"].to_list()

    return bondedAtoms

#####################################################################
def calculate_distance(row, atom_x, atom_y, atom_z):
    return np.sqrt((row['X'] - atom_x) ** 2 +
                   (row['Y'] - atom_y) ** 2 +
                   (row['Z'] - atom_z) ** 2)

#####################################################################


def find_capping_pdbs():
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
    configYaml = "/home/esp/scriptDevelopment/drFrankenstein/02_NMH_outputs/drFrankenstein.yaml"
    with open(configYaml, "r") as yamlFile:
        config = yaml.safe_load(yamlFile)

    capping_protocol(config)