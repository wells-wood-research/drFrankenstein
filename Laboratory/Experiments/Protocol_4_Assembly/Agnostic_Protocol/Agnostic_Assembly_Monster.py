import os
from os import path as p

import yaml

import numpy as np
from scipy.spatial.distance import pdist, squareform
from collections import defaultdict
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler
import pandas as pd
import networkx as nx
from subprocess import call , PIPE
from pdbUtils import pdbUtils


from OperatingTools import file_parsers
from . import Agnostic_Assembly_Assistant
from .Agnostic_Assembly_Assistant import (
    average_dict_values,
    get_compliance_matrix,
    measure_angle,
    measure_bond_length,
    _extract_and_scale_features,
    _find_optimal_clusters,
    find_non_bonded_lookup
)

## CLEAN CODE ##
class DirectoryPath:
    pass
class FilePath:
    pass

def calculate_angleStrengths(complianceMatrix: np.ndarray, angles: list, pdbDf: pd.DataFrame) -> dict:
    """
    Calculates the relaxed force constants (angle strengths) for a list of angles.
    
    Parameters:
    hessianMatrix (np.ndarray): 3Nx3N Cartesian Hessian matrix (atomic units).
    angles (list of tuples): List of (i, j, k) indices where j is the central atom.
    pdbDf (pd.DataFrame): DataFrame containing 'X', 'Y', 'Z' coordinates.
    
    Returns:
    dict: Mapping of angle tuple (i, j, k) to its AMBER-compatible force constant in kcal/(mol * rad^2).
    """

    # 2. Extract coordinates
    coords = pdbDf[['X', 'Y', 'Z']].to_numpy()
    N = len(coords)
    
    angleStrengths = {}
    
    # Conversion factor from atomic units to kcal/mol
    AU_TO_KCAL_PER_MOL_PER_ANG = 2240.874995
    
    for i, j, k in angles:
        # Vectors from the central atom j to outer atoms i and k
        vec_ji = coords[i] - coords[j]
        vec_jk = coords[k] - coords[j]
        
        norm_ji = np.linalg.norm(vec_ji)
        norm_jk = np.linalg.norm(vec_jk)
        
        # Skip if atoms overlap (shouldn't happen in a valid geometry)
        if norm_ji == 0 or norm_jk == 0:
            continue
            
        # Unit vectors
        u_ji = vec_ji / norm_ji
        u_jk = vec_jk / norm_jk
        
        # Calculate the angle using dot product
        cos_theta = np.dot(u_ji, u_jk)
        # Clip to [-1.0, 1.0] to prevent floating point issues with arccos
        cos_theta = np.clip(cos_theta, -1.0, 1.0)
        sin_theta = np.sqrt(1.0 - cos_theta**2)
        
        # If the angle is too close to 0 or 180 degrees, the B-vector diverges (collinear atoms)
        if sin_theta < 1e-5:
            angleStrengths[(i, j, k)] = 0.0
            continue
        
        # Construct the Wilson B-vector for angle bending (length 3N)
        b = np.zeros(3 * N)
        
        # Gradients of the internal angle w.r.t coordinates
        b_i = (cos_theta * u_ji - u_jk) / (norm_ji * sin_theta)
        b_k = (cos_theta * u_jk - u_ji) / (norm_jk * sin_theta)
        b_j = -b_i - b_k  # The central atom moves in the opposite direction
        
        b[3*i : 3*i+3] = b_i
        b[3*j : 3*j+3] = b_j
        b[3*k : 3*k+3] = b_k
        
        # Calculate compliance constant for this angle
        c_val = np.dot(b.T, np.dot(complianceMatrix, b))
        
        # Relaxed force constant
        if c_val > 0:
            # Multiply by 0.5 for AMBER compatibility: E = K(theta - theta_eq)^2
            k_relaxed = (1.0 / c_val) * AU_TO_KCAL_PER_MOL_PER_ANG * 0.5
        else:
            k_relaxed = 0.0
            
        angleStrengths[(i, j, k)] = k_relaxed
        
    return angleStrengths

def calculate_bondStrengths(complianceMatrix: np.ndarray, bonds: list, pdbDf: pd.DataFrame) -> dict:
    """
    Calculates the relaxed force constants (bond strengths) for a list of bonds.
    
    Parameters:
    hessianMatrix (np.ndarray): 3Nx3N Cartesian Hessian matrix (atomic units).
    bonds (list of tuples): List of (atom_i, atom_j) indices.
    pdbDf (pd.DataFrame): DataFrame containing 'X', 'Y', 'Z' coordinates.
    
    Returns:
    dict: Mapping of bond tuple (i, j) to its relaxed force constant in kcal/mol/A^2.
    """

    
    # 2. Extract coordinates to build the geometric vectors
    coords = pdbDf[['X', 'Y', 'Z']].to_numpy()
    N = len(coords)
    
    bondStrengths = {}
    
    # Conversion factor from atomic units (Hartree/Bohr^2) to mdyn/Angstrom
    AU_TO_MDYN_ANG = 15.5689

    AU_TO_KCAL_PER_MOL_PER_ANG = 2240.874995
    
    for i, j in bonds:
        # Calculate the unit vector pointing from Atom i to Atom j
        vec = coords[j] - coords[i]
        norm = np.linalg.norm(vec)
        
        if norm == 0:
            continue
            
        u = vec / norm
        
        # Construct the Wilson B-vector for a bond stretch (length 3N)
        b = np.zeros(3 * N)
        
        # The derivative of the bond length w.r.t coordinates
        b[3*i : 3*i+3] = -u  # Atom i moves away
        b[3*j : 3*j+3] =  u  # Atom j moves away
        
        # Calculate compliance constant for this specific internal coordinate
        # c = b^T * C * b
        c_val = np.dot(b.T, np.dot(complianceMatrix, b))
        
        # Relaxed force constant is the inverse of the compliance constant
        if c_val > 0:
            k_relaxed = (1.0 / c_val) * AU_TO_KCAL_PER_MOL_PER_ANG * 0.5 # additional 1/2 for amber compatability
        else:
            # If c_val is negative, it implies an imaginary frequency (transition state) 
            # or unconverged geometry along this mode.
            k_relaxed = 0.0 
            
        bondStrengths[(i, j)] = k_relaxed
        
    return bondStrengths

def extract_topology_lists(adjacencyMatrix: np.ndarray, angleMatrix: np.ndarray):
    """
    Converts Adjacency and Angle matrices into lists of unique integer tuples.
    
    Returns:
    bonds: list of (i, j) tuples representing bonds.
    angles: list of (i, j, k) tuples representing angles (j is the vertex).
    """

    bond_i, bond_j = np.where(np.triu(adjacencyMatrix, k=1) == 1)
    
    # Convert numpy int arrays to standard Python ints inside tuples
    bonds =[(int(i), int(j)) for i, j in zip(bond_i, bond_j)]
    
    # Find all coordinates where the 3D matrix is 1
    ang_i, ang_j, ang_k = np.where(angleMatrix == 1)
    
    # Create a boolean mask to keep only the angles where i < k
    # This filters out the duplicate reverse angles (k-j-i)
    unique_mask = ang_i < ang_k
    
    # Apply the mask and zip into tuples
    angles =[
        (int(i), int(j), int(k)) 
        for i, j, k in zip(ang_i[unique_mask], ang_j[unique_mask], ang_k[unique_mask])
    ]
    
    return bonds, angles
def construct_connectivity_matrices(df, cutoff=1.6):
    """
    Convert a PDB DataFrame directly to an Adjacency Matrix and an Angle Matrix.
    Assumes the DataFrame is already in the correct atom order.
    
    Parameters:
    df (pd.DataFrame): Must contain 'X', 'Y', 'Z' columns.
    cutoff (float): Distance threshold for a bond in Angstroms.
    
    Returns:
    A (np.ndarray): N x N Adjacency Matrix (1 for bond, 0 for no bond).
    Theta (np.ndarray): N x N x N Angle Matrix in radians.
    """
    # 1. Extract coordinates as an N x 3 numpy array
    coords = df[['X', 'Y', 'Z']].to_numpy()
    N = len(coords)
    
    # 2. Create the Adjacency Matrix
    # pdist calculates pairwise distances, squareform turns it into an NxN matrix
    dist_matrix = squareform(pdist(coords))
    
    # A bond exists if distance is less than cutoff
    adjacencyMatrix = (dist_matrix < cutoff).astype(int)
    
    # Ensure atoms aren't bonded to themselves (zero out the diagonal)
    np.fill_diagonal(adjacencyMatrix, 0)
    
    angleMatrix = adjacencyMatrix[:,:,None] * adjacencyMatrix[None,:,:]

    indexes = np.arange(N)

    angleMatrix[indexes,:,indexes] = 0
                        
    return adjacencyMatrix, angleMatrix

def parse_orca_hessian(hessianFile: FilePath):
    """
    Parses an ORCA .hess file and returns the Cartesian Hessian 
    as a 2D numpy array.
    """
    with open(hessianFile, 'r') as f:
        lines = f.readlines()
        
    in_hessian = False
    dim = 0
    hessianMatrix = None
    col_indices =[]
    
    for i, line in enumerate(lines):
        line = line.strip()
        
        # Skip empty lines
        if not line:
            continue
            
        # Detect the start of the Hessian block
        if line == "$hessian":
            in_hessian = True
            # The line immediately after "$hessian" is the dimension (e.g., 66)
            dim = int(lines[i+1].strip())
            hessianMatrix = np.zeros((dim, dim))
            continue
            
        if in_hessian:
            # Stop if we hit the next block (e.g., $vibrational_frequencies)
            if line.startswith("$") and line != "$hessian":
                break
                
            # Skip the dimension line (since we already processed it)
            if line == str(dim):
                continue
                
            parts = line.split()
            
            # Differentiate between column header lines and data lines.
            # Data lines contain floats (which will have a '.' or 'E' in them).
            if not any('E' in p or '.' in p for p in parts):
                # This is a column index header (e.g., "0 1 2 3 4")
                col_indices =[int(p) for p in parts]
            else:
                # This is a row data line (e.g., "0  8.242E-01 -2.091E-01 ...")
                row_idx = int(parts[0])
                values = [float(p) for p in parts[1:]]
                
                # Assign values to their proper place in the NxN matrix
                for j, val in enumerate(values):
                    hessianMatrix[row_idx, col_indices[j]] = val
                    
    return hessianMatrix



def get_dihedral_indices(adjacencyMatrix):
    # Find all bonded pairs (central bond J-K)
    
    Js, Ks = np.where(adjacencyMatrix == 1)
    
    dihedrals =[]
    for j, k in zip(Js, Ks):
        # Find all I bonded to J
        Is = np.where(adjacencyMatrix[j] == 1)[0]
        # Find all L bonded to K
        Ls = np.where(adjacencyMatrix[k] == 1)[0]
        
        for i in Is:
            if i == k: continue # Prevent I-J-I
            for l in Ls:
                if l == j or l == i: continue # Prevent J-K-J and I-J-K-I
                
                # We found a valid I-J-K-L path!
                dihedrals.append((i, j, k, l))
                
    return dihedrals



def get_improper_indices(adjacencyMatrix):
    """
    Finds improper angle groups where exactly three atoms are bonded to a central atom.
    
    Parameters:
    adjacencyMatrix (np.ndarray): N x N Adjacency Matrix (1 for bond, 0 for no bond).
    
    Returns:
    list[tuple[int, int, int, int]]: A list of (central_atom, i, j, k) indices.
    """
    impropers =[]
    
    # 1. Sum along the rows to get the number of bonds for each atom
    degrees = np.sum(adjacencyMatrix, axis=1)
    
    # 2. Find all central atoms that have exactly 3 bonds
    central_atoms = np.where(degrees == 3)[0]
    
    # 3. For each central atom, get the indices of its neighbors
    for central in central_atoms:
        bonded_atoms = np.where(adjacencyMatrix[central] == 1)[0]
        
        # Unpack the 3 surrounding atoms
        i, j, k = bonded_atoms
        
        # Convert numpy int to standard Python int for the final tuple
        # NOTE: central as 3rd element
        impropers.append((int(i), int(j), int(central), int(k)))
        
    return impropers

def get_molecule_connectivity(config:dict) -> tuple[list[tuple[int, int]], list[tuple[int, int, int]]]:
    cappedPdb = config["runtimeInfo"]["madeByCapping"]["cappedPdb"]
    cappedDf = pdbUtils.pdb2df(cappedPdb)
    adjacencyMatrix , angleMatrix = construct_connectivity_matrices(cappedDf)
    ## get a list of bonds and angles from the matrices
    bonds, angles = extract_topology_lists(adjacencyMatrix, angleMatrix)

    dihedrals = get_dihedral_indices(adjacencyMatrix)
    impropers = get_improper_indices(adjacencyMatrix)

    return bonds, angles, dihedrals, impropers



def process_hessian(hessianFile: FilePath, xyzFile: FilePath,
                    bonds: list[tuple[int, int]], angles: list[tuple[int, int, int]],
                                                     config:dict) -> tuple[dict[tuple[int, int], float],        ## bonds strengths
                                                                            dict[tuple[int, int], float],       ## bond lengths
                                                                            dict[tuple[int, int, int], float],  ## angle strengths
                                                                            dict[tuple[int, int, int], float]]: ## angle lengths
    ## extract hessian from orca .hess file
    hessianMatrix = parse_orca_hessian(hessianFile)
    ## invert hessian to get compliance matrix
    complianceMatrix = get_compliance_matrix(hessianMatrix)

    ## construct bond and angle matrices from pdb 
    ## TODO: replace XYZ from opt_freq calculation
    cappedPdb = config["runtimeInfo"]["madeByCapping"]["cappedPdb"]
    pdbDf = pdbUtils.pdb2df(cappedPdb)
    xyzDf = file_parsers.xyz2df(xyzFile)
    pdbDf.loc[:, ["X", "Y", "Z"]] = xyzDf[["x", "y", "z"]].values    



    ## get bond and angle strengths
    kbValues = calculate_bondStrengths(complianceMatrix, bonds, pdbDf)
    kThetaValues = calculate_angleStrengths(complianceMatrix, angles, pdbDf)

    r0values = measure_bond_length(pdbDf, bonds)
    theta0values = measure_angle(pdbDf, angles)

    return kbValues, r0values, kThetaValues,  theta0values


def construct_atom_features(kbValues, r0values, kThetaValues, theta0values, config:dict) -> tuple[pd.DataFrame, dict]:
    pdbDf = pdbUtils.pdb2df(config["runtimeInfo"]["madeByCapping"]["cappedPdb"])

    atomFeaturesDict = pdbDf.loc[:,["ATOM_ID", "ATOM_NAME", "ELEMENT"]].copy().T.to_dict()

    
    for bondedTerm, keyName in zip([kbValues, r0values], ["kB_Values", "r0_Values"]):
        for (i, j) in bondedTerm:
            if not keyName in atomFeaturesDict[i]:
                atomFeaturesDict[i][keyName]  = []
            if not keyName in atomFeaturesDict[j]:
                atomFeaturesDict[j][keyName]  = []
            atomFeaturesDict[i][keyName].append(bondedTerm[(i, j)])
            atomFeaturesDict[j][keyName].append(bondedTerm[(i, j)])

    for angleTerm, keyName in zip([kThetaValues, theta0values], ["kTheta_Values", "theta0_Values"]):
        for (i, j, k) in angleTerm:
            if not keyName in atomFeaturesDict[i]:
                atomFeaturesDict[i][keyName]  = []
            if not keyName in atomFeaturesDict[j]:
                atomFeaturesDict[j][keyName]  = []
            if not keyName in atomFeaturesDict[k]:
                atomFeaturesDict[k][keyName]  = []
            atomFeaturesDict[i][keyName].append(angleTerm[(i, j, k)])
            atomFeaturesDict[j][keyName].append(angleTerm[(i, j, k)])
            atomFeaturesDict[k][keyName].append(angleTerm[(i, j, k)])

    dictForClustering = atomFeaturesDict.copy()

    for atomZeroIndex, values in dictForClustering.items():
        dictForClustering[atomZeroIndex]["kB_Mean"] = np.mean(dictForClustering[atomZeroIndex]["kB_Values"])
        dictForClustering[atomZeroIndex]["r0_Mean"] = np.mean(dictForClustering[atomZeroIndex]["r0_Values"])
        dictForClustering[atomZeroIndex]["kTheta_Mean"] = np.mean(dictForClustering[atomZeroIndex]["kTheta_Values"])
        dictForClustering[atomZeroIndex]["theta0_Mean"] = np.mean(dictForClustering[atomZeroIndex]["theta0_Values"])
        dictForClustering[atomZeroIndex]["nBonds"] = len(dictForClustering[atomZeroIndex]["kB_Values"])
    

    ## remove >1 dimentional entries
    dictForClustering = {idx: {k: v for k, v in atom.items()
                               if not k.endswith("_Values")} 
                               for idx, atom in atomFeaturesDict.items()}

    atomFeaturesDf = pd.DataFrame(dictForClustering).T


    chargesCsv = config["runtimeInfo"]["madeByCharges"]["chargesCsv"]
    chargesDf = pd.read_csv(chargesCsv)
    atomFeaturesDf["Charge"] = chargesDf["Charge"].values

    return atomFeaturesDf, atomFeaturesDict

def group_params_by_atom_type(atomTypesDf, kbValues, r0values, kThetaValues, theta0values, dihedrals, impropers):

    atomTypes = atomTypesDf["ATOM_TYPE"].to_dict()

    bondedParams = {}
    for (atomZeroIndexes, kB), r0 in zip(kbValues.items(), r0values.values()):
        atomTypeI = atomTypes[atomZeroIndexes[0]]
        atomTypeJ = atomTypes[atomZeroIndexes[1]]
        if not (atomTypeI, atomTypeJ) in bondedParams:
            bondedParams[(atomTypeI, atomTypeJ)] = []
        bondedParams[(atomTypeI, atomTypeJ)].append({"kB": kB, "r0": r0})

    meanBondParams = {}
    for key, value in bondedParams.items():
        meanBondParams[key] = {"kB": np.mean([v["kB"] for v in value]), "r0": np.mean([v["r0"] for v in value])}

    ##
    angleParams = {}
    for (atomZeroIndexes, kTheta), theta0 in zip(kThetaValues.items(), theta0values.values()):
        atomTypeI = atomTypes[atomZeroIndexes[0]]
        atomTypeJ = atomTypes[atomZeroIndexes[1]]
        atomTypeK = atomTypes[atomZeroIndexes[2]]
        if not (atomTypeI, atomTypeJ, atomTypeK) in angleParams:
            angleParams[(atomTypeI, atomTypeJ, atomTypeK)] = []
        angleParams[(atomTypeI, atomTypeJ, atomTypeK)].append({"kTheta": kTheta, "theta0": theta0})

    meanAngleParams = {}
    for key, value in angleParams.items():
        meanAngleParams[key] = {"kTheta": np.mean([v["kTheta"] for v in value]), "theta0": np.mean([v["theta0"] for v in value])}

    ##
    groupedDihedrals = set()
    for dihedral in dihedrals:  
        groupedDihedrals.add((atomTypes[dihedral[0]], atomTypes[dihedral[1]], atomTypes[dihedral[2]], atomTypes[dihedral[3]]))

    groupedImpropers = set()
    for improper in impropers:
        groupedImpropers.add((atomTypes[improper[0]], atomTypes[improper[1]], atomTypes[improper[2]], atomTypes[improper[3]]))

    return meanBondParams, meanAngleParams, groupedDihedrals, groupedImpropers


def reassign_capping_types(atomTypesDf, config) -> tuple[dict, pd.DataFrame]:
    amberDefaultMap = {
                    "N_N": "N",
                    "H_N": "H",
                    "C_N": "CX", 
                    "H1_N": "H3",
                    "H2_N": "H3",
                    "H2_N": "H3",
                    "C_C": "C",
                    "O_C": "O",
                    "CM_C": "CX", 
                    "H1_C": "H3",
                    "H2_C": "H3",
                    "H3_C": "H3"
     }
    atomTypesDf["ATOM_TYPE"] = atomTypesDf.apply(lambda row: amberDefaultMap.get(row["ATOM_NAME"], row["ATOM_TYPE"]), axis=1)

    return atomTypesDf

def reassign_backbone_types(atomTypesDf, config) -> tuple[dict, pd.DataFrame]:
    amberDefaultMap = {
                    "N_N": "N",
                    "H_N": "H",
                    "C_N": "CX", 
                    "H1_N": "H3",
                    "H2_N": "H3",
                    "H2_N": "H3",
                    "C_C": "C",
                    "O_C": "O",
                    "CM_C": "CX", 
                    "H1_C": "H3",
                    "H2_C": "H3",
                    "H3_C": "H3"
     }
    atomTypesDf["ATOM_TYPE"] = atomTypesDf.apply(lambda row: amberDefaultMap.get(row["ATOM_NAME"], row["ATOM_TYPE"]), axis=1)

    return atomTypesDf




def assign_atom_types_by_clustering(atomFeaturesDf, maxClusters=10):
    """
    Main function. Groups DataFrame by ELEMENT, determines the optimal
    number of clusters for each element, and assigns a cluster label.
    """
    # Work on a copy to avoid SettingWithCopyWarning
    atomTypesDf = atomFeaturesDf.copy()
    
    # Initialize the new columns
    atomTypesDf['CLUSTER_LABEL'] = -1
    atomTypesDf['ATOM_TYPE'] = ""
    
    colsToDrop =['ATOM_ID', 'ATOM_NAME', 'ELEMENT', 'CLUSTER_LABEL', 'ATOM_TYPE']
    
    atomTypeSuffixes = pd.Series([f"{letter}" for letter in "zyxwvutsrqponmlkjihgfedcba"])

    for element, elementDf in atomTypesDf.groupby("ELEMENT"):
            
            # 1. Preprocess and scale features
            featuresScaled = _extract_and_scale_features(elementDf, colsToDrop)
            
            # 2. Find best K and get labels (assuming labels are 0, 1, 2...)
            labels = _find_optimal_clusters(featuresScaled, kMax=maxClusters)
            
            # 3. Assign labels back to the DataFrame corresponding to their original index
            atomTypesDf.loc[elementDf.index, 'CLUSTER_LABEL'] = labels
            
            # Create a Series of the labels aligned with the current element's indices
            cluster_series = pd.Series(labels, index=elementDf.index)
            
            # Map the integer cluster labels to their corresponding suffix 
            # (e.g., 0 -> 'v1', 1 -> 'v2')
            mapped_suffixes = cluster_series.map(atomTypeSuffixes)
            
            # Concatenate the Element name and the mapped suffix
            atomTypesDf.loc[elementDf.index, 'ATOM_TYPE'] = element.lower() + mapped_suffixes

    
    return  atomTypesDf


def  assign_non_bonded_by_analogy(atomTypesDf):
    nonBondedLookup, absentNonBondedLookup = find_non_bonded_lookup()

    with open(nonBondedLookup, "r") as f:
        nonBondedLookup = yaml.safe_load(f)

    with open(absentNonBondedLookup, "r") as f:
        absentNonBondedLookup = yaml.safe_load(f)

    nonBondedParams = {}
    for atomType, atomTypeDf, in atomTypesDf.groupby("ATOM_TYPE"):
        meanCharge = atomTypeDf["Charge"].mean()
        element = atomTypeDf["ELEMENT"].iloc[0]
        ## use gaff2 lookup (by-analogy will be for the correct element)
        elementLookup = nonBondedLookup.get(element, None)
        if elementLookup is None:   ## if element not in gaff2 lookup, use the absent lookup. This will be by-analogy to the closest element in the gaff2 lookup based on charge
            elementLookup = absentNonBondedLookup.get(element, None)
            print(f"WARNING: Non-Bonded Parameters for element {element} not found in GAFF2. Non-Bonded parameters will be assigned by-analogy to the closest available element")
        if elementLookup is None:
            raise ValueError(f"No non-bonded parameters found for element {element}\n Please check your input PDB file.")
        elementLookupDf = pd.DataFrame.from_dict(elementLookup).T
        elementLookupDf["Delta_Charge"] = elementLookupDf["CHARGE"] - meanCharge
        minDeltaCharge = elementLookupDf["Delta_Charge"].min()
        closestTypesDf = elementLookupDf[elementLookupDf["Delta_Charge"] == minDeltaCharge]
        meanRadius = closestTypesDf["RADIUS"].mean()
        meanWellDepth = closestTypesDf["WELL_DEPTH"].mean()

        nonBondedParams[atomType] = {
            "RADIUS": meanRadius,
            "WELL_DEPTH": meanWellDepth
        }

    return nonBondedParams
  

def construct_frcmod(atomTypeMasses: dict,
                      bondParamsByType: dict, angleParamsByType: dict,
                        groupedDihedrals: list[tuple[str, str, str, str]], groupedImpropers: list[tuple[str]],
                          nonBondedParams: dict,
                          config: dict) -> dict:
    assemblyDir = config["runtimeInfo"]["madeByAssembly"]["assemblyDir"]
    moleculeName = config["moleculeInfo"]["moleculeName"]

    frcmodFile = f"{moleculeName}_assembled.frcmod"
    frcmodFile = p.join(assemblyDir, frcmodFile)


    with open(frcmodFile, "w") as f:
        f.write("CREATED BY drFRANKENSTEIN\n\n")
        f.write("MASS\n")
        for atomType, mass in atomTypeMasses.items():
            f.write(f"{atomType}\t{mass:.3f}\n")
        f.write("\n")
        f.write("BOND\n")
        for (atomI, atomJ), bondData in bondParamsByType.items():
            (atomI, atomJ) = (atom.ljust(2) for atom in (atomI, atomJ))
            kB = bondData["kB"]
            r0 = bondData["r0"]
            f.write(f"{atomI}-{atomJ}\t{kB:.2f}\t{r0:.3f}\n")
        f.write("\n")
        f.write("ANGLE\n")
        for (atomI, atomJ, atomK), angleData in angleParamsByType.items():
            (atomI, atomJ, atomK) = (atom.ljust(2) for atom in (atomI, atomJ, atomK))
            kTheta = angleData["kTheta"]
            theta0 = angleData["theta0"]
            f.write(f"{atomI}-{atomJ}-{atomK}\t{kTheta:.2f}\t{theta0:.3f}\n")
        f.write("\n")
        f.write("DIHE\n")
        (divfactor, kTorsion, phase, periodicity) = ("1", "0.000", "180.000", "1.000")
        for (atomI, atomJ, atomK, atomL) in groupedDihedrals:
            (atomI, atomJ, atomK, atomL) = (atom.ljust(2) for atom in (atomI, atomJ, atomK, atomL))
            f.write(f"{atomI}-{atomJ}-{atomK}-{atomL}\t{divfactor}\t{kTorsion}\t{phase}\t{periodicity}\t!DUMMY INPUT\n")
        f.write("\n")
        f.write("IMPROPER\n")
        for (atomI, atomJ, atomK, atomL) in groupedImpropers:
            (atomI, atomJ, atomK, atomL) = (atom.ljust(2) for atom in (atomI, atomJ, atomK, atomL))
            f.write(f"{atomI}-{atomJ}-{atomK}-{atomL}\t\t1.1\t180.000\t2.0\t!DEFAULT VALUE\n")
        f.write("\n")
        f.write("NONB\n")
        for atomType, nonBondedData in nonBondedParams.items():
            radius = nonBondedData["RADIUS"]
            wellDepth = nonBondedData["WELL_DEPTH"]
            f.write(f"{atomType}\t{radius:.3f}\t{wellDepth:.3f}\n")
        f.write("\n")
            
    config["runtimeInfo"]["madeByAssembly"]["assembledFrcmod"] = frcmodFile


    return config
############################################################
def construct_mol2(atomTypesDf: pd.DataFrame, config: dict):

    cappedPdb = config["runtimeInfo"]["madeByCapping"]["cappedPdb"]
    assemblyDir = config["runtimeInfo"]["madeByAssembly"]["assemblyDir"] 
    moleculeName = config["moleculeInfo"]["moleculeName"]
    chargesCsv = config["runtimeInfo"]["madeByCharges"]["chargesCsv"]
    chargesDf = pd.read_csv(chargesCsv, index_col="Unnamed: 0")

    chargeFittingProtocol = config["chargeFittingInfo"]["chargeFittingProtocol"]
    chargeQMMethod = config["chargeFittingInfo"]["singlePointMethod"]

    tmpMol2 = p.join(assemblyDir, "tmp.mol2")

    obabelCall = ["obabel", cappedPdb,  "-O", tmpMol2]
    call(obabelCall, stdout=PIPE, stderr=PIPE)
    
    mol2File = p.join(assemblyDir, f"{moleculeName}.mol2")

    with open(tmpMol2, "r") as inFile, open(mol2File, "w") as outFile:
        section = None
        
        for line in inFile:
            
            # --- Detect Block Headers ---
            if line.startswith("@<TRIPOS>MOLECULE"):
                section = "MOLECULE"
                outFile.write(line)
                
                # 1. Replace default Obabel molecule name with ours
                next(inFile) 
                outFile.write(f"{moleculeName}\n")
                
                # 2. Fix the header counts to declare exactly 1 substructure
                counts_line = next(inFile)
                counts = counts_line.split()
                if len(counts) >= 3:
                    counts[2] = "1" # 3rd integer represents the number of substructures
                outFile.write(" " + " ".join(counts) + "\n")
                continue
                
            elif line.startswith("@<TRIPOS>ATOM"):
                section = "ATOM"
                outFile.write(line)
                continue
                
            elif line.startswith("@<TRIPOS>BOND"):
                section = "BOND"
                outFile.write(line)
                continue
                
            elif line.startswith("@<TRIPOS>SUBSTRUCTURE"):
                section = "SUBSTRUCTURE"
                # We skip OpenBabel's substructure block entirely
                continue

            # --- Process Block Contents ---
            if section == "ATOM":
                if not line.strip():
                    continue
                atomData = line.split()
                atomName = atomData[1]
                
                # Assign forcefield atom type
                atomType = atomTypesDf[atomTypesDf["ATOM_NAME"] == atomName]["ATOM_TYPE"].iloc[0]
                atomData[5] = atomType
                
                # Force all atoms into a single molecular residue for tleap
                atomData[6] = "1"     # Residue Number
                atomData[7] = moleculeName   # Residue Name
                
                # Assign partial charges
                atomCharge = chargesDf[chargesDf["ATOM_NAME"] == atomName]["Charge"].iloc[0]
                atomData[8] = f"{atomCharge:.4f}"
                
                outFile.write("\t".join(atomData) + "\n")
                
            elif section == "BOND":
                if not line.strip():
                    continue
                bondData = line.split()
                
                # Convert TRIPOS 'am' (amide) bond orders to standard '1' 
                if len(bondData) >= 4 and bondData[3].lower() == "am":
                    bondData[3] = "1"
                    
                outFile.write("\t".join(bondData) + "\n")
                
            elif section == "MOLECULE":
                # Inject the custom QM method name instead of GASTEIGER
                if "GASTEIGER" in line or "NO_CHARGES" in line:
                    outFile.write(f"{chargeFittingProtocol} [{chargeQMMethod}]\n")
                else:
                    outFile.write(line)
                    
            elif section == "SUBSTRUCTURE":
                # Do nothing; wait for the end of the file
                pass
                
            else:
                outFile.write(line)

        # --- Append Safe Substructure Block ---
        # AMBER's tleap needs this precise formatting to allocate the 1-4 bonding memory correctly
        outFile.write("@<TRIPOS>SUBSTRUCTURE\n")
        outFile.write("     1 MOL         1 TEMP              0 ****  ****    0 ROOT\n")

    config["runtimeInfo"]["madeByAssembly"]["cappedMol2"] = mol2File

    return config


def get_atom_type_masses(atomTypesDf: pd.DataFrame) -> dict:
    massLookup = Agnostic_Assembly_Assistant.init_mass_lookup()
    atomTypeMasses = {}
    for atomType, groupDf in atomTypesDf.groupby("ATOM_TYPE"):
        element = groupDf["ELEMENT"].iloc[0]
        mass = massLookup[element]
        atomTypeMasses[atomType] = mass
    return atomTypeMasses