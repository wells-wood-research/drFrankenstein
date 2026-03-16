
import os
from os import path as p

import numpy as np
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler
import pandas as pd
from collections import defaultdict


## CLEAN CODE ##
class DirectoryPath:
    pass
class FilePath:
    pass

def average_dict_values(list_of_strength_dicts):
    """
    Calculates the average strength for each unique key across a list of dictionaries.
    """
    sums = defaultdict(float)
    counts = defaultdict(int)

    for strength_dict in list_of_strength_dicts:
        for key, strength in strength_dict.items():
            sums[key] += strength
            counts[key] += 1
    
    averages = {key: sums[key] / counts[key] for key in sums}
    return averages


def get_compliance_matrix(hessianMatrix: np.ndarray, tol: float = 1e-4) -> np.ndarray:
    """
    Calculates the compliance matrix by inverting the Hessian, 
    explicitly removing near-zero eigenvalues (translations and rotations)
    that would otherwise blow up to massive values.
    """
    eigenvalues, eigenvectors = np.linalg.eigh(hessianMatrix)
    
    inv_eigenvalues = np.zeros_like(eigenvalues)
    
    for i, eig in enumerate(eigenvalues):
        if abs(eig) > tol:
            inv_eigenvalues[i] = 1.0 / eig
            
    complianceMatrix = eigenvectors @ np.diag(inv_eigenvalues) @ eigenvectors.T
    
    return complianceMatrix


def measure_angle(pdbDf: pd.DataFrame, angles: list[tuple[int, int, int]]) -> dict[tuple[int, int, int], float]:
    theta0values = {}

    for i, j, k in angles:
        coord_i = pdbDf.loc[i, ["X", "Y", "Z"]].values
        coord_j = pdbDf.loc[j, ["X", "Y", "Z"]].values
        coord_k = pdbDf.loc[k, ["X", "Y", "Z"]].values
        
        vec1 = coord_i - coord_j
        vec2 = coord_k - coord_j
        
        cosine_angle = np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2))
        cosine_angle = np.clip(cosine_angle, -1.0, 1.0)
        
        rad_angle = np.arccos(cosine_angle)
        theta0values[(i, j, k)] = np.degrees(rad_angle)

    return theta0values


def measure_bond_length(pdbDf: pd.DataFrame, bonds: list[tuple[int, int]]) -> float:
    r0values = {}

    for i, j in bonds:
        vec = pdbDf.loc[[i, j], ["X", "Y", "Z"]].values[0] - pdbDf.loc[[i, j], ["X", "Y", "Z"]].values[1]
        r0values[(i, j)] = np.linalg.norm(vec)

    return r0values


def _extract_and_scale_features(df, cols_to_ignore):
    """
    Removes ignored columns and standardizes the numerical features.
    """
    features_df = df.drop(columns=cols_to_ignore, errors='ignore')

    scaler = StandardScaler()
    scaled_features = scaler.fit_transform(features_df)

    return scaled_features


def _find_optimal_clusters(featuresScaled, kMax=10):
    """
    Iterates through possible values of K to find the best silhouette score.
    Returns the cluster labels for the optimal K.
    """
    nAtoms = featuresScaled.shape[0]
    
    if nAtoms < 3:
        return np.zeros(nAtoms, dtype=int)
        
    bestK = 2
    bestSilhouetteScore = -1.0
    bestLabels = np.zeros(nAtoms, dtype=int)
    
    kLimit = min(kMax, nAtoms - 1)
    if kLimit < 2:
        return bestLabels
        
    for k in range(2, kLimit + 1):
        kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
        labels = kmeans.fit_predict(featuresScaled)
        
        if len(set(labels)) > 1:
            score = silhouette_score(featuresScaled, labels)
            if score > bestSilhouetteScore:
                bestSilhouetteScore = score
                bestK = k
                bestLabels = labels
                
    return bestLabels


def find_non_bonded_lookup() -> FilePath:
    thisDir = p.dirname(p.abspath(__file__))
    labDir = p.dirname(p.dirname(thisDir))
    nonBondedLookup = p.join(labDir, "Ingredients", "NonBonded_LookUp.yaml")
    return nonBondedLookup
