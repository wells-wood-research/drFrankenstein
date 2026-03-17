
import os
from os import path as p

import numpy as np
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler
import pandas as pd
from collections import defaultdict
from subprocess import call, PIPE

## CLEAN CODE ##
class DirectoryPath:
    pass
class FilePath:
    pass


def run_tleap_to_make_params(config: dict) -> dict:

    """
    Uses TLEAP to create a prmtop inpcrd pair

    """
    inMol2: FilePath = config["runtimeInfo"]["madeByAssembly"]["moleculeMol2"]
    molFrcmod: FilePath = config["runtimeInfo"]["madeByAssembly"]["moleculeFrcmod"]
    outDir: DirectoryPath = config["runtimeInfo"]["madeByAssembly"]["assemblyDir"]
    moleculeName: str = config["moleculeInfo"]["moleculeName"]

    prmtop: FilePath = p.join(outDir, f"{moleculeName}_cappped.prmtop")
    inpcrd: FilePath = p.join(outDir, f"{moleculeName}_capped.inpcrd")

    tleapInput: FilePath = p.join(outDir, f"leap.in")
    with open(tleapInput, "w") as f:
        f.write("source leaprc.gaff2\n")
        f.write(f"mol  = loadmol2 {inMol2} \n")
        f.write(f"loadamberparams {molFrcmod} \n") # use frcmod previously made
        f.write(f"saveamberparm mol {prmtop} {inpcrd}  \n")
        f.write("quit")

    tleapOutput: FilePath = p.join(outDir, f"tleap.out")

    tleapCommand: list = ["tleap", "-f", tleapInput, ">", tleapOutput]

    call(tleapCommand)

    config["runtimeInfo"]["madeByAssembly"]["assembledPrmtop"] = prmtop

    return config


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
def init_mass_lookup() -> dict:
    massLookup = {
    'H': 1.0080, 'He': 4.0026, 'Li': 6.9400, 'Be': 9.0122, 'B': 10.8100, 'C': 12.0110, 'N': 14.0070, 'O': 15.9990,
    'F': 18.9984, 'Ne': 20.1797, 'Na': 22.9898, 'Mg': 24.3050, 'Al': 26.9815, 'Si': 28.0850, 'P': 30.9738, 'S': 32.0600,
    'Cl': 35.4500, 'Ar': 39.9480, 'K': 39.0983, 'Ca': 40.0780, 'Sc': 44.9559, 'Ti': 47.8670, 'V': 50.9415, 'Cr': 51.9961,
    'Mn': 54.9380, 'Fe': 55.8450, 'Co': 58.9332, 'Ni': 58.6934, 'Cu': 63.5460, 'Zn': 65.3800, 'Ga': 69.7230, 'Ge': 72.6300,
    'As': 74.9216, 'Se': 78.9710, 'Br': 79.9040, 'Kr': 83.7980, 'Rb': 85.4678, 'Sr': 87.6200, 'Y': 88.9058, 'Zr': 91.2240,
    'Nb': 92.9064, 'Mo': 95.9500, 'Tc': 98.0000, 'Ru': 101.0700, 'Rh': 102.9055, 'Pd': 106.4200, 'Ag': 107.8682, 'Cd': 112.4140,
    'In': 114.8180, 'Sn': 118.7100, 'Sb': 121.7600, 'Te': 127.6000, 'I': 126.9045, 'Xe': 131.2930, 'Cs': 132.9055, 'Ba': 137.3270,
    'La': 138.9055, 'Ce': 140.1160, 'Pr': 140.9077, 'Nd': 144.2420, 'Pm': 145.0000, 'Sm': 150.3600, 'Eu': 151.9640, 'Gd': 157.2500,
    'Tb': 158.9253, 'Dy': 162.5000, 'Ho': 164.9303, 'Er': 167.2590, 'Tm': 168.9342, 'Yb': 173.0540, 'Lu': 174.9668, 'Hf': 178.4900,
    'Ta': 180.9479, 'W': 183.8400, 'Re': 186.2070, 'Os': 190.2300, 'Ir': 192.2170, 'Pt': 195.0840, 'Au': 196.9665, 'Hg': 200.5920,
    'Tl': 204.3800, 'Pb': 207.2000, 'Bi': 208.9804, 'Po': 209.0000, 'At': 210.0000, 'Rn': 222.0000, 'Fr': 223.0000, 'Ra': 226.0000,
    'Ac': 227.0000, 'Th': 232.0377, 'Pa': 231.0359, 'U': 238.0289, 'Np': 237.0000, 'Pu': 244.0000, 'Am': 243.0000, 'Cm': 247.0000,
    'Bk': 247.0000, 'Cf': 251.0000, 'Es': 252.0000, 'Fm': 257.0000, 'Md': 258.0000, 'No': 259.0000, 'Lr': 262.0000, 'Rf': 267.0000,
    'Db': 270.0000, 'Sg': 269.0000, 'Bh': 270.0000, 'Hs': 270.0000, 'Mt': 278.0000, 'Ds': 281.0000, 'Rg': 282.0000, 'Cn': 285.0000,
    'Nh': 286.0000, 'Fl': 289.0000, 'Mc': 290.0000, 'Lv': 293.0000, 'Ts': 294.0000, 'Og': 294.0000,
}
    return massLookup