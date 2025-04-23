import os
from os import path as p
from pdbUtils import pdbUtils
import pandas as pd
import numpy as np
import mdtraj as md
## CLEAN CODE ##
class FilePath:
    pass
class DirectoryPath:
    pass

##TODO: incorporate into final creation step??
def get_donor_acceptors(config: dict, debug:bool = False) -> dict:
    """
    Looks though optimised solvated conformers for hydrogen bonds donors and acceptors
    Creates a dict of DONOR and ACCEPTOR atom pairs,
    where the first entry is the atom participating in teh Hydrogen Bond, 
    and the second is the (a) heavy atom covalently bound to the first

    Args:
        config (dict): dictionary containing all run information
        debug (bool): Whether to run in debug mode (toggles multiprocessing)

    Returns:
        config (dict): updated config
    """

    ## unpack config ##
    chargeDir = config["runtimeInfo"]["madeByCharges"]["chargeDir"]
    moleculeName = config["moleculeInfo"]["moleculeName"]

    hydrogenBondDir = p.join(chargeDir, "05_donor_acceptors")
    os.makedirs(hydrogenBondDir, exist_ok=True)
    config["runtimeInfo"]["madeByCharges"]["hydrogenBondDir"] = hydrogenBondDir
    solvatedOptimisedXyzs = config["runtimeInfo"]["madeByCharges"]["solvatedOptimisedXyzs"]

    config, solvatedDf = pad_pdb_with_waters(config)

    allDonors = []
    allAcceptors = []

    i = 0
    for solvatedOptXyz in solvatedOptimisedXyzs:
        i += 1
        atoms, coords = (solvatedOptXyz)
        hBonds = detect_hydrogen_bonds(atoms, coords)
        donors, acceptors = construct_donor_acceptors(hBonds, solvatedDf, moleculeName)
        allDonors.append(donors)
        allAcceptors.append(acceptors)
    print(allDonors)
    print(allAcceptors)
    uniqueDonors = list(set(allDonors))
    uniqueAcceptors = list(set(allAcceptors))

    donorAcceptors = {"DONORS": uniqueDonors,
                      "ACCEPTORS": uniqueAcceptors}
    
    config["runtimeInfo"]["madeByCharges"]["donorAcceptors"] = donorAcceptors

    return config
    




def load_xyz(file_path):
    """Load an XYZ file and return atom symbols and coordinates."""
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    num_atoms = int(lines[0].strip())
    atoms = []
    coords = []
    
    for line in lines[2:2+num_atoms]:
        parts = line.strip().split()
        if len(parts) >= 4:
            atoms.append(parts[0])
            coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    
    return atoms, np.array(coords)

def calculate_distance(coord1, coord2):
    """Calculate Euclidean distance between two points."""
    return np.linalg.norm(coord1 - coord2)

def calculate_angle(coord1, coord2, coord3):
    """Calculate angle (in degrees) between three points (coord2 is vertex)."""
    vec1 = coord1 - coord2
    vec2 = coord3 - coord2
    
    norm1 = np.linalg.norm(vec1)
    norm2 = np.linalg.norm(vec2)
    
    # Check for zero-length vectors
    if norm1 == 0 or norm2 == 0:
        return 0.0  # Default to 0 degrees for undefined cases
    
    vec1_norm = vec1 / norm1
    vec2_norm = vec2 / norm2
    cos_theta = np.clip(np.dot(vec1_norm, vec2_norm), -1.0, 1.0)
    return np.degrees(np.arccos(cos_theta))

def detect_hydrogen_bonds(atoms, coords, distance_cutoff=3.5, angle_cutoff=120.0):
    """
    Detect hydrogen bonds based on geometric criteria.
    Criteria: D-H...A where D is donor (O, N), A is acceptor (O, N),
    H...A distance < distance_cutoff (Angstroms), D-H...A angle > angle_cutoff (degrees).
    Returns list of (donor_idx, hydrogen_idx, acceptor_idx) tuples.
    """
    hBonds = []
    donorAcceptorTypes = ["N", "O", "F"]  # Exclude C and N
    
    # Find potential donors (O, N) and hydrogens
    for i, (atom_i, coord_i) in enumerate(zip(atoms, coords)):
        if atom_i  in donorAcceptorTypes:  # Potential donor
            # Look for bonded hydrogens (simple distance-based check)
            for j, (atom_j, coord_j) in enumerate(zip(atoms, coords)):
                if atom_j == 'H':
                    dist_dh = calculate_distance(coord_i, coord_j)
                    if 0.8 < dist_dh < 1.2:  # Typical D-H bond length
                        # Found D-H pair, now look for acceptors
                        for k, (atom_k, coord_k) in enumerate(zip(atoms, coords)):
                            if atom_k  in donorAcceptorTypes and k != i:  # Potential acceptor
                                dist_ha = calculate_distance(coord_j, coord_k)
                                if dist_ha < distance_cutoff:  # H...A distance check
                                    angle_dha = calculate_angle(coord_i, coord_j, coord_k)
                                    if angle_dha > angle_cutoff:  # D-H...A angle check
                                        hBonds.append((i, j, k))

    return hBonds

def construct_donor_acceptors(hBonds: list, pdbDf: pd.DataFrame, moleculeName: str) -> dict:
    donors = []
    acceptors = []
    for hBond in hBonds:
        donorDf = pdbDf.iloc[hBond[0]]
        hydrogenDf = pdbDf.iloc[hBond[1]]
        acceptorDf = pdbDf.iloc[hBond[2]]
        if donorDf["RES_NAME"] == moleculeName:
            donorPair = (hydrogenDf["ATOM_NAME"], donorDf["ATOM_NAME"])
            donors.extend(donorPair)
        if acceptorDf["RES_NAME"] == moleculeName:
            bondedAtoms = find_bonded_atoms(pdbDf, acceptorDf["ATOM_NAME"])
            ## exclude tri-valent (or more!) Nitrogen atoms
            if acceptorDf["ELEMENT"] == "N" and len(bondedAtoms) > 2:
                print(bondedAtoms)
                continue
            bondedAtom = choose_bonded_atom_for_acceptor(bondedAtoms, pdbDf, moleculeName)
            acceptorPair = (acceptorDf["ATOM_NAME"], bondedAtom)
            acceptors.extend(acceptorPair)

    donors = list(set(donors))
    acceptors = list(set(acceptors))
            
    return donors, acceptors
    



def choose_bonded_atom_for_acceptor(bondedAtoms, pdbDf, moleculeName):
    bondedAtoms = [atom for atom in bondedAtoms if not atom.startswith("H")]
    ## remove 
    bondedAtomDf = pdbDf[pdbDf["ATOM_NAME"].isin(bondedAtoms)]
    bondedAtomDf = bondedAtomDf[bondedAtomDf["RES_NAME"] == moleculeName]
    chosenBondedAtom = bondedAtomDf["ATOM_NAME"].to_list()[0]
    return chosenBondedAtom
    



def pad_pdb_with_waters(config: dict) -> FilePath:
    """
    Takes original capped PDB file and adds water molecules with dummy coords
    
    Args:
        pdbFile (FilePath): an unsolvated PDB file
        nWaters (int): The number of water molecules to add
    """
    ## unpack config
    cappedPdb = config["runtimeInfo"]["madeByCapping"]["cappedPdb"]
    nWatersAdded = config["runtimeInfo"]["madeByCharges"]["nWaters"]

    pdbDf = pdbUtils.pdb2df(cappedPdb)
    maxResId = pdbDf["RES_ID"].max()

    dummyWaterDf = pd.DataFrame([
        {'ATOM': 'HETATM', 'ATOM_ID': i, 'ATOM_NAME': name, 'RES_NAME': 'HOH', 
         'CHAIN_ID': 'A', 'RES_ID': 1, 'X': 0.0, 'Y': 0.0, 'Z': 0.0, 
         'OCCUPANCY': 1.0, 'BETAFACTOR': 0.0, 'ELEMENT': elem}
        for i, (name, elem) in enumerate([('OH2', 'O'), ('H1', 'H'), ('H2', 'H')], 1)
    ])

    dfsToConcat = [pdbDf]
    for i in range(0, nWatersAdded):
        waterResId = maxResId + i + 1
        tmpDf = dummyWaterDf.copy()
        tmpDf["RES_ID"] = waterResId
        dfsToConcat.append(tmpDf)
    
    solvatedDf = pd.concat(dfsToConcat, axis=0)
    solvatedDf["ATOM_ID"] = np.arange(0,len(solvatedDf["ATOM_ID"])) + 1

    solvatedPdb = p.join(config["runtimeInfo"]["madeByCharges"]["hydrogenBondDir"], "solvated.pdb")
    pdbUtils.df2pdb(solvatedDf, solvatedPdb)

    config["runtimeInfo"]["madeByCharges"]["solvatedPdb"] = solvatedPdb

    return config, solvatedDf

def calculate_distance_apply(row: pd.Series,
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
    tmpDf['distance_to_atom'] = tmpDf.apply(calculate_distance_apply,
                                             axis=1, atomX=atomX, atomY=atomY, atomZ=atomZ)
    ## find all atoms within distance cutoff
    bondedDf = tmpDf[tmpDf['distance_to_atom'] < distanceCutoff]
    ## convert to list
    bondedAtoms = bondedDf["ATOM_NAME"].to_list()
    ## remove atom of interest
    bondedAtoms = [atom for atom in bondedAtoms if atom != atomName]
    return bondedAtoms


def main():
    config = {}
    config["runtimeInfo"] = {}
    config["runtimeInfo"]["madeByCharges"] = {}

    config["runtimeInfo"]["madeByCharges"]["chargeDir"] =    "/home/esp/scriptDevelopment/drFrankenstein/05_AIB_outputs/04_charge_calculations"
    qmmmDir = "/home/esp/scriptDevelopment/drFrankenstein/05_AIB_outputs/04_charge_calculations/02_QMMM_optimisations"
    xyzFiles = [p.join(qmmmDir,runName,runName+".xyz") for runName in os.listdir(qmmmDir)]

    config["runtimeInfo"]["madeByCharges"]["solvatedOptimisedXyzs"] = xyzFiles

    config["runtimeInfo"]["madeByCapping"] = {}
    config["runtimeInfo"]["madeByCapping"]["cappedPdb"] = "/home/esp/scriptDevelopment/drFrankenstein/05_AIB_outputs/01_termini_capping/AIB_capped.pdb"

    config["runtimeInfo"]["madeByCharges"]["nWaters"] = 44

    config["moleculeInfo"] = {}
    config["moleculeInfo"]["moleculeName"] = "AIB"
    get_donor_acceptors(config)

if __name__ == "__main__":
    main()