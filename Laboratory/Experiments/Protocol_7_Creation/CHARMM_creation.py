import os
from os import path as p
from pdbUtils import pdbUtils
import pandas as pd
import numpy as np
from shutil import copy
## MULTIPROCESSING AND LOADING BAR LIBRARIES ##
import multiprocessing as mp
from tqdm import tqdm

## drFrankenstein LIBRARIES ##e
from OperatingTools import file_parsers
from OperatingTools import symmetry_tool
from ..Protocol_1_Capping import Capping_Assistant
## CLEAN CODE ##
class FilePath:
    pass
class DirectoryPath:
    pass

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

def copy_final_prm(config: dict) -> None:
    """
    Copies final PRM file into the final creation directory

    Args:
        config (dict): contains all run information

    Returns:
        None
    """
    moleculePrm = config["runtimeInfo"]["madeByStitching"]["moleculePrm"]
    finalCreationDir = config["runtimeInfo"]["madeByCreator"]["finalCreationDir"]
    moleculeName = config["moleculeInfo"]["moleculeName"]
    finalPrm = p.join(finalCreationDir, f"{moleculeName}.prm")

    copy(moleculePrm, finalPrm)


def create_final_rtf(config):
    ## unpack config
    finalCreationDir = config["runtimeInfo"]["madeByCreator"]["finalCreationDir"]
    moleculeRtf = config["runtimeInfo"]["madeByStitching"]["moleculeRtf"]

    
    moleculeName = config["moleculeInfo"]["moleculeName"]
    nTerminialAtoms = config["moleculeInfo"]["nTermini"]
    cTerminalAtoms = config["moleculeInfo"]["cTermini"]
    donorAcceptors = config["runtimeInfo"]["madeByCreator"]["donorAcceptors"]
    finalRtf = p.join(finalCreationDir, f"{moleculeName}.rtf")

    cappingAtomNames = ["CN", "NN", "HNN1", "HCN1", "HCN2", "HCN3", "CC1", "OC", "CC2", "HC1", "HC2", "HC3"]

    terminalSectionWritten = False
    with open(moleculeRtf, "r") as inRtf, open(finalRtf, "w") as outRtf:
        for line in inRtf.readlines():
            if line.startswith("END"): 
                continue
            if any(atomName in line for atomName in cappingAtomNames):
                continue
            if line.startswith("BOND") and not terminalSectionWritten:
                for atomName in cTerminalAtoms:
                    outRtf.write(f"BOND {atomName} +N\n")
                for atomName in nTerminialAtoms:
                    outRtf.write(f"BOND {atomName} -C\n")
                terminalSectionWritten = True
            outRtf.write(line)

        if config["torsionScanInfo"]["preserveBackboneTorsions"]:
            cmapTerms = create_cmap_terms(config)
            for cmapTerm in cmapTerms:
                outRtf.write(f"{cmapTerm}\n")
        
        for donorPair in donorAcceptors["DONORS"]:
            outRtf.write(f"DONOR {donorPair[0]} {donorPair[1]}\n")
        for acceptorPair in donorAcceptors["ACCEPTORS"]:
            outRtf.write(f"ACCEPTOR {acceptorPair[0]} {acceptorPair[1]}\n")
        outRtf.write("END\n")
            
        ##TODO: CMAP for Phi/Psi torsion angles if possible


def create_cmap_terms(config: dict) -> dict:
    backboneAliases = config["moleculeInfo"]["backboneAliases"]
    cmapTerms = []
    for nAtom, caAtom, cAtom in zip(backboneAliases["N"], backboneAliases["CA"], backboneAliases["C"]):
        cmapTerms.append("CMAP -C "+nAtom+" "+caAtom+" "+cAtom+ " N CA C +N")
    return cmapTerms
    



def get_donor_acceptors(config: dict, debug:bool = True) -> dict:
    """
    Looks though optimised solvated conformers for hydrogen bonds donors and acceptors
    Creates a dict of DONOR and ACCEPTOR atom pairs,
    where the first entry is the atom participating in the Hydrogen Bond, 
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
    print(solvatedOptimisedXyzs)
    config, solvatedDf = pad_pdb_with_waters(config)

    allDonors = []
    allAcceptors = []

    ## run serial
    if debug:
        for solvatedOptXyz in solvatedOptimisedXyzs:
            donors, acceptors = detect_hydrogen_bonds_construct_donor_acceptors(solvatedOptXyz, solvatedDf, moleculeName)
            allDonors.extend(donors)
            allAcceptors.extend(acceptors)
    ## run in parallel
    else:
        purpleText = "\033[35m"
        resetTextColor = "\033[0m"
        ## set up progress bar
        tqdmBarOptions = {
            "desc": f"{purpleText}Assigning Donor / Acceptor Pairs {resetTextColor}",
            "ascii": "-🗲→",    
            "colour": "yellow",
            "unit":  "conformer",
            "dynamic_ncols": True,
        }
        ## construct an arglist to pass to multiprocessing
        argsList = [(solvatedOptXyz, solvatedDf, moleculeName) for solvatedOptXyz in solvatedOptimisedXyzs]
        nCpus = min(len(argsList), config["miscInfo"]["availableCpus"])
        with mp.Pool(processes=nCpus) as pool:
            results = pool.starmap(detect_hydrogen_bonds_construct_donor_acceptors, argsList)
        for donors, acceptors in results:
            allDonors.extend(donors)
            allAcceptors.extend(acceptors)


    uniqueDonors = list(set(allDonors))
    uniqueAcceptors = list(set(allAcceptors))

    donorAcceptors = {"DONORS": uniqueDonors,
                      "ACCEPTORS": uniqueAcceptors}
    
    print(donorAcceptors)
        
    config["runtimeInfo"]["madeByCreator"]["donorAcceptors"] = donorAcceptors

    return config
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

def detect_hydrogen_bonds_construct_donor_acceptors(solvatedOptXyz: FilePath,
                                                     solvatedDf: pd.DataFrame,
                                                       moleculeName: str) -> tuple:
    xyzDf = file_parsers.xyz2df(solvatedOptXyz)
    hBonds = detect_hydrogen_bonds(xyzDf)
    donors, acceptors = construct_donor_acceptors(hBonds, solvatedDf, moleculeName)

    return donors, acceptors



def construct_donor_acceptors(hBonds: list, pdbDf: pd.DataFrame, moleculeName: str) -> dict:
    print(moleculeName)
    exit()
    donors = []
    acceptors = []
    for hBond in hBonds:
        donorDf = pdbDf.iloc[hBond[0]]
        hydrogenDf = pdbDf.iloc[hBond[1]]
        acceptorDf = pdbDf.iloc[hBond[2]]
        print(donorDf["RES_NAME"], acceptorDf["RES_NAME"])
        if donorDf["RES_NAME"] == moleculeName:
            donorPair = (hydrogenDf["ATOM_NAME"], donorDf["ATOM_NAME"])
            donors.append(donorPair)
        if acceptorDf["RES_NAME"] == moleculeName:
            bondedAtoms = Capping_Assistant.find_bonded_atoms(pdbDf, acceptorDf["ATOM_NAME"])
            ## exclude tri-valent (or more!) Nitrogen atoms
            if acceptorDf["ELEMENT"] == "N" and len(bondedAtoms) > 2:
                continue
            bondedAtom = choose_bonded_atom_for_acceptor(bondedAtoms, pdbDf, moleculeName)
            acceptorPair = (acceptorDf["ATOM_NAME"], bondedAtom)
            acceptors.append(acceptorPair)

    donors = list(set(donors))
    acceptors = list(set(acceptors))
            
    return donors, acceptors


# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def choose_bonded_atom_for_acceptor(bondedAtoms, pdbDf, moleculeName):
    bondedAtoms = [atom for atom in bondedAtoms if not atom.startswith("H")]
    ## remove 
    bondedAtomDf = pdbDf[pdbDf["ATOM_NAME"].isin(bondedAtoms)]
    bondedAtomDf = bondedAtomDf[bondedAtomDf["RES_NAME"] == moleculeName]
    chosenBondedAtom = bondedAtomDf["ATOM_NAME"].to_list()[0]
    return chosenBondedAtom


# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def detect_hydrogen_bonds(xyzDf: pd.DataFrame, hydrogenBondDistanceCutoff: float = 3.5, hydrogenBondAngleCutoff: float = 120.0):
    """
    Detect hydrogen bonds based on geometric criteria.
    Criteria: D-H...A where D is donor (O, N), A is acceptor (O, N),
    H...A distance < distance_cutoff (Angstroms), D-H...A angle > angle_cutoff (degrees).
    Returns list of (donor_idx, hydrogen_idx, acceptor_idx) tuples.
    """
    hBonds = []
    donorAcceptorTypes = ["N", "O", "F"]  # Exclude C and N
    
    # Find potential donors (O, N) and hydrogens
    for indexAtomA, rowAtomA in  xyzDf.iterrows():
        coordsA = rowAtomA[["x", "y", "z"]].values
        if rowAtomA["element"]  in donorAcceptorTypes:  # Potential donor
            # Look for bonded hydrogens (simple distance-based check)
            for indexAtomB, rowAtomB in  xyzDf.iterrows():
                coordsB = rowAtomB[["x", "y", "z"]].values
                if rowAtomB["element"] == 'H':
                    protonBondLength = calculate_distance(coordsA, coordsB)
                    if 0.8 < protonBondLength < 1.2:  # Typical D-H bond length
                        # Found D-H pair, now look for acceptors
                        for indexAtomC, rowAtomC in  xyzDf.iterrows():
                            coordsC = rowAtomC[["x", "y", "z"]].values
                            if rowAtomC["element"]  in donorAcceptorTypes and indexAtomC != indexAtomA:  # Potential acceptor
                                hydrogenBondLength = calculate_distance(coordsB, coordsC)
                                if hydrogenBondLength < hydrogenBondDistanceCutoff:  # H...A distance check
                                    angle_dha = calculate_angle(coordsA, coordsB, coordsC)
                                    if angle_dha > hydrogenBondAngleCutoff:  # D-H...A angle check
                                        hBonds.append((indexAtomA, indexAtomB, indexAtomC))

    return hBonds




# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def calculate_distance(coord1, coord2):
    """Calculate Euclidean distance between two points."""
    return np.linalg.norm(coord1 - coord2)
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
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
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
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
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
