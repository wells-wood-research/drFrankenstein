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

from typing import List, Tuple, Optional, Dict # Added for type hinting

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


def create_final_rtf(config: dict) -> None:
    """
    Creates the final RTF (Residue Topology File) for the uncapped molecule,
    adding BONDs for termini and CMAP terms if specified.

    Args:
        config: Configuration dictionary containing paths and molecule information.
    """
    ## unpack config
    finalCreationDir: DirectoryPath = config["runtimeInfo"]["madeByCreator"]["finalCreationDir"] # type: ignore
    moleculeRtf: FilePath = config["runtimeInfo"]["madeByStitching"]["moleculeRtf"] # type: ignore
    
    moleculeName: str = config["moleculeInfo"]["moleculeName"]
    nTerminalAtoms: List[str] = config["moleculeInfo"]["nTermini"]
    cTerminalAtoms: List[str] = config["moleculeInfo"]["cTermini"]
    donorAcceptors: Dict[str, List[Tuple[str, str]]] = config["runtimeInfo"]["madeByCreator"]["donorAcceptors"]
    finalRtf: FilePath = p.join(finalCreationDir, f"{moleculeName}.rtf") # type: ignore

    cappingAtomNames = ["CN", "NN", "HNN1", "HCN1", "HCN2", "HCN3", "CC1", "OC", "CC2", "HC1", "HC2", "HC3"]

    terminalSectionWritten = False
    with open(moleculeRtf, "r") as inRtf, open(finalRtf, "w") as outRtf: # type: ignore
        for line in inRtf.readlines():
            if line.startswith("END"): 
                continue
            if any(atomName in line for atomName in cappingAtomNames):
                continue
            if line.startswith("BOND") and not terminalSectionWritten:
                for atomName in cTerminalAtoms:
                    outRtf.write(f"BOND {atomName} +N\n")
                for atomName in nTerminalAtoms:
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


def create_cmap_terms(config: dict) -> List[str]:
    """Creates CMAP term strings based on backbone atom aliases."""
    backboneAliases: Dict[str, List[str]] = config["moleculeInfo"]["backboneAliases"]
    cmapTerms: List[str] = []
    # Ensure all lists have the same length before zipping
    min_len = min(len(backboneAliases.get("N", [])), len(backboneAliases.get("CA", [])), len(backboneAliases.get("C", [])))
    
    for i in range(min_len):
        nAtom = backboneAliases["N"][i]
        caAtom = backboneAliases["CA"][i]
        cAtom = backboneAliases["C"][i]
        cmapTerms.append(f"CMAP -C {nAtom} {caAtom} {cAtom} N CA C +N")
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
    solvatedOptimisedXyzs: List[FilePath] = config["runtimeInfo"]["madeByCharges"]["solvatedOptimisedXyzs"]
    config, solvatedDf = pad_pdb_with_waters(config)

    allDonors: List[Tuple[str, str]] = []
    allAcceptors: List[Tuple[str, str]] = []

    ## run serial
    if debug:
        for solvatedOptXyz_item in solvatedOptimisedXyzs: # Renamed loop variable
            donors, acceptors = detect_hydrogen_bonds_construct_donor_acceptors(solvated_opt_xyz=solvatedOptXyz_item, solvated_df=solvatedDf, molecule_name=moleculeName)
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
        argsList = [(solvatedOptXyz, solvatedDf, moleculeName) for solvatedOptXyz in solvatedOptimisedXyzs] # solvatedOptXyz is fine here as it's local to list comprehension
        nCpus = min(len(argsList), config["miscInfo"]["availableCpus"])
        with mp.Pool(processes=nCpus) as pool:
            results = pool.starmap(detect_hydrogen_bonds_construct_donor_acceptors, argsList)
        for donors_list, acceptors_list in results: # Renamed to avoid conflict
            allDonors.extend(donors_list)
            allAcceptors.extend(acceptors_list)


    uniqueDonors = list(set(allDonors))
    uniqueAcceptors = list(set(allAcceptors))

    donorAcceptorsDict: Dict[str, List[Tuple[str, str]]] = {"DONORS": uniqueDonors, # Renamed for clarity
                                                           "ACCEPTORS": uniqueAcceptors}
    
        
    config["runtimeInfo"]["madeByCreator"]["donorAcceptors"] = donorAcceptorsDict

    return config
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

def detect_hydrogen_bonds_construct_donor_acceptors(solvated_opt_xyz: FilePath,
                                                     solvated_df: pd.DataFrame,
                                                       molecule_name: str) -> Tuple[List[Tuple[str, str]], List[Tuple[str, str]]]:
    xyzDf = file_parsers.xyz2df(solvated_opt_xyz)
    hBonds = detect_hydrogen_bonds(xyz_df=xyzDf)
    donors, acceptors = construct_donor_acceptors(h_bonds=hBonds, pdb_df=solvated_df, molecule_name=molecule_name)

    return donors, acceptors



def construct_donor_acceptors(h_bonds: list, pdb_df: pd.DataFrame, molecule_name: str) -> Tuple[List[Tuple[str, str]], List[Tuple[str, str]]]:
    donors: List[Tuple[str,str]] = []
    acceptors: List[Tuple[str,str]] = []
    for hBond_tuple in h_bonds: # Renamed loop variable
        donorDf = pdb_df.iloc[hBond_tuple[0]]
        hydrogenDf = pdb_df.iloc[hBond_tuple[1]]
        acceptorDf = pdb_df.iloc[hBond_tuple[2]]
        if donorDf["RES_NAME"] == molecule_name:
            donorPair = (str(hydrogenDf["ATOM_NAME"]), str(donorDf["ATOM_NAME"]))
            donors.append(donorPair)
        if acceptorDf["RES_NAME"] == molecule_name:
            bondedAtoms = Capping_Assistant.find_bonded_atoms(molDf=pdb_df, atomName=str(acceptorDf["ATOM_NAME"])) # type: ignore
            ## exclude tri-valent (or more!) Nitrogen atoms
            if acceptorDf["ELEMENT"] == "N" and len(bondedAtoms) > 2:
                continue
            bondedAtom = choose_bonded_atom_for_acceptor(bonded_atoms=bondedAtoms, pdb_df=pdb_df, molecule_name=molecule_name)
            acceptorPair = (str(acceptorDf["ATOM_NAME"]), bondedAtom)
            acceptors.append(acceptorPair)

    donors = list(set(donors))
    acceptors = list(set(acceptors))
            
    return donors, acceptors


# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def choose_bonded_atom_for_acceptor(bonded_atoms: List[str], pdb_df: pd.DataFrame, molecule_name: str) -> str:
    """Selects a non-hydrogen bonded atom for an acceptor from a list of bonded atoms."""
    bondedAtomsFiltered = [atom for atom in bonded_atoms if not atom.startswith("H")] # Renamed
    ## remove 
    bondedAtomDf = pdb_df[pdb_df["ATOM_NAME"].isin(bondedAtomsFiltered)]
    bondedAtomDf = bondedAtomDf[bondedAtomDf["RES_NAME"] == molecule_name]
    chosenBondedAtom = bondedAtomDf["ATOM_NAME"].to_list()[0]
    return str(chosenBondedAtom)


# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def detect_hydrogen_bonds(xyz_df: pd.DataFrame, hydrogen_bond_distance_cutoff: float = 3.5, hydrogen_bond_angle_cutoff: float = 120.0) -> List[Tuple[int, int, int]]:
    """
    Detect hydrogen bonds based on geometric criteria.
    Criteria: D-H...A where D is donor (O, N), A is acceptor (O, N),
    H...A distance < distance_cutoff (Angstroms), D-H...A angle > angle_cutoff (degrees).
    Returns list of (donor_idx, hydrogen_idx, acceptor_idx) tuples.
    """
    hBonds: List[Tuple[int, int, int]] = []
    donorAcceptorTypes = ["N", "O", "F"]  # Exclude C and N
    
    # Find potential donors (O, N) and hydrogens
    for indexAtomA, rowAtomA in  xyz_df.iterrows():
        coordsA: np.ndarray = rowAtomA[["x", "y", "z"]].values
        if rowAtomA["element"]  in donorAcceptorTypes:  # Potential donor
            # Look for bonded hydrogens (simple distance-based check)
            for indexAtomB, rowAtomB in  xyz_df.iterrows():
                coordsB: np.ndarray = rowAtomB[["x", "y", "z"]].values
                if rowAtomB["element"] == 'H':
                    protonBondLength = calculate_distance(coord1=coordsA, coord2=coordsB)
                    if 0.8 < protonBondLength < 1.2:  # Typical D-H bond length
                        # Found D-H pair, now look for acceptors
                        for indexAtomC, rowAtomC in  xyz_df.iterrows():
                            coordsC: np.ndarray = rowAtomC[["x", "y", "z"]].values
                            if rowAtomC["element"]  in donorAcceptorTypes and indexAtomC != indexAtomA:  # Potential acceptor
                                hydrogenBondLength = calculate_distance(coord1=coordsB, coord2=coordsC)
                                if hydrogenBondLength < hydrogen_bond_distance_cutoff:  # H...A distance check
                                    angleDha = calculate_angle(coord1=coordsA, coord2=coordsB, coord3=coordsC) # Renamed
                                    if angleDha > hydrogen_bond_angle_cutoff:  # D-H...A angle check
                                        hBonds.append((int(indexAtomA), int(indexAtomB), int(indexAtomC))) # Ensure int

    return hBonds




# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def calculate_distance(coord1: np.ndarray, coord2: np.ndarray) -> float:
    """Calculate Euclidean distance between two points."""
    return float(np.linalg.norm(coord1 - coord2))
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def calculate_angle(coord1: np.ndarray, coord2: np.ndarray, coord3: np.ndarray) -> float:
    """Calculate angle (in degrees) between three points (coord2 is vertex)."""
    vec1: np.ndarray = coord1 - coord2
    vec2: np.ndarray = coord3 - coord2
    
    norm1: float = np.linalg.norm(vec1)
    norm2: float = np.linalg.norm(vec2)
    
    # Check for zero-length vectors
    if norm1 == 0 or norm2 == 0:
        return 0.0  # Default to 0 degrees for undefined cases
    
    vec1Norm: np.ndarray = vec1 / norm1 # Renamed
    vec2Norm: np.ndarray = vec2 / norm2 # Renamed
    cosTheta: float = np.clip(np.dot(vec1Norm, vec2Norm), -1.0, 1.0) # Renamed
    return float(np.degrees(np.arccos(cosTheta)))
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def pad_pdb_with_waters(config: dict) -> Tuple[dict, pd.DataFrame]:
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
