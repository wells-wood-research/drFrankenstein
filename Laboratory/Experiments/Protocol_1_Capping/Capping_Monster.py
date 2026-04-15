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

## drFRANKENSTIEN LIBRARIES ##
from . import Capping_Assistant
from OperatingTools import drOrca
from OperatingTools import file_parsers

def optimise_capped_structures(cappedPdb: FilePath, config:dict) -> FilePath:
    """
    Runs a quick ORCA XTB2 optimisation of the capped structure
    This avoids odd distances between the capping groups and the rest of the molecule

    Args:
        cappedPdb (FilePath): path to the capped PDB file
        config (dict): drFrankenstein config
    Returns:
        optimisedPdb (FilePath): path to the capped PDB file
    """

    ## unpack config
    cappingDir = config["runtimeInfo"]["madeByCapping"]["cappingDir"]
    optDir = p.join(cappingDir, "geometry_optimisation")
    moleculeName = config["moleculeInfo"]["moleculeName"]
    os.makedirs(optDir, exist_ok=True)

    ## convert PDB to XYZ
    xyzFile = p.join(optDir, f"{moleculeName}_capped.xyz")
    file_parsers.pdb2xyz(cappedPdb, xyzFile)


    ## make an ORCA input file for optimisation
    optOrcaInput: FilePath = drOrca.make_orca_input_for_opt(inputXyz=xyzFile,
                                                   outDir = optDir,
                                                   moleculeInfo=config["moleculeInfo"],
                                                   qmMethod="XTB2",
                                                    solvationMethod="ALPB(water)",
                                                    geomOptions="MaxIter 100")
    optOrcaOutput: FilePath = p.join(optDir, "orca_opt.out")
    drOrca.run_orca(optOrcaInput, optOrcaOutput, config)

    optimisedXyz = p.join(optDir, "orca_opt.xyz")
    optimisedPdb = p.join(optDir, f"{moleculeName}_capped_opt.pdb")
    file_parsers.update_pdb_coords(cappedPdb, optimisedXyz, optimisedPdb)

    return optimisedPdb

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
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
    nTerminalAtoms = config["moleculeInfo"]["backboneAliases"]["N"]
    cTerminalAtoms = config["moleculeInfo"]["backboneAliases"]["C"]
    ## get the proton (or equivalent) atoms bound to N-Terminus
    allAtomsToDelete = []
    for nTerminalAtom in nTerminalAtoms:
        bondedAtoms = Capping_Assistant.find_bonded_atoms(molDf, nTerminalAtom)
        nTerminalProtons = Capping_Assistant.decide_atom_to_delete_N_termini(bondedAtoms, config)
        allAtomsToDelete.extend(nTerminalProtons)

    ## get the oxygen (or equivalent) atoms bound to C-Terminus, and anything bound to that oxygen
    for cTerminalAtom in cTerminalAtoms:
        bondedAtoms = Capping_Assistant.find_bonded_atoms(molDf, cTerminalAtom)
        cTerminalOxygen = Capping_Assistant.decide_atom_to_delete_C_termini(bondedAtoms)
        if cTerminalOxygen != "None":
            allAtomsToDelete.append(cTerminalOxygen)
        if not cTerminalOxygen == "None":
            bondedAtoms = Capping_Assistant.find_bonded_atoms(molDf, cTerminalOxygen)
            cTerminalProton = Capping_Assistant.decide_atom_to_delete_C_terminal_proton(bondedAtoms)
            if cTerminalProton != "None":
                allAtomsToDelete.append(cTerminalProton)

    molDf = molDf[~molDf["ATOM_NAME"].isin(allAtomsToDelete)]

    return molDf

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def place_cc(cappedDf: pd.DataFrame,
               nTerminalAtom: str,
                 aceDf: pd.DataFrame) -> pd.DataFrame:
    """
    Places the C_C capping atom in the Acetyl capping group
    
    Args:
        cappedDf (pd.DataFrame): DataFrame of capped molecule
        nTerminalAtom (str): name of N-Terminus atom
        aceDf (pd.DataFrame): DataFrame of Acetyl capping group
    Returns:
        aceDf (pd.DataFrame): DataFrame with C_C capping atom placed
    """
    
    ## get coordinates for N-Terminus
    nTerminalCoords = cappedDf[cappedDf["ATOM_NAME"] == nTerminalAtom].loc[:, ["X", "Y", "Z"]].values[0]
    ##  find atoms bound to N-Terminus
    bondedAtoms = Capping_Assistant.find_bonded_atoms(cappedDf, nTerminalAtom)
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
    ## place the C_C capping atom
    aceDf = Capping_Assistant.place_capping_atom_by_average_of_vectors(capDf = aceDf,
                                                      vectorAtomCoordsA = hCoords,
                                                       vectorAtomCoordsB = nTerminalCoords,
                                                        vectorAtomCoordsC = caCoords,
                                                          vectorAtomCoordsD =  nTerminalCoords,
                                                             bondedAtomCoords = nTerminalCoords,
                                                              atomNameToPlace = "C_C",
                                                               bondLength = 1.4)

    return aceDf
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲𗲗🗲🗲
def place_oc(cappedDf: pd.DataFrame,
              nTerminalAtom: str,
                aceDf: pd.DataFrame) -> pd.DataFrame:
    """
    Places the O_C capping atom in the Acetyl capping group
    
    Args:
        cappedDf (pd.DataFrame): DataFrame of capped molecule
        nTerminalAtom (str): name of N-Terminus atom
        aceDf (pd.DataFrame): DataFrame of Acetyl capping group
    Returns:
        aceDf (pd.DataFrame): DataFrame with O_C capping atom placed
    """

    ## get coordinates for N-Terminus
    nTerminalCoords = cappedDf[cappedDf["ATOM_NAME"] == nTerminalAtom].loc[:, ["X", "Y", "Z"]].values[0]
    ## find atoms bound to N-Terminus
    bondedAtoms = Capping_Assistant.find_bonded_atoms(cappedDf, nTerminalAtom)
    ## find coordinates for proton bound to N-Terminus
    hAtomNames = [atom for atom in bondedAtoms if atom.startswith("H")]
    if len(hAtomNames)  == 0:
        hAtomName = [atom for atom in bondedAtoms if atom.startswith("C") and not atom == "CA"][0]
    else:
        hAtomName = hAtomNames[0]
    hCoords = cappedDf[cappedDf["ATOM_NAME"] == hAtomName].loc[:, ["X", "Y", "Z"]].values[0]
    ## get coordinates for CC1 in the Acetyl capping group
    ccCoords = aceDf[aceDf["ATOM_NAME"] == "C_C"].loc[:, ["X", "Y", "Z"]].values[0]
    ## place the OC capping atom
    aceDf = Capping_Assistant.place_capping_atom_by_vector(capDf=aceDf,
                                        vectorAtomCoordsA = hCoords,
                                        vectorAtomCoordsB = nTerminalCoords,
                                        bondedAtomCoords = ccCoords,
                                        atomNameToPlace = "O_C",
                                        bondLength = 1.2)

    return aceDf
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def place_c2c(cappedDf: pd.DataFrame,
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
    ccCoords = aceDf[aceDf["ATOM_NAME"] == "C_C"].loc[:, ["X", "Y", "Z"]].values[0]
    ocCoords = aceDf[aceDf["ATOM_NAME"] == "O_C"].loc[:, ["X", "Y", "Z"]].values[0]
    ## place the CC2 capping atom
    aceDf = Capping_Assistant.place_capping_atom_by_average_of_vectors(capDf = aceDf,
                                                      vectorAtomCoordsA = nTerminalCoords,
                                                       vectorAtomCoordsB = ccCoords,
                                                        vectorAtomCoordsC = ocCoords,
                                                         vectorAtomCoordsD = ccCoords, 
                                                          bondedAtomCoords = ccCoords,
                                                           atomNameToPlace = "C2_C",
                                                            bondLength= 1.45)
    return aceDf

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def place_nn(cappedDf: pd.DataFrame,
              cTerminalAtom: str,
                nmeDf: pd.DataFrame) -> pd.DataFrame:
    """
    Places the N_N capping atom in the N-Methyl capping group
    
    Args:
        cappedDf (pd.DataFrame): DataFrame of capped molecule
        cTerminalAtom (str): name of C-Terminus atom
        nmeDf (pd.DataFrame): DataFrame of N-Methyl capping group
    Returns:
        nmeDf (pd.DataFrame): DataFrame with NN capping atom placed
    """    
    # find atoms bonded to C-Terminus
    bondedAtoms = Capping_Assistant.find_bonded_atoms(cappedDf, cTerminalAtom)
    ## get coordinates for CA and O atoms, or best equivalent
    caEquivalentName = [atomName for atomName in bondedAtoms if atomName.startswith("C")][0]
    oEquivalentName = [atomName for atomName in bondedAtoms if atomName.startswith("O")][0]

    cTerminalCoords = cappedDf[cappedDf["ATOM_NAME"] == cTerminalAtom].loc[:, ["X", "Y", "Z"]].values[0]
    caCoords = cappedDf[cappedDf["ATOM_NAME"] == caEquivalentName].loc[:, ["X", "Y", "Z"]].values[0]
    oCoords = cappedDf[cappedDf["ATOM_NAME"] == oEquivalentName].loc[:, ["X", "Y", "Z"]].values[0]
    ## place NN
    nmeDf = Capping_Assistant.place_capping_atom_by_average_of_vectors(capDf=nmeDf,
                                                     vectorAtomCoordsA=oCoords,
                                                      vectorAtomCoordsB=cTerminalCoords,
                                                       vectorAtomCoordsC=caCoords,
                                                        vectorAtomCoordsD=cTerminalCoords,
                                                         bondedAtomCoords=cTerminalCoords,
                                                          atomNameToPlace="N_N",
                                                           bondLength=1.4)

    
    return nmeDf
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def place_hn(cappedDf: pd.DataFrame,
                cTerminalAtom: str,
                  nmeDf: pd.DataFrame) -> pd.DataFrame:
    """
    Places the H_N capping atom in the N-Methyl capping group
    
    Args:
        cappedDf (pd.DataFrame): DataFrame of capped molecule
        cTerminalAtom (str): name of C-Terminus atom
        nmeDf (pd.DataFrame): DataFrame of N-Methyl capping group
    Returns:
        nmeDf (pd.DataFrame): DataFrame with H_N capping atom placed
    """
    print(nmeDf)
    # Get coordinates from protein
    bondedAtoms = Capping_Assistant.find_bonded_atoms(cappedDf, cTerminalAtom)

    oEquivalentName = [atomName for atomName in bondedAtoms if atomName.startswith("O")][0]
    
    cTerminalCoords = cappedDf[cappedDf["ATOM_NAME"] == cTerminalAtom].loc[:, ["X", "Y", "Z"]].values[0]
    oCoords = cappedDf[cappedDf["ATOM_NAME"] == oEquivalentName].loc[:, ["X", "Y", "Z"]].values[0]
    
    # Get original NN coordinates
    nnCoords = nmeDf[nmeDf["ATOM_NAME"] == "N_N"].loc[:, ["X", "Y", "Z"]].values[0]

    nmeDf = Capping_Assistant.place_capping_atom_by_vector(capDf=nmeDf,
                                        vectorAtomCoordsA = oCoords,
                                        vectorAtomCoordsB = cTerminalCoords,
                                        bondedAtomCoords = nnCoords,
                                        atomNameToPlace = "H_N",
                                        bondLength = 1.0)

    return nmeDf
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def place_cn(cappedDf: pd.DataFrame,
              cTerminalAtom: str,
                nmeDf: pd.DataFrame) -> pd.DataFrame:
    """
    Places the C_N capping atom in the N-Methyl capping group
    
    Args:
        cappedDf (pd.DataFrame): DataFrame of capped molecule
        cTerminalAtom (str): name of C-Terminus atom
        nmeDf (pd.DataFrame): DataFrame of N-Methyl capping group
    Returns:
        nmeDf (pd.DataFrame): DataFrame with CN capping atom placed
    """
    ## get coordinates for C-Terminus, NN and HNN atoms
    cTerminalCoords = cappedDf[cappedDf["ATOM_NAME"] == cTerminalAtom].loc[:, ["X", "Y", "Z"]].values[0]
    nnCoords = nmeDf[nmeDf["ATOM_NAME"] == "N_N"].loc[:, ["X", "Y", "Z"]].values[0]
    hnCoords = nmeDf[nmeDf["ATOM_NAME"] == "H_N"].loc[:, ["X", "Y", "Z"]].values[0]
    ## place CN
    nmeDf = Capping_Assistant.place_capping_atom_by_average_of_vectors(capDf = nmeDf,
                                             vectorAtomCoordsA = cTerminalCoords,
                                               vectorAtomCoordsB = nnCoords,
                                                 vectorAtomCoordsC = hnCoords,
                                                   vectorAtomCoordsD = nnCoords,
                                                     bondedAtomCoords = nnCoords,
                                                       atomNameToPlace = "C_N",
                                                         bondLength = 1.4)
    
    return nmeDf
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
