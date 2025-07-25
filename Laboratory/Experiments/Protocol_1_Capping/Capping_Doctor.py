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
from . import Capping_Monster
from . import Capping_Assistant
from OperatingTools import Timer, cleaner
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
@Timer.time_function("Termini Capping", "TERMINI_CAPPING")
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
    inputDir = config["pathInfo"]["inputDir"]
    outputDir = config["pathInfo"]["outputDir"]
    moleculeName = config["moleculeInfo"]["moleculeName"]
    molPdb = p.join(inputDir, f"{moleculeName}.pdb")
    backboneAliases = config["moleculeInfo"].get("backboneAliases", None)
    ## Skip for ligands ##
    if not backboneAliases:
        optimedPdb = Capping_Monster.optimise_capped_structures(molPdb, config)
        config["runtimeInfo"]["madeByCapping"]["cappedPdb"] = optimedPdb
        ## update config, this lets drMD know that capping is complete and lets it skip this step next time
        config["checkpointInfo"]["cappingComplete"] = True
        return config

    ## create an entry in runtimeInfo for capping
    config["runtimeInfo"]["madeByCapping"] = {}

    ## make a dir for capping, update config ##
    cappingDir = p.join(outputDir, "01_termini_capping")
    os.makedirs(cappingDir, exist_ok=True)
    config["runtimeInfo"]["madeByCapping"]["cappingDir"] = cappingDir

    ## load pdb file into a dataframe
    molDf = pdbUtils.pdb2df(molPdb)
    ## remove any extra atoms bound to the termini to make room for capping groups
    molDf = Capping_Monster.trim_termini(molDf, config)
    ## get capping PDBs from the drCapper directory
    nmePdb, acePdb = Capping_Assistant.find_capping_pdbs()
    ## add N-Methyl (NCH3) caps to each C-Terminus
    nmeCappedDf = add_nmethyl_caps(molDf, nmePdb, config)
    ## add Acetyl (COCH3) caps to each N-Terminus
    aceCappedDf = add_acetyl_caps(nmeCappedDf, acePdb, config)
    ## edit the ATOM_ID column to make it 1 -> N
    cappedDf = Capping_Assistant.reorder_atom_ids(aceCappedDf)
    ## write to PDB file, update config
    cappedPdb = p.join(cappingDir, f"{moleculeName}_capped.pdb")
    pdbUtils.df2pdb(cappedDf, cappedPdb)

    optimedPdb = Capping_Monster.optimise_capped_structures(cappedPdb, config)

    config["runtimeInfo"]["madeByCapping"]["cappedPdb"] = optimedPdb
    ## update config, this lets drMD know that capping is complete and lets it skip this step next time
    config["checkpointInfo"]["cappingComplete"] = True

    cleaner.clean_capping(config)


    return config


#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def add_nmethyl_caps(molDf: pd.DataFrame,
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
    cTerminalAtoms = config["moleculeInfo"]["backboneAliases"]["C"]
    nmeDf = pdbUtils.pdb2df(nmePdb)
    originalNmeDf = nmeDf.copy()

    dfsToConcat = [molDf]

    maxResId = molDf["RES_ID"].max()
    for i, cTerminalAtom in enumerate(cTerminalAtoms):
        tmpNmeDf = Capping_Monster.place_nn(cappedDf, cTerminalAtom, nmeDf)
        tmpNmeDf = Capping_Monster.place_hnn1(cappedDf, cTerminalAtom, tmpNmeDf)
        tmpNmeDf = Capping_Monster.place_cn(cappedDf, cTerminalAtom, tmpNmeDf)
        tmpNmeDf = Capping_Assistant.transform_whole(originalNmeDf, tmpNmeDf, atomNames=["NN", "HNN1", "CN"])
        tmpNmeDf["RES_ID"] = maxResId + i + 1
        dfsToConcat.append(tmpNmeDf)

    nmeCappedDf = pd.concat(dfsToConcat, axis=0)
    return nmeCappedDf
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
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
    nTerminalAtoms = config["moleculeInfo"]["backboneAliases"]["N"]
    aceDf = pdbUtils.pdb2df(acePdb)

    dfsToConcat = [molDf]
    maxResId = molDf["RES_ID"].max()
    for i, nTerminalAtom in enumerate(nTerminalAtoms):
        tmpAceDf = Capping_Monster.place_cc1(cappedDf, nTerminalAtom, aceDf)
        tmpAceDf = Capping_Monster.place_oc(cappedDf, nTerminalAtom, aceDf)
        tmpAceDf = Capping_Monster.place_cc2(cappedDf, nTerminalAtom, aceDf)
        aceDf = pdbUtils.pdb2df(acePdb)
        tmpAceDf = Capping_Assistant.transform_whole(aceDf, tmpAceDf, atomNames=["CC1", "OC", "CC2"])
        tmpAceDf["RES_ID"] = maxResId + i + 1
        dfsToConcat.append(tmpAceDf)

    aceCappedDf = pd.concat(dfsToConcat, axis=0)
    return aceCappedDf
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
if __name__ == "__main__":
    raise NotImplementedError