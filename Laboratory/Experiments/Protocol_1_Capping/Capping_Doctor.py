import os
from os import path as p
import numpy as np
import pandas as pd
from pdbUtils import pdbUtils

class FilePath:
    pass

class DirectoryPath:
    pass
from typing import List, Tuple
from . import Capping_Monster
from . import Capping_Assistant
from OperatingTools import Timer, cleaner

@Timer.time_function('Termini Capping', 'TERMINI_CAPPING')
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
    inputDir = config['pathInfo']['inputDir']
    outputDir = config['pathInfo']['outputDir']
    moleculeName = config['moleculeInfo']['moleculeName']
    molPdb = p.join(inputDir, f'{moleculeName}.pdb')
    config['runtimeInfo']['madeByCapping'] = {}
    cappingDir = p.join(outputDir, '01_termini_capping')
    os.makedirs(cappingDir, exist_ok=True)
    config['runtimeInfo']['madeByCapping']['cappingDir'] = cappingDir
    molDf = pdbUtils.pdb2df(molPdb)
    molDf = Capping_Monster.trim_termini(molDf, config)
    (nmePdb, acePdb) = Capping_Assistant.find_capping_pdbs()
    nmeCappedDf = add_nMethyl_caps(molDf, nmePdb, config)
    aceCappedDf = add_acetyl_caps(nmeCappedDf, acePdb, config)
    cappedDf = Capping_Assistant.reorder_atom_ids(aceCappedDf)
    cappedPdb = p.join(cappingDir, f'{moleculeName}_capped.pdb')
    pdbUtils.df2pdb(cappedDf, cappedPdb)
    optimedPdb = Capping_Monster.optimise_capped_structures(cappedPdb, config)
    config['runtimeInfo']['madeByCapping']['cappedPdb'] = optimedPdb
    config['checkpointInfo']['cappingComplete'] = True
    cleaner.clean_capping(config)
    return config

def add_n_methyl_caps(molDf: pd.DataFrame, nmePdb: FilePath, config: dict) -> pd.DataFrame:
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
    cTerminalAtoms = config['moleculeInfo']['cTermini']
    nmeDf = pdbUtils.pdb2df(nmePdb)
    orginalNmeDf = nmeDf.copy()
    dfsToConcat = [molDf]
    maxResId = molDf['RES_ID'].max()
    for (i, cTerminalAtom) in enumerate(cTerminalAtoms):
        tmpNmeDf = Capping_Monster.place_NN(cappedDf, cTerminalAtom, nmeDf)
        tmpNmeDf = Capping_Monster.place_HNN1(cappedDf, cTerminalAtom, tmpNmeDf)
        tmpNmeDf = Capping_Monster.place_CN(cappedDf, cTerminalAtom, tmpNmeDf)
        tmpNmeDf = Capping_Assistant.transform_whole(orginalNmeDf, tmpNmeDf, atomNames=['NN', 'HNN1', 'CN'])
        tmpNmeDf['RES_ID'] = maxResId + i + 1
        dfsToConcat.append(tmpNmeDf)
    nmeCappedDf = pd.concat(dfsToConcat, axis=0)
    return nmeCappedDf

def add_acetyl_caps(molDf: pd.DataFrame, acePdb: FilePath, config: dict) -> pd.DataFrame:
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
    nTerminalAtoms = config['moleculeInfo']['nTermini']
    aceDf = pdbUtils.pdb2df(acePdb)
    dfsToConcat = [molDf]
    maxResId = molDf['RES_ID'].max()
    for (i, nTerminalAtom) in enumerate(nTerminalAtoms):
        tmpAceDf = Capping_Monster.place_CC1(cappedDf, nTerminalAtom, aceDf)
        tmpAceDf = Capping_Monster.place_OC(cappedDf, nTerminalAtom, aceDf)
        tmpAceDf = Capping_Monster.place_CC2(cappedDf, nTerminalAtom, aceDf)
        aceDf = pdbUtils.pdb2df(acePdb)
        tmpAceDf = Capping_Assistant.transform_whole(aceDf, tmpAceDf, atomNames=['CC1', 'OC', 'CC2'])
        tmpAceDf['RES_ID'] = maxResId + i + 1
        dfsToConcat.append(tmpAceDf)
    aceCappedDf = pd.concat(dfsToConcat, axis=0)
    return aceCappedDf
if __name__ == '__main__':
    raise NotImplementedError