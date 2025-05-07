## BASIC LIBRARIES ##
from os import path as p
import pandas as pd
import numpy as np

## PARMED LIBRARIES ##
import parmed
from parmed.charmm import CharmmParameterSet

## drFRANKENSTEIN LIBRARIES ##
from .. import Stitching_Assistant
from . import CHARMM_helper_functions

## CLEAN CODE CLASSES ##
from typing import Tuple
class FilePath:
    pass
class DirectoryPath:
    pass


# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def get_MM_torsion_energies(config: dict, torsionTag: str, debug: bool = False) -> Tuple[dict, dict]:    
    """
    Gets MM[torsion] energy for each torsion we have scanned
    This is done by extracting torsion parameters from FRCMOD file
    Then reconstructing energy as a sum of the cosine functions described

    Args:
        config (dict): config containing all run information
        torsionTag (str): tag for torsion we are interested in

    Returns:
        mmTorsionEnergies (dict): energies 
        mmCosineComponents (dict): cosine components
    """
    if debug:
        print(f"Getting MM torsion energies for torsion {torsionTag}")
    ## get torsion parameters from PRM file
    mmTorsionParameters = extract_torsion_parameters_from_prm(config, torsionTag)
    ## reconstruct torsion energies from parameters
    mmTorsionEnergies, mmCosineComponents = construct_MM_torsion_energies(mmTorsionParameters)

    return  mmTorsionEnergies, mmCosineComponents

def extract_torsion_parameters_from_prm(config: dict, torsionTag: str) -> dict:
    """
    Extracts torsion parameters from PRM file

    Args:
        config (dict): config containing all run information
        torsionTag (str): tag for torsion we are interested in

    Returns:
        torsionParameters (dict): torsion parameters
    """
    ## unpack config
    moleculePrm = config["runtimeInfo"]["madeByStitching"]["moleculePrm"]
    moleculeRtf = config["runtimeInfo"]["madeByStitching"]["moleculeRtf"]
    rotatableDihedrals = config["runtimeInfo"]["madeByTwisting"]["rotatableDihedrals"]
    
    ## get atom types from previously made dict
    targetAtomTypes = tuple(rotatableDihedrals[torsionTag]["ATOM_TYPES"])

    ## load PRM and RTF in to Parmed
    parmedPrm = CharmmParameterSet(moleculePrm, moleculeRtf)

    ## find old dihedral types
    for prmDihedralTypes in parmedPrm.dihedral_types:
        if prmDihedralTypes == targetAtomTypes or prmDihedralTypes == targetAtomTypes[::-1]:
            diheralParams = parmedPrm.dihedral_types[prmDihedralTypes]
            break


    torsionParameters = []
    for cosineTerm in diheralParams:
        torsionParameters.append({
            "k": cosineTerm.phi_k,
            "period": cosineTerm.per,
            "phase": cosineTerm.phase
        })

    return torsionParameters


def construct_MM_torsion_energies(mmTorsionParameters) -> dict:
    """
    Constructs MM energies from mmTorsionParameters using:

        E(torsion) = K  * (1 + cos(periodicity * (Angle - Phase)))
    https://ambermd.org/FileFormats.php#frcmod

    Args:
        mmTorsionParameters (dict): dict containing torsion parameters for each torsion

    Returns:
        mmTorsionEnergies (dict): dict containing MM energies for each torsion
    """

    ## init angle 
    angle = np.radians(np.arange(0, 360, 10, dtype=float))
    ## init empty array
    mmTorsionEnergy = np.zeros_like(angle)
    ## loop through terms for torsion parameter
    mmCosineComponents = {}
    for parameter in mmTorsionParameters:
        ## extract params from dict
        potentialConstant = float(parameter["k"])
        periodicityNumber = abs(float(parameter["period"]))
        phase = np.radians(float(parameter["phase"]))
        ## construct cosine component
        cosineComponent: np.array = potentialConstant  * (1 + np.cos(periodicityNumber * angle - phase)) 
        ## add to torsion energy
        mmTorsionEnergy += cosineComponent
        mmCosineComponents[periodicityNumber] = cosineComponent
        
    return mmTorsionEnergy, mmCosineComponents
