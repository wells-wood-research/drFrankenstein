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
def get_MM_torsion_energies(config: dict, torsionTag: str, paramFile: FilePath, debug: bool = False) -> Tuple[dict, dict]:    
    """Compute MM torsion energies and cosine components for a CHARMM scan."""
    if debug:
        print(f"Getting MM torsion energies for torsion {torsionTag}")
    ## get torsion parameters from PRM file
    mmTorsionParameters = extract_torsion_parameters_from_prm(paramFile, torsionTag, config)
    ## reconstruct torsion energies from parameters
    mmTorsionEnergies, mmCosineComponents = construct_MM_torsion_energies(mmTorsionParameters)

    return  mmTorsionEnergies, mmCosineComponents

def extract_torsion_parameters_from_prm(moleculePrm: FilePath, torsionTag: str, config: dict) -> dict:
    """Extract torsion parameters from the CHARMM PRM/RTF pair."""
    ## unpack config
    moleculeRtf = config["runtimeInfo"]["madeByStitching"]["moleculeRtf"]
    torsionsToScan = config["runtimeInfo"]["madeByTwisting"]["torsionsToScan"]
    
    ## get atom types from previously made dict
    targetAtomTypes = tuple(torsionsToScan[torsionTag]["ATOM_TYPES"])

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
    """Construct CHARMM torsion energies from parameter dictionaries."""

    ## init angle 
    angle = np.radians(np.arange(0, 360, 10, dtype=float))
    ## init empty array
    mmTorsionEnergy = np.zeros_like(angle)
    ## loop through terms for torsion parameter
    mmCosineComponents = {}
    for termIndex, parameter in enumerate(mmTorsionParameters, start=1):
        ## extract params from dict
        potentialConstant = float(parameter["k"])
        periodicityNumber = abs(float(parameter["period"]))
        phase = np.radians(float(parameter["phase"]))
        ## construct cosine component
        cosineComponent: np.array = potentialConstant  * (1 + np.cos(periodicityNumber * angle - phase)) 
        ## add to torsion energy
        mmTorsionEnergy += cosineComponent
        componentKey = f"term_{termIndex}_n{periodicityNumber:g}"
        mmCosineComponents[componentKey] = cosineComponent
        
    return mmTorsionEnergy, mmCosineComponents
