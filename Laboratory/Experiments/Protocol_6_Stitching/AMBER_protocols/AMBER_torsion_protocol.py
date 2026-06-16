## BASIC LIBRARIES ##
from os import path as p
import pandas as pd
import parmed
## drFRANKENSTEIN LIBRARIES ##
from . import AMBER_helper_functions

## CLEAN CODE CLASSES ##
from typing import Tuple
class FilePath:
    pass
class DirectoryPath:
    pass

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def get_MM_torsion_energies(config: dict, torsionTag: str, moleculeFrcmod) -> Tuple[dict, dict]:
    """Compute MM torsion energies and cosine components for one torsion."""

    ## get torsion parameters from FRCMOD file
    mmTorsionParameters = extract_torsion_parameters_from_frcmod(moleculeFrcmod, torsionTag, config)

    ## reconstruct torsion energies from parameters
    mmTorsionEnergies, mmCosineComponents = AMBER_helper_functions.construct_MM_torsion_energies(mmTorsionParameters)

    return  mmTorsionEnergies, mmCosineComponents

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def extract_torsion_parameters_from_prmtop(config: dict, torsionTag: str) -> dict:
    """Extract torsion parameters from the PRMTOP structure."""

    molPrmtop = config["runtimeInfo"]["madeByStitching"]["moleculePrmtop"]
    parmedPrmtop  = parmed.load_file(molPrmtop, structure=True)

    torsionsToScan = config["runtimeInfo"]["madeByTwisting"]["torsionsToScan"]
    atomNames = tuple(torsionsToScan[torsionTag]["ATOM_NAMES"])

    torsionParameters = []
    for dihedral in parmedPrmtop.dihedrals:
        paramAtoms = (dihedral.atom1.name, dihedral.atom2.name, dihedral.atom3.name, dihedral.atom4.name)
        if paramAtoms == atomNames or paramAtoms[::-1] == atomNames:
            params = parse_torsion_params(dihedral.type)
            torsionParameters.append(params)


    return torsionParameters
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def extract_torsion_parameters_from_frcmod( moleculeFrcmod: FilePath,  torsionTag: str, config: dict,) -> dict:
    """Extract torsion parameters from the FRCMOD file."""
    torsionsToScan = config["runtimeInfo"]["madeByTwisting"]["torsionsToScan"]

    molParams = parmed.load_file(moleculeFrcmod)

    atomTypes = tuple(torsionsToScan[torsionTag]["ATOM_TYPES"])
    torsionParameters = []
    # Replace DIHEDRAL parameters
    for dihedralKey in molParams.dihedral_types:
        if dihedralKey == atomTypes or dihedralKey == atomTypes[::-1]:
            for cosineTerm in molParams.dihedral_types[dihedralKey]:
                params = parse_torsion_params(cosineTerm)
                if not params in torsionParameters:
                    torsionParameters.append(params)
    return torsionParameters
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def parse_torsion_params(dihedralType: parmed.topologyobjects.Dihedral) -> dict:
    """Convert a ParmEd dihedral object into a torsion parameter dict."""
    params = {
        "k": dihedralType.phi_k,
        "periodicity": dihedralType.per,
        "phase": dihedralType.phase,
    }
    return params
    
