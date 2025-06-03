from os import path as p
import pandas as pd
import parmed
from . import AMBER_helper_functions
from typing import Tuple

class FilePath:
    pass

class DirectoryPath:
    pass

def get_mm_torsion_energies(config: dict, torsionTag: str) -> Tuple[dict, dict]:
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
    mmTorsionParameters = extract_torsion_parameters_from_frcmod(config, torsionTag)
    (mmTorsionEnergies, mmCosineComponents) = AMBER_helper_functions.construct_MM_torsion_energies(mmTorsionParameters)
    return (mmTorsionEnergies, mmCosineComponents)

def extract_torsion_parameters_from_prmtop(config: dict, torsionTag: str) -> dict:
    """
    Reads through a frcmod file 
    Finds torsion parameters for each torsion that we have scanned
    Returns a dict with the torsion tag as the key and the torsion parameters as the value

    Args:
        molFrcmod (FilePath): frcmod file
        atomTypeMap (dict): dict mapping atom names to atom types
        config (dict): config dict

    Returns:
        mmTorsionParameters (dict): dict with the torsion tag as the key and the torsion parameters as the value
    """
    molPrmtop = config['runtimeInfo']['madeByStitching']['moleculePrmtop']
    parmedPrmtop = parmed.load_file(molPrmtop, structure=True)
    rotatableDihedrals = config['runtimeInfo']['madeByTwisting']['rotatableDihedrals']
    atomNames = tuple(rotatableDihedrals[torsionTag]['ATOM_NAMES'])
    torsionParameters = []
    for dihedral in parmedPrmtop.dihedrals:
        paramAtoms = (dihedral.atom1.name, dihedral.atom2.name, dihedral.atom3.name, dihedral.atom4.name)
        if paramAtoms == atomNames or paramAtoms[::-1] == atomNames:
            params = parse_torsion_params(dihedral.type)
            torsionParameters.append(params)
    return torsionParameters

def extract_torsion_parameters_from_frcmod(config: dict, torsionTag: str) -> dict:
    """
    Reads through a frcmod file 
    Finds torsion parameters for each torsion that we have scanned
    Returns a dict with the torsion tag as the key and the torsion parameters as the value

    Args:
        molFrcmod (FilePath): frcmod file
        atomTypeMap (dict): dict mapping atom names to atom types
        config (dict): config dict

    Returns:
        mmTorsionParameters (dict): dict with the torsion tag as the key and the torsion parameters as the value
    """
    rotatableDihedrals = config['runtimeInfo']['madeByTwisting']['rotatableDihedrals']
    molFrcmod = config['runtimeInfo']['madeByStitching']['moleculeFrcmod']
    molParams = parmed.load_file(molFrcmod)
    atomTypes = tuple(rotatableDihedrals[torsionTag]['ATOM_TYPES'])
    torsionParameters = []
    for dihedralKey in molParams.dihedral_types:
        if dihedralKey == atomTypes or dihedralKey == atomTypes[::-1]:
            for cosineTerm in molParams.dihedral_types[dihedralKey]:
                params = parse_torsion_params(cosineTerm)
                if not params in torsionParameters:
                    torsionParameters.append(params)
    return torsionParameters

def parse_torsion_params(dihedralType: parmed.topologyobjects.Dihedral) -> dict:
    params = {'k': dihedralType.phi_k, 'periodicity': dihedralType.per, 'phase': dihedralType.phase}
    return params