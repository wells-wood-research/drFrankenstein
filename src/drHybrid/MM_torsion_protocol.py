## BASIC LIBRARIES ##
from os import path as p
import pandas as pd

## drFRANKENSTEIN LIBRARIES ##
from drHybrid import shared_utils

## CLEAN CODE CLASSES ##
from typing import Tuple
class FilePath:
    pass
class DirectoryPath:
    pass
#######################
def get_generate_initial_frcmod(config: dict) -> dict:
    """
    Gets MM(torsion) energy for each torsion we have scanned
    This is done by extracting torsion parameters from FRCMOD file 
    And reconstructing energy as a sum of the cosine functions described

    Args:
        config (dict): config containing all run information

    Returns:
        mmTorsionEnergies (dict): energies
        config (dict): updated config dict 
    """

    ## get dir to write files to
    mmTorsionDir: DirectoryPath = config["pathInfo"]["mmTorsionCalculationDir"]
    ## get capped pdb file as input
    cappedPdb: FilePath = config["moleculeInfo"]["cappedPdb"]
    ## get molecule name from filepath
    moleculeName: str = "_".join(p.basename(config['moleculeInfo']['moleculePdb']).split(".")[0].split("_"))

    ## convert capped PDB to MOL2
    cappedMol2: FilePath = p.join(mmTorsionDir, f"{moleculeName}_capped.mol2")
    shared_utils.pdb2mol2(cappedPdb, cappedMol2, mmTorsionDir)
    config["pathInfo"]["cappedMol2"] = cappedMol2

    ## read calculated partial charges
    chargesCsv = config["chargeFittingInfo"]["chargesCsv"]
    chargesDf = pd.read_csv(chargesCsv, index_col="Unnamed: 0")
    chargesDf["Charge"] = chargesDf["Charge"].round(4)

    ## paste charges from prior QM calculations into MOL2 file
    chargesMol2 = p.join(config["pathInfo"]["mmTorsionCalculationDir"], f"{moleculeName}_charges.mol2")
    shared_utils.edit_mol2_partial_charges(cappedMol2, chargesDf, chargesMol2)

    ## fix atom types in MOL2 file
    fixedAtomMol2 = p.join(mmTorsionDir, f"{moleculeName}_fixed_atoms.mol2")
    shared_utils.edit_mo2_atom_types(chargesMol2, fixedAtomMol2)
    config["pathInfo"]["finalMol2"] = fixedAtomMol2

    ## get a map of atomName -> atomType
    atomTypeMap = shared_utils.create_atom_type_map(fixedAtomMol2)
    config["moleculeInfo"]["atomTypeMap"] = atomTypeMap
    ## create a FRCMOD file from MOL2 file
    molFrcmod = shared_utils.create_frcmod_file(fixedAtomMol2, moleculeName, config)
    ## update config 
    config["pathInfo"]["moleculeFrcmod"] = molFrcmod
    return config

#######################
def get_MM_torsion_energies(config: dict, torsionTag: str) -> Tuple[dict, dict]:    
    ## get torsion parameters from FRCMOD file
    mmTorsionParameters = shared_utils.extract_torsion_parameters(config, torsionTag)

    ## reconstruct torsion energies from parameters
    mmTorsionEnergies, mmCosineComponents = shared_utils.construct_MM_torsion_energies(mmTorsionParameters)

    return  mmTorsionEnergies, mmCosineComponents

#######################
