## BASIC IMPORTS ##
import os
from os import path as p
from subprocess import call, PIPE

## CLEAN CODE CLASSES ##
class FilePath:
    pass
class DirectoryPath:
    pass
from typing import List, Tuple

## drFRANKESTEN IMPORTS ##
from OperatingTools import drOrca
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def conformer_generation_protocol(config: dict) -> dict:
    """
    Runs conformer generation protocol using ORCA's GOAT program

    Args:
        config (dict): dictionary containing all information

    Returns:
        config (dict): updated config
    """
    print("--> GENERATING CONFORMERS WITH GOAT") #TODO: add progress bar and better SPLASH
    ## make new dirs and add to config ##
    config = sort_out_directories(config)

    ## unpack config ##
    cappedPdb = config["moleculeInfo"]["cappedPdb"]
    moleculeName = config["moleculeInfo"]["moleculeName"]
    conformerDir = config["pathInfo"]["conformerDir"]

    ## convert capped PDB file to XYZ format
    cappedXyz = p.join(conformerDir, f"{moleculeName}_capped.xyz")
    pdb_to_xyz(cappedPdb, cappedXyz)

    ## create an ORCA input file for GOAT
    goatOrcaInput = drOrca.write_goat_input(conformerDir, cappedXyz, config)
    goatOrcaOutput = p.join(conformerDir, f"GOAT_orca.out")
    if not p.isfile(goatOrcaOutput):
        drOrca.run_orca(goatOrcaInput, goatOrcaOutput, config)

    ## Split output XYZ file into multiple separate conformers
    conformerXyzs = split_conformers(conformerDir, moleculeName)
    config["pathInfo"]["conformerXyzs"] = conformerXyzs

    ## delete unnecessary files
    clean_up(conformerDir, config)

    ## update checkpoint flags
    config["checkpointInfo"]["conformersComplete"] = True

    return config
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def  sort_out_directories(config: dict) -> dict:
    """
    Creates a directory structure for conformers and updates config

    Args:
        config (dict): dictionary containing all information

    Returns:
        config (dict): updated config
    """
    outputDir = config["pathInfo"]["outputDir"]
    conformerDir = p.join(outputDir, "02_GOAT_conformers")
    os.makedirs(conformerDir, exist_ok=True) 
    config["pathInfo"]["conformerDir"] = conformerDir

    return config   

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

def split_conformers(conformerDir: DirectoryPath, moleculeName: str) -> List[FilePath]:
    """
    Separates the GOAT final ensemble into individual conformers in the XYZ format

    Args:
        conformerDir (DirectoryPath): directory containing conformers
        moleculeName (str): name of molecule

    Returns:
        conformerXyzs (List[FilePath]): list of conformer XYZ files
    
    """
    ## split GOAT output into individual conformers
    goatFinalXyz = p.join(conformerDir, "GOAT_orca.finalensemble.xyz")
    if not p.isfile(goatFinalXyz):
        raise FileNotFoundError(f"GOAT final ensemble not found in {conformerDir}")
    os.chmod(goatFinalXyz, 0o755)
    conformerXyzFilePattern = p.join(conformerDir, f"{moleculeName}_conformer_.xyz")
    obabelCommand = ["obabel", goatFinalXyz, "-O", conformerXyzFilePattern, "-m"]
    call(obabelCommand, stdout=PIPE, stderr=PIPE)

    conformerXyzs = [p.join(conformerDir, f) for f in os.listdir(conformerDir) if f.startswith(f"{moleculeName}_conformer_")]

    return conformerXyzs

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

def pdb_to_xyz(inPdb: FilePath, outXyz: FilePath) -> None:
    """
    Uses obabel to convert PDB to XYZ

    Args:
        inPdb (FilePath): input file in PDB format
        outXyz (FilePath): output file in XYZ format
    """
    obabelCommand = ["obabel", "-i", "pdb", inPdb, "-o", "xyz", "-O", outXyz]
    call(obabelCommand, stdout=PIPE, stderr=PIPE)
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

def clean_up(conformerDir: DirectoryPath, config: dict) -> None:
    """
    Removes unnecessary files from conformer directory

    Args:
        conformerDir (DirectoryPath): directory containing conformers
        config (dict): dictionary containing all information
    """
    if config["cleanUpLevel"] in ["basic", "full"]:
        filesToRemove = [p.join(conformerDir, f) for f in os.listdir(conformerDir) if f.startswith("GOAT")]
        for f in filesToRemove:
            os.remove(f)
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

if __name__ == "__main__":
    raise NotImplementedError