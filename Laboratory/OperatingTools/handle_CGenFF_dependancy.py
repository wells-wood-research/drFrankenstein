import os
from os import path as p

from OperatingTools import drSplash
from OperatingTools import file_parsers
from typing import Dict, Any # Added for type hinting

# Placeholder classes (extend if needed)
class FilePath:
    pass
class DirectoryPath:
    pass

def handle_cgenff_dependency(config: dict) -> dict: # Renamed function
    """
    Handles the CGenFF dependency by checking for a local CGenFF executable
    or an existing stream file. If neither is found, it guides the user
    to generate the stream file using the CGenFF webserver.

    Args:
        config: Configuration dictionary.

    Returns:
        The potentially modified configuration dictionary.
    """
    if config["pathInfo"]["cgenffExe"] is None: # Changed == None to is None
        inputDir: DirectoryPath = config["pathInfo"]["inputDir"] # type: ignore
        moleculeName: str = config["moleculeInfo"]["moleculeName"]
        moleculeStr: FilePath = p.join(inputDir, f"{moleculeName}_capped.str") # type: ignore
        if p.isfile(moleculeStr):
            config["runtimeInfo"]["madeByCapping"]["cappedMoleculeStr"] = moleculeStr
            return config
        else:

            cappedPdb: FilePath = config["runtimeInfo"]["madeByCapping"]["cappedPdb"] # type: ignore
            cappingDir: DirectoryPath = config["runtimeInfo"]["madeByCapping"]["cappingDir"] # type: ignore
            cappedMol2: FilePath = p.join(cappingDir, f"{moleculeName}_capped.mol2") # type: ignore
            config["runtimeInfo"]["madeByCapping"]["cappedMol2"] = cappedMol2
            file_parsers.pdb2mol2(pdb_file=cappedPdb, mol2_file=cappedMol2)
            rename_molecule_mol2(mol2_file=cappedMol2, molecule_name=moleculeName)

            drSplash.show_need_cgenff_str(capped_mol2=cappedMol2)
            exit(0)
    else:  
        return config


def rename_molecule_mol2(mol2_file: FilePath, molecule_name: str) -> None:
    """
    Renames the molecule within a MOL2 file. 
    The first line after @<TRIPOS>MOLECULE is replaced with molecule_name.

    Args:
        mol2_file: Path to the MOL2 file to be modified.
        molecule_name: The new name for the molecule.
    """
    tmpMol2: FilePath = p.join(p.dirname(mol2_file), f"tmp.mol2") # type: ignore
    with open(mol2_file, "r") as inMol2, open(tmpMol2, "w") as outMol2: # type: ignore
        moleculeNextLine = False
        for line in inMol2:
            if line.startswith("@<TRIPOS>MOLECULE"):
                outMol2.write(line)
                moleculeNextLine = True
            elif moleculeNextLine:  
                outMol2.write(f"{molecule_name}\n")
                moleculeNextLine = False
            else:
                outMol2.write(line)
    os.replace(tmpMol2, mol2_file) # type: ignore