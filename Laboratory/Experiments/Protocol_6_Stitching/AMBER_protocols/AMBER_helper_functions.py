## BASIC LIBRARIES ##
import os
from os import path as p
from subprocess import call, PIPE
import pandas as pd
from shutil import move, copy
import numpy as np
import parmed
## drFRANKENSTEIN LIBRARIES ##
from .. import Stitching_Assistant

## CLEAN CODE CLASSES ##
from typing import Tuple
class FilePath:
    pass
class DirectoryPath:
    pass
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def construct_MM_torsion_energies(mmTorsionParameters) -> dict:
    """
    Constructs MM energies from mmTorsionParameters using:

        E(torsion) = (K / Mult) * (1 + cos(periodicity * Angle - Phase))
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
        # inverseDivisionFactor = float(parameter["multiplicity"])
        periodicityNumber = abs(float(parameter["periodicity"]))
        phase = np.radians(float(parameter["phase"]))
        ## construct cosine component
        cosineComponent: np.array = (potentialConstant) * (1 + np.cos(periodicityNumber * angle - phase)) 
        ## add to torsion energy
        mmTorsionEnergy += cosineComponent
        mmCosineComponents[periodicityNumber] = cosineComponent
        
    return mmTorsionEnergy, mmCosineComponents

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def update_frcmod(moleculeFrcmod: FilePath,
                  config: dict,
                  torsionTag: str,
                  torsionParamDf: pd.DataFrame,
                  shuffleIndex: int) -> dict:
    """
    Uses parmed to update torsion parameters in a frcmod file based on a DataFrame.
    """
    torsionsToScan = config["runtimeInfo"]["madeByTwisting"]["torsionsToScan"]
    moleculeParameterDir = config["runtimeInfo"]["madeByStitching"]["moleculeParameterDir"]
    moleculeName = config["moleculeInfo"]["moleculeName"]

    targetAtomTypes = tuple(torsionsToScan[torsionTag]["ATOM_TYPES"])

    # Prepare new dihedral types from DataFrame
    new_dihedral_types = [
        parmed.DihedralType(
            phi_k=row['Amplitude'],
            per=int(row['Period']), # Ensure periodicity is integer
            phase=row['Phase'],
            scee=1.2, # Default AMBER scaling
            scnb=2.0  # Default AMBER scaling
        )
        for _, row in torsionParamDf.iterrows()
    ]

    # Load frcmod as parameter set
    parmedFrcmod = parmed.load_file(moleculeFrcmod)
   # 3. Find and replace the matching dihedral entry
    keysToDelete = []
    for frcmodAtomNames in parmedFrcmod.dihedral_types:
        if frcmodAtomNames == targetAtomTypes or frcmodAtomNames == targetAtomTypes[::-1]:
            keysToDelete.append(frcmodAtomNames)
    
    for key_to_delete in keysToDelete:
        del parmedFrcmod.dihedral_types[key_to_delete]

   
    parmedFrcmod.dihedral_types[targetAtomTypes] = new_dihedral_types

    # Define output path and write the updated frcmod file
    outputFrcmod = p.join(moleculeParameterDir, f"{moleculeName}_{shuffleIndex}.frcmod")
    parmedFrcmod.write(outputFrcmod, style='frcmod')


    return outputFrcmod


def edit_mol2_partial_charges(config: dict) -> None:
    """
    Gets partial charges stored in a dataframe and pastes them into the
    charges column of a MOL2 file

    Args:
        config (dict): config with all run information

    Returns:
        None [path to updated mol2 file already in config]
    
    """

    ## unpack config ##
    moleculeMol2 = config["runtimeInfo"]["madeByStitching"]["moleculeMol2"]
    chargesCsv = config["runtimeInfo"]["madeByCharges"]["chargesCsv"]

    outMol2 = p.splitext(moleculeMol2)[0] + "_tmp.mol2"
    chargesDf = pd.read_csv(chargesCsv)

    ## init atom index counter
    atomIndex: int = 0
    ## open inMol2 for reading and outMol2 for writing
    with open(moleculeMol2, 'r') as inMol2, open(outMol2, "w") as writeMol2:
        ## read through inMol2 until we get to the atom section
        mol2Lines = inMol2.readlines()
        isAtomLine = False
        for line in mol2Lines:
            if line.strip() == "@<TRIPOS>ATOM":
                isAtomLine = True
            elif line.strip() == "@<TRIPOS>BOND":
                isAtomLine = False
            elif isAtomLine:
                ## increment atom index
                atomIndex += 1
                ## get atom charge for this index
                atomCharge = chargesDf.loc[chargesDf['atomIndex'] == atomIndex, 'Charge'].values[0]
                ## Format to 4 decimal places
                atomCharge = f"{atomCharge:.4f}"
                ## sort out spaces for the case of negative charges
                if not atomCharge.startswith("-"):
                    atomCharge = " "+atomCharge
                newLine = line[:-10]  + atomCharge
                line = newLine+"\n"
            ## write to outMol2
            writeMol2.write(line)

    ## overwrite inMol2 with outMol2    
    move(outMol2, moleculeMol2)
    return None
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

def copy_assembled_parameters(config: dict) -> dict:
    """
    Copies the FRCMOD created during the assembly stage
    to the stitching directory

    Args:
        config (dict): contains all information needed for run

    Returns:
        config (dict): updated config

    """
    ## unpack config
    assembledFrcmod = config["runtimeInfo"]["madeByAssembly"]["assembledFrcmod"]
    assembledPrmtop = config["runtimeInfo"]["madeByAssembly"]["assembledPrmtop"]
    cappedMol2 = config["runtimeInfo"]["madeByAssembly"]["cappedMol2"]
    moleculeParameterDir = config["runtimeInfo"]["madeByStitching"]["moleculeParameterDir"]
    moleculeName = config["moleculeInfo"]["moleculeName"]


    destFrcmod = p.join(moleculeParameterDir, f"{moleculeName}_capped.frcmod")
    copy(assembledFrcmod, destFrcmod)

    destPrmtop = p.join(moleculeParameterDir, f"{moleculeName}_capped.prmtop")
    copy(assembledPrmtop, destPrmtop)

    destMol2 = p.join(moleculeParameterDir, f"{moleculeName}_capped.mol2")
    copy(cappedMol2, destMol2)

    config["runtimeInfo"]["madeByStitching"]["moleculeFrcmod"] = destFrcmod
    config["runtimeInfo"]["madeByStitching"]["moleculePrmtop"] = destPrmtop
    config["runtimeInfo"]["madeByStitching"]["moleculeMol2"] = destMol2
    return config

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def run_tleap_to_make_params(moleculeFrcmod: FilePath, config: dict) -> None:
    """
    Uses TLEAP to create a prmtop inpcrd pair

    Args:
        inMol2 (FilePath): input MOL2 file
        molFrcmod (FilePath): input FRCMOD file
        outDir (DirectoryPath): output directory
        index (str): identifier for output files

    Returns: 
        prmtop (FilePath): topology file for AMBER
        impcrd (FilePath): coordinate file for AMBER
    """
    ## unpack config
    inMol2 = config["runtimeInfo"]["madeByStitching"]["moleculeMol2"]
    outDir = config["runtimeInfo"]["madeByStitching"]["moleculeParameterDir"]
    moleculeName = config["moleculeInfo"]["moleculeName"]


    prmtop: FilePath = p.join(outDir, f"{moleculeName}_capped.prmtop")
    inpcrd: FilePath = p.join(outDir, f"{moleculeName}_capped.inpcrd")

    tleapInput: FilePath = p.join(outDir, f"leap_{moleculeName}.in")
    with open(tleapInput, "w") as f:
        f.write("source leaprc.gaff2\n")
        f.write(f"mol  = loadmol2 {inMol2} \n")
        f.write(f"loadamberparams {moleculeFrcmod} \n") # use frcmod previously made
        # this frcmod will need to be updated after each torsion fit, so the next batch of mol2 are parameterised with the updated file
        f.write(f"saveamberparm mol {prmtop} {inpcrd}  \n")
        f.write("quit")

    tleapOutput: FilePath = p.join(outDir, f"tleap.out")

    tleapCommand: list = ["tleap", "-f", tleapInput, ">", tleapOutput]

    call(tleapCommand, stdout=PIPE)

    return prmtop