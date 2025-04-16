## BASIC LIBRARIES ##
import os
from os import path as p
from subprocess import call, PIPE
import pandas as pd
from shutil import move

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
        inverseDivisionFactor = float(parameter["multiplicity"])
        periodicityNumber = abs(float(parameter["periodicity"]))
        phase = np.radians(float(parameter["phase"]))
        ## construct cosine component
        cosineComponent: np.array = (potentialConstant / inverseDivisionFactor) * (1 + np.cos(periodicityNumber * angle - phase)) 
        ## add to torsion energy
        mmTorsionEnergy += cosineComponent
        mmCosineComponents[periodicityNumber] = cosineComponent
        
    return mmTorsionEnergy, mmCosineComponents

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def update_frcmod(config: dict,
                   torsionTag: str,
                     torsionParamDf: pd.DataFrame) -> dict:
    """
    Updates an AMBER frcmod file with a torsion parameter block
    
    Args:
        config (dict): the drFrankenstein config containing all run information
        torsionTag (str): the torsion tag for the torsion we are updating
        torsionParamDf (pd.DataFrame): the torsion parameters for the torsion we are updating
    Returns:    
        config (dict): updated config
    """
    ## unpack config ##
    inFrcmod = config["runtimeInfo"]["madeByStitching"]["moleculeFrcmod"]
    tmpFrcmod = p.splitext(inFrcmod)[0] + "_tmp.frcmod"
    atomTypeMap = config["runtimeInfo"]["madeByStitching"]["atomTypeMap"]
    
    ## construct torsion identifier forwards and reversed for finding line in frcmod ##
    atomNames = torsionTag.split("-")
    atomTypes = [atomTypeMap[name] for name in atomNames]
    atomTypesReversed = atomTypes[::-1]
    torsionIdentifier = "-".join([el + ' ' if len(el) == 1 else el for el in atomTypes])
    torsionIdentifierReversed = "-".join([el + ' ' if len(el) == 1 else el for el in atomTypesReversed])

    ## init a bool to determine whether to write lines to frcmod
    writeLines = False
    ## init new torsion block as an empty string
    newTorsionBlock = ""

    ##TODO: multiplicity from symmetry
    multiplicity = 1

    ## read through torsionParamDf and write new torsion block to add to frcmod
    for _, row in torsionParamDf.iloc[:-1].iterrows():
        amplitude = f"{row['Amplitude']:.3f}"
        phase = f"{row['Phase']:.3f}"
        period = f"{-row['Period']:.3f}"
        newTorsionBlock += (f"{torsionIdentifier:<14}{multiplicity:<5}"
                            f"{amplitude:>5}{phase:>14}{period:>16}"
                            f"\t\tMADE BY drFRANKENSTEIN\n")
    ## write last row TODO: (why is this separate????)
    lastRow = torsionParamDf.iloc[-1]
    amplitude = f"{lastRow['Amplitude']:.3f}"
    phase = f"{lastRow['Phase']:.3f}"
    period = f"{lastRow['Period']:.3f}"
    newTorsionBlock += (f"{torsionIdentifier:<14}{multiplicity:<5}"
                        f"{amplitude:>5}{phase:>14}{period:>16}"
                        f"\t\tMADE BY drFRANKENSTEIN\n")
    
    ## init a bool to determine whether to write lines to frcmod
    writeLines = False
    ## read through frcmod and write to tmp frcmod
    with open(inFrcmod, "r") as f, open(tmpFrcmod, "w") as tmp:
        for line in f:
            ## dont copy params for this torsion
            if  line.startswith(torsionIdentifier) or line.startswith(torsionIdentifierReversed):
                continue
            elif line.startswith("DIHE"):
                writeLines = True
                tmp.write(line)
            elif writeLines:
                tmp.write(newTorsionBlock)
                writeLines = False
                tmp.write(line)
            else:
                tmp.write(line)
    ## owerwrite frcmod
    move(tmpFrcmod, inFrcmod)
    return config
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
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
    mmTorsionDir: DirectoryPath = config["runtimeInfo"]["madeByStitching"]["mmTorsionCalculationDir"]
    ## get capped pdb file as input
    cappedPdb: FilePath = config["runtimeInfo"]["madeByCapping"]["cappedPdb"]
    ## get molecule name from filepath
    moleculeName: str = config["moleculeInfo"]["moleculeName"]

    ## convert capped PDB to MOL2
    cappedMol2: FilePath = p.join(mmTorsionDir, f"{moleculeName}_capped.mol2")
    Stitching_Assistant.pdb2mol2(cappedPdb, cappedMol2, mmTorsionDir)
    config["runtimeInfo"]["madeByStitching"]["cappedMol2"] = cappedMol2

    ## read calculated partial charges
    chargesCsv = config["runtimeInfo"]["madeByCharges"]["chargesCsv"]
    chargesDf = pd.read_csv(chargesCsv, index_col="Unnamed: 0")
    chargesDf["Charge"] = chargesDf["Charge"].round(4)

    ## paste charges from prior QM calculations into MOL2 file
    chargesMol2 = p.join(config["runtimeInfo"]["madeByStitching"]["mmTorsionCalculationDir"], f"{moleculeName}_charges.mol2")
    Stitching_Assistant.edit_mol2_partial_charges(cappedMol2, chargesDf, chargesMol2)

    ## fix atom types in MOL2 file
    fixedAtomMol2 = p.join(mmTorsionDir, f"{moleculeName}_fixed_atoms.mol2")
    Stitching_Assistant.edit_mo2_atom_types(chargesMol2, fixedAtomMol2)
    config["runtimeInfo"]["madeByStitching"]["finalMol2"] = fixedAtomMol2

    ## get a map of atomName -> atomType
    atomTypeMap = Stitching_Assistant.create_atom_type_map(fixedAtomMol2)
    config["runtimeInfo"]["madeByStitching"]["atomTypeMap"] = atomTypeMap
    ## create a FRCMOD file from MOL2 file
    molFrcmod = Stitching_Assistant.create_frcmod_file(fixedAtomMol2, moleculeName, config)
    ## update config 
    config["runtimeInfo"]["madeByStitching"]["moleculeFrcmod"] = molFrcmod
    return config

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

def make_prmtop_inpcrd(trajXyz, fittingRoundDir, cappedPdb, chargesDf, moleculeFrcmod, debug=False):

    trajIndex = trajXyz.split(".")[1]

    os.makedirs(fittingRoundDir, exist_ok=True)
    os.chdir(fittingRoundDir)

    trajPdb = p.join(fittingRoundDir, f"orca_{trajIndex}.pdb")
    Stitching_Assistant.update_pdb_coords(cappedPdb, trajXyz, trajPdb)

    trajMol2 = p.join(fittingRoundDir, f"orca_{trajIndex}.mol2")
    Stitching_Assistant.pdb2mol2(trajPdb, trajMol2, fittingRoundDir)

    chargedMol2 = p.join(fittingRoundDir, f"orca_{trajIndex}_charged.mol2")
    Stitching_Assistant.edit_mol2_partial_charges(trajMol2, chargesDf, chargedMol2)

    renamedMol2 = p.join(fittingRoundDir, f"orca_{trajIndex}_charged_renamed.mol2")
    Stitching_Assistant.edit_mo2_atom_types(chargedMol2, renamedMol2)
    prmtop, inpcrd = run_tleap_to_make_params(renamedMol2, moleculeFrcmod, fittingRoundDir, trajIndex)

    return prmtop, inpcrd
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def parse_frcmod(molFrcmod: FilePath) -> dict:
    """
    Reads through a FRCMOD file and returns a dict with the contents of the file

    Args:
        molFrcmod (FilePath): frcmod file

    Returns:
        parsedFrcmod (dict): dict with the contents of the frcmod file
    
    """
    ## init a bunch of bools
    readingMass: bool       = False
    readingBonds: bool      = False
    readingAngles: bool     = False
    readingDihedrals: bool  = False
    readingImpropers: bool  = False
    readingNonbonded: bool  = False

    ## init empty parsedFrcmod dict
    parsedFrcmod = {"ATOMS": [], "BONDS": [], "ANGLES": [], "DIHEDRALS": [], "IMPROPERS": [], "NONBONDED": []}
    ## open molFrcmod for reading
    with open(molFrcmod, 'r') as f:
        ## loop through lines
        for line in f:
            ## skip empty lines
            if line.strip() == "":
                continue
            ## use bools to determine which section of the frcmod file we are in
            elif line.startswith("MASS"):
                readingMass = True
            elif line.startswith("BOND"):
                readingBonds = True
                readingMass = False
            elif line.startswith("ANGLE"):
                readingAngles = True
                readingBonds = False
            elif line.startswith("DIHE"):
                readingDihedrals = True
                readingAngles = False
            elif line.startswith("IMPROPER"):
                readingImpropers = True
                readingDihedrals = False
            elif line.startswith("NONBON"):
                readingNonbonded = True
                readingImpropers = False
            ## if we are in a data section, parse the line
            else:
                ## process mass data
                if readingMass:
                    lineData = line.split()
                    lineParsed = {"atomType": lineData[0], "mass": lineData[1], "vdw-radius": lineData[2]}
                    parsedFrcmod["ATOMS"].append(lineParsed)
                ## process bond data
                elif readingBonds:
                    atomData = get_frcmod_atom_data(line, [[0, 2], [4, 6]])
                    paramData = "".join(line[6:]).split()[0:2]
                    lineParsed = {"atoms": atomData, "r0": paramData[0], "k": paramData[1]}
                    parsedFrcmod["BONDS"].append(lineParsed)
                ## process angle data
                elif readingAngles:
                    atomData = get_frcmod_atom_data(line, [[0, 2], [3, 5], [6, 8]])
                    paramData = "".join(line[9:]).split()[0:2]
                    lineParsed = {"atoms": atomData, "theta0": paramData[0], "k": paramData[1]}
                    parsedFrcmod["ANGLES"].append(lineParsed)
                ## process dihedral data
                elif readingDihedrals:
                    atomData = get_frcmod_atom_data(line, [[0, 2], [3, 5], [6, 8], [9, 11]])
                    paramData = "".join(line[12:]).split()[0:4]
                    lineParsed = {"atoms": atomData, "multiplicity": paramData[0], "k": paramData[1], "phase": paramData[2], "periodicity": paramData[3]}  
                    parsedFrcmod["DIHEDRALS"].append(lineParsed)  
                ## process improper data
                elif readingImpropers:
                    atomData = get_frcmod_atom_data(line, [[0, 2], [3, 5], [6, 8], [9, 11]])
                    paramData = "".join(line[12:]).split()[0:3]
                    lineParsed = {"atoms": atomData, "k": paramData[0], "phi0": paramData[1], "periodicity": paramData[2]}
                    parsedFrcmod["IMPROPERS"].append(lineParsed)
                ## process nonbonded data
                elif readingNonbonded:
                    atomData = get_frcmod_atom_data(line, [[0, 6]])
                    paramData = "".join(line[6:]).split()[0:2]
                    lineParsed = {"atoms": atomData, "vdw-radius": paramData[0], "well-depth": paramData[1]}
                    parsedFrcmod["NONBONDED"].append(lineParsed)
    return parsedFrcmod

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def get_frcmod_atom_data(line: str, indexes: list) -> list:
    """
    Small function for getting atom names from FRCMOD file lines

    Args:
        line (str): line from a FRCMOD file
        indexes (List[int,int]): location of atom types in FRCMOD line

    Returns:
        atomData (List[str]): list of atom types 
    """

    return [line[start:end].strip() for start, end in indexes]

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def run_tleap_to_make_params(inMol2: FilePath,
                            molFrcmod: FilePath,
                              outDir: DirectoryPath,
                                index: str):
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

    prmtop: FilePath = p.join(outDir, f"orca_{index}.prmtop")
    inpcrd: FilePath = p.join(outDir, f"orca_{index}.inpcrd")

    tleapInput: FilePath = p.join(outDir, f"leap_{index}.in")
    with open(tleapInput, "w") as f:
        f.write("source leaprc.gaff2\n")
        f.write(f"mol  = loadmol2 {inMol2} \n")
        f.write(f"loadamberparams {molFrcmod} \n") # use frcmod previously made
        # this frcmod will need to be updated after each torsion fit, so the next batch of mol2 are parameterised with the updated file
        f.write(f"saveamberparm mol {prmtop} {inpcrd}  \n")
        f.write("quit")

    tleapOutput: FilePath = p.join(outDir, f"tleap_{index}.out")

    tleapCommand: list = ["tleap", "-f", tleapInput, ">", tleapOutput]

    call(tleapCommand, stdout=PIPE)

    return prmtop, inpcrd