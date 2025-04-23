## BASIC LIBRARIES ##
from os import path as p
import pandas as pd
import numpy as np
from pdbUtils import pdbUtils
from shutil import move

## OPENMM LIBRARIES
import openmm.app as app
import openmm as openmm
import  openmm.unit  as unit

## CHARMM LIBRARIES ##
from psfgen import PsfGen

## CLEAN CODE CLASSES ##
class FilePath:
    pass
class DirectoryPath:
    pass


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


# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲##
def update_prm(config: dict, 
               torsionTag: str,
               torsionParamDf: pd.DataFrame) -> dict:
    """
    Updates a CHARMM prm file with a torsion parameter block
    
    Args:
        config (dict): the drFrankenstein config containing all run information
        torsionTag (str): the torsion tag for the torsion we are updating
        torsionParamDf (pd.DataFrame): the torsion parameters for the torsion we are updating
    Returns:    
        config (dict): updated config
    """

    ## unpack config
    inPrm  = config["runtimeInfo"]["madeByStitching"]["moleculePrm"]
    tmpPrm = p.splitext(inPrm)[0] + "_tmp.prm"
    atomTypeMap = config["runtimeInfo"]["madeByStitching"]["atomTypeMap"]

    ## construct torsion identifier forwards and reversed for finding line in frcmod ##
    atomNames = torsionTag.split("-")
    atomTypes = [atomTypeMap[name] for name in atomNames]
    atomTypesReversed = atomTypes[::-1]
    torsionIdentifier = f"{atomTypes[0]:<7}{atomTypes[1]:<7}{atomTypes[2]:<7}{atomTypes[3]:<7}"
    torsionIdentifierReversed = f"{atomTypes[3]:<7}{atomTypes[2]:<7}{atomTypes[1]:<7}{atomTypes[0]:<7}"

    ## init a bool to determine whether to write lines to prm file
    writeLines = False
    ## init new torsion block as an empty string
    newTorsionBlock = ""


    ## read through torsionParamDf and write new torsion block to add to prm
    for _, row in torsionParamDf.iterrows():
        amplitude = f"{row['Amplitude']:.4f}"
        phase = f"{row['Phase']:.2f}"
        period = f"{int(row['Period'])}"
        newTorsionBlock += (f"{torsionIdentifier:<32}{amplitude:<8}"
                            f"{period:<4}{phase:>6}"
                            f" !\t\tMADE BY drFRANKENSTEIN\n")
    ## init a bool to determine whether to write lines to frcmod
    writeLines = False
    ## read through frcmod and write to tmp frcmod
    with open(inPrm, "r") as f, open(tmpPrm, "w") as tmp:
        for line in f:
            ## dont copy params for this torsion
            if  line.startswith(torsionIdentifier) or line.startswith(torsionIdentifierReversed):
                continue
            elif line.startswith("DIHEDRALS"):
                writeLines = True
                tmp.write(line)
            elif writeLines:
                tmp.write(newTorsionBlock)
                writeLines = False
                tmp.write(line)
            else:
                tmp.write(line)
    ## owerwrite frcmod
    move(tmpPrm, inPrm)
    return config


# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲##
def parse_prm(moleculePrm: FilePath) -> dict:
    """
    Reads through a CHARMM 
    """
    ## init a bunch of bools
    readingBonds: bool      = False
    readingAngles: bool     = False
    readingDihedrals: bool  = False
    readingImpropers: bool  = False

    ## init empty parsedPrm dict
    parsedPrm = {"BONDS": [], "ANGLES": [], "DIHEDRALS": [], "IMPROPERS": []}

    with open(moleculePrm, "r") as f:
        for line in f.readlines():
            if line.strip() == "":
                continue
            if line.startswith("BONDS"):
                readingBonds = True
            elif line.startswith("ANGLES"):
                readingBonds = False
                readingAngles = True
            elif line.startswith("DIHEDRALS"):
                readingAngles = False
                readingDihedrals = True
            elif line.startswith("IMPROPERS"):
                readingDihedrals = False
                readingImpropers = True
            elif line.startswith("END"):
                readingImpropers = False
            else:
                lineData = line.split()
                if "!" in line:
                    comment = "!" + line.split("!")[1].strip()
                else: comment = ""
                if readingBonds:
                    lineParsed = {"atoms": [lineData[0], lineData[1]],
                                    "k": lineData[2], "r0": lineData[3],
                                    "comment": comment}
                    parsedPrm["BONDS"].append(lineParsed)
                elif readingAngles:
                    lineParsed = {"atoms": [lineData[0], lineData[1], lineData[2]],
                                    "k": lineData[3], "theta0": lineData[3],
                                    "comment": comment}
                    parsedPrm["ANGLES"].append(lineParsed)
                elif readingDihedrals:
                    lineParsed = {"atoms": [lineData[0], lineData[1], lineData[2], lineData[3]],
                                    "k": lineData[4], "period": lineData[5], "phase": lineData[6],
                                    "comment": comment}
                    parsedPrm["DIHEDRALS"].append(lineParsed)
                elif readingImpropers:
                    lineParsed = {"atoms": [lineData[0], lineData[1], lineData[2], lineData[3]],
                                    "k": lineData[4], "period": lineData[5], "phase": lineData[6],
                                    "comment": comment}
                    parsedPrm["IMPROPERS"].append(lineParsed)
    return parsedPrm


def write_prm(parsedPrm: dict, outPrm: FilePath) -> None:
    header = """
read param card flex append
* Parameters generated by analogy by
* CHARMM General Force Field (CGenFF) program version 4.0
* And modified by drFRANKENSTEIN
*

! Penalties lower than 10 indicate the analogy is fair; penalties between 10
! and 50 mean some basic validation is recommended; penalties higher than
! 50 indicate poor analogy and mandate extensive validation/optimization.\n\n
"""
    bondsSection = ""
    for param in parsedPrm["BONDS"]:
        bondsSection += f"{param['atoms'][0]:<7}{param['atoms'][1]:<8}{param['k']:<11}{param['r0']:<7} {param['comment']}\n"
    anglesSection = ""
    for param in parsedPrm["ANGLES"]:
        anglesSection += f"{param['atoms'][0]:<7}{param['atoms'][1]:<7}{param['atoms'][2]:<7}"
        anglesSection += f"{param['k']:>7}{param['theta0']:>10} {param['comment']}\n"
    dihedralSection = ""
    for param in parsedPrm["DIHEDRALS"]:
        dihedralSection += f"{param['atoms'][0]:<7}{param['atoms'][1]:<7}{param['atoms'][2]:<7}{param['atoms'][3]:<7}"
        dihedralSection += f"{param['k']:>10}{param['period']:>3}{param['phase']:>9} {param['comment']}\n"
    improperSection = ""
    for param in parsedPrm["IMPROPERS"]:
        improperSection += f"{param['atoms'][0]:<7}{param['atoms'][1]:<7}{param['atoms'][2]:<7}{param['atoms'][3]:<7}"
        improperSection += f"{param['k']:>10}{param['period']:>3}{param['phase']:>9} {param['comment']}\n"

    with open(outPrm, "w") as f:
        f.write(header)
        f.write("BONDS\n")
        f.write(bondsSection+ "\n")
        f.write("ANGLES\n")
        f.write(anglesSection+ "\n")
        f.write("DIHEDRALS\n")
        f.write(dihedralSection+ "\n")
        f.write("IMPROPERS\n")
        f.write(improperSection+ "\n")
        f.write("END\nRETURN")





# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲##
def create_initial_CHARMM_parameters(config: dict) -> dict:
    """
    Given a STR file (made by CGenFF)
    1. Splits STR -> PRM + RTF
    2. Edits RTF Charges 
    2. Uses PSFGen to create a PSF file

    Args: 
        config (dict): contains all run info
    Returns:
        config (dict): updated config 
    
    """
    ## find cgenff RTF and PRM files
    config = find_cgenff_params(config)
    ## split STR file into RTF and PRM for OpenMM
    config = split_charmm_str(config)
    ## update charges in RTF file 
    edit_rtf_charges(config)
    ## set capping groups to CHARMM defaults (NH1, H, C, O etc)
    ##TODO: WIP
    # set_capping_types_to_charmm_defaults(config)
    ## create a PSF file from PRM and RTF file
    config = make_charmm_psf(config)
    ## create an atom type map
    config = create_atom_type_map(config)

    return config

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲##
def set_capping_types_to_charmm_defaults(config: dict) -> None:
    """
    Sets the capping groups to CHARMM defaults for backbone atoms

    Args:
        config (dict): contains all run info
    Returns:
        None
    """
    mapCgenffTypeToCharmmType = set_capping_types_to_charmm_defaults_rtf(config)
    set_capping_types_to_charmm_defaults_prm(config, mapCgenffTypeToCharmmType)

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲##
def set_capping_types_to_charmm_defaults_prm(config: dict, mapCgenffTypeToCharmmType: dict) -> None:
    moleculePrm = config["runtimeInfo"]["madeByStitching"]["moleculePrm"]
    mmTorsionCalculationDir = config["runtimeInfo"]["madeByStitching"]["mmTorsionCalculationDir"]

    tmpPrm = p.join(mmTorsionCalculationDir, "tmp.prm")
    
    parsedPrm = parse_prm(moleculePrm)
    for cgenffType, charmmType in mapCgenffTypeToCharmmType.items():
        for paramEntry in ["BONDS", "ANGLES", "DIHEDRALS", "IMPROPERS"]:
            for param in parsedPrm[paramEntry]:
                try:
                    param["atoms"][param["atoms"].index(cgenffType)] = charmmType
                except ValueError:
                    pass
    write_prm(parsedPrm, tmpPrm)
    move(tmpPrm, moleculePrm)

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲##
def set_capping_types_to_charmm_defaults_rtf(config: dict) -> dict:
    moleculeRtf = config["runtimeInfo"]["madeByStitching"]["moleculeRtf"]
    mmTorsionCalculationDir = config["runtimeInfo"]["madeByStitching"]["mmTorsionCalculationDir"]

    mapAtomNameToDefaultCharmmType = {
        "NN": "NH1",        ## N-Me Nitrogen -> BB Nitrogen
        "CN": "CT1",        ## N-Me Carbon -> BB CA
        "HNN1": "H",        ## N-Me Hydrogen -> BB N-H Hydrogen
        "OC": "O",          ## ACE oxygen -> BB amide carbon
        "CC1": "C",         ## ACE amide carbon -> BB amide carbon
        "CC2": "CT1"        ## ACE methyl carbon -> BB CA
    }


    mapCgenffTypeToCharmmType = {}
        
    tmpRtf = p.join(mmTorsionCalculationDir, "tmp.rtf")
    with open(moleculeRtf, "r") as inRtf, open(tmpRtf, "w") as outRtf:
        for line in inRtf.readlines():
            if line.startswith("ATOM"):
                lineData = line.split()
                if lineData[1] in mapAtomNameToDefaultCharmmType:
                    mapCgenffTypeToCharmmType[lineData[2]] = mapAtomNameToDefaultCharmmType[lineData[1]]
                    lineData[2] = mapAtomNameToDefaultCharmmType[lineData[1]]
                line = "\t". join(lineData) + "\n"
            outRtf.write(line)

    move(tmpRtf, moleculeRtf)

    return mapCgenffTypeToCharmmType
                
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲##
def create_atom_type_map(config: dict) -> dict:
    """
    Creates a map of atom types for CHARMM given a RTF file

    Args: 
        config (dict): contains all run info
    Returns:
        config (dict): updated config
    """
    moleculeRtf = config["runtimeInfo"]["madeByStitching"]["moleculeRtf"]

    atomTypeMap = {}
    with open(moleculeRtf, "r") as inRtf:
        for line in inRtf.readlines():
            if line.startswith("ATOM"):
                lineData = line.split()
                atomTypeMap[lineData[1]] = lineData[2]

    config["runtimeInfo"]["madeByStitching"]["atomTypeMap"] = atomTypeMap

    return config

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲##
def edit_rtf_charges(config: dict) -> dict:
    """
    Edits an RTF file to change the charges to those calculated by drFrankenstein
    Args: 
        config (dict): contains all run info
    Returns:
        config (dict): updated config 
    """
    moleculeRtf = config["runtimeInfo"]["madeByStitching"]["moleculeRtf"]
    chargesCsv = config["runtimeInfo"]["madeByCharges"]["chargesCsv"]
    mmTotalCalculationDir = config["runtimeInfo"]["madeByStitching"]["mmTotalCalculationDir"]

    chargesDf = pd.read_csv(chargesCsv, index_col="Unnamed: 0")
    atomIndex = 0
    tmpRtf = p.join(mmTotalCalculationDir, "update_charges.rtf")
    with open(moleculeRtf, "r") as inRtf, open(tmpRtf, "w") as outRtf:
        for line in inRtf.readlines():
            if line.startswith("ATOM"):
                lineData = line.split()
                lineData[3] = round(chargesDf.loc[atomIndex, "Charge"],3).astype(str)
                line = "\t". join(lineData) + "\n"
                atomIndex += 1
            outRtf.write(line)

    move(tmpRtf, moleculeRtf)

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲##
def split_charmm_str(config: dict) -> dict:
    """
    Splits a CHARMM STR (stream) file into RTF and PRM files

    Args: 
        config (dict): contains all run info
    Returns:
        config (dict): updated config 
    """

    moleculeName = config["moleculeInfo"]["moleculeName"]
    cappedStr = config["runtimeInfo"]["madeByCapping"]["cappedMoleculeStr"]
    mmTorsionCalculationDir = config["runtimeInfo"]["madeByStitching"]["mmTorsionCalculationDir"]


    moleculePrm = p.join(mmTorsionCalculationDir, f"{moleculeName}.prm")
    moleculeRtf = p.join(mmTorsionCalculationDir, f"{moleculeName}.rtf")

    with open (cappedStr, "r") as stream, open(moleculePrm, "w") as prm, open(moleculeRtf, "w") as rtf:
        writeRtf = True
        writePrm = False
        for line in stream:
            if line.startswith("read param card flex append"):
                writeRtf = False
                writePrm = True
            if writeRtf:
                rtf.write(line)
            if writePrm:
                prm.write(line)


    config["runtimeInfo"]["madeByStitching"]["moleculeRtf"] = moleculeRtf
    config["runtimeInfo"]["madeByStitching"]["moleculePrm"] = moleculePrm
    return config
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲##
def make_charmm_psf(config: dict) -> dict:
    """
    Uses PSFGen to generate a PSF file from RTF and PRM files

    Args: 
        config (dict): contains all run info
    Returns:
        config (dict): updated config 
    """
    ## unpack molecule
    moleculeRtf = config["runtimeInfo"]["madeByStitching"]["moleculeRtf"]
    cgenffRtf = config["runtimeInfo"]["madeByStitching"]["cgenffRtf"]
    cappedPdb = config["runtimeInfo"]["madeByCapping"]["cappedPdb"]
    mmTorsionCalculationDir = config["runtimeInfo"]["madeByStitching"]["mmTorsionCalculationDir"]
    moleculeName = config["moleculeInfo"]["moleculeName"]

    ## make a PDB file with only one RES_ID and RES_NAME value for all atoms
    cappedDf = pdbUtils.pdb2df(cappedPdb)
    cappedDf["RES_ID"] = "1"
    cappedDf["RES_NAME"] = moleculeName
    tmpPdb = p.join(mmTorsionCalculationDir, f"united_capped_{moleculeName}.pdb")
    pdbUtils.df2pdb(cappedDf, tmpPdb)


    gen = PsfGen(output="/dev/null")  # Suppress output since there's too much
    gen.read_topology(moleculeRtf)
    gen.read_topology(cgenffRtf)
    # Read TIP3P water
    gen.add_segment(segid="A", pdbfile=tmpPdb)
    gen.read_coords(segid="A", filename=tmpPdb)


    psfFile = p.join(mmTorsionCalculationDir, f"{moleculeName}_capped.psf")
    gen.write_psf(psfFile)
    config["runtimeInfo"]["madeByStitching"]["moleculePsf"] = psfFile
    return config

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲##
def find_cgenff_params(config: dict) -> dict:
    """
    Looks inside the drFrankenstein dir to find CGenFF PRM and RTF files

    Args: 
        config (dict): contains all run info
    Returns:
        config (dict): updated config 
    """

    stitchingSrcDir =  p.dirname(p.abspath(__file__))
    cgenffPrm = p.join(stitchingSrcDir, "par_all36_cgenff.prm")
    cgenffRtf = p.join(stitchingSrcDir, "top_all36_cgenff.rtf")

    if not p.isfile(cgenffPrm):
        raise FileNotFoundError(f"Cannot find CGenFF PRM file at {cgenffPrm}")
    if not p.isfile(cgenffRtf):
        raise FileNotFoundError(f"Cannot find CGenFF RTF file at {cgenffRtf}")


    config["runtimeInfo"]["madeByStitching"]["cgenffPrm"] = cgenffPrm
    config["runtimeInfo"]["madeByStitching"]["cgenffRtf"] = cgenffRtf

    return config