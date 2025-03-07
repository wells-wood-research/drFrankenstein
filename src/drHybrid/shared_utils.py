## BASIC LIBRARIES ##
import os
from os import path as p
from subprocess import call, PIPE, run, STDOUT
import pandas as pd
import re
import numpy as np
from pdbUtils import pdbUtils
from shutil import move
import random

## OPENMM LIBRARIES
import openmm.app as app
import openmm as openmm
import  openmm.unit  as unit
## CLEAN CODE CLASSES ##
class FilePath:
    pass
class DirectoryPath:
    pass

def shuffle_torsion_tags(torsionTags: list) -> list:
    shuffledTorsionTags = torsionTags
    random.shuffle(shuffledTorsionTags)

    return shuffledTorsionTags

################################################################
def update_frcmod(config, torsionTag, torsionParamDf):
    inFrcmod = config["pathInfo"]["moleculeFrcmod"]
    tmpFrcmod = p.splitext(inFrcmod)[0] + "_tmp.frcmod"

    ##TODO: sort out multiplicities
    multiplicity = 1


    atomTypeMap = config["moleculeInfo"]["atomTypeMap"]
    atomNames = torsionTag.split("-")
    atomTypes = [atomTypeMap[name] for name in atomNames]
    atomTypesReversed = atomTypes[::-1]
    torsionIdentifier = "-".join([el + ' ' if len(el) == 1 else el for el in atomTypes])
    torsionIdentifierReversed = "-".join([el + ' ' if len(el) == 1 else el for el in atomTypesReversed])

    writeLines = False

    newTorsionBlock = ""
    for _, row in torsionParamDf.iloc[:-1].iterrows():
        amplitude = f"{row['Amplitude']:.3f}"
        phase = f"{row['Phase']:.3f}"
        period = f"{-row['Period']:.3f}"
        newTorsionBlock += (f"{torsionIdentifier:<14}{multiplicity:<5}"
                            f"{amplitude:>5}{phase:>14}{period:>16}"
                            f"\t\tMADE BY drFRANKENSTEIN\n")
    lastRow = torsionParamDf.iloc[-1]
    amplitude = f"{lastRow['Amplitude']:.3f}"
    phase = f"{lastRow['Phase']:.3f}"
    period = f"{lastRow['Period']:.3f}"
    newTorsionBlock += (f"{torsionIdentifier:<14}{multiplicity:<5}"
                        f"{amplitude:>5}{phase:>14}{period:>16}"
                        f"\t\tMADE BY drFRANKENSTEIN\n")

    writeLines = False
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
    move(tmpFrcmod, inFrcmod)
    return config
################################################################
def pdb2mol2(inPdb: FilePath,
              outMol2: FilePath,
                workingDir: DirectoryPath) -> None:
    """
    Uses antechamber to convert pdb to mol2
    Writes some unwanted temporary files
    TODO: clean these up

    Args:
        inPdb (FilePath): input file in PDB format
        outPdb (FilePath): output file in MOL2 format
        workingDir (DirectoryPath): working directory (vital for cleanup)

    Returns:
        None (outMol2 has already been defined!)
    
    """
    os.chdir(workingDir)

    ## set RES_ID to 1 for all atoms to keep antechamber happy
    pdbDf = pdbUtils.pdb2df(inPdb)
    pdbDf["RES_ID"] = 1
    tmpPdb = p.join(workingDir, "tmp.pdb")
    pdbUtils.df2pdb(pdbDf, tmpPdb)
    ## get index, set path for antechamber to write outputs
    index: str = p.basename(inPdb).split("_")[1].split(".")[0]
    antechamberOut: FilePath = p.join(workingDir, f"antechamber_{index}.out")
    ## run antechamber to create MOL2 file from PDB
    antechamberCommand: list = [
        "antechamber", "-i", tmpPdb, "-fi", "pdb", "-o", outMol2,
        "-fo", "mol2", "-at", "gaff2", "-rn", "MOL", "-s", "2",
    ]
    with open(antechamberOut, 'w') as outfile:
        run(antechamberCommand, stdout=outfile, stderr=STDOUT)
################################################################
def edit_mol2_partial_charges(inMol2: FilePath,
                               chargesDf: pd.DataFrame,
                                 outMol2: FilePath) -> None:
    """
    Gets partial charges stored in a dataframe and pastes them into the
    charges column of a MOL2 file

    Args:
        inMol2 (FilePath): input MOL2 file
        chargesDf (pd.DataFrame): dataframe of partial charges
        outMol2 (FilePath): output MOL2 file

    Returns:
        None (outMol2 is already defined!)
    
    """
    ## init atom index counter
    atomIndex: int = 0
    ## open inMol2 for reading and outMol2 for writing
    with open(inMol2, 'r') as inMol2, open(outMol2, "w") as writeMol2:
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

################################################################
def edit_mo2_atom_types(inMol2: FilePath,
                         outMol2: FilePath) -> None:
    
    """
    Edits the atom types in a MOL2 file to make compatable with gaff2
    This may need to be added to in the case of new weird atom types
    from antechamber

    Args:
        inMol2 (FilePath): input MOL2 file
        outMol2 (FilePath): output MOL2 file

    Returns:
        None (outMol2 is already defined!)
    """


    ## open inMol2 for reading and outMol2 for writing
    with open(inMol2, 'r') as inMol2, open(outMol2, "w") as outMol2:
        ## loop through inMol2 until we get to the atom section
        mol2Lines = inMol2.readlines()
        isAtomLine = False
        for line in mol2Lines:
            if line.strip() == "@<TRIPOS>ATOM":
                isAtomLine = True
            elif line.strip() == "@<TRIPOS>BOND":
                isAtomLine = False
            elif isAtomLine:
                ## get atom name
                atomData = line.split()
                atomName = atomData[1]
                ## set alpha-carbon to CX rather than c3
                if atomName == "CA":
                    line = re.sub(r'(\s)c3(\s)', r'\1CX\2', line)
                ## other substitutions
                line = re.sub(r'(\s)os(\s)', r'\1o \2', line)
                line = re.sub(r'(\s)n7(\s)', r'\1ns\2', line)
                line = re.sub(r'(\s)h1(\s)', r'\1hc\2', line)
            outMol2.write(line)

################################################################
def create_atom_type_map(inMol2: FilePath) -> dict:
    """
    Creates a dict that maps atom names to atom types
    using a MOL2 file

    Args:
        inMol2 (FilePath): input MOL2 file

    Returns:
        atomTypeMape (dict): dict mapping atom names to atom types
    
    """
    ## init empty dict
    atomTypeMap = {}
    ## open inMol2 for reading
    with open(inMol2, 'r') as inMol2:
        ## loop through inMol2 until we get to the atoms section
        mol2Lines = inMol2.readlines()
        isAtomLine = False
        for line in mol2Lines:
            if line.strip() == "@<TRIPOS>ATOM":
                isAtomLine = True
            elif line.strip() == "@<TRIPOS>BOND":
                isAtomLine = False
            elif isAtomLine:
                ## get atom name and type, add to dict
                atomData = line.split()
                atomTypeMap[atomData[1]] = atomData[5]
    return atomTypeMap

################################################################
def create_frcmod_file(chargesMol2: FilePath,
                        moleculeName: str,
                          config: dict) -> FilePath:
    """
    uses parmchk2 to create a frcmod file from a mol2 file

    Args:
        chargesMol2 (FilePath): mol2 file of charges
        moleculeName (str): name of molecule
        config (dict): config dict

    Returns:
        molFrcmod (FilePath): path to frcmod file
    
    """
    ## get path to gaff2.dat fom config file
    gaff2Dat: FilePath = config["pathInfo"]["gaff2Dat"]
    ## init molFrcmod output file
    molFrcmod = p.join(config["pathInfo"]["mmTorsionCalculationDir"], f"{moleculeName}.frcmod")
    ## run run parmchk2
    parmchk2Command = ["parmchk2", "-i", chargesMol2, "-f", "mol2", "-o", molFrcmod, "-a", "-Y", "-p", gaff2Dat]
    call(parmchk2Command)

    return molFrcmod

################################################################
def extract_torsion_parameters(config: dict, torsionTag: str) -> dict:
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

    molFrcmod = config["pathInfo"]["moleculeFrcmod"]
    parsedFrcmod: dict = parse_frcmod(molFrcmod)

    atomTypeMap = config["moleculeInfo"]["atomTypeMap"]

    ## get torsion atom names
    torsionAtoms: list = torsionTag.split("-")
    ## get torsion atom types
    torsionAtomTypes: list = [atomTypeMap[atom] for atom in torsionAtoms]
    ## find torsion parameters for these atom types (account for reverse torsion order)
    torsionAtomTypes_reversed = torsionAtomTypes[::-1]
    lookForAtoms = [torsionAtomTypes, torsionAtomTypes_reversed]
    try:
        torsionParameters = [entry for entry in parsedFrcmod["DIHEDRALS"] if entry["atoms"] in lookForAtoms]

    except:
        raise Exception(f"Torsion {torsionTag} not found in frcmod file")
    ## add torsion parameters to dict

    return torsionParameters
################################################################
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
################################################################
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
################################################################
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

    return mmTorsionEnergy

################################################################

def sort_out_directories(config) -> dict:
    """
    Makes directories for the parameter fitting stage
    Adds these paths to the "pathInfo" section of config

    Args: 
        config (dict): contains all information needed for run
    
    Returns:
        config (dict): updated config with pathInfo of new dirs
    
    """
    outputDir: DirectoryPath = config["pathInfo"]["outputDir"]

    parameterFittingTopDir: DirectoryPath = p.join(outputDir, "05_parameter_fitting")
    os.makedirs(parameterFittingTopDir, exist_ok=True)
    config["pathInfo"]["parameterFittingTopDir"] = parameterFittingTopDir

    mmTorsionCalculationDir: DirectoryPath = p.join(parameterFittingTopDir, "mm_torsion_energies")
    os.makedirs(mmTorsionCalculationDir, exist_ok=True)
    config["pathInfo"]["mmTorsionCalculationDir"] = mmTorsionCalculationDir

    mmTotalCalculationDir: DirectoryPath = p.join(parameterFittingTopDir, "mm_total_energies")
    os.makedirs(mmTotalCalculationDir, exist_ok=True)
    config["pathInfo"]["mmTotalCalculationDir"] = mmTotalCalculationDir

    return config
################################################################
def get_completed_torsion_scan_dirs(config: dict, torsionTag:str):
    """
    Looks through ORCA scan directories 
    Works out whether scan completed
    Returns a list of completed scan dirs

    Args:
        config (dict): contains all information needed for run
        torsionTag (str): identifier for torsion (format "W-X-Y-Z")

    Returns:
        completedTorsionScanDirs (list): list of completed scan directories
    
    """
    ## init empty list to store scan dirs in
    completedTorsionScanDirs: list = []

    ## find dir where ORCA torsion scans were performed
    torsionTopDir: DirectoryPath = config["pathInfo"]["torsionTopDir"]

    torsionDir: DirectoryPath = p.join(torsionTopDir, f"torsion_{torsionTag}")

    ## loop through batches of scans
    for conformerName in os.listdir(torsionDir):
        if not "conformer" in conformerName:
            continue
        ## find conformer dir
        conformerDir: DirectoryPath = p.join(torsionDir, conformerName)
        for calculationName in os.listdir(conformerDir):
            ## get calculation dir (forwards and backwards)
            calculationDir: DirectoryPath = p.join(conformerDir, calculationName)
            if  calculationName.endswith("forwards") or calculationName.endswith("backwards"):
                ## scan is completed if orca_trj.xyz exists - add to list
                if p.isfile(p.join(calculationDir, "orca_trj.xyz")):
                    completedTorsionScanDirs.append(calculationDir)
    return completedTorsionScanDirs
################################################################
def get_scan_angles_from_orca_inp(scanDir: DirectoryPath):
    """
    Read through an ORCA input file and find start and end angles of a dihedral scan
    Creates a range between these

    Args:
        scanDir (DirectoryPath): location of orca.inp

    Returns:
        angleRange (list): list of angles between start and end angles
    """
    ## find orca input
    orcaInp = p.join(scanDir, "orca.inp")
    if not p.isfile(orcaInp):
        raise FileNotFoundError(f"orca.inp not found in {scanDir}")
    ## open orca input for reading
    with open(orcaInp, "r") as f:
        for line in f:
            ## find dihedral scan, get start, end angle
            if line.startswith("D"):
                lineData = line.split()
                startAngle = float(lineData[6][:-1])
                endAngle = float(lineData[7][:-1])
                break
    ## create range
    if startAngle < endAngle:
        stepSize = 10
    else:
        stepSize = -10
    angleRange = np.arange(startAngle, endAngle + stepSize, stepSize)
    
    return angleRange
################################################################
def rescale_angles_0_360(angle: np.array) -> np.array:
    """
    puts angles on a -180 to +180 scale
    """
    angle = angle % 360  # First, reduce the angle to the 0-360 range
    # if angle > 180:
    #     angle -= 360  # Shift angles greater than 180 to the negative side
    return angle
################################################################
def merge_energy_dfs(energyDfs: list) -> pd.DataFrame:
    """
    merges energy DataFrames on the Angle column

    Args:
        energyDfs (list): list of DataFrames with Angle column
    
    Returns:
        mergedDf (pd.DataFrame): merged DataFrame
    """
    ##TODO: is this efficient?
    energyDfs = [df for df in energyDfs if not df is None]

    mergedDf = energyDfs[0][['Angle', 'Energy']].rename(
        columns={'Energy': 'Energy_0'}
    )

    for i, df in enumerate(energyDfs[1:], start=1):
        mergedDf = mergedDf.merge(
            df[['Angle', 'Energy']],
            on='Angle',
            how='outer',
            suffixes=('', f'_df{i}')
        ).rename(columns={'Energy': f'Energy_{i}'})
    mergedDf = mergedDf.groupby('Angle', as_index=False).first()

    return mergedDf
################################################################
def update_pdb_coords(inPdb: FilePath, xyzFile: FilePath, outPdb: FilePath) -> None:
    """
    updates PDB file with XYZ coords
    
    Args:
        inPdb (FilePath): input PDB file
        xyzFile (FilePath): input XYZ file
        outPdb (FilePath): output PDB file

    Returns:
        None (outPdb already defined!)
    """

    inDf = pdbUtils.pdb2df(inPdb)
    xyzDf = xyz2df(xyzFile)

    inDf["X"] = xyzDf["x"]
    inDf["Y"] = xyzDf["y"]
    inDf["Z"] = xyzDf["z"]

    pdbUtils.df2pdb(inDf, outPdb)
################################################################
def xyz2df(xyzFile: FilePath) -> pd.DataFrame:
    """
    converts an xyz file to a pd.DataFrame

    Args:
        xyzFile (FilePath): input XYZ file

    Returns:
        xyzDf (pd.DataFrame): dataframe containing index, element and coords
    """

    xyzData = []
    atomIndex = 0
    with open(xyzFile, "r") as f:
        lines = f.readlines()[2:]
        for line in lines:
            atomIndex +=1
            lineData = line.split()
            xyzData.append({"index": atomIndex,
                            "element": lineData[0],
                              "x": float(lineData[1]), 
                              "y": float(lineData[2]),
                                "z": float(lineData[3])})
    return pd.DataFrame(xyzData)
################################################################

def make_prmtop_and_inpcrd(inMol2: FilePath,
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

################################################################
def run_mm_singlepoints(trajPdbs: list, moleculePrmtop: FilePath, moleculeInpcrd: FilePath) -> float:
    """
    Runs a singlepoint energy calculation at the MM level 
    using OpenMM

    Args:
        prmtop (FilePath): topology file for AMBER
        impcrd (FilePath): coordinate file for AMBER

    Returns:
        singlePointEnergy (float): energy of prmtop // inpcrd
    """
    # Load Amber files and create system
    prmtop: app.Topology = app.AmberPrmtopFile(moleculePrmtop)
    inpcrd: app.InpcrdFile = app.AmberInpcrdFile(moleculeInpcrd)

    # Create the system.
    system: openmm.System = prmtop.createSystem(nonbondedMethod=app.NoCutoff,
                                                nonbondedCutoff=1 * unit.nanometer,
                                                constraints=None)

    integrator = openmm.LangevinIntegrator(300, 1/unit.picosecond,  0.0005*unit.picoseconds)
    platform = openmm.Platform.getPlatformByName('CPU')

    simulation = app.Simulation(prmtop.topology, system, integrator, platform)
    ## set coordinates of simulation 
    singlePointEnergies = []
    for trajPdb in trajPdbs:
        pdbFile = app.PDBFile(trajPdb)
        simulation.context.setPositions(pdbFile.positions)
        state: openmm.State = simulation.context.getState(getPositions=True, getEnergy=True)
        singlePointEnergy = state.getPotentialEnergy() / unit.kilocalories_per_mole

        trajIndex = trajPdb.split(".")[0].split("_")[1]
        singlePointEnergies.append((trajIndex,singlePointEnergy))

    return singlePointEnergies
################################################################