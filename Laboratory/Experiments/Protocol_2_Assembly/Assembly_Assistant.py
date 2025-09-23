import os
from os import path as p
import warnings
from subprocess import call, run, PIPE, STDOUT
from pdbUtils import pdbUtils
import parmed as pmd
from parmed.charmm import CharmmParameterSet, CharmmPsfFile
from parmed.exceptions import ParameterWarning

# Suppress ParameterWarning
warnings.filterwarnings('ignore', category=ParameterWarning)

## drFRANKENSTEIN LIBRARIES ##
from OperatingTools import file_parsers

# Placeholder classes (extend if needed)
class FilePath:
    pass

class DirectoryPath:
    pass

def run_tleap_to_make_params(inMol2: FilePath,
                            molFrcmod: FilePath,
                              outDir: DirectoryPath,
                                moleculeName: str):
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

    prmtop: FilePath = p.join(outDir, f"{moleculeName}_cappped.prmtop")
    inpcrd: FilePath = p.join(outDir, f"{moleculeName}_capped.inpcrd")

    tleapInput: FilePath = p.join(outDir, f"leap.in")
    with open(tleapInput, "w") as f:
        f.write("source leaprc.gaff2\n")
        f.write(f"mol  = loadmol2 {inMol2} \n")
        f.write(f"loadamberparams {molFrcmod} \n") # use frcmod previously made
        f.write(f"saveamberparm mol {prmtop} {inpcrd}  \n")
        f.write("quit")

    tleapOutput: FilePath = p.join(outDir, f"tleap.out")

    tleapCommand: list = ["tleap", "-f", tleapInput, ">", tleapOutput]

    call(tleapCommand, stdout=PIPE)

    return prmtop


def find_default_amber_parameters(amberHome: DirectoryPath) -> tuple[FilePath, FilePath]:
    """
    Finds gaff2.dat and parm19.dat from AMBERHOME

    Args:
        amberHome (DirectoryPath): path to AMBERHOME

    Returns:
        gaff2Dat (FilePath): path to gaff2.dat
        parm19Dat (FilePath): path to parm19.dat
    
    """

    gaff2Dat = p.join(amberHome, "dat", "leap", "parm", "gaff2.dat")
    parm19Dat = p.join(amberHome, "dat", "leap", "parm", "parm19.dat")

    if not p.exists(gaff2Dat):
        raise FileNotFoundError(f"Could not find {gaff2Dat}")
    if not p.exists(parm19Dat):
        raise FileNotFoundError(f"Could not find {parm19Dat}")

    return gaff2Dat, parm19Dat

def create_frcmod_file(mol2File: FilePath,
                        frcmodFile: FilePath,
                          gaff2Dat: FilePath) -> FilePath:
    """
    uses parmchk2 to create a frcmod file from a mol2 file

    Args:
        chargesMol2 (FilePath): mol2 file of charges
        moleculeName (str): name of molecule
        config (dict): config dict

    Returns:
        molFrcmod (FilePath): path to frcmod file
    
    """

    ## run run parmchk2
    parmchk2Command = ["parmchk2",
                        "-i", mol2File,
                          "-f", "mol2",
                            "-o", frcmodFile,
                              "-a", "Y",
                                "-p", gaff2Dat]

    try:
        call(parmchk2Command)
    except Exception as e:
        raise(e)

    return None


def pdb2mol2(inPdb: FilePath,
              outMol2: FilePath,
                workingDir: DirectoryPath,
                config: dict) -> None:
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


    if p.exists(outMol2):
        return None
    os.chdir(workingDir)

    ## unpack config
    netCharge = config["moleculeInfo"]["charge"]

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
        "-c", "bcc", "-nc", str(netCharge)
    ]
    with open(antechamberOut, 'w') as outfile:
        run(antechamberCommand, stdout=outfile, stderr=STDOUT)

    ## clean up temporary files
    os.remove(tmpPdb)
    filesToRemove = [f for f in os.listdir(workingDir) if f.startswith("ANTECHAMBER")]
    for f in filesToRemove:
        os.remove(p.join(workingDir, f))
    return None


def  update_rtf_types(inRtf: FilePath, nameToDesiredType: dict,  config: dict) -> dict:
    """
    Update RTF file with updated atom types

    Args:
        inRtf (FilePath): input RTF file
        outRtf (FilePath): output RTF file
        config (dict): configuration dictionary

    Returns:
        config (dict): updated configuration dictionary
    """
    assemblyDir = config["runtimeInfo"]["madeByAssembly"]["assemblyDir"]
    moleculeName = config["moleculeInfo"]["moleculeName"]
    parsedRtf = file_parsers.parse_rtf(inRtf)

    for atom in parsedRtf["residues"][moleculeName]["atoms"]:
        if not atom["name"] in nameToDesiredType:
            continue
        atom["type"] = nameToDesiredType[atom["name"]]

    outRtf = p.join(assemblyDir, f"{moleculeName}_assembled.rtf")
    file_parsers.write_rtf(parsedRtf, outRtf)


    return outRtf

def save_modified_prm_file(parmedPsf: CharmmPsfFile, config: dict) -> dict:
    """
    Save modified RTF files.

    Args:
        parmedPsf (CharmmPsfFile): ParmEd CharmmPsfFile object with updated atom types.
        config (dict): Configuration dictionary.

    Returns:
        config (dict): Updated configuration dictionary with file paths for RTF, PRM, and PSF files.
    """
    ## unpack config 
    assemblyDir = config["runtimeInfo"]["madeByAssembly"]["assemblyDir"]
    moleculeName = config["moleculeInfo"]["moleculeName"]

    outputParams = CharmmParameterSet.from_structure(parmedPsf)
    outPrm = p.join(assemblyDir, f"{moleculeName}_assembled.prm")
    outPsf = p.join(assemblyDir, f"{moleculeName}_assembled.psf")

    outputParams.write(par=outPrm)
    parmedPsf.save(outPsf, overwrite=True)

    return outPrm, outPsf

    # config["runtimeInfo"]["madeByAssembly"]["assembledPsf"] = outPsf 

    
def assign_missing_impropers(missingParams: CharmmParameterSet,
                        parmedPsf: CharmmPsfFile,
                          nameToCgenffType: dict,
                            completeParameterSet: CharmmParameterSet) -> CharmmParameterSet:
        
    """
    Check for missing improper parameters and assign from CGenFF parameters.

    Args:
        missingParams: (CharmmParameterSet) object to store missing parameters.
        parmedPsf: (CharmmPsfFile) object with updated atom types.
        nameToCgenffType: (dict) mapping atom names to original CGenFF types.
        completeParameterSet: (CharmmParameterSet) object with all parameters.

    Returns:
        missingParams: (CharmmParameterSet) object with added missing improper parameters.
    """
    for improper in parmedPsf.impropers:
        atomNames=[improper.atom1.name, improper.atom2.name, improper.atom3.name, improper.atom4.name]
        newTypes = [improper.atom1.type, improper.atom2.type, improper.atom3.type, improper.atom4.type]
        origTypes = [nameToCgenffType.get(name, atomType) for name, atomType in zip(atomNames, newTypes)]

        improperKey = tuple(newTypes)
        if improperKey not in completeParameterSet.improper_types and improperKey not in missingParams.improper_types:
            origImproperKey = tuple(origTypes)
            if origImproperKey in completeParameterSet.improper_types:
                improperTypes = completeParameterSet.improper_types[origImproperKey]
                missingParams.improper_types[improperKey] = improperTypes
            else:
                raise ValueError(f"Warning: No CGenFF improper parameters for \n" \
                    f"CGenFF: {origTypes}, CHARMM: {newTypes} NAME {atomNames}")

    return missingParams
def assign_missing_dihedrals(missingParams: CharmmParameterSet,
                        parmedPsf: CharmmPsfFile,
                          nameToCgenffType: dict,
                            completeParameterSet: CharmmParameterSet) -> CharmmParameterSet:
    # Check for missing dihedral parameters
    """
    Check for missing dihedral parameters and assign from CGenFF parameters.

    Args:
        missingParams: (CharmmParameterSet) object to store missing parameters.
        parmedPsf: (CharmmPsfFile) object with updated atom types.
        nameToCgenffType: (dict) mapping atom names to original CGenFF types.
        completeParameterSet: (CharmmParameterSet) object with all parameters.

    Returns:
        missingParams: (CharmmParameterSet) object with added missing dihedral parameters.
    """
    for dihedral in parmedPsf.dihedrals:
        atomNames=[dihedral.atom1.name, dihedral.atom2.name, dihedral.atom3.name, dihedral.atom4.name]
        newTypes = [dihedral.atom1.type, dihedral.atom2.type, dihedral.atom3.type, dihedral.atom4.type]
        origTypes = [nameToCgenffType.get(name, atomType) for name, atomType in zip(atomNames, newTypes)]

        dihedralKey = tuple(newTypes)
        if dihedralKey not in completeParameterSet.dihedral_types and dihedralKey not in missingParams.dihedral_types:
            origDihedralKey = tuple(origTypes)
            origDihedralKeyRev = tuple(reversed(origDihedralKey))
            if origDihedralKey in completeParameterSet.dihedral_types:
                dihedralTypes = completeParameterSet.dihedral_types[origDihedralKey]
                missingParams.dihedral_types[dihedralKey] = dihedralTypes
            elif origDihedralKeyRev in completeParameterSet.dihedral_types:
                dihedralTypes = completeParameterSet.dihedral_types[origDihedralKeyRev]
                missingParams.dihedral_types[dihedralKey] = dihedralTypes
            else:
                raise ValueError(f"Warning: No CGenFF dihedral parameters for \n" \
                    f"CGenFF: {origTypes}, CHARMM: {newTypes} NAME {atomNames}")

    return missingParams
def assign_missing_angles(missingParams: CharmmParameterSet,
                        parmedPsf: CharmmPsfFile,
                          nameToCgenffType: dict,
                            completeParameterSet: CharmmParameterSet) -> CharmmParameterSet:
    # Check for missing angle parameters
    """
    Check for missing angle parameters and assign from CGenFF parameters.

    Args:
        missingParams: (CharmmParameterSet) object to store missing parameters.
        parmedPsf: (CharmmPsfFile) object with updated atom types.
        nameToCgenffType: (dict) mapping atom names to original CGenFF types.
        completeParameterSet: (CharmmParameterSet) object with all parameters.

    Returns:
        missingParams: (CharmmParameterSet) object with added missing angle parameters.
    """
    for angle in parmedPsf.angles:
        atomNames=[angle.atom1.name, angle.atom2.name, angle.atom3.name]
        newTypes = [angle.atom1.type, angle.atom2.type, angle.atom3.type]
        origTypes = [nameToCgenffType.get(name, atomType) for name, atomType in zip(atomNames, newTypes)] 
        angleKey = tuple(newTypes)
        if angleKey not in completeParameterSet.angle_types and angleKey not in missingParams.angle_types:
            origAngleKey = tuple(origTypes)
            if origAngleKey in completeParameterSet.angle_types:
                angleType = completeParameterSet.angle_types[origAngleKey]
                missingParams.angle_types[angleKey] = angleType
            elif origAngleKey[::-1] in completeParameterSet.angle_types:
                angleType = completeParameterSet.angle_types[origAngleKey[::-1]]
                missingParams.angle_types[angleKey] = angleType
            else:
                raise ValueError(f"Warning: No CGenFF angle parameters for \n" \
                                 f"CGenFF: {origTypes}, CHARMM: {newTypes} NAME {atomNames}")
            
    return missingParams
def assign_missing_bonds(missingParams: CharmmParameterSet,
                        parmedPsf: CharmmPsfFile,
                          nameToCgenffType: dict,
                            completeParameterSet: CharmmParameterSet) -> CharmmParameterSet:
    """
    Check for missing bond parameters and assign from CGenFF parameters.

    Args:
        missingParams: (CharmmParameterSet) object to store missing parameters.
        parmedPsf: (CharmmPsfFile) object with updated atom types.
        nameToCgenffType: (dict) mapping atom names to original CGenFF types.
        completeParameterSet: (CharmmParameterSet) object with all parameters.

    Returns:
        missingParams: (CharmmParameterSet) object with added missing bond parameters.
    """
    for bond in parmedPsf.bonds:
        atomNames = [bond.atom1.name, bond.atom2.name]
        newTypes = [bond.atom1.type, bond.atom2.type]
        origTypes = [nameToCgenffType.get(name, atomType) for name, atomType in zip(atomNames, newTypes)] 
        bondKey = tuple(sorted(newTypes))  # Bonds are order-independent
        if bondKey not in completeParameterSet.bond_types and bondKey not in missingParams.bond_types:
            # Map to original CGenFF types
            origBondKey = tuple(sorted(origTypes))
            if origBondKey in completeParameterSet.bond_types:
                bondType = completeParameterSet.bond_types[origBondKey]
                missingParams.bond_types[bondKey] = bondType
            else:
                raise ValueError(f"Warning: No CGenFF bond parameters for \n" \
                                 f"CGenFF: {origTypes}, CHARMM: {newTypes} NAME {atomNames}")
    return missingParams


def update_psf_atom_types(parmedPsf: CharmmPsfFile, nameToDesiredType: dict) -> CharmmPsfFile:
    """
    Update atom types in the PSF structure based on a mapping.

    Args:
        parmedPsf: ParmEd CharmmPsfFile object.
        nameToDesiredType: Dictionary mapping atom names to new CHARMM36m types.

    Returns:
        parmedPsf: ParmEd CharmmPsfFile object with updated atom types

    """
    for atom in parmedPsf.atoms:
        if atom.name in nameToDesiredType:
            oldType = atom.type
            atom.type = nameToDesiredType[atom.name]
    return parmedPsf

def get_cgenff_atom_types(parmedPsf: CharmmPsfFile, nameToDesiredType: dict) -> dict:
    """
    Get the original CGenFF atom types for atoms to be retyped.

    Args:
        psf: ParmEd CharmmPsfFile object with CGenFF types.
        nameToDesiredType: Dictionary mapping atom names to new CHARMM36m types.

    Returns:
        Dictionary mapping atom names to their original CGenFF types.
    """
    cgenffTypes = {}
    for atom in parmedPsf.atoms:
        if atom.name in nameToDesiredType:
            cgenffTypes[atom.name] = atom.type
    return cgenffTypes


def load_psf_with_params(psfFile:FilePath, params:tuple) -> CharmmPsfFile:
    """
    Reads a PSF file into Parmed and loads parameters

    Args:
        psfFile (FilePath): Path to the PSF file to read
        params (tuple): Tuple of parameter paths 
    Returns:
        psf (CharmmPsfFile): Parmed PSF object
    
    """

    parmedPsf = CharmmPsfFile(psfFile)
    originalCgenffParams = CharmmParameterSet(*params)
    parmedPsf.load_parameters(originalCgenffParams)

    return parmedPsf


def find_default_charmm_parameters() -> dict:
    """
    Uses relative locations to find default CHARMM parameters

    Args:
        None [Uses location of this file as reference]
    Returns:
        charmmDefaultParams (dict): contains paths to all default CHARMM parameters
    """
    assemblySrcDir = p.dirname(p.abspath(__file__))
    labDir = p.dirname(p.dirname(assemblySrcDir))
    ingredientDir = p.join(labDir, "Ingredients")
    charmmDir = p.join(ingredientDir, "CHARMM")
    charmmDownloadDir = [p.join(charmmDir, dir) for dir in os.listdir(charmmDir) 
                 if p.isdir(p.join(charmmDir, dir))][0]
    topparDir = p.join(charmmDownloadDir, "toppar")
    


    charmmRtf = p.join(topparDir, "top_all36_prot.rtf")
    charmmPrm = p.join(topparDir, "par_all36m_prot.prm")

    cgenffRtf = p.join(topparDir , "top_all36_cgenff.rtf")
    cgenffPrm = p.join(topparDir , "par_all36_cgenff.prm")

    filesNotFound = []
    for file in [charmmRtf, charmmPrm, cgenffRtf, cgenffPrm]:
        if not p.isfile(file):
            filesNotFound.append(file)

    if len(filesNotFound) > 0:
        raise FileNotFoundError(f"Cannot find CHARMM parameters at {filesNotFound}")

    charmmDefaultParams = {
        "charmmRtf": charmmRtf,
        "charmmPrm": charmmPrm,
        "cgenffRtf": cgenffRtf,
        "cgenffPrm": cgenffPrm
    }

    return charmmDefaultParams
