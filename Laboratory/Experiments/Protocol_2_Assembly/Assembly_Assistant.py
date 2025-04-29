import os
from os import path as p
import warnings
import parmed as pmd
from parmed.charmm import CharmmParameterSet, CharmmPsfFile
from parmed.exceptions import ParameterWarning
from psfgen import PsfGen

# Suppress ParameterWarning
warnings.filterwarnings('ignore', category=ParameterWarning)


from . import Assembly_Monster
# Placeholder classes (extend if needed)
class FilePath:
    pass

class DirectoryPath:
    pass



def save_modified_parameter_files(parmedPsf: CharmmPsfFile, config: dict) -> dict:
    """
    Save modified RTF, PRM, and PSF files.

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
    outRtf = p.join(assemblyDir, f"{moleculeName}_assembled.rtf")
    outPrm = p.join(assemblyDir, f"{moleculeName}_assembled.prm")
    outPsf = p.join(assemblyDir, f"{moleculeName}_assembled.psf")

    outputParams.write(top=outRtf, par=outPrm)
    parmedPsf.save(outPsf, overwrite=True)

    config["runtimeInfo"]["madeByAssembly"]["assembledRtf"] = outRtf 
    config["runtimeInfo"]["madeByAssembly"]["assembledPrm"] = outPrm
    config["runtimeInfo"]["madeByAssembly"]["assembledPsf"] = outPsf 

    return config

def assign_missing_dihedrals(missingParams: CharmmParameterSet,
                        parmedPsf: CharmmPsfFile,
                          nameToCgenffType: dict,
                            completeParameterSet: CharmmParameterSet) -> CharmmParameterSet:
    # Check for missing dihedral parameters
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
                raise ValueError(f"Warning: No CGenFF angle parameters for \n" \
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
