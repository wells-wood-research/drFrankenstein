import os
from os import path as p
import warnings
import parmed
from parmed.charmm import CharmmParameterSet, CharmmPsfFile
from parmed.exceptions import ParameterWarning
from psfgen import PsfGen
from pdbUtils import pdbUtils 

# Suppress ParameterWarning
warnings.filterwarnings('ignore', category=ParameterWarning)

from . import Assembly_Assistant

from typing import Tuple
class FilePath:
    pass
class DirectoryPath:
    pass

def replace_parameters(moleculeFrcmod: FilePath,
                        parameterDataset: FilePath) -> None:
    """
    Reads though a FRCMOD file and replaces parameters with those
    from another FRCMOD / DAT file
    Overwrites original file

    Args:
        moleculeFrcmod (FilePath): path to original FRCMOD file
        paraneterDataset (FilePath): path to FRCMOD / DAT file
    Returns:
        None
    """

    ## load parameter files into parmed
    molParams = parmed.load_file(moleculeFrcmod)
    ffParams = parmed.amber.AmberParameterSet(parameterDataset)

    # Replace BOND parameters
    for bondKey in molParams.bond_types:
        if bondKey in ffParams.bond_types:
            bondKey = ffParams.bond_types[bondKey]

    # Replace ANGLE parameters
    for angleKey in molParams.angle_types:
        if angleKey in ffParams.angle_types:
            angleKey = ffParams.angle_types[angleKey]

    # Replace DIHEDRAL parameters
    for dihedralKey in molParams.dihedral_types:
        wildCardDihedralKey = ("X", dihedralKey[1], dihedralKey[2], "X")
        if dihedralKey in ffParams.dihedral_types:
            dihedralKey = ffParams.dihedral_types[dihedralKey]
        elif wildCardDihedralKey in ffParams.dihedral_types:
            dihedralKey = ffParams.dihedral_types[wildCardDihedralKey]
            
    ## hard-coded changes to PHI and PSI using "CX" wildcard ##
    ## manual change of PHI
    molParams.dihedral_types[("C", "N", "CA", "C")] = ffParams.dihedral_types[("C", "N", "CX", "C")] 
    ## manual change of PSI
    molParams.dihedral_types[("N", "CA", "C", "N")] = ffParams.dihedral_types[("N", "CX", "C", "N")] 


    for improperKey in molParams.improper_types:
        if improperKey in ffParams.improper_types:
            improperKey = ffParams.improper_types[improperKey]
    
    molParams.write(moleculeFrcmod, style="frcmod")
    return None

def change_backbone_types_amber(mol2File: FilePath, config: dict) -> None:
    """
    Uses a map to change backbone atom types from auto-assigned gaff2
    to protein-style amber types

    Args:
        mol2File (FilePath): path to mol2 file
        config (dict): config dict

    Returns:
        None

    """

    # Define CHARMM36m atom type mapping
    amberDefaultMap = {
                    "N" : "N",
                    "HN": "H",
                    "CA": "CT",
                    "HA": "H1",
                    "C" : "C",
                    "O" : "O",
                    "NN": "N",
                    "HNN1": "H",
                    "CN": "CT", 
                    "HCN1": "H3",
                    "HCN2": "H3",
                    "HCN3": "H3",
                    "CC1": "C",
                    "OC": "O",
                    "CC2": "CT", 
                    "HC1": "H3",
                    "HC2": "H3",
                    "HC3": "H3"
    }

    parmedStructure = parmed.load_file(mol2File)
    for atom in parmedStructure.atoms:
        if atom.name in amberDefaultMap:
            atom.type = amberDefaultMap[atom.name]

    parmedStructure.save(mol2File, overwrite=True)





def create_missing_prm(parmedPsf: CharmmPsfFile,
                       cappedRtf: FilePath,
                       cappedPrm: FilePath,
                        charmmDefaultParams: dict,
                        nameToCgenffType: dict,
                            config: dict) -> FilePath:
    """
    Create a PRM file with parameters for missing terms after atom type reassignment.

    Args:
        parmedPsf: (CharmmPsfFile) object with updated atom types.
        cappedRtf: (FilePath) to the capped RTF file.
        cappedPrm: (FilePath) to the capped PRM file.
        charmmDefaultParams: (dict) with default CHARMM parameters.
        nameToCgenffType: (dict) mapping atom names to original CGenFF types.
        config: (dict) with configuration information.

    Returns:
        missingPrm (FilePath): PRM file containing missing parameters.
    """
    ## unpack config
    moleculeName = config["moleculeInfo"]["moleculeName"]
    assemblyDir = config["runtimeInfo"]["madeByAssembly"]["assemblyDir"]
    ## make a param set using default CHARMM parameters, CGenFF and by-analogy params 
    completeParameterSet = CharmmParameterSet(cappedRtf, cappedPrm, *charmmDefaultParams.values())

    ## init an empty CharmmParameterSet object
    missingParams = CharmmParameterSet()

    missingParams = Assembly_Assistant.assign_missing_bonds(missingParams,
                                                             parmedPsf,
                                                               nameToCgenffType,
                                                                 completeParameterSet)
                
    missingParams = Assembly_Assistant.assign_missing_angles(missingParams,
                                                             parmedPsf,
                                                               nameToCgenffType,
                                                                 completeParameterSet)
    missingParams = Assembly_Assistant.assign_missing_dihedrals(missingParams,
                                                             parmedPsf,
                                                               nameToCgenffType,
                                                                 completeParameterSet)



    missingPrm = p.join(assemblyDir, f"{moleculeName}_missing.prm")
    # Write missing parameters to PRM file
    with open(missingPrm, "w") as f:
        if missingParams.bond_types or missingParams.angle_types or missingParams.dihedral_types:
            missingParams.write(par=f.name)
        else:
            pass

    return missingPrm

def set_backbone_types_psf(parmedPsf: CharmmPsfFile,
                            config: dict) -> Tuple[CharmmPsfFile, dict]:
    """
    Modifies BACKBONE atom's ATOM_TYPE in PSF file.
    Replaces exotic CGenFF ATOM_TYPES with standard CHARMM backbone ATOM_TYPES.
    eg. [NG2S1 -> NH1 (N), CG311 -> CT1 (CA)]

    Args:
        parmedPsf (CharmmPsfFile): ParmEd CharmmPsfFile object with CGenFF types.
        config (dict): A dictionary containing configuration information.

    Returns:
        parmedPsf (CharmmPsfFile): ParmEd CharmmPsfFile object with CHARMM36m types.
        nameToCgenffType (dict): Dictionary mapping atom names to their original CGenFF types.

    """
    # Define CHARMM36m atom type mapping
    ## TODO: Automate using config input
    nameToDesiredType_mol = {
                    "N" : "NH1",
                    "HN": "H",
                    "CA": "CT1",
                    "HA": "HB1",
                    "C" : "C",
                    "O" : "O" }

    nameToDesiredType_bb = {                    
                    "NN": "NH1",
                    "HNN1": "H",
                    "CN": "CT1", 
                    "HCN1": "HB1",
                    "HCN2": "HB1",
                    "HCN3": "HB1",
                    "CC1": "C",
                    "OC": "O",
                    "CC2": "CT1", 
                    "HC1": "HB1",
                    "HC2": "HB1",
                    "HC3": "HB1"}
    
    nameToDesiredType = {**nameToDesiredType_bb, **nameToDesiredType_mol}
    nameToCgenffType = Assembly_Assistant.get_cgenff_atom_types(parmedPsf, nameToDesiredType)

    parmedPsf = Assembly_Assistant.update_psf_atom_types(parmedPsf, nameToDesiredType)

    return parmedPsf, nameToCgenffType, nameToDesiredType



def split_charmm_str(config: dict) -> tuple[str, str]:
    """
    Split a CHARMM stream file into RTF and PRM files.

    Args:
        config (dict): A dictionary containing configuration information.

    Returns:
        Tuple of paths to the generated RTF and PRM files.
    """
    ## unpack config 
    outDir = config["runtimeInfo"]["madeByAssembly"]["assemblyDir"]
    moleculeName = config["moleculeInfo"]["moleculeName"]
    strFile = config["runtimeInfo"]["madeByCapping"]["cappedMoleculeStr"] 

    rtfFile = p.join(outDir, f"{moleculeName}_capped.rtf")
    prmFile = p.join(outDir, f"{moleculeName}_capped.prm")

    if not p.exists(strFile):
        raise FileNotFoundError(f"STR file not found: {strFile}")

    with open(strFile, "r") as stream, open(rtfFile, "w") as rtf, open(prmFile, "w") as prm:
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

    return rtfFile, prmFile


def make_charmm_psf(moleculeRtf: FilePath, cgenffRtf: FilePath, config: dict) -> str:
    """
    Generate a PSF file using PsfGen from RTF and PDB files.

    Args:
        moleculeRtf (FilePath): Residue Topology File for capped molecule
        cgenffRtf (FilePath): Residue Topology File for CGenFF
        config (dict) contains all run info
    Returns:
        Path to the generated PSF file.
    """
    ## unpack config
    cappedPdb = config["runtimeInfo"]["madeByCapping"]["cappedPdb"]
    moleculeName = config["moleculeInfo"]["moleculeName"]
    outDir = config["runtimeInfo"]["madeByAssembly"]["assemblyDir"]

    # Standardize PDB with single residue ID and name
    cappedDf = pdbUtils.pdb2df(cappedPdb)
    cappedDf["RES_ID"] = "1"
    cappedDf["RES_NAME"] = moleculeName
    tmpPdb = p.join(outDir, f"united_capped_{moleculeName}.pdb")
    pdbUtils.df2pdb(cappedDf, tmpPdb)

    # Generate PSF using PsfGen
    gen = PsfGen(output="/dev/null")
    gen.read_topology(moleculeRtf)
    gen.read_topology(cgenffRtf)
    gen.add_segment(segid="A", pdbfile=tmpPdb)
    gen.read_coords(segid="A", filename=tmpPdb)

    psfFile = p.join(outDir, f"{moleculeName}_capped.psf")
    gen.write_psf(psfFile)

    return psfFile