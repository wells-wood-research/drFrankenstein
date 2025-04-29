import os
from os import path as p
import warnings
import parmed as pmd
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

    return parmedPsf, nameToCgenffType



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