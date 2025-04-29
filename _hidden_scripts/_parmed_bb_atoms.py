import os
from os import path as p
import warnings
import parmed as pmd
from parmed.charmm import CharmmParameterSet, CharmmPsfFile
from parmed.exceptions import ParameterWarning
from psfgen import PsfGen

# Suppress ParameterWarning
warnings.filterwarnings('ignore', category=ParameterWarning)

# Placeholder classes (extend if needed)
class FilePath:
    pass

class DirectoryPath:
    pass

def get_dummy_inputs():
    """Return dummy input parameters for the script."""
    return {
        "strFile": "/home/esp/scriptDevelopment/drFrankenstein/Inputs/PHG_capped.str",
        "cappedPdb": "/home/esp/scriptDevelopment/drFrankenstein/06_PHG_outputs/01_termini_capping/geometry_optimisation/PHG_capped_opt.pdb",
        "moleculeName": "PHG",
        "outDir": "/home/esp/scriptDevelopment/drFrankenstein/_hidden_scripts/BB_types",
        "cgenffRtf": "/home/esp/scriptDevelopment/drFrankenstein/Laboratory/Experiments/Protocol_5_Stitching/CHARMM_protocols/top_all36_cgenff.rtf",
        "cgenffPrm": "/home/esp/scriptDevelopment/drFrankenstein/Laboratory/Experiments/Protocol_5_Stitching/CHARMM_protocols/par_all36_cgenff.prm",
        "charmmPrm": "/home/esp/scriptDevelopment/drFrankenstein/_hidden_scripts/BB_types/par_all36m_prot.prm",
        "charmmRtf": "/home/esp/scriptDevelopment/drFrankenstein/_hidden_scripts/BB_types/top_all36_prot.rtf",
        "chargeGroups": {
            "cTerminal": {"atoms": ["C", "O"], "charge": 0},
            "nTerminal": {"atoms": ["CA", "N"], "charge": 0},
            "CB": {"atoms": ["CB"], "charge": 0},
            "ring": {"atoms": ["CG1", "CG2", "CD1", "CD2", "CE"], "charge": 0}
        }
    }

def split_charmm_str(str_file: str, molecule_name: str, out_dir: str) -> tuple[str, str]:
    """
    Split a CHARMM stream file into RTF and PRM files.

    Args:
        str_file: Path to the input STR file.
        molecule_name: Name of the molecule for output file naming.
        out_dir: Output directory for RTF and PRM files.

    Returns:
        Tuple of paths to the generated RTF and PRM files.
    """
    moleculeRtf = p.join(out_dir, f"{molecule_name}.rtf")
    moleculePrm = p.join(out_dir, f"{molecule_name}.prm")

    if not p.exists(str_file):
        raise FileNotFoundError(f"STR file not found: {str_file}")

    with open(str_file, "r") as stream, open(moleculeRtf, "w") as rtf, open(moleculePrm, "w") as prm:
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

    return moleculeRtf, moleculePrm

def make_charmm_psf(molecule_rtf: str, cgenff_rtf: str, capped_pdb: str, out_dir: str, molecule_name: str) -> str:
    """
    Generate a PSF file using PsfGen from RTF and PDB files.

    Args:
        molecule_rtf: Path to the molecule's RTF file.
        cgenff_rtf: Path to the CGenFF RTF file.
        capped_pdb: Path to the input PDB file.
        out_dir: Output directory for the PSF file.
        molecule_name: Name of the molecule for output file naming.

    Returns:
        Path to the generated PSF file.
    """
    from pdbUtils.pdbUtils import pdb2df, df2pdb  # Assumed available

    # Standardize PDB with single residue ID and name
    cappedDf = pdb2df(capped_pdb)
    cappedDf["RES_ID"] = "1"
    cappedDf["RES_NAME"] = molecule_name
    tmpPdb = p.join(out_dir, f"united_capped_{molecule_name}.pdb")
    df2pdb(cappedDf, tmpPdb)

    # Generate PSF using PsfGen
    gen = PsfGen(output="/dev/null")
    gen.read_topology(molecule_rtf)
    gen.read_topology(cgenff_rtf)
    gen.add_segment(segid="A", pdbfile=tmpPdb)
    gen.read_coords(segid="A", filename=tmpPdb)

    psfFile = p.join(out_dir, f"{molecule_name}_capped.psf")
    gen.write_psf(psfFile)

    return psfFile

def get_cgenff_atom_types(psf: CharmmPsfFile, atom_type_map: dict) -> dict:
    """
    Get the original CGenFF atom types for atoms to be retyped.

    Args:
        psf: ParmEd CharmmPsfFile object with CGenFF types.
        atom_type_map: Dictionary mapping atom names to new CHARMM36m types.

    Returns:
        Dictionary mapping atom names to their original CGenFF types.
    """
    cgenffTypes = {}
    for atom in psf.atoms:
        if atom.name in atom_type_map:
            cgenffTypes[atom.name] = atom.type
    return cgenffTypes

def create_missing_prm(psf: CharmmPsfFile,
                        cgenff_params: CharmmParameterSet,
                          nameToOrigType: dict,
                          charmmDefaultMap: dict,
                            out_dir: str, moleculeName: str) -> str:
    
    """
    Create a PRM file with parameters for missing terms after atom type reassignment.

    Args:
        psf: ParmEd CharmmPsfFile object with updated atom types.
        cgenff_params: CharmmParameterSet with original CGenFF parameters.
        type_map: Dictionary mapping old CGenFF types to new CHARMM36m types.
        out_dir: Output directory for the PRM file.
        file_prefix: Prefix for the output PRM file name.

    Returns:
        Path to the generated PRM file.
    """
    print("#"*72)

    missingParams = CharmmParameterSet()



    # Check for missing bond parameters
    for bond in psf.bonds:
        atomNames = [bond.atom1.name, bond.atom2.name]
        newTypes = [bond.atom1.type, bond.atom2.type]
        origTypes = [nameToOrigType.get(name, atomType) for name, atomType in zip(atomNames, newTypes)] 
        print(atomNames, newTypes, origTypes)

        bondKey = tuple(sorted(newTypes))  # Bonds are order-independent
        if bondKey not in cgenff_params.bond_types and bondKey not in missingParams.bond_types:
            # Map to original CGenFF types
            origBondKey = tuple(sorted(origTypes))
            if origBondKey in cgenff_params.bond_types:
                bondType = cgenff_params.bond_types[origBondKey]
                missingParams.bond_types[bondKey] = bondType
                print(f"Added missing bond: {newTypes} (from {origTypes})")
            else:
                print(f"Warning: No CGenFF bond parameters for {origTypes}")
                exit()
    print("#"*72)

    # Check for missing angle parameters
    for angle in psf.angles:
        atomNames=[angle.atom1.name, angle.atom2.name, angle.atom3.name]
        newTypes = [angle.atom1.type, angle.atom2.type, angle.atom3.type]
        origTypes = [nameToOrigType.get(name, atomType) for name, atomType in zip(atomNames, newTypes)] 
        print(origTypes, newTypes,atomNames)
        angleKey = tuple(newTypes)
        if angleKey not in cgenff_params.angle_types and angleKey not in missingParams.angle_types:
            origAngleKey = tuple(origTypes)
            if origAngleKey in cgenff_params.angle_types:
                angleType = cgenff_params.angle_types[origAngleKey]
                missingParams.angle_types[angleKey] = angleType
                print(f"Added missing angle: {newTypes} (from {origTypes})")
            elif origAngleKey[::-1] in cgenff_params.angle_types:
                angleType = cgenff_params.angle_types[origAngleKey[::-1]]
                missingParams.angle_types[angleKey] = angleType
                print(f"Added missing angle: {newTypes} (from {origTypes})")
            else:
                print(f"Warning: No CGenFF angle parameters for {origTypes} -> {newTypes} -> {atomNames}")
                exit()
    print("#"*72)


    # Check for missing dihedral parameters
    for dihedral in psf.dihedrals:
        atomNames=[dihedral.atom1.name, dihedral.atom2.name, dihedral.atom3.name, dihedral.atom4.name]
        newTypes = [dihedral.atom1.type, dihedral.atom2.type, dihedral.atom3.type, dihedral.atom4.type]
        origTypes = [nameToOrigType.get(name, atomType) for name, atomType in zip(atomNames, newTypes)]

        dihedralKey = tuple(newTypes)
        if dihedralKey not in cgenff_params.dihedral_types and dihedralKey not in missingParams.dihedral_types:
            origDihedralKey = tuple(origTypes)
            origDihedralKeyRev = tuple(reversed(origDihedralKey))
            if origDihedralKey in cgenff_params.dihedral_types:
                dihedralTypes = cgenff_params.dihedral_types[origDihedralKey]
                missingParams.dihedral_types[dihedralKey] = dihedralTypes
                print(f"Added missing dihedral: {newTypes} (from {origTypes})")
            elif origDihedralKeyRev in cgenff_params.dihedral_types:
                dihedralTypes = cgenff_params.dihedral_types[origDihedralKeyRev]
                missingParams.dihedral_types[dihedralKey] = dihedralTypes
                print(f"Added missing dihedral REV: {newTypes} (from {origTypes})")
            else:
                print(f"Warning: No CGenFF dihedral parameters for {origTypes}: {newTypes} ")
                exit()
    print("#"*72)
    missingPrm = p.join(out_dir, f"{moleculeName}_missing.prm")

    # Write missing parameters to PRM file
    with open(missingPrm, "w") as f:
        f.write("* Missing parameters for retyped atoms\n")
        if missingParams.bond_types or missingParams.angle_types or missingParams.dihedral_types:
            missingParams.write(par=f.name)
        else:
            print("No missing parameters found; creating empty PRM file")

    print("#"*72)

    return missingPrm

def load_psf_and_parameters(psf_file: str, cgenff_rtf: str, cgenff_prm: str, molecule_rtf: str, molecule_prm: str, charmm_rtf: str, charmm_prm: str, missing_prm: str) -> tuple[CharmmPsfFile, CharmmParameterSet]:
    """
    Load PSF file and parameter sets into ParmEd.

    Args:
        psf_file: Path to the PSF file.
        cgenff_rtf: Path to the CGenFF RTF file.
        cgenff_prm: Path to the CGenFF PRM file.
        molecule_rtf: Path to the molecule's RTF file.
        molecule_prm: Path to the molecule PRM file.
        charmm_rtf: Path to the CHARMM36m RTF file.
        charmm_prm: Path to the CHARMM36m PRM file.
        missing_prm: Path to the PRM file with missing parameters.

    Returns:
        Tuple of the loaded PSF object and combined parameter set.
    """
    if not p.exists(psf_file):
        raise FileNotFoundError(f"PSF file not found: {psf_file}")

    # Load PSF
    psf = CharmmPsfFile(psf_file)

    # Load CGenFF parameters
    cgenffParams = CharmmParameterSet(cgenff_rtf, cgenff_prm, molecule_rtf, molecule_prm)
    psf.load_parameters(cgenffParams)

    # Load all parameters (CGenFF + CHARMM36m + missing)
    allParams = CharmmParameterSet(cgenff_rtf, cgenff_prm, molecule_rtf, molecule_prm, charmm_rtf, charmm_prm, missing_prm)
    return psf, allParams

def update_atom_types(psf: CharmmPsfFile, atom_type_map: dict) -> None:
    """
    Update atom types in the PSF structure based on a mapping.

    Args:
        psf: ParmEd CharmmPsfFile object.
        atom_type_map: Dictionary mapping atom names to new CHARMM36m types.
    """
    for atom in psf.atoms:
        if atom.name in atom_type_map:
            oldType = atom.type
            atom.type = atom_type_map[atom.name]
            print(f"Changed atom {atom.name} type from {oldType} to {atom.type}")

def verify_atom_types(psf: CharmmPsfFile, atom_type_map: dict) -> None:
    """
    Verify updated atom types and print their charges.

    Args:
        psf: ParmEd CharmmPsfFile object.
        atom_type_map: Dictionary of atom names to verify.
    """
    for atom in psf.atoms:
        if atom.name in atom_type_map:
            print(f"Atom {atom.name}: Type={atom.type}, Charge={atom.charge}")


def save_modified_files(psf: CharmmPsfFile, out_dir: str, molecule_name: str) -> None:
    """
    Save modified RTF, PRM, and PSF files.

    Args:
        psf: ParmEd CharmmPsfFile object with updated atom types.
        out_dir: Output directory for the files.
        molecule_name: Name of the molecule for output file naming.
    """
    outputParams = CharmmParameterSet.from_structure(psf)
    outRtf = p.join(out_dir, f"{molecule_name}_modified.rtf")
    outPrm = p.join(out_dir, f"{molecule_name}_modified.prm")
    outPsf = p.join(out_dir, f"{molecule_name}_modified.psf")

    outputParams.write(top=outRtf, par=outPrm)
    psf.save(outPsf, overwrite=True)

def main():
    """Main function to orchestrate the workflow."""

    # Load inputs
    inputs = get_dummy_inputs()

    # Create output directory
    os.makedirs(inputs["outDir"], exist_ok=True)

    # Split STR file into RTF and PRM
    moleculeRtf, moleculePrm = split_charmm_str(inputs["strFile"], inputs["moleculeName"], inputs["outDir"])

    # Generate PSF file
    moleculePsf = make_charmm_psf(moleculeRtf, inputs["cgenffRtf"], inputs["cappedPdb"], inputs["outDir"], inputs["moleculeName"])

    # Load PSF with CGenFF parameters to get original types
    psf = CharmmPsfFile(moleculePsf)
    cgenffParams = CharmmParameterSet(inputs["cgenffRtf"], inputs["cgenffPrm"], moleculeRtf, moleculePrm)
    psf.load_parameters(cgenffParams)

    allParams = CharmmParameterSet(inputs["cgenffRtf"], inputs["cgenffPrm"], moleculeRtf, moleculePrm, inputs["charmmRtf"], inputs["charmmPrm"])


    # Define CHARMM36m atom type mapping
    charmmDefaultMap = {
                    "N" : "NH1",
                    "HN": "H",
                    "CA": "CT1",
                    "HA": "HB1",
                    "C" : "C",
                    "O" : "O",
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
                    "HC3": "HB1"
    }

    # Get original CGenFF types for atoms to be retyped
    nameToOrigType = get_cgenff_atom_types(psf, charmmDefaultMap)


    # Update atom types to identify missing parameters
    update_atom_types(psf, charmmDefaultMap)

    # Create PRM file for missing parameters
    missingPrm = create_missing_prm(psf, allParams, nameToOrigType, charmmDefaultMap,  inputs["outDir"], inputs["moleculeName"]) 
    print(f"Created missing PRM file: {missingPrm}")


    # Load PSF and all parameters, including missing PRM
    psf, allParams = load_psf_and_parameters(
        moleculePsf, inputs["cgenffRtf"], inputs["cgenffPrm"], moleculeRtf, moleculePrm, inputs["charmmRtf"], inputs["charmmPrm"], missingPrm
    )
    print(f"Loaded PSF and parameters")

    # Update atom types again (since PSF was reloaded)
    update_atom_types(psf, charmmDefaultMap)
    psf.load_parameters(allParams)
    verify_atom_types(psf, charmmDefaultMap)

    # Save modified files
    save_modified_files(psf, inputs["outDir"], inputs["moleculeName"])
    print(f"Saved modified files in {inputs['outDir']}")

if __name__ == "__main__":
    main()