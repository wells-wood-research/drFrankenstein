import os
from os import path as p
from subprocess import call
## CLEAN CODE ##
class FilePath:
    pass
class DirectoryPath:
    pass
from OperatingTools import file_parsers
from typing import List, Dict, Optional, Any # Added for type hinting

########################################
def run_orca(orca_input: FilePath, orca_output: FilePath, config: dict) -> None:
    """
    Runs an ORCA calculation.

    Args:
        orca_input: Path to the ORCA input file.
        orca_output: Path to write the ORCA output.
        config: Configuration dictionary containing the path to the ORCA executable.
    """
    orcaExe: str = config["pathInfo"]["orcaExe"]
    orcaCommand: List[str] = [orcaExe, orca_input] # type: ignore
    with open(orca_output, 'w') as output_file: # type: ignore
        try:
            call(orcaCommand, stdout=output_file, stderr=output_file)
        except Exception as e:
            raise(e)
########################################
def make_orca_input_qmmm_opt(input_xyz: FilePath,
                                out_dir: DirectoryPath,
                                molecule_info: dict,
                                qm_method: str,
                                qm_atoms: str,
                                parameter_file: FilePath,
                                n_cpus: int) -> FilePath:
    """
    Creates an ORCA input file for QM/MM geometry optimization.

    Args:
        input_xyz: Path to the input XYZ file.
        out_dir: Directory to write the output input file.
        molecule_info: Dictionary containing charge and multiplicity.
        qm_method: QM method string.
        qm_atoms: String specifying QM atoms for ORCA.
        parameter_file: Path to the ORCAFF parameter file.
        n_cpus: Number of CPUs for the calculation.

    Returns:
        Path to the generated ORCA input file.
    """
    ## write QM/MM opt orca input file
    charge: int = molecule_info["charge"]
    multiplicity: int = molecule_info["multiplicity"]
    qmmmOptOrcaInput: FilePath = p.join(out_dir, f"QMMM_orca_opt.inp") # type: ignore
    with open(qmmmOptOrcaInput, "w") as f:
        f.write(" # --------------------------------- #\n")
        f.write(" #  QM/MM Optimisation               #\n")
        f.write(" # --------------------------------- #\n")  
        ## simpleinput line for method and QMMM flag      
        f.write(f"! QMMM {qm_method} Opt\n")
        f.write(f"%pal nprocs {n_cpus}\nend\n")
        ## qmmmm options
        f.write("%qmmm\n")
        f.write(f"{' '*4}QMAtoms {qm_atoms} end\n")
        f.write(f"{' '*4}ORCAFFFilename \"{parameter_file}\"\n")
        f.write("Rigid_MM_Water TRUE\n")
        f.write("end\n")
        f.write(f"*xyzfile {charge} {multiplicity} {input_xyz}\n")
    return qmmmOptOrcaInput



########################################

def make_orca_input_qmmm_singlepoint(input_xyz: FilePath,
                                out_dir: DirectoryPath,
                                molecule_info: dict,
                                qm_method: str,
                                qm_atoms: str,
                                parameter_file: FilePath) -> FilePath:
    """
    Creates an ORCA input file for QM/MM single point energy calculation.

    Args:
        input_xyz: Path to the input XYZ file.
        out_dir: Directory to write the output input file.
        molecule_info: Dictionary containing charge and multiplicity.
        qm_method: QM method string.
        qm_atoms: String specifying QM atoms for ORCA.
        parameter_file: Path to the ORCAFF parameter file.

    Returns:
        Path to the generated ORCA input file.
    """
    ## write QM/MM single-point orca input file
    charge: int = molecule_info["charge"]
    multiplicity: int = molecule_info["multiplicity"]
    qmmmSinglepointOrcaInput: FilePath = p.join(out_dir, f"QMMM_orca_sp.inp") # type: ignore
    with open(qmmmSinglepointOrcaInput, "w") as f:
        f.write(" # --------------------------------- #\n")
        f.write(" #  QM/MM Single-Point Energy        #\n")
        f.write(" # --------------------------------- #\n")  
        ## simpleinput line for method and QMMM flag      
        f.write(f"! QMMM {qm_method} SP\n")
        ## qmmmm options
        f.write("%qmmm\n")
        f.write(f"{' '*4}QMAtoms {qm_atoms} end\n")
        f.write(f"{' '*4}ORCAFFFilename \"{parameter_file}\"\n")
        f.write("end\n")
        f.write(f"*xyzfile {charge} {multiplicity} {input_xyz}\n")
    return qmmmSinglepointOrcaInput

########################################
def make_orca_input_for_solvator(input_xyz: FilePath,
                                  out_dir: DirectoryPath,
                                    molecule_info: dict,
                                      qm_method: str,
                                        solvation_method: str,
                                          n_waters: int) -> FilePath:
    """
    Creates an ORCA input file to run the SOLVATOR protocol.
    """
    ## write GOAT orca input file
    charge: int = molecule_info["charge"]
    multiplicity: int = molecule_info["multiplicity"]
    solvatorOrcaInput: FilePath = p.join(out_dir, f"SOLVATOR_orca.inp") # type: ignore
    with open(solvatorOrcaInput, "w") as f:
        f.write(" # --------------------------------- #\n")
        f.write(" #  SOLVATOR solvent addition        #\n")
        f.write(" # --------------------------------- #\n")
        ## method and solvation
        f.write(f"!{qm_method} {solvation_method}\n")
        ## solvator params
        f.write(f"%SOLVATOR\n")
        f.write(f"{' '*4}NSOLV {n_waters}\n")
        f.write(f"{' '*4}CLUSTERMODE STOCHASTIC\n")
        f.write("END\n")
        ## geom, charge, multiplicity
        f.write(f"*xyzfile {charge} {multiplicity} {input_xyz}\n")
    return solvatorOrcaInput


########################################



def write_goat_input(conformer_dir: DirectoryPath, capped_xyz: FilePath, config: dict) -> FilePath:
    """
    Writes an ORCA input file for GOAT conformer generation.

    Args:
        conformer_dir: Directory to write the GOAT input file.
        capped_xyz: Path to the capped XYZ file for the molecule.
        config: Configuration dictionary containing molecule info and other settings.

    Returns:
        Path to the generated ORCA GOAT input file.
    """
    ## write GOAT orca input file
    charge: int = config["moleculeInfo"]["charge"]
    multiplicity: int = config["moleculeInfo"]["multiplicity"]
    goatOrcaInput: FilePath = p.join(conformer_dir, f"GOAT_orca.inp") # type: ignore
    with open(goatOrcaInput, "w") as f:
        f.write(" # --------------------------------- #\n")
        f.write(" #  GOAT conformation generation     #\n")
        f.write(" # --------------------------------- #\n")
        f.write("!XTB2 GOAT\n")
        f.write("%PAL NPROCS 16 END\n")  ##TODO: work out how to use more cores or have a stable way of inputting this
        f.write("%GOAT\n")
        f.write("\tMAXITERMULT 1\n")            ## only do one round (save some time)
        f.write("\tFREEZEAMIDES TRUE\n")           ## freeze amides
        f.write("\tFREEZECISTRANS TRUE\n ")         ## freeze cis-trans bonds
        f.write("END\n")
        f.write(f"*xyzfile {charge} {multiplicity} {capped_xyz}\n")
    return goatOrcaInput

########################################
def make_orca_input_for_opt(input_xyz: FilePath, 
                             out_dir: DirectoryPath, 
                             molecule_info: dict, 
                             qm_method: str, 
                             solvation_method: Optional[str]) -> FilePath:
    """Creates an ORCA input file for geometry optimization."""
    ## unpack moleculeInfo
    charge: int = molecule_info["charge"]
    multiplicity: int = molecule_info["multiplicity"]
    ## create orca input file
    orcaInputFile: FilePath = p.join(out_dir, "orca_opt.inp") # type: ignore
    with open(orcaInputFile, "w") as f:
        f.write(" # --------------------------------- #\n")
        f.write(" #  Geometry Optimisation            #\n")
        f.write(" # --------------------------------- #\n")
        ## METHOD
        if solvation_method is None:
            f.write(f"! {qm_method} Opt\n")
        else:
            f.write(f"! {qm_method} {solvation_method} Opt\n") 
        ## GEOMETRY
        f.write(f"*xyzfile {charge} {multiplicity} {input_xyz}\n\n")

    return orcaInputFile

###################################################################################
def  make_orca_input_for_scan(input_xyz: FilePath,
                               out_dir: DirectoryPath,
                                 molecule_info: dict,
                                   qm_method: str,
                                     solvation_method: Optional[str],
                                       torsion_indexes: List[int],
                                         scan_angles: str) -> FilePath:
    """Creates an ORCA input file for a torsion scan."""
    ## unpack moleculeInfo
    charge: int = molecule_info["charge"]
    multiplicity: int = molecule_info["multiplicity"]
    ## create orca input file
    orcaInputFile: FilePath = p.join(out_dir, "orca_scan.inp") # type: ignore
    with open(orcaInputFile, "w") as f:
        f.write(" # --------------------------------- #\n")
        f.write(" #  Torsion Scan                     #\n")
        f.write(" # --------------------------------- #\n")
        ## METHOD
        if solvation_method is None:
            f.write(f"! {qm_method} Opt\n")
        else:
            f.write(f"! {qm_method} {solvation_method} Opt\n") 
        ## SCAN 
        if not scan_angles is None: # scan_angles can be an empty string
            f.write("%geom Scan\n")
            torsionText: str = f"D {' '.join(map(str, torsion_indexes))} = {scan_angles}\n"
            f.write(torsionText)
            f.write("end\n")
            f.write("end\n") # This was present in original, assuming it's for %geom Scan block

        ## GEOMETRY
        f.write(f"*xyzfile {charge} {multiplicity} {input_xyz}\n\n")

    return orcaInputFile
###################################################################################
def  make_orca_input_for_singlepoint(input_xyz: FilePath, 
                                      out_dir: DirectoryPath, 
                                      molecule_info: dict, 
                                      qm_method: str, 
                                      solvation_method: Optional[str]) -> FilePath:
    """Creates an ORCA input file for a single point energy calculation."""
    ## unpack moleculeInfo
    charge: int = molecule_info["charge"]
    multiplicity: int = molecule_info["multiplicity"]
    ## create orca input file
    orcaInputFile: FilePath = p.join(out_dir, "orca_sp.inp") # type: ignore
    with open(orcaInputFile, "w") as f:
        f.write(" # --------------------------------- #\n")
        f.write(" #  Single Point Calculation         #\n")
        f.write(" # --------------------------------- #\n")
        ## METHOD
        if solvation_method is None:
            f.write(f"! {qm_method} SP\n")
        else:
            f.write(f"! {qm_method} {solvation_method} SP\n") 
        ## GEOMETRY
        f.write(f"*xyzfile {charge} {multiplicity} {input_xyz}\n\n")

    return orcaInputFile


###################################################################################
def did_orca_finish_normallly(orca_out: FilePath) -> bool:
    """Checks if an ORCA calculation terminated normally by reading its output file."""
    with open(orcaOut, "r") as f:
        lines = f.readlines()
        for line in reversed(lines):
            if "****ORCA TERMINATED NORMALLY****" in line:
                return True
            elif "ORCA finished with an error in the energy calculation" in line:
                return False
    return False
###################################################################################