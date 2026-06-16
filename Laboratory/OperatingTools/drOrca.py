import os
from os import path as p
from subprocess import call
## CLEAN CODE ##
class FilePath:
    pass
class DirectoryPath:
    pass
from OperatingTools import file_parsers
from OperatingTools import drLogger
########################################
def run_orca(orcaInput: FilePath, orcaOutput: FilePath, config: dict) -> None:
    """Run ORCA with the configured executable and capture stdout/stderr to a file."""
    orcaExe = config["pathInfo"]["orcaExe"]
    orcaCommand = [orcaExe, orcaInput]
    logger = drLogger.get_logger()
    with open(orcaOutput, 'w') as outputFile:
        try:
            result = call(orcaCommand, stdout=outputFile, stderr=outputFile)
            if logger:
                logger.log_subprocess_call(orcaCommand, result)
        except Exception as e:
            if logger:
                logger.log_subprocess_error(orcaCommand, e)
            raise(e)
########################################

def make_orca_input_for_constrained_opt(inputXyz: FilePath,
                             outDir: DirectoryPath,
                               moleculeInfo: dict,
                                 qmMethod: str,
                                   solvationMethod: str,
                                     geomOptions:str = None,
                                     constrainIndexes: list = None) -> FilePath:
    """Write an ORCA constrained optimisation input file."""
    ## unpack moleculeInfo
    charge = moleculeInfo["charge"]
    multiplicity = moleculeInfo["multiplicity"]
    ## create orca input file
    orcaInputFile = p.join(outDir, "orca_opt.inp")
    with open(orcaInputFile, "w") as f:
        f.write(" # --------------------------------- #\n")
        f.write(" #  Geometry Optimisation            #\n")
        f.write(" # --------------------------------- #\n")
        ## METHOD
        if solvationMethod is None:
            f.write(f"! {qmMethod} Opt\n")
        else:
            f.write(f"! {qmMethod} {solvationMethod} Opt\n") 

        if not geomOptions is None:
            f.write(f"%geom\n{geomOptions}\nend\n")
        ## GEOMETRY
        f.write(f"*xyzfile {charge} {multiplicity} {inputXyz}\n\n")

    return orcaInputFile

########################################
def make_orca_input_qmmm_opt(inputXyz: FilePath,
                                outDir: DirectoryPath,
                                moleculeInfo: dict,
                                qmMethod: str,
                                qmAtoms: str,
                                parameterFile: FilePath,
                                nCpus: int) -> FilePath:
    """Write an ORCA QM/MM optimisation input file."""

    ## write QM/MM opt orca input file
    charge = moleculeInfo["charge"]
    multiplicity = moleculeInfo["multiplicity"]
    qmmmOptOrcaInput = p.join(outDir, f"QMMM_orca_opt.inp")
    with open(qmmmOptOrcaInput, "w") as f:
        f.write(" # --------------------------------- #\n")
        f.write(" #  QM/MM Optimisation               #\n")
        f.write(" # --------------------------------- #\n")  
        ## simpleinput line for method and QMMM flag      
        f.write(f"! QMMM {qmMethod} L-Opt\n")
        f.write(f"%pal nprocs {nCpus}\nend\n")
        ## qmmmm options
        f.write("%qmmm\n")
        f.write(f"{' '*4}QMAtoms {qmAtoms} end\n")
        f.write(f"{' '*4}ORCAFFFilename \"{parameterFile}\"\n")
        f.write("Rigid_MM_Water TRUE\n")
        f.write("end\n")
        f.write("%geom\n  maxIter 500\nend\n")
        f.write(f"*xyzfile {charge} {multiplicity} {inputXyz}\n")
    return qmmmOptOrcaInput





########################################

def make_orca_input_qmmm_singlepoint(inputXyz: FilePath,
                                outDir: DirectoryPath,
                                moleculeInfo: dict,
                                qmMethod: str,
                                qmAtoms: str,
                                parameterFile: FilePath,
                                nCpus: int | None = None) -> FilePath:
    """Write an ORCA QM/MM single-point input file."""
    ## write QM/MM single-point orca input file
    charge = moleculeInfo["charge"]
    multiplicity = moleculeInfo["multiplicity"]
    qmmmSinglepointOrcaInput = p.join(outDir, f"QMMM_orca_sp.inp")
    with open(qmmmSinglepointOrcaInput, "w") as f:
        f.write(" # --------------------------------- #\n")
        f.write(" #  QM/MM Single-Point Energy        #\n")
        f.write(" # --------------------------------- #\n")  
        ## simpleinput line for method and QMMM flag      
        f.write(f"! QMMM {qmMethod} SP\n")
        if nCpus is not None:
            f.write(f"%pal nprocs {nCpus}\nend\n")
        ## qmmmm options
        f.write("%qmmm\n")
        f.write(f"{' '*4}QMAtoms {qmAtoms} end\n")
        f.write(f"{' '*4}ORCAFFFilename \"{parameterFile}\"\n")
        f.write("end\n")
        f.write(f"*xyzfile {charge} {multiplicity} {inputXyz}\n")
    return qmmmSinglepointOrcaInput

########################################
def make_orca_input_for_solvator(inputXyz: FilePath,
                                  outDir: DirectoryPath,
                                    moleculeInfo: dict,
                                      qmMethod: str,
                                        solvationMethod: str,
                                          nWaters: int,
                                            nCpus: int | None = None) -> FilePath:
    """Write an ORCA SOLVATOR input file."""
    ## write GOAT orca input file
    charge = moleculeInfo["charge"]
    multiplicity = moleculeInfo["multiplicity"]
    solvatorOrcaInput = p.join(outDir, f"SOLVATOR_orca.inp")
    with open(solvatorOrcaInput, "w") as f:
        f.write(" # --------------------------------- #\n")
        f.write(" #  SOLVATOR solvent addition        #\n")
        f.write(" # --------------------------------- #\n")
        ## method and solvation
        f.write(f"!{qmMethod} {solvationMethod}\n")
        if nCpus is not None:
            f.write(f"%pal nprocs {nCpus}\nend\n")
        ## solvator params
        f.write(f"%SOLVATOR\n")
        f.write(f"{' '*4}NSOLV {nWaters}\n")
        f.write(f"{' '*4}CLUSTERMODE STOCHASTIC\n")
        f.write("END\n")
        ## geom, charge, multiplicity
        f.write(f"*xyzfile {charge} {multiplicity} {inputXyz}\n")
    return solvatorOrcaInput


########################################



def write_goat_input(conformerDir: DirectoryPath, cappedXyz: FilePath, config: dict) -> FilePath:
    """Write an ORCA GOAT conformer-generation input file."""
    ## write GOAT orca input file
    charge = config["moleculeInfo"]["charge"]
    multiplicity = config["moleculeInfo"]["multiplicity"]
    availableCpus = config["miscInfo"]["availableCpus"]

    conformerGenerationInfo = config["conformerGenerationInfo"]
    goatMode = conformerGenerationInfo["goatMode"] ## GOAT-ENTROPY or GOAT
    energyCutoff = conformerGenerationInfo["energyCutoff"] ## energy cutoff for accepting conformers (kcal/mol)
    goatMethod = conformerGenerationInfo["goatMethod"] ## XTB2 or GFN-FF

    goatOrcaInput = p.join(conformerDir, f"GOAT_orca.inp")

    useCpus = min(availableCpus, 16) ## limit to 16 cpus for GOAT to avoid overloading the system
    with open(goatOrcaInput, "w") as f:
        f.write(" # --------------------------------- #\n")
        f.write(" #  GOAT conformation generation     #\n")
        f.write(" # --------------------------------- #\n")
        f.write(f"!{goatMethod} {goatMode}\n")
        f.write(f"%PAL NPROCS {useCpus} END\n")  
        f.write("%GOAT\n")
        f.write("\tMAXITERMULT 1\n")            ## only do one round (save some time)
        f.write("\tFREEZEAMIDES TRUE\n")         ## freeze amides
        f.write("\tFREEZECISTRANS TRUE\n ")      ## freeze cis-trans bonds
        f.write((f"\tMAXEN {str(energyCutoff)}\n"))   ## only accept conformers within {energyCutoff} kcal/mol of the lowest energy
        f.write("END\n")
        ## constrain bond lengths to prevent dissociation during conformer generation
        f.write("%geom Constraints\n")
        f.write("{ B * * C }\n")
        f.write("end\n")
        f.write("end\n")
        f.write(f"*xyzfile {charge} {multiplicity} {cappedXyz}\n")
    return goatOrcaInput

########################################
def make_orca_input_for_opt(inputXyz: FilePath,
                             outDir: DirectoryPath,
                               moleculeInfo: dict,
                                 qmMethod: str,
                                   solvationMethod: str,
                                     geomOptions:str = None,
                                       nCpus: int | None = None) -> FilePath:
    """Write an ORCA geometry optimisation input file."""
    ## unpack moleculeInfo
    charge = moleculeInfo["charge"]
    multiplicity = moleculeInfo["multiplicity"]
    ## create orca input file
    orcaInputFile = p.join(outDir, "orca_opt.inp")
    with open(orcaInputFile, "w") as f:
        f.write(" # --------------------------------- #\n")
        f.write(" #  Geometry Optimisation            #\n")
        f.write(" # --------------------------------- #\n")
        ## METHOD
        if solvationMethod is None:
            f.write(f"! {qmMethod} Opt\n")
        else:
            f.write(f"! {qmMethod} {solvationMethod} Opt\n") 
        if nCpus is not None:
            f.write(f"%pal nprocs {nCpus}\nend\n")

        if not geomOptions is None:
            f.write(f"%geom\n{geomOptions}\nend\n")
        ## GEOMETRY
        f.write(f"*xyzfile {charge} {multiplicity} {inputXyz}\n\n")

    return orcaInputFile

###################################################################################
def  make_orca_input_for_scan(inputXyz,
                               outDir,
                                 moleculeInfo,
                                   qmMethod,
                                     solvationMethod,
                                       torsionIndexes,
                                         scanAngles,
                                           nCpus: int | None = None):
    ## unpack moleculeInfo
    charge = moleculeInfo["charge"]
    multiplicity = moleculeInfo["multiplicity"]
    ## create orca input file
    orcaInputFile = p.join(outDir, "orca_scan.inp")
    with open(orcaInputFile, "w") as f:
        f.write(" # --------------------------------- #\n")
        f.write(" #  Torsion Scan                     #\n")
        f.write(" # --------------------------------- #\n")
        ## METHOD
        if solvationMethod is None:
            f.write(f"! {qmMethod} Opt\n")
        else:
            f.write(f"! {qmMethod} {solvationMethod} Opt\n") 
        if nCpus is not None:
            f.write(f"%pal nprocs {nCpus}\nend\n")
        ## SCAN 
        if not scanAngles is None:
            f.write("%geom Scan\n")
            torsionText: str = f"D {' '.join(map(str, torsionIndexes))} = {scanAngles}\n"
            f.write(torsionText)
            f.write("end\n")
            f.write("end\n")

        ## GEOMETRY
        f.write(f"*xyzfile {charge} {multiplicity} {inputXyz}\n\n")

    return orcaInputFile
###################################################################################
def make_orca_input_for_sp(inputXyz: FilePath, outDir: DirectoryPath, moleculeInfo: dict, qmMethod: str, solvationMethod: str) -> FilePath:
    ## unpack moleculeInfo
    charge = moleculeInfo["charge"]
    multiplicity = moleculeInfo["multiplicity"]
    ## create orca input file
    orcaInputFile = p.join(outDir, "orca_sp.inp")
    with open(orcaInputFile, "w") as f:
        f.write(" # --------------------------------- #\n")
        f.write(" #  Single Point Calculation         #\n")
        f.write(" # --------------------------------- #\n")
        ## METHOD
        if solvationMethod is None:
            f.write(f"! {qmMethod} SP\n")
        else:
            f.write(f"! {qmMethod} {solvationMethod} SP\n") 
        ## GEOMETRY
        f.write(f"*xyzfile {charge} {multiplicity} {inputXyz}\n\n")

    return orcaInputFile
###################################################################################



def make_orca_input_for_opt_freq(inputXyz: FilePath, outDir: DirectoryPath, moleculeInfo: dict, qmMethod: str, solvationMethod: str, nCpus: int | None = None) -> FilePath:
    ## unpack moleculeInfo
    charge = moleculeInfo["charge"]
    multiplicity = moleculeInfo["multiplicity"]
    ## create orca input file
    orcaInputFile = p.join(outDir, "orca_opt_freq.inp")
    with open(orcaInputFile, "w") as f:
        f.write(" # --------------------------------- #\n")
        f.write(" #  Frequency Calculation            #\n")
        f.write(" # --------------------------------- #\n")
        ## METHOD
        if solvationMethod is None:
            f.write(f"! {qmMethod} Opt Freq\n")
        else:
            f.write(f"! {qmMethod} {solvationMethod} Opt Freq\n") 
        if nCpus is not None:
            f.write(f"%pal nprocs {nCpus}\nend\n")
        ## GEOMETRY
        f.write(f"*xyzfile {charge} {multiplicity} {inputXyz}\n\n")

    return orcaInputFile


def make_orca_input_for_singlepoint(inputXyz: FilePath, outDir: DirectoryPath, moleculeInfo: dict, qmMethod: str, solvationMethod: str, nCpus: int | None = None) -> FilePath:
    ## unpack moleculeInfo
    charge = moleculeInfo["charge"]
    multiplicity = moleculeInfo["multiplicity"]
    ## create orca input file
    orcaInputFile = p.join(outDir, "orca_sp.inp")
    with open(orcaInputFile, "w") as f:
        f.write(" # --------------------------------- #\n")
        f.write(" #  Single Point Calculation         #\n")
        f.write(" # --------------------------------- #\n")
        ## METHOD
        if solvationMethod is None:
            f.write(f"! {qmMethod} SP\n")
        else:
            f.write(f"! {qmMethod} {solvationMethod} SP\n") 
        if nCpus is not None:
            f.write(f"%pal nprocs {nCpus}\nend\n")
        ## GEOMETRY
        f.write(f"*xyzfile {charge} {multiplicity} {inputXyz}\n\n")

    return orcaInputFile



def make_orca_input_for_resp(inputXyz: FilePath, outDir: DirectoryPath, moleculeInfo: dict, qmMethod: str, solvationMethod: str, nCpus: int | None = None) -> FilePath:
    ## unpack moleculeInfo
    charge = moleculeInfo["charge"]
    multiplicity = moleculeInfo["multiplicity"]
    ## create orca input file
    orcaInputFile = p.join(outDir, "orca_sp.inp")
    with open(orcaInputFile, "w") as f:
        f.write(" # ------------------------------------ #\n")
        f.write(" # Restrained EletcroStatic Potential   #\n")
        f.write(" # ------------------------------------ #\n")
        ## METHOD
        if solvationMethod is None:
            f.write(f"! {qmMethod} RESP\n")
        else:
            f.write(f"! {qmMethod} {solvationMethod} SP\n") 
        if nCpus is not None:
            f.write(f"%pal nprocs {nCpus}\nend\n")
        ## GEOMETRY
        f.write(f"*xyzfile {charge} {multiplicity} {inputXyz}\n\n")

    return orcaInputFile

###################################################################################
def did_orca_finish_normallly(orcaOut: FilePath) -> bool:
    with open(orcaOut, "r") as f:
        lines = f.readlines()
        for line in reversed(lines):
            if "****ORCA TERMINATED NORMALLY****" in line:
                return True
            elif "ORCA finished with an error in the energy calculation" in line:
                return False
    return False
###################################################################################
