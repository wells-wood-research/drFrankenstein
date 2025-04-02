import os
from os import path as p
from subprocess import call


########################################
def run_orca(orcaInput, orcaOutput, config):
    orcaExe = config["pathInfo"]["orcaExe"]
    orcaCommand = [orcaExe, orcaInput]
    with open(orcaOutput, 'w') as output_file:
        try:
            call(orcaCommand, stdout=output_file, stderr=output_file)
        except Exception as e:
            raise(e)

########################################
def write_goat_input(conformerDir, cappedXyz, config):
    ## write GOAT orca input file
    charge = config["moleculeInfo"]["charge"]
    multiplicity = config["moleculeInfo"]["multiplicity"]
    goatOrcaInput = p.join(conformerDir, f"GOAT_orca.inp")
    with open(goatOrcaInput, "w") as f:
        f.write("!XTB2 GOAT\n")
        f.write("%PAL NPROCS 16 END\n")  ##TODO: work out how to use more cores or have a stable way of inputting this
        f.write("%GOAT\n")
        f.write("\tMAXITERMULT 1\n")            ## only do one round (save some time)
        f.write("\tFREEZEAMIDES TRUE\n")           ## freeze amides
        f.write("\tFREEZECISTRANS TRUE\n ")         ## freeze cis-trans bonds
        f.write("END\n")
        f.write(f"*xyzfile {charge} {multiplicity} {cappedXyz}\n")
    return goatOrcaInput

########################################
def make_orca_input_for_opt(inputXyz, outDir, moleculeInfo, qmMethod, solvationMethod):
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
                                         scanAngles):
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
def  make_orca_input_for_singlepoint(inputXyz, outDir, moleculeInfo, qmMethod, solvationMethod):
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
def did_orca_finish_normallly(orcaOut):
    with open(orcaOut, "r") as f:
        lines = f.readlines()
        for line in reversed(lines):
            if "****ORCA TERMINATED NORMALLY****" in line or "* finished run" in line:
                return True
    return False
###################################################################################