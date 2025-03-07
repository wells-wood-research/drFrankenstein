## BASIC IMPORTS ##
import os
from os import path as p
from subprocess import call, PIPE
import yaml

def  conformer_generation_protocol(config):
    print("--> GENERATING CONFORMERS WITH GOAT")

    config = sort_out_directories(config)
    cappedPdb = config["moleculeInfo"]["cappedPdb"]
    moleculeName = config["moleculeInfo"]["moleculeName"]
    conformerDir = config["pathInfo"]["conformerDir"]
    cappedXyz = p.join(conformerDir, f"{moleculeName}_capped.xyz")
    pdb_to_xyz(cappedPdb, cappedXyz)
    goatOrcaInput = write_goat_input(conformerDir, cappedXyz, config)

    run_goat_conformer_generation(goatOrcaInput, conformerDir, config)

    conformerXyzs = split_conformers(conformerDir, moleculeName)
    config["pathInfo"]["conformerXyzs"] = conformerXyzs

    clean_up(conformerDir)

    config["checkpointInfo"]["conformersComplete"] = True
    return config

def  sort_out_directories(config):
    ## sort out directories
    outputDir = config["pathInfo"]["outputDir"]
    conformerDir = p.join(outputDir, "02_GOAT_conformers")
    os.makedirs(conformerDir, exist_ok=True) 
    config["pathInfo"]["conformerDir"] = conformerDir

    return config   

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

def run_goat_conformer_generation(goatOrcaInput, conformerDir, config):
    ## run GOAT
    os.chdir(conformerDir)
    orcaExe = config["pathInfo"]["orcaExe"]
    goatOrcaOutput = p.join(conformerDir, f"GOAT_orca.out")
    with open(goatOrcaOutput, 'w') as goatOrcaOutputFile:
            goatOrcaCommand = [orcaExe, goatOrcaInput]
            call(goatOrcaCommand, stdout=goatOrcaOutputFile, 
                stderr=goatOrcaOutputFile)
            
def split_conformers(conformerDir, moleculeName):
    ## split GOAT output into individual conformers
    goatFinalXyz = p.join(conformerDir, "GOAT_orca.finalensemble.xyz")
    os.chmod(goatFinalXyz, 0o755)
    conformerXyzFilePattern = p.join(conformerDir, f"{moleculeName}_conformer_.xyz")
    obabelCommand = ["obabel", goatFinalXyz, "-O", conformerXyzFilePattern, "-m"]
    call(obabelCommand, stdout=PIPE, stderr=PIPE)

    conformerXyzs = [p.join(conformerDir, f) for f in os.listdir(conformerDir) if f.startswith(f"{moleculeName}_conformer_")]
    return conformerXyzs


def pdb_to_xyz(inPdb, outXyz):
    obabelCommand = ["obabel", "-i", "pdb", inPdb, "-o", "xyz", "-O", outXyz]
    call(obabelCommand, stdout=PIPE, stderr=PIPE)

def clean_up(conformerDir):
    filesToRemove = [p.join(conformerDir, f) for f in os.listdir(conformerDir) if f.startswith("GOAT")]
    for f in filesToRemove:
        os.remove(f)

if __name__ == "__main__":
    configYaml = "/home/esp/scriptDevelopment/drFrankenstein/NMH_outputs/drFrankenstein.yaml"
    with open(configYaml, "r") as yamlFile:
        config = yaml.safe_load(yamlFile)
    conformer_generation_protocol(config)