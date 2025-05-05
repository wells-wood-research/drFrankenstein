
import os
from os import path as p
from subprocess import run, STDOUT, call, PIPE
from pdbUtils import pdbUtils
import parmed
import openmm 
from openmm import app, unit
class FilePath:
    pass

class DirectoryPath:
    pass


"""
1. unpack config
2. make a dir for parameters
3. make a mol2 with antechamber
4. make a frcmod file with tleap
5. map backbone atom types from gaff2 to amber
6. modify frcmod file to have all parameters from gaff2 relative to molecule
7. change types in frcmod file to match amber backbone types

"""


def dummy_inputs():
    config = {}
    config["moleculeInfo"] = {}
    config["moleculeInfo"]["moleculeName"] = "PHG"
    config["pathInfo"] = {}
    config["pathInfo"]["outputDir"] = "/home/esp/scriptDevelopment/drFrankenstein/_hidden_scripts/BB_types/AMBER"
    config["pathInfo"]["gaff2Dat"] = "/home/esp/anaconda3/envs/Igor/dat/leap/parm/gaff2.dat"
    config["pathInfo"]["parm19Dat"] = "/home/esp/anaconda3/envs/Igor/dat/leap/parm/parm19.dat"

    config["runtimeInfo"] = {}
    config["runtimeInfo"]["madeByCapping"] = {}
    config["runtimeInfo"]["madeByCapping"]["cappedPdb"] = "/home/esp/scriptDevelopment/drFrankenstein/06_PHG_outputs/01_termini_capping/geometry_optimisation/PHG_capped_opt.pdb"
    return config


def run_tleap_to_make_params(inMol2: FilePath,
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
        f.write(f"saveamberparm mol {prmtop} {inpcrd}  \n")
        f.write("quit")

    tleapOutput: FilePath = p.join(outDir, f"tleap_{index}.out")

    tleapCommand: list = ["tleap", "-f", tleapInput, ">", tleapOutput]

    call(tleapCommand, stdout=PIPE)

    return prmtop, inpcrd

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
    if p.exists(outMol2):
        return None
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
        "-c", "bcc"
    ]
    with open(antechamberOut, 'w') as outfile:
        run(antechamberCommand, stdout=outfile, stderr=STDOUT)

    ## clean up temporary files
    os.remove(tmpPdb)
    filesTtoRemove = [f for f in os.listdir(workingDir) if f.startswith("ANTECHAMBER")]
    for f in filesTtoRemove:
        os.remove(p.join(workingDir, f))
    return None

def change_backbone_types_amber(mol2File):


    # Define CHARMM36m atom type mapping
    amberDefaultMAp = {
                    "N" : "N",
                    "HN": "H",
                    "CA": "CT",
                    "HA": "H1",
                    "C" : "C",
                    "O" : "O",
                    "NN": "N",
                    "HNN1": "H",
                    "CN": "CT", 
                    "HCN1": "H1",
                    "HCN2": "H1",
                    "HCN3": "H1",
                    "CC1": "C",
                    "OC": "O",
                    "CC2": "CT", 
                    "HC1": "H1",
                    "HC2": "H1",
                    "HC3": "H1"
    }

    parmedStructure = parmed.load_file(mol2File)
    for atom in parmedStructure.atoms:
        if atom.name in amberDefaultMAp:
            atom.type = amberDefaultMAp[atom.name]

    parmedStructure.save(mol2File, overwrite=True)


def create_frcmod_file(mol2File: FilePath,
                        frcmodFile: FilePath,
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

    ## run run parmchk2
    parmchk2Command = ["parmchk2",
                        "-i", mol2File,
                          "-f", "mol2",
                            "-o", frcmodFile,
                              "-a", "Y",
                                "-p", gaff2Dat]

    call(parmchk2Command)
    
    return None

def replace_parameters(moleculeFrcmod: FilePath,
                        paraneterDataset: FilePath) -> None:
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
    ffParams = parmed.amber.AmberParameterSet(paraneterDataset)

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

def main():
    
    config = dummy_inputs()
    ## unpack config
    outDir = config["pathInfo"]["outputDir"]
    moleculeName = config["moleculeInfo"]["moleculeName"]
    ## make a dir, add back to config
    assemblyDir = p.join(outDir, "02_assembly")
    os.makedirs(assemblyDir, exist_ok= True)
    config["runtimeInfo"]["assemblyDir"] = assemblyDir
    ## create a mol2 file with antechamber
    moleculeMol2 = p.join(assemblyDir, f"{moleculeName}_capped.mol2")
    pdb2mol2(config["runtimeInfo"]["madeByCapping"]["cappedPdb"], moleculeMol2, assemblyDir)
    ## reset backbone atom types to AMBER defaults
    change_backbone_types_amber(moleculeMol2)
    ## create frcmod file
    moleculeFrcmod = p.join(assemblyDir, f"{moleculeName}_capped.frcmod")
    create_frcmod_file(moleculeMol2, moleculeFrcmod, config)
    ## replace gaff2 params with amber19sb parameters
    replace_parameters(moleculeMol2, moleculeFrcmod, config["pathInfo"]["parm19Dat"])
    ## TESTING - NOT IN REAL SCRIPT
    prmtop, inpcrd = run_tleap_to_make_params(moleculeMol2, moleculeFrcmod, assemblyDir, "molecule")

    energies = run_mm_singlepoints([config["runtimeInfo"]["madeByCapping"]["cappedPdb"]], prmtop, inpcrd)
    print(energies)



if __name__ == "__main__":
    main()


    