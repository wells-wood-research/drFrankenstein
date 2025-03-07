from rdkit import Chem
from rdkit.Chem import AllChem
from subprocess import call, DEVNULL 
from os import path as p
import os
import glob
from shutil import rmtree, copy
from pdbUtils import pdbUtils
import pandas as pd
import numpy as np
import pexpect
import yaml
import re

from tqdm import tqdm
## drFrankenstein LIBRARIES ##
import drTwist

## MULTIPROCESSING AND LOADING BAR LIBRARIES ##c
from mpire import WorkerPool
from mpire.utils import make_single_arguments


class FilePath:
    pass
class DirectoryPath:
    pass

"""
1. Generate conformers
2. load each conformer molecule from pdb
3. identify capping groups
4. set up a ORCA job that does:
    - Geom optimisation with XTB
    - Single point with HF [ must output molden file]
5. use MultiWFN  to calculate charges
    7 (charge fitting)
    18 (RESP)
    5 (equivalence)
    6 (charge restraint)
        1 (yes)
        PATH/TO/charge_group_constraints.txt
    -1 (load conformer list and weights)
        PATH/TO/conformer_list.txt
    2 (Start one-stage ESP fitting calculation with constraints)

6. set up while loop to check for convergence
7. save converged charges  to file
"""

def dummy_inputs():

    config = {

        "moleculeInfo" : {
        "charge": 0,
        "multiplicity": 1,
        "moleculePdb": "/home/esp/scriptDevelopment/drFrankenstein/ALA_outputs/capped_amino_acids/ALA_capped.pdb"
        },

        "pathInfo" : {
        "outputDir" : "/home/esp/scriptDevelopment/drFrankenstein/ALA_outputs",
        "orcaExe": "/home/esp/bin/orca_6_0_1_linux_x86-64_shared_openmpi416/orca",
        "multiWfnDir":  "/home/esp/bin/Multiwfn_3.8_dev_bin_Linux_noGUI/"
    },

        "chargeFittingInfo" : {
        "optMethod": "! HF 6-31G(d)",
        # "singlePointMethod": "! revPBE def2-TZVP D3BJ SP CPCM(water)"
        "singlePointMethod": "! HF 6-31G(d)",
        "nCoresPerCalculation": 8,
        "nConformers": 4
    }

    }

    return config

###########################################################################
###########################################################################
def charge_protocol(config, debug=False):
    ## unpack config 
    moleculeInfo = config["moleculeInfo"]
    chargeFittingInfo = config["chargeFittingInfo"]

    config = set_up_directories(config)

    pathInfo = config["pathInfo"]

    config = get_charge_group_indexes(moleculeInfo["cappedPdb"], config)
    chargeConstraintsTxt = generate_charge_constraints_file(config)

    pathInfo.update({"chargeConstraints": chargeConstraintsTxt})

    run_qm_calculations_for_charge_fitting(pathInfo,  moleculeInfo, chargeFittingInfo, config, debug)
    conformerListTxt = generate_conformer_list_file(pathInfo["orcaCalculationsDir"], pathInfo["chargeFittingDir"])
    pathInfo.update({"conformerListForChargeFitting": conformerListTxt})

    rawMultiWfnOutputs = run_charge_fitting(config)
    chargesDf = parse_multiwfn_output(rawMultiWfnOutputs)
    chargesCsv = p.join(pathInfo["chargeFittingDir"], "charges.csv")
    chargesDf.to_csv(chargesCsv)


    config["checkpointInfo"]["chargesComplete"] = True
    return config

###########################################################################
###########################################################################
def set_up_directories(config) -> dict:
    ## get pathInfo dict
    pathInfo: dict = config["pathInfo"]

    outputDir: DirectoryPath = pathInfo["outputDir"]
    ## charge dir - acts as a topDir for all charge-related processes
    chargeDir: DirectoryPath = p.join(outputDir, "04_charge_calculations")
    os.makedirs(chargeDir, exist_ok=True)


    ## directory to run ORCA QM calculations
    orcaCalculationsDir: DirectoryPath = p.join(chargeDir, "orca_calculations")
    os.makedirs(orcaCalculationsDir, exist_ok=True)
    ## directory to run MultiWFN for charge fitting
    chargeFittingDir: DirectoryPath = p.join(chargeDir, "charge_fitting")
    os.makedirs(chargeFittingDir, exist_ok=True)

    ## update config
    config["pathInfo"].update({
        "chargeDir": chargeDir,
        "chargeFittingDir": chargeFittingDir,
        "orcaCalculationsDir": orcaCalculationsDir
    })

    return config
###########################################################################
def parse_multiwfn_output(rawMultiWfnOutputs):
    """
    Parses the output from Multiwfn to extract atomic charges.

    Args:
        rawMultiwfnOutputs (str): Path to the file containing the raw Multiwfn output.

    Returns:
        pd.DataFrame: A DataFrame containing columns 'atomIndex', 'atomElement', and 'Charge', 
                      representing the index, element, and charge of each atom, respectively.

    This function reads the Multiwfn output file, identifies the section after the 
    "Successfully converged!" line, and extracts atomic charges until the line containing 
    "Sum of charges:". The extracted data is structured into a DataFrame for further analysis.
    """

    with open(rawMultiWfnOutputs, 'r') as file:
        lines = file.readlines()

    # Find the index of the line after "Successfully converged!"
    start_index = next(i for i, line in enumerate(lines) 
                       if "Successfully converged!" in line) + 3

    # Extract the charge data lines
    chargeDataLines = []
    for line in lines[start_index:]:
        if "Sum of charges:" in line.strip():
            break
        chargeDataLines.append(line.strip())

    # Parse the charge data into a DataFrame
    data = []
    for line in chargeDataLines:
        parts = line.split(")")
        charge = float(parts[1].strip())
        parts = parts[0].split("(")
        atomIndex = int(parts[0])
        atomElement = parts[1].strip()

        data.append({'atomIndex': atomIndex, 'atomElement': atomElement, 'Charge': charge})

    chargesDf = pd.DataFrame(data)
    return chargesDf

###########################################################################
def run_charge_fitting(config):
    ## Unpack pathInfo
    pathInfo = config["pathInfo"]
    multiWfnDir = pathInfo["multiWfnDir"]
    chargeFittingDir = pathInfo["chargeFittingDir"]
    chargeConstraintsTxt = pathInfo["chargeConstraints"]
    conformerListTxt = pathInfo["conformerListForChargeFitting"]

    nConformers = config["chargeFittingInfo"]["nConformers"]

    ## Change dir to MultiWFN build dir so it can find settings.ini
    os.chdir(multiWfnDir)

    ## Get MultiWFN binary
    multiWfnExe = p.join(multiWfnDir, "Multiwfn_noGUI")
    if not p.isfile(multiWfnExe):
        raise FileNotFoundError("Multiwfn_noGUI not found")

    ## Get a molden file for the input command
    moldenFile = glob.glob(p.join(chargeFittingDir, "*.molden.input"))[0]

    ## Define output file
    outputFile = p.join(chargeFittingDir, "MultiWfn_raw_outputs.txt")

    ## Open a file object for continuous logging
    with open(outputFile, 'w') as logFile:
        # Clear file and write a starting message
        logFile.write("Starting charge fitting process...\n")
        logFile.flush()  # Ensure it’s written immediately

        ## Spawn the child process with logfile enabled
        chargeFittingCommand = f"{multiWfnExe} {moldenFile}"
        child = pexpect.spawn(chargeFittingCommand, encoding='utf-8', logfile=logFile)


        ####### GET TO RESP MAIN MENU #######
        child.expect(r'.*300.*\n?\r?')
        child.sendline("7")

        child.expect(r'.*20.*\n?\r?')
        child.sendline("18")

        ###### LOAD CONFORMERS #######
        child.expect(r'.*11.*\n?\r?')
        child.sendline("-1")

        child.expect(r'.*Input.*\n?\r?')
        child.sendline(conformerListTxt)

        ####### LOAD CHARGE CONSTRAINTS #######
        child.expect(r'.*11.*\n?\r?')
        child.sendline("6")

        child.expect(r'.*1.*\n?\r?')
        child.sendline("1")

        child.expect(r'.*Input.*\n?\r?')
        child.sendline(chargeConstraintsTxt)

        ####### RUN CALCULATIONS #######
        child.expect(r'.*11.*\n?\r?')
        child.sendline("2")



        # Monitor the output for progress updates

        # Initialize tqdm progress bar
        tqdmBarOptions = {
        "desc": f"Performing Charge Fitting",
        "ascii": "-ϟ",  
        "colour": "yellow",
        "unit":  "calculations",
        "dynamic_ncols": True
            }
        totalTickProgress = nConformers * 100.00
        progress_bar = tqdm(total=totalTickProgress, **tqdmBarOptions)
        try:
            while True:
                # Expect either a progress line or the end of the process (EOF)
                index = child.expect([r"Progress: \[.*?\]\s+(\d+\.\d+) %", 
                                      r" 11 Choose ESP type, current: Nuclear \+ Electronic[\n\r]?", pexpect.TIMEOUT], timeout=5)
                if index == 0:  # Progress line matched
                    progress_bar.update(1)  # Update the tqdm bar by 1 tick
                
                elif index == 1:  # (process finished)
                    break
                
                elif index == 2:  # Timeout (no output for 5 seconds)
                    continue

        except pexpect.ExceptionPexpect as e:
            raise e

        # Close the child process
        child.close()

    return outputFile

###########################################################################
def generate_conformer_list_file(orcaCalculationsDir, chargeFittingDir):
    singlePointData = {}
    for conformerName in os.listdir(orcaCalculationsDir):
        singlePointData[conformerName] = {}
        calculationDir = p.join(orcaCalculationsDir, conformerName)
        if not p.isdir:
            continue
        singlePointOutFile = p.join(calculationDir, "orca_single_point.out")
        singlePointEnergy = find_final_single_point_energy(singlePointOutFile)
        singlePointData[conformerName]["Energy"] = singlePointEnergy
        singlePointData[conformerName]["Path"] = p.join(chargeFittingDir, f"{conformerName}.molden.input")


    singlePointDf =pd.DataFrame.from_dict(singlePointData, orient='index')

    singlePointDf = calculate_boltzmann_probabilities(singlePointDf)

    conformerListTxt = p.join(chargeFittingDir, "conformer_list.txt")
    with open(conformerListTxt, "w") as f:
        for _, row in singlePointDf.iterrows():
            f.write(f"{row['Path']} {str(row['Probability'])}\n")

    return conformerListTxt

###########################################################################
def calculate_boltzmann_probabilities(df, temperature=298.15):
    boltzmannConstant = 0.0083144621  # kJ/(mol*K)
    energies = df['Energy'].values
    minEnergy = np.min(energies)
    deltaEnergies = energies - minEnergy
    boltzmannFactors = np.exp(-deltaEnergies / 
                              (boltzmannConstant * temperature))
    probabilities = boltzmannFactors / np.sum(boltzmannFactors)
    df['Probability'] = probabilities
    return df    

###########################################################################
def find_final_single_point_energy(outFilePath):
    with open(outFilePath, 'r') as file:
        for line in file:
            if "FINAL SINGLE POINT ENERGY" in line:
                return float(line.split()[-1])



###########################################################################
def run_qm_calculations_for_charge_fitting(pathInfo,  moleculeInfo, chargeFittingInfo, config,  debug = False):
    ## unpack molecule info to get input pdb
    
    conformerXyzs = config["pathInfo"]["conformerXyzs"]

    sampledConformerXyzs = conformerXyzs[:config["chargeFittingInfo"]["nConformers"]]


    argsList = [(conformerXyz, pathInfo, chargeFittingInfo,  moleculeInfo) for conformerXyz in sampledConformerXyzs]
    
    
    tqdmBarOptions = {
        "desc": f"Running Charge Calculations",
        "ascii": "-ϟ",  
        "colour": "yellow",
        "unit":  "scan",
        "dynamic_ncols": True
    }
    
    if debug:
        for arg in argsList:
            process_conformer(arg)

    else:
        with WorkerPool(n_jobs = chargeFittingInfo["nConformers"]) as pool:
            pool.map(process_conformer,
                    make_single_arguments(argsList),
                    progress_bar=True,
                    iterable_len = len(argsList),
                    progress_bar_options=tqdmBarOptions)


###########################################################################
def process_conformer(args):
    conformerXyz, pathInfo, chargeFittingInfo, moleculeInfo,  = args
    conformerName = p.basename(conformerXyz).split(".")[0]
    
    conformerQmDir = p.join(pathInfo["orcaCalculationsDir"], conformerName)
    os.makedirs(conformerQmDir, exist_ok=True)

    orcaOptInput = generate_orca_input_for_optimisation(conformerXyz,
                                                            conformerQmDir,
                                                            moleculeInfo["charge"], 
                                                            moleculeInfo["multiplicity"],
                                                            chargeFittingInfo["optMethod"],
                                                            chargeFittingInfo["optSolvationMethod"],
                                                            chargeFittingInfo["nCoresPerCalculation"])
    
    run_orca(pathInfo["orcaExe"], orcaOptInput, p.join(conformerQmDir, "orca_geom_opt.out"))

    # orcaSolvateInput = generate_orca_input_for_solvation(conformerXyz,
    #                                                         conformerQmDir,
    #                                                         moleculeInfo["charge"], 
    #                                                         moleculeInfo["multiplicity"],
    #                                                         nWaters=10)
    # run_orca(pathInfo["orcaExe"], orcaSolvateInput, p.join(conformerQmDir, "orca_solvation.out"))


    optXyz = p.join(conformerQmDir, "orca_geom_opt.xyz")
    orcaSinglePointInput = generate_orca_input_for_single_point(optXyz,
                                                                conformerQmDir,
                                                                moleculeInfo["charge"],
                                                                moleculeInfo["multiplicity"],
                                                            chargeFittingInfo["singlePointMethod"],
                                                            chargeFittingInfo["singlePointSolvationMethod"],
                                                            chargeFittingInfo["nCoresPerCalculation"])
    
    run_orca(pathInfo["orcaExe"], orcaSinglePointInput, p.join(conformerQmDir, "orca_single_point.out"))


    singlePointName = p.splitext(p.basename(orcaSinglePointInput))[0]
    call(["orca_2mkl", p.join(conformerQmDir,singlePointName), "-molden"], stdout=DEVNULL)
    singlePointMolden = p.join(conformerQmDir, "orca_single_point.molden.input")

    destMolden = p.join(pathInfo["chargeFittingDir"], f"{conformerName}.molden.input")

    copy(singlePointMolden, destMolden)
###########################################################################
def generate_orca_input_for_solvation(conformerXyz, conformerQmDir, charge, multiplicity, nWaters):
    ## create orca input file
    orcaInputFile = p.join(conformerQmDir, "orca_solvation.inp")
    with open(orcaInputFile, "w") as f:

        f.write("# --------------------------------- #\n")
        f.write("#  SOLVATION WITH XTB2              #\n")
        f.write("# --------------------------------- #\n")
        ## METHOD
        f.write("! XTB2 ALPB(WATER)\n")
        ## SOLVATOR KEYWORDS
        f.write("%SOLVATOR\n")
        f.write(f"\tNSOLV\t{nWaters}\n")
        f.write("\tCLUSTERMODE STOCHASTIC\n")
        f.write("END\n")
        ## GEOMETRY
        f.write(f"*xyz {charge} {multiplicity} {conformerXyz}\n")
        f.write("END\n")


    return orcaInputFile

##########################################################
def generate_charge_constraints_file(config):
    chargeFittingDir = config["pathInfo"]["chargeFittingDir"]

    chargeGroups: dict = config["moleculeInfo"]["chargeGroups"]

    chargeConstraintsTxt = p.join(chargeFittingDir, "charge_group_constraints.txt")
    with open(chargeConstraintsTxt, "w") as f:
        for _, chargeGroupData in chargeGroups.items():
            indexes = chargeGroupData["indexes"]
            if len(indexes) == 0:
                continue
            charge = chargeGroupData["charge"]
            formattedIndexes = ",".join(map(str, indexes))
            f.write(f"{formattedIndexes} {str(charge)}\n")
          

    return chargeConstraintsTxt
##########################################################
def run_orca(orcaExe, orcaInput, orcaOutput):
    orcaCommand = [orcaExe, orcaInput]
    with open(orcaOutput, 'w') as output_file:
        try:
            call(orcaCommand, stdout=output_file, stderr=output_file)
        except Exception as e:
            pass
##########################################################
def enforce_carboxyl_bonding(mol):
    # Define the initial and target SMARTS patterns
    amideSmarts = '[C]([O])([N])([C])'

    # Find the substructure match
    amidePattern = Chem.MolFromSmarts(amideSmarts)
    matches = mol.GetSubstructMatches(amidePattern)
    
    if not matches:
        raise ValueError("No matching substructure found.")
    
    # Modify the first match found
    emol = Chem.EditableMol(mol)
    for match in matches:
        # Get the indices of the atoms in the match
        carbon_idx = match[0]
        oxygen_idx = match[1]        
        # Remove the single bond and add a double bond
        emol.RemoveBond(carbon_idx, oxygen_idx)
        emol.AddBond(carbon_idx, oxygen_idx, Chem.BondType.DOUBLE)
    mol = emol.GetMol()

    # Return the modified molecule
    return emol.GetMol()   
##########################################################
def get_charge_group_indexes(pdbFile, config) -> dict:

    pdbDf = pdbUtils.pdb2df(pdbFile)

    chargeGroups = config["moleculeInfo"]["chargeGroups"]
    overallCharge = config["moleculeInfo"]["charge"]

    ## remove capping groups
    uncappedDf = pdbDf[~pdbDf["RES_NAME"].isin(["NME", "ACE"])]

    ## deal with user-defined charge groups
    userDefinedAtoms = []
    userDefinedCharge = 0
    for chargeGroupName, chargeGroupData in chargeGroups.items():
        chargeGroupAtoms = chargeGroupData["atoms"]
        chargeGroupIndexes = uncappedDf[uncappedDf["ATOM_NAME"].isin(chargeGroupAtoms)]["ATOM_ID"].to_list()
        chargeGroups[chargeGroupName]["indexes"] = chargeGroupIndexes
        userDefinedCharge += chargeGroupData["charge"]
        userDefinedAtoms.extend(chargeGroupAtoms)

    ## deal with left-over atoms, these will all go in one group
    leftOverDf = uncappedDf[~uncappedDf["ATOM_NAME"].isin(userDefinedAtoms)]
    leftOverAtoms = leftOverDf["ATOM_NAME"].to_list()
    leftOverIndexes = leftOverDf["ATOM_ID"].to_list()
    leftOverCharge = overallCharge - userDefinedCharge 
    ## deal with left-over atoms, these will all go in one group
    leftOverDf = uncappedDf[~uncappedDf["ATOM_NAME"].isin(userDefinedAtoms)]
    leftOverAtoms = leftOverDf["ATOM_NAME"].to_list()
    leftOverIndexes = leftOverDf["ATOM_ID"].to_list()
    leftOverCharge = overallCharge - userDefinedCharge 
    
    chargeGroups["left-over-atoms"] = {
        "atoms" : leftOverAtoms,
        "charge" : leftOverCharge,
        "indexes" : leftOverIndexes
    }
    ## deal with Terminal Caps 
    for capResName in ["NME", "ACE"]:
        capDf = pdbDf[pdbDf["RES_NAME"]==capResName]
        for capIndex, capResDf in capDf.groupby("RES_ID"):
            chargeGroupAtoms = capResDf["ATOM_NAME"].to_list()
            chargeGroupIndexes = capResDf["ATOM_ID"].to_list()
            chargeGroupName = f"{capResName}_cap_{capIndex}"
            chargeGroups[chargeGroupName] = {
                "atoms": chargeGroupAtoms,
                "charge": 0,
                "indexes": chargeGroupIndexes
            }

    config["moleculeInfo"]["chargeGroups"] = chargeGroups

    return config
##########################################################
def generate_orca_input_for_single_point(optXyz, conformerChargeDir, charge, multiplicity, qmMethod, solvationMethod, nCores):
    ## create orca input file
    orcaInputFile = p.join(conformerChargeDir, "orca_single_point.inp")
    with open(orcaInputFile, "w") as f:

        f.write("# --------------------------------- #\n")
        f.write("#  Single Point Calculation         #\n")
        f.write("# --------------------------------- #\n")
        ## METHOD
        f.write(f"! {qmMethod} {solvationMethod}\n")
        f.write(f"! pal{str(nCores)}\n")
        f.write("\n")
        ## GEOMETRY
        f.write(f"*xyzfile {charge} {multiplicity} {optXyz}\n")

    return orcaInputFile
##########################################################
def generate_orca_input_for_optimisation(conformerXyz, conformerChargeDir, charge, multiplicity, qmMethod, solvationMethod, nCores):
    ## create orca input file
    orcaInputFile = p.join(conformerChargeDir, "orca_geom_opt.inp")
    with open(orcaInputFile, "w") as f:

        f.write(" # --------------------------------- #\n")
        f.write(" #  Geometry Optimisation            #\n")
        f.write(" # --------------------------------- #\n")
        ## METHOD
        f.write(f"! {qmMethod} {solvationMethod} OPT\n")
        # f.write(f"! pal{str(nCores)}\n")
        f.write("\n")

        ## GEOMETRY
        f.write(f"*xyzfile {charge} {multiplicity} {conformerXyz}\n\n")

    return orcaInputFile

##########################################################
##########################################################



if __name__ == "__main__":
    configYaml = "/home/esp/scriptDevelopment/drFrankenstein/02_NMH_outputs/drFrankenstein.yaml"
    with open(configYaml, "r") as yamlFile:
        config = yaml.safe_load(yamlFile)
    charge_protocol(config, debug=False)





