from rdkit import Chem
from rdkit.Chem import AllChem
import subprocess
from subprocess import call, DEVNULL 
from os import path as p
import os
import glob
from shutil import rmtree, copy
from pdbUtils import pdbUtils
import pandas as pd
import numpy as np
import time
import pexpect

## drFrankenstein LIBRARIES ##
import drTwist
import drInputs

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

    chargeGroupIndexes, nAtoms = get_charge_group_indexes(moleculeInfo["cappedPdb"])
    chargeConstraintsTxt = generate_charge_constraints_file(chargeGroupIndexes, moleculeInfo["charge"], pathInfo["chargeFittingDir"])

    pathInfo.update({"chargeConstraints": chargeConstraintsTxt})

    run_qm_calculations_for_charge_fitting(pathInfo,  moleculeInfo, chargeFittingInfo)
    conformerListTxt = generate_conformer_list_file(pathInfo["orcaCalculationsDir"], pathInfo["chargeFittingDir"])
    pathInfo.update({"conformerListForChargeFitting": conformerListTxt})

    rawMultiWfnOutputs = run_charge_fitting(pathInfo)
    chargesDf = parse_multiwfn_output(rawMultiWfnOutputs)
    chargesCsv = p.join(pathInfo["chargeFittingDir"], "charges.csv")
    chargesDf.to_csv(chargesCsv)

###########################################################################
###########################################################################
def set_up_directories(config) -> dict:
    ## get pathInfo dict
    pathInfo: dict = config["pathInfo"]

    outputDir: DirectoryPath = pathInfo["outputDir"]
    ## charge dir - acts as a topDir for all charge-related processes
    chargeDir: DirectoryPath = p.join(outputDir, "03_charge_calculations")
    os.makedirs(chargeDir, exist_ok=True)

    ## directory to store conformer PDB files
    ## NOTE these do not interact with torsion scans at all!
    conformerDir: DirectoryPath = p.join(chargeDir, "conformers")
    os.makedirs(conformerDir, exist_ok=True)
    ## directory to run ORCA QM calculations
    orcaCalculationsDir: DirectoryPath = p.join(chargeDir, "orca_calculations")
    os.makedirs(orcaCalculationsDir, exist_ok=True)
    ## directory to run MultiWFN for charge fitting
    chargeFittingDir: DirectoryPath = p.join(chargeDir, "charge_fitting")
    os.makedirs(chargeFittingDir, exist_ok=True)

    ## update config
    config["pathInfo"].update({
        "chargeDir": chargeDir,
        "conformerDir": conformerDir,
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
def run_charge_fitting(pathInfo):
    ## unpack pathInfo
    multiWfnDir = pathInfo["multiWfnDir"]
    chargeFittingDir = pathInfo["chargeFittingDir"]
    chargeConstraintsTxt = pathInfo["chargeConstraints"]
    conformerListTxt = pathInfo["conformerListForChargeFitting"]

    ## change dir to MultiWFN build dir so it can find settings.ini
    os.chdir(multiWfnDir)

    ## get MultoiWFN binary
    multiWfnExe = p.join(multiWfnDir, "Multiwfn_noGUI")
    if not p.isfile(multiWfnExe):
        raise(FileNotFoundError, "Multiwfn_noGUI not found")

    ## get a molden file for the input command
    moldenFile = glob.glob(p.join(chargeFittingDir,"*.molden.input"))[0]

    ## create a child with pexpect using a command
    chargeFittingCommand = f"{multiWfnExe} {moldenFile}"
    child = pexpect.spawn(chargeFittingCommand, encoding='utf-8')

    # Capture all output until the expected prompt
    output = ""

    ####### GET TO RESP MAIN MENU #######

    ## look for main menu, input "7" for charge calculations
    child.expect(r'.*300.*\n?\r?')
    child.sendline("7")

    ## look for charge calculation menu, input "18" for RESP calculations
    child.expect(r'.*20.*\n?\r?')
    child.sendline("18")

    ###### LOAD CONFORMERS #######

    ## look for RESP main menu, input "-1" to load conformer file
    child.expect(r'.*11.*\n?\r?')
    child.sendline("-1")

    ## look for directory input prompt, input directory
    child.expect(r'.*Input.*\n?\r?')
    child.sendline(conformerListTxt)

    ####### LOAD CHARGE CONSTRAINTS #######
    ## look for RESP main menu, input "6" to load charge constraint file
    child.expect(r'.*11.*\n?\r?')
    child.sendline("6")

    ## look for conformation, input "1" to proceed
    child.expect(r'.*1.*\n?\r?')
    child.sendline("1")

    ## look for directory input prompt, input directory
    child.expect(r'.*Input.*\n?\r?')
    child.sendline(chargeConstraintsTxt)

    ####### RUN CALCULATIONS #######

    ## look for RESP main menu, input "2" to run charge calculations
    child.expect(r'.*11.*\n?\r?')
    child.sendline("2")
    output += child.after
    ## look for Note at the end of calculation, input "q" to quit
    child.expect(r'.*Note: Because.*\n?\r?', timeout=None)
    child.sendline("q")
    output += child.after
    outputFile = p.join(chargeFittingDir, "MultiWfn_raw_outputs.txt")

    # Write the output to the file
    with open(outputFile, 'w') as file:
        file.write(output)

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
def run_qm_calculations_for_charge_fitting(pathInfo,  moleculeInfo, chargeFittingInfo, debug = False):
    ## unpack molecule info to get input pdb
    inputPdb = moleculeInfo["cappedPdb"]
    ## make a pre-set number of conformers
    conformerPdbs = drTwist.gen_conformers(inputPdb, pathInfo["chargeFittingDir"], 1, chargeFittingInfo["nConformers"])
    argsList = [(conformerPdb, pathInfo, chargeFittingInfo,  moleculeInfo) for conformerPdb in conformerPdbs]
    
    
    tqdmBarOptions = {
        "desc": f"Running Charge Calculations",
        "ascii": "-🗲",  # Use the character for the bar
        "colour": "yellow",
        "unit":  "Calc",
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
    conformerPdb, pathInfo, chargeFittingInfo, moleculeInfo,  = args
    conformerName = p.basename(conformerPdb).split(".")[0]
    
    conformerQmDir = p.join(pathInfo["orcaCalculationsDir"], conformerName)
    os.makedirs(conformerQmDir, exist_ok=True)

    orcaOptInput = generate_orca_input_for_optimisation(conformerPdb,
                                                            conformerQmDir,
                                                            moleculeInfo["charge"], 
                                                            moleculeInfo["multiplicity"],
                                                            chargeFittingInfo["optMethod"],
                                                            chargeFittingInfo["optSolvationMethod"],
                                                            chargeFittingInfo["nCoresPerCalculation"])
    
    run_orca(pathInfo["orcaExe"], orcaOptInput, p.join(conformerQmDir, "orca_geom_opt.out"))

    # orcaSolvateInput = generate_orca_input_for_solvation(conformerPdb,
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
def generate_orca_input_for_solvation(conformerPdb, conformerQmDir, charge, multiplicity, nWaters):
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
        conformerDf = pdbUtils.pdb2df(conformerPdb)
        f.write(f"*xyz {charge} {multiplicity}\n\n")
        for index, row in conformerDf.iterrows():
            geomLine = f"{row['ELEMENT']} {row['X']} {row['Y']} {row['Z']}\n"
            f.write(geomLine)
        f.write("*")


    return orcaInputFile

##########################################################
def generate_charge_constraints_file(chargeGroupIndexes, charge, chargeFittingDir):
    chargeConstraintsTxt = p.join(chargeFittingDir, "charge_group_constraints.txt")
    with open(chargeConstraintsTxt, "w") as f:
        for chargeGroup in ["N-Methyl", "Acyl", "Backbone"]:
            formattedIndexes = ",".join(map(str, chargeGroupIndexes[chargeGroup]))
            f.write(f"{formattedIndexes} 0\n")
          
        formattedIndexes = ",".join(map(str, chargeGroupIndexes["Sidechain"]))
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
def get_charge_group_indexes(pdbFile) -> dict:

    ##TODO: Implement fragmentation in sidechain C-C bonds so we can constrain charges nicely!

    rdkitMol = Chem.MolFromPDBFile(pdbFile, removeHs = False)
    rdkitMol = enforce_carboxyl_bonding(rdkitMol)

    backBoneAndCapsSmarts = "[C]([H])([H])([H])[N]([H])[C](=[O])[C@@]([H])[N]([H])[C](=[O])[C]([H])([H])([H])"
    # Find the substructure match
    backBoneAndCapsPattern = Chem.MolFromSmarts(backBoneAndCapsSmarts)
    backBoneAndCapsMatches = rdkitMol.GetSubstructMatches(backBoneAndCapsPattern)

    ## get caps and backbone using matching 
    backBoneAndCapsIndexes = [list(match) for match in backBoneAndCapsMatches]
    nMethylCapIndexes = backBoneAndCapsIndexes[0][0:6]
    acylCapIndexes = backBoneAndCapsIndexes[0][-6:]
    backBoneIndexes = backBoneAndCapsIndexes[0][6:-6]


    ## get sidechain indexes 
    allAtomIndexes = set(range(rdkitMol.GetNumAtoms()))
    sideChainIndexes = [index for index in allAtomIndexes if not index in backBoneAndCapsIndexes[0]]
    chargeGroupIndexes = {
        "N-Methyl": nMethylCapIndexes,
        "Acyl": acylCapIndexes,
        "Backbone": backBoneIndexes,
        "Sidechain": sideChainIndexes
    }

    nAtoms = rdkitMol.GetNumAtoms()
    return chargeGroupIndexes, nAtoms
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
        ## GEOMETRY
        f.write(f"*xyzfile {charge} {multiplicity} {optXyz}\n\n")

    return orcaInputFile
##########################################################
def generate_orca_input_for_optimisation(conformerPdb, conformerChargeDir, charge, multiplicity, qmMethod, solvationMethod, nCores):
    ## create orca input file
    orcaInputFile = p.join(conformerChargeDir, "orca_geom_opt.inp")
    with open(orcaInputFile, "w") as f:

        f.write(" # --------------------------------- #\n")
        f.write(" #  Geometry Optimisation            #\n")
        f.write(" # --------------------------------- #\n")
        ## METHOD
        f.write(f"! {qmMethod} {solvationMethod} OPT\n")
        f.write(f"! pal{str(nCores)}\n")
        ## GEOMETRY
        conformerDf = pdbUtils.pdb2df(conformerPdb)
        f.write(f"*xyz {charge} {multiplicity}\n\n")
        for index, row in conformerDf.iterrows():
            geomLine = f"{row['ELEMENT']} {row['X']} {row['Y']} {row['Z']}\n"
            f.write(geomLine)
        f.write("*")

    return orcaInputFile

##########################################################
##########################################################



if __name__ == "__main__":
    main()





