## BASIC IMPORTS ##
import os
from os import path as p
from pdbUtils import pdbUtils
from subprocess import call
import pandas as pd
from functools import reduce
import matplotlib.pyplot as plt

## CLEAN CODE ##
class FilePath:
    pass
class DirectoryPath:
    pass

from typing import List, Tuple
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def dummy_inputs() -> dict:
    config: dict = {
        "pathInfo": {
            "inputDir": "/home/esp/scriptDevelopment/drFrankenstein/inputs",
            "outputDir": "/home/esp/scriptDevelopment/drFrankenstein/outputs"
        },
        "moleculeInfo": {
            "charge": 0,
            "multiplicity": 1
        } 
    }
    return config
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def main():
    config = dummy_inputs()
    ## unpack config
    inputDir = config["pathInfo"]["inputDir"]
    outputDir = config["pathInfo"]["outputDir"]

    moleculeInfo = config["moleculeInfo"]
    ## make top level output dir
    os.makedirs(outputDir,exist_ok=True)
    conformerPdbs = [p.join(inputDir,file) for file in os.listdir(inputDir) if file.endswith(".pdb")]
    conformerIds = [f"conformer_{file.split('.')[0].split('_')[-1]}" for file in os.listdir(inputDir) if file.endswith(".pdb")]

    conformers: dict = {conformerId: conformerPdb for conformerId, conformerPdb in zip(conformerIds, conformerPdbs)}

    torsionIndexes = identify_rotatable_bonds(conformerPdbs[0])
    for torsionIndex in torsionIndexes:
        torsionId = "-".join(map(str, torsionIndex)) 
        torsionTopDir = p.join(outputDir, f"torsion_{torsionId}" )
        os.makedirs(torsionTopDir, exist_ok=True)
        scanDfsToMerge = []
        for conformerId, conformerPdb in conformers.items():
            
            ## do a forwards scan
            torsionConformerDir_forwards: DirectoryPath = p.join(torsionTopDir, f"{conformerId}_forwards")
            os.makedirs(torsionConformerDir_forwards, exist_ok=True)
            scanDf_forwards = do_the_twist(conformerPdb, torsionIndex, torsionConformerDir_forwards, moleculeInfo, "0, 360, 10")
            scanDf_forwards.columns = ["Angle", f"{conformerId}_f"]
            scanDfsToMerge.append(scanDf_forwards)

            ## do a backwards scan
            torsionConformerDir_backwards: DirectoryPath = p.join(torsionTopDir, f"{conformerId}_backwards")
            os.makedirs(torsionConformerDir_backwards, exist_ok=True)
            scanDf_backwards = do_the_twist(conformerPdb, torsionIndex, torsionConformerDir_backwards, moleculeInfo, "360, 0, 10")
            scanDf_backwards.columns = ["Angle", f"{conformerId}_b"]
            scanDfsToMerge.append(scanDf_backwards)

        torsionScanDf = reduce(lambda left, right: pd.merge(left, right, on="Angle"), scanDfsToMerge)
        plot_torsion_scan(torsionScanDf, torsionTopDir)
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def do_the_twist(conformerPdb, torsionIndex, conformerTorsionDir, moleculeInfo, scanAngles):
    os.chdir(conformerTorsionDir)
    charge = moleculeInfo["charge"]
    multiplicity = moleculeInfo["multiplicity"]

    # extract xyz and element from input pdb
    orcaInput: FilePath = generate_orca_input(conformerPdb, torsionIndex, conformerTorsionDir, charge, multiplicity, scanAngles)
    orcaOutput: FilePath = p.join(conformerTorsionDir, "orca_scan.out")

    if not p.isfile(orcaOutput):
        run_orca(orcaInput, orcaOutput)

    scanDf: pd.DataFrame = read_scan_data(conformerTorsionDir)
    return scanDf
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def plot_torsion_scan(torsionScanDf, torsionTopDir):
    # Plotting
    
    
    # Plotting
    plt.figure(figsize=(10, 6))

    # Iterate over each column except 'Angle'
    for column in torsionScanDf.columns:
        if column != 'Angle':
            plt.plot(torsionScanDf['Angle'], torsionScanDf[column], label=column, marker='o')

    # Adding labels and title
    plt.xlabel('Angle')
    plt.ylabel('Energy')
    plt.title('Conformer Energies vs. Angle')
    plt.legend()

    plt.savefig(p.join(torsionTopDir, "torsion_scan.png"))
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def identify_rotatable_bonds(inputPdb) -> List[Tuple[int,int,int,int]]:

    dummyTorsions = [(8, 7, 10, 11), (2, 10, 9, 8)]
    return dummyTorsions

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def read_scan_data(conformerTorsionDir):
    scanDat: FilePath = p.join(conformerTorsionDir, "orca_scan.relaxscanact.dat")
    scanDf: pd.DataFrame = pd.read_csv(scanDat, delim_whitespace=True, header=None)
    return scanDf
    
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def run_orca(orcaInput, orcaOutput):
    orcaCommand = ["orca", orcaInput]
    with open(orcaOutput, 'w') as output_file:
        call(orcaCommand, stdout=output_file)
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def generate_orca_input(conformerPdb, torsionIndex, conformerTorsionDir, charge, multiplicity, scanAngles):
    orcaInputFile = p.join(conformerTorsionDir, "orca_scan.inp")
    with open(orcaInputFile, "w") as f:
        f.write("! XTB2 Opt\n")
        f.write("%geom Scan\n")
        torsionText: str = f"D {' '.join(map(str, torsionIndex))} = {scanAngles}\n"
        f.write(torsionText)
        f.write("end\n")
        f.write("end\n")
        f.write(f"*xyz {charge} {multiplicity}\n\n")
        
        conformerDf = pdbUtils.pdb2df(conformerPdb)
        for index, row in conformerDf.iterrows():
            geomLine = f"{row['ELEMENT']} {row['X']} {row['Y']} {row['Z']}\n"
            f.write(geomLine)
        f.write("*")

    return orcaInputFile
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

if __name__ == "__main__":
    main()