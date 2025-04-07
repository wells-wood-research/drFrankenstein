import os
from os import path as p
import pandas as pd
from subprocess import call, PIPE
from shutil import copy


################################################################################

def create_the_monster(config):
    ## create a runtimeInfo 
    config["runtimeInfo"]["madeByCreator"] = {}

    ## make a dir
    outputDir = config["pathInfo"]["outputDir"]
    finalCreationDir = p.join(outputDir, "06_final_creation")
    os.makedirs(finalCreationDir, exist_ok=True)
    config["runtimeInfo"]["madeByCreator"]["finalCreationDir"] = finalCreationDir

    ##TODO: this could be done during the capping step instead
    cappingAtomIds = get_capping_atom_ids(config)

    create_final_lib_and_mol2(cappingAtomIds, config)
    copy_final_frcmod(config)

################################################################################
def copy_final_frcmod(config):
    finalCreationDir = config["runtimeInfo"]["madeByCreator"]["finalCreationDir"]
    mmTorsionCalculationDir = config["runtimeInfo"]["madeByStitching"]["mmTorsionCalculationDir"]
    moleculeName = config["moleculeInfo"]["moleculeName"]

    frcmodFile = p.join(mmTorsionCalculationDir, f"{moleculeName}.frcmod")
    finalFrcmod = p.join(finalCreationDir, f"{moleculeName}.frcmod")

    copy(frcmodFile, finalFrcmod)
    
################################################################################
def get_capping_atom_ids(config):
    cappedMol2 = config["runtimeInfo"]["madeByStitching"]["finalMol2"]

    cappingHeteroAtomNames = ["NN", "CN", "CC1", "OC", "CC2"]
    atomDf, bondDf  = parse_mol2(cappedMol2)
    cappingHeteroAtomIds = atomDf[atomDf["ATOM_NAME"].isin(cappingHeteroAtomNames)]["ATOM_ID"].to_list()
    cappingProtonIds = get_capping_proton_ids(cappingHeteroAtomNames, atomDf, bondDf)

    cappingAtomIds = cappingHeteroAtomIds + cappingProtonIds

    return cappingAtomIds
################################################################################

def create_final_lib_and_mol2(cappingAtomIds, config):
    cappedMol2 = config["runtimeInfo"]["madeByStitching"]["finalMol2"]

    finalCreationDir = config["runtimeInfo"]["madeByCreator"]["finalCreationDir"]
    moleculeName = config["moleculeInfo"]["moleculeName"]
    finalMol2 = p.join(finalCreationDir, f"{moleculeName}.mol2")
    finalLib= p.join(finalCreationDir, f"{moleculeName}.lib")

    finalTleapInput = p.join(finalCreationDir, "tleap.in")
    
    with open(finalTleapInput, "w") as f:
        f.write("source leaprc.gaff2 \n")
        f.write(f"{moleculeName}  = loadmol2  {cappedMol2}\n")
        f.write(f"set {moleculeName} head {moleculeName}.1.N \n")
        f.write(f"set {moleculeName} tail {moleculeName}.1.C \n")
        for atomId in cappingAtomIds:
            f.write(f"remove {moleculeName} {moleculeName}.1.{atomId}\n")
        f.write(f"check {moleculeName} \n")
        f.write(f"saveoff  {moleculeName} {finalLib} \n")
        f.write(f"saveMol2 {moleculeName} {finalMol2} 1\n")
        f.write("quit")

    tleapOutput = p.join(finalCreationDir, f"tleap.out")

    tleapCommand: list = ["tleap", "-f", finalTleapInput, ">", tleapOutput]

    os.chdir(finalCreationDir)
    call(tleapCommand, stdout=PIPE)
################################################################################




def get_capping_proton_ids(cappingHeteroAtomNames, atomDf, bondDf):
    cappingProtonIds = []
    for cappingHeteroAtomName in cappingHeteroAtomNames:
        heteroAtomId = atomDf[atomDf["ATOM_NAME"]==cappingHeteroAtomName]["ATOM_ID"].to_list()[0]

        bondedToHeteroAtom_A = bondDf[(bondDf["ATOM_A_ID"] == heteroAtomId)]["ATOM_B_ID"].to_list()
        bondedToHeteroAtom_B = bondDf[(bondDf["ATOM_B_ID"] == heteroAtomId)]["ATOM_A_ID"].to_list()

        bondedToHeteroAtomIds = bondedToHeteroAtom_A + bondedToHeteroAtom_B
        protonIds = atomDf[atomDf["ATOM_ID"].isin(bondedToHeteroAtomIds) & atomDf["ATOM_TYPE"].str.startswith("h")]["ATOM_ID"].values
        cappingProtonIds.extend(protonIds)

    return cappingProtonIds

################################################################################

def parse_mol2(mol2File):
    atomData = []
    bondData = []
    readingAtoms = False
    readingBonds = False

    with open(mol2File, "r") as mol2:
        for line in mol2:
            if line.strip() == "":
                continue
            if line.startswith("@<TRIPOS>ATOM"):
                readingAtoms=True
                continue
            if line.startswith("@<TRIPOS>BOND"):
                readingBonds=True
                readingAtoms=False
                continue
            if line.startswith("@<TRIPOS>SUBSTRUCTURE"):
                break
            if readingAtoms:
                atomData.append(line.split())
            elif readingBonds:
                bondData.append(line.split())


        
    atomDataColumns = ["ATOM_ID", "ATOM_NAME", "X", "Y", "Z", "ATOM_TYPE", "RES_ID", "RES_NAME", "CHARGE"]
    atomDf = pd.DataFrame(atomData, columns=atomDataColumns)
    bondDataColumns = ["BOND_ID", "ATOM_A_ID", "ATOM_B_ID", "BOND_ORDER"]
    bondDf = pd.DataFrame(bondData, columns=bondDataColumns)

    return atomDf, bondDf
################################################################################

if __name__ == "__main__":
    raise NotImplementedError