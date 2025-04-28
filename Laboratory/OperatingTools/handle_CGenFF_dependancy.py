import os
from os import path as p

from OperatingTools import drSplash
from OperatingTools import file_parsers

def handle_CGenFF_dependancy(config):
    if config["pathInfo"]["cgenffExe"] == None:
        inputDir = config["pathInfo"]["inputDir"]
        moleculeName = config["moleculeInfo"]["moleculeName"]
        moleculeStr = p.join(inputDir, f"{moleculeName}_capped.str")
        if p.isfile(moleculeStr):
            print("WONK")
            config["runtimeInfo"]["madeByCapping"]["cappedMoleculeStr"] = moleculeStr
            return config
        else:

            cappedPdb = config["runtimeInfo"]["madeByCapping"]["cappedPdb"]
            cappingDir = config["runtimeInfo"]["madeByCapping"]["cappingDir"]
            cappedMol2 = p.join(cappingDir, f"{moleculeName}_capped.mol2")
            config["runtimeInfo"]["madeByCapping"]["cappedMol2"] = cappedMol2
            file_parsers.pdb2mol2(cappedPdb, cappedMol2)
            rename_molecule_mol2(cappedMol2, moleculeName)

            drSplash.show_need_cgenff_str(cappedMol2)
            exit()
    else:  
        print("SPLONK")
        return config


def rename_molecule_mol2(mol2File, moleculeName):
    tmpMol2 = p.join(p.dirname(mol2File), f"tmp.mol2")
    with open(mol2File, "r") as inMol2, open(tmpMol2, "w") as outMol2:
        moleculeNextLine = False
        for line in inMol2:
            if line.startswith("@<TRIPOS>MOLECULE"):
                outMol2.write(line)
                moleculeNextLine = True
            elif moleculeNextLine:  
                outMol2.write(f"{moleculeName}\n")
                moleculeNextLine = False
            else:
                outMol2.write(line)
    os.replace(tmpMol2, mol2File)