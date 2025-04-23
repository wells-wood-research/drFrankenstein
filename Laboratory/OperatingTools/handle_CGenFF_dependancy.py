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
            drSplash.show_need_cgenff_str(cappedMol2)
            exit()
    else:  
        print("SPLONK")
        return config
