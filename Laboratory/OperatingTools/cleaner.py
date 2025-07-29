import os
from os import path as p
from shutil import rmtree
from zipfile import ZipFile, ZIP_DEFLATED
## CLEAN CODE ##
class FilePath:
    pass
class DirectoryPath:
    pass


#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def clean_capping(config: dict) -> None:
    """
    Cleans up the capping directory
    Depends on cleanUpLevel in config

    Args:
        config (dict): dictionary containing all information
    Returns:
        None
    """
    cleanUpLevel = config["miscInfo"]["cleanUpLevel"]
 
    if cleanUpLevel == 0:
        return
    elif cleanUpLevel == 1:
        return
    elif cleanUpLevel in [2, 3]:
        cappingDir = config["runtimeInfo"]["madeByCapping"]["cappingDir"]
        moleculeName = config["moleculeInfo"]["moleculeName"]
        optDir = p.join(cappingDir, "geometry_optimisation")
        if p.isdir(optDir):
            cappedPdb = config["runtimeInfo"]["madeByCapping"]["cappedPdb"]
        keepFiles = [p.join(optDir, file) for file in os.listdir(optDir) if file.endswith(ext) for ext in [".inp", ".out", ".pdb"]]
        for file in os.listdir(optDir):
            if file in keepFiles:
                continue
            os.remove(p.join(optDir, file))
        wrongCappedPdb = p.join(cappingDir, f"{moleculeName}_capped.pdb")
        if p.isfile(wrongCappedPdb):
            os.remove(wrongCappedPdb)


#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def clean_assembly(config: dict) -> None:
    """
    Cleans up the assembly directory
    Depends on cleanUpLevel in config

    Args:
        config (dict): dictionary containing all information
    Returns:
        None
    """
    cleanUpLevel = config["miscInfo"]["cleanUpLevel"]
 
    if cleanUpLevel == 0:
        return
    elif cleanUpLevel == 1:
        return
    elif cleanUpLevel in [2, 3]:
        assemblyDir = config["runtimeInfo"]["madeByAssembly"]["assemblyDir"]
        filesToKeep = [file for file in os.listdir(assemblyDir) if "assembled" in file]
        for file in os.listdir(assemblyDir):
            if file in filesToKeep: 
                continue
            os.remove(p.join(assemblyDir, file))

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def clean_wriggle(config: dict) -> None:
    """
    Cleans up the wriggle directory
    Depends on cleanUpLevel in config

    Args:
        config (dict): dictionary containing all information
    Returns:
        None
    """
    cleanUpLevel = config["miscInfo"]["cleanUpLevel"]
 
    if cleanUpLevel == 0:
        return
    elif cleanUpLevel in [1, 2, 3]:
        conformerDir = config["runtimeInfo"]["madeByConformers"]["conformerDir"]
        moleculeName = config["moleculeInfo"]["moleculeName"]
        filesToKeep = [file for file in os.listdir(conformerDir) 
                       if file.startswith(f"{moleculeName}_conformer_") and file.endswith(".xyz")]
        filesToKeep.extend([file for file in os.listdir(conformerDir) if p.splitext(file)[1] in [".inp", ".out"]])
        for file in os.listdir(conformerDir):
            if file in filesToKeep: 
                continue
            os.remove(p.join(conformerDir, file))


#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def clean_up_opt_dir(optDir: DirectoryPath, config: dict):
    """
    For an opt dir:
    vital = orca_opt.xyz (for use in forwards)
    nice2have = orca_opt.out (for debugging)
    
    OPTIONS: 0, 1, 2, 3
 
    """
    cleanUpLevel = config["miscInfo"]["cleanUpLevel"]

    if cleanUpLevel == 0:
        return
    elif cleanUpLevel in [1, 2, 3]:
        keepFiles = ["orca_opt.xyz", "ORCA_FINISHED_NORMALLY", "ORCA_CRASHED"]
        for file in os.listdir(optDir):
            if file in keepFiles:
                continue
            elif cleanUpLevel in ["basic", "full"]:
                os.remove(p.join(optDir, file))                                                 
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def clean_up_scan_dir(scanDir: DirectoryPath, config: dict) -> None:
    """
    For a scan dir
    vital = orca_scan_XYZ.xyz, orca_scan.relaxscanact.dat, orca_scan.out
    """
    cleanUpLevel = config["miscInfo"]["cleanUpLevel"]

    if cleanUpLevel == 0:
        return
    elif cleanUpLevel in [1, 2, 3]:
        keepFiles = [file for file in os.listdir(scanDir) if  file.endswith(".xyz")]
        keepFiles.extend(["orca_scan.relaxscanact.dat",  ## need for getting energies
                        "orca_scan.inp",               ## need for getting angles
                        "ORCA_FINISHED_NORMALLY",      ## need as a positive flag for charge fitting later
                        "ORCA_CRASHED"])                ## need as a negative flag for charge fitting later

        for file in os.listdir(scanDir):
            if file in keepFiles:
                continue
            os.remove(p.join(scanDir, file))
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

def clean_up_singlepoint_dir(spDir: DirectoryPath, config: dict, keepGbw: bool = False) -> None:
    """
    Cleans up orca singlepoint calculation dir

    """

    cleanUpLevel = config["miscInfo"]["cleanUpLevel"]

    if cleanUpLevel == 0:
        return
    elif cleanUpLevel == 1:
        keepFiles = ["orca_sp.inp", "orca_sp.out", "ORCA_FINISHED_NORMALLY", "ORCA_CRASHED"]
        keepFiles.extend([file for file in os.listdir(spDir) if file.endswith(".molden.input")])
        keepFiles.extend([file for file in os.listdir(spDir) if file.endswith(".out")])

        if keepGbw:
            keepFiles.extend([file for file in os.listdir(spDir) if file.endswith(".gbw")])
        for file in os.listdir(spDir):
            if file in keepFiles:
                continue
            os.remove(p.join(spDir, file))
    elif cleanUpLevel in [2, 3]:
        rmtree(spDir)

def  clean_up_solvator_dir(solvatorDir: DirectoryPath, config: dict) -> None:
    """
    Cleans up an orca SOLVATOR calculation Dir
    """

    cleanUpLevel = config["miscInfo"]["cleanUpLevel"]

    if cleanUpLevel == 0:
        keepFiles = [file for file in os.listdir(solvatorDir)]
    elif cleanUpLevel == 1:
        keepFiles = [file for file in os.listdir(solvatorDir) if file.endswith(".inp") 
                     or file.endswith(".out") or file.endswith(".xyz")]
    elif cleanUpLevel in [2, 3]:
        keepFiles = [file for file in os.listdir(solvatorDir) if file.endswith(".xyz")]   

    for file in os.listdir(solvatorDir):
        if file in keepFiles:
            continue
        os.remove(p.join(solvatorDir, file))

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def clean_up_charges(config: dict) -> None:
    """
    Cleans up the charges directory
    Depends on cleanUpLevel in config

    Args:
        config (dict): dictionary containing all information
    Returns:
        None
    """
    cleanUpLevel = config["miscInfo"]["cleanUpLevel"]
 
    if cleanUpLevel == 0:
        return
    elif cleanUpLevel == 1:
        return
    elif cleanUpLevel in [2, 3]:
        chargesDir = config["runtimeInfo"]["madeByCharges"]["chargesDir"]
        chargeFittingDir = config["runtimeInfo"]["madeByCharges"]["chargeFittingDir"]
        filesToKeep = [file for file in os.listdir(chargeFittingDir) if file.endswith(".csv")]
        for file in os.listdir(chargesDir):
            if file in filesToKeep:
                continue
            os.remove(p.join(chargesDir, file))

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def clean_up_stitching(config: dict) -> None:
    """
    Cleans up parameter fitting directory
    Depends on cleanUpLevel in config

    Args:
        config (dict): dictionary containing all information
    Returns:
        None
    """
    cleanUpLevel = config["miscInfo"]["cleanUpLevel"]
    if cleanUpLevel == 0:
        return
    ## unpack config
    qmmmParameterFittingDir = config["runtimeInfo"]["madeByStitching"]["qmmmParameterFittingDir"]
    moleculeParameterDir = config["runtimeInfo"]["madeByStitching"]["moleculeParameterDir"]
    nShuffles = config["parameterFittingInfo"]["nShuffles"]

    if cleanUpLevel in [1, 2, 3]:
        print("\n")

        ## DEAL WITH EXTRA PNG FILES ##
        for torsionTag in os.listdir(qmmmParameterFittingDir):
            dirPath = p.join(qmmmParameterFittingDir, torsionTag)
            if not p.isdir(dirPath):
                continue
            ## ZIP PNG FILES ##
            pngsToZip = [p.join(qmmmParameterFittingDir, torsionTag, file) for file in os.listdir(p.join(qmmmParameterFittingDir, torsionTag)) if file.endswith(".png")]
            if cleanUpLevel == 1:
                zip_files(pngsToZip, p.join(qmmmParameterFittingDir, torsionTag, "PNG_FILES.zip"))
            ## delete extra PNG files, except the last and MAE ##
            pngsToDelete = [file for file in pngsToZip if not file.endswith(f"_{nShuffles}.png") and not file.endswith("mean_average_error.png")]
            for png in pngsToDelete:
                os.remove(p.join(qmmmParameterFittingDir, torsionTag, png))
        ## zip PRM / FRCMOD files
        paramFiles = [p.join(moleculeParameterDir, file) for file in os.listdir(moleculeParameterDir) 
                      if file.endswith(".prm") or file.endswith(".frcmod")]
        ## remove the last PRM / FRCMOD file
        paramFiles = [file for file in paramFiles if not p.splitext(file)[0].endswith(f"_{nShuffles}")]
                
        if cleanUpLevel == 1:
            zip_files(paramFiles, p.join(moleculeParameterDir, "PARAM_FILES.zip"))
        for file in paramFiles:
            os.remove(p.join(moleculeParameterDir, file))
    ########################

 
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def zip_files(files: list[FilePath], outZip:FilePath):
    """
    Zips a list of files into a single zip file
    Args:
        files (list): list of files to zip
        outZip (str): path to output zip file
    Returns:
        None
    """
    with ZipFile(outZip, "w", ZIP_DEFLATED) as zip:
        for file in files:
            zip.write(file, p.basename(file))
   
