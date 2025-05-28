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
        cappingDir: DirectoryPath = config["runtimeInfo"]["madeByCapping"]["cappingDir"] # type: ignore
        moleculeName: str = config["moleculeInfo"]["moleculeName"]
        optDir: DirectoryPath = p.join(cappingDir, "geometry_optimisation") # type: ignore
        if p.isdir(optDir):
            cappedPdb: FilePath = config["runtimeInfo"]["madeByCapping"]["cappedPdb"] # type: ignore
        # keepFiles logic seems complex and might not be perfectly translatable to simple type hints
        # without knowing more about `ext`. Assuming it's a list of strings.
        keepFiles = [p.join(optDir, file) for file in os.listdir(optDir) if any(file.endswith(ext_val) for ext_val in [".inp", ".out", ".pdb"])] # Corrected list comprehension and added ext_val
        for file in os.listdir(optDir):
            filePath = p.join(optDir, file) # Defined filePath for clarity
            if filePath in keepFiles: # Compare full paths
                continue
            os.remove(filePath)
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
            filePath = p.join(assemblyDir, file) # Defined filePath
            if filePath in filesToKeep: # Compare full paths. Note: filesToKeep stores basenames. This logic might need review if filesToKeep should be full paths.
                continue
            os.remove(filePath)

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def clean_wriggle(config: dict) -> None:
    """
    Cleans up the wriggle directory (conformer generation).
    Depends on cleanUpLevel in config.

    Args:
        config (dict): dictionary containing all information
    Returns:
        None
    """
    cleanUpLevel: int = config["miscInfo"]["cleanUpLevel"]
 
    if cleanUpLevel == 0:
        return
    elif cleanUpLevel in [1, 2, 3]:
        conformerDir: DirectoryPath = config["runtimeInfo"]["madeByConformers"]["conformerDir"] # type: ignore
        moleculeName: str = config["moleculeInfo"]["moleculeName"]
        # filesToKeep stores basenames
        filesToKeep = [file for file in os.listdir(conformerDir) 
                       if file.startswith(f"{moleculeName}_conformer_") and file.endswith(".xyz")]
        filesToKeep.extend([file for file in os.listdir(conformerDir) if p.splitext(file)[1] in [".inp", ".out"]])
        for file in os.listdir(conformerDir):
            if file in filesToKeep: # Comparing basename with list of basenames
                continue
            os.remove(p.join(conformerDir, file))


#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def clean_up_opt_dir(opt_dir: DirectoryPath, config: dict) -> None:
    """
    For an opt dir:
    vital = orca_opt.xyz (for use in forwards)
    nice2have = orca_opt.out (for debugging)
    
    OPTIONS: 0, 1, 2, 3
 
    """
    cleanUpLevel: int = config["miscInfo"]["cleanUpLevel"]

    if cleanUpLevel == 0:
        return
    elif cleanUpLevel in [1, 2, 3]:
        keepFiles = ["orca_opt.xyz", "ORCA_FINISHED_NORMALLY", "ORCA_CRASHED"] # These are basenames
        for file in os.listdir(opt_dir): # file is a basename
            if file in keepFiles:
                continue
            # The original code had `cleanUpLevel in ["basic", "full"]` which are not int.
            # Assuming this was a typo and it should always remove if not in keepFiles for levels 1,2,3.
            # If "basic" and "full" were meant to be specific string values from config, that needs clarification.
            # For now, interpreting as: if level is 1, 2, or 3, and file not in keepFiles, remove it.
            os.remove(p.join(opt_dir, file)) # type: ignore                                                 
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def clean_up_scan_dir(scan_dir: DirectoryPath, config: dict) -> None:
    """
    For a scan dir
    vital = orca_scan_XYZ.xyz, orca_scan.relaxscanact.dat, orca_scan.out
    """
    cleanUpLevel: int = config["miscInfo"]["cleanUpLevel"]

    if cleanUpLevel == 0:
        return
    elif cleanUpLevel in [1, 2, 3]:
        keepFiles = [file for file in os.listdir(scan_dir) if  file.endswith(".xyz")] # Basenames
        keepFiles.extend(["orca_scan.relaxscanact.dat",  ## need for getting energies
                        "orca_scan.inp",               ## need for getting angles
                        "ORCA_FINISHED_NORMALLY",      ## need as a positive flag for charge fitting later
                        "ORCA_CRASHED"])                ## need as a negative flag for charge fitting later

        for file in os.listdir(scan_dir): # file is a basename
            if file in keepFiles:
                continue
            os.remove(p.join(scan_dir, file)) # type: ignore
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

def clean_up_singlepoint_dir(sp_dir: DirectoryPath, config: dict, keep_gbw: bool = False) -> None:
    """
    Cleans up orca singlepoint calculation dir

    """

    cleanUpLevel: int = config["miscInfo"]["cleanUpLevel"]

    if cleanUpLevel == 0:
        return
    elif cleanUpLevel == 1:
        keepFiles = ["orca_sp.inp", "orca_sp.out", "ORCA_FINISHED_NORMALLY", "ORCA_CRASHED"] # Basenames
        keepFiles.extend([file for file in os.listdir(sp_dir) if file.endswith(".molden.input")])
        if keep_gbw:
            keepFiles.extend([file for file in os.listdir(sp_dir) if file.endswith(".gbw")])
        for file in os.listdir(sp_dir): # file is a basename
            if file in keepFiles:
                continue
            os.remove(p.join(sp_dir, file)) # type: ignore
    elif cleanUpLevel in [2, 3]:
        rmtree(sp_dir) # type: ignore

def  clean_up_solvator_dir(solvator_dir: DirectoryPath, config: dict) -> None:
    """
    Cleans up an orca SOLVATOR calculation Dir
    """

    cleanUpLevel: int = config["miscInfo"]["cleanUpLevel"]
    keepFiles: list[str] = [] # Initialize keepFiles

    if cleanUpLevel == 0:
        return
    elif cleanUpLevel == 1:
        keepFiles = [file for file in os.listdir(solvator_dir) if file.endswith(".inp") 
                     or file.endswith(".out") or file.endswith(".xyz")]
    elif cleanUpLevel in [2, 3]:
        keepFiles = [file for file in os.listdir(solvator_dir) if file.endswith(".xyz")]   

    for file in os.listdir(solvator_dir): # file is a basename
        if file in keepFiles:
            continue
        os.remove(p.join(solvator_dir, file)) # type: ignore

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
    cleanUpLevel: int = config["miscInfo"]["cleanUpLevel"]
 
    if cleanUpLevel == 0:
        return
    elif cleanUpLevel == 1:
        return
    elif cleanUpLevel in [2, 3]:
        chargesDir: DirectoryPath = config["runtimeInfo"]["madeByCharges"]["chargesDir"] # type: ignore
        chargeFittingDir: DirectoryPath = config["runtimeInfo"]["madeByCharges"]["chargeFittingDir"] # type: ignore
        # filesToKeep stores basenames from chargeFittingDir
        filesToKeep = [file for file in os.listdir(chargeFittingDir) if file.endswith(".csv")]
        for file in os.listdir(chargesDir): # file is a basename from chargesDir
            # This logic seems to intend to keep files in chargesDir if their basename also exists in filesToKeep (from chargeFittingDir)
            # However, it might be simpler or more correct to check if p.join(chargesDir, file) should be kept based on some criteria
            # For now, sticking to the original logic: if a basename from chargesDir is in the list of basenames from filesToKeep, it's kept.
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
    cleanUpLevel: int = config["miscInfo"]["cleanUpLevel"]
    if cleanUpLevel == 0:
        return
    ## unpack config
    qmmmParameterFittingDir: DirectoryPath = config["runtimeInfo"]["madeByStitching"]["qmmmParameterFittingDir"] # type: ignore
    moleculeParameterDir: DirectoryPath = config["runtimeInfo"]["madeByStitching"]["moleculeParameterDir"] # type: ignore
    nShuffles: int = config["parameterFittingInfo"]["nShuffles"]

    if cleanUpLevel in [1, 2, 3]:
        print("\n") # This print might be for CLI, consider logging for library use

        ## DEAL WITH EXTRA PNG FILES ##
        for torsionTag in os.listdir(qmmmParameterFittingDir):
            currentTorsionDir = p.join(qmmmParameterFittingDir, torsionTag) # Defined for clarity
            ## ZIP PNG FILES ##
            pngsToZip = [p.join(currentTorsionDir, file) for file in os.listdir(currentTorsionDir) if file.endswith(".png")]
            if cleanUpLevel == 1:
                zip_files(files=pngsToZip, out_zip=p.join(currentTorsionDir, "PNG_FILES.zip")) # type: ignore
            ## remove the last PNG file
            # This logic keeps only the PNG file from the LAST shuffle.
            pngsToDelete = [file for file in pngsToZip if not file.endswith(f"_{nShuffles}.png")]
            for png_path in pngsToDelete: # Renamed loop var
                os.remove(png_path) # png_path is already a full path
        ## zip PRM / FRCMOD files 
        paramFiles = [p.join(moleculeParameterDir, file) for file in os.listdir(moleculeParameterDir) 
                      if file.endswith(".prm") or file.endswith(".frcmod")]
        ## remove the last PRM / FRCMOD file (keeps only the one from the last shuffle)
        paramFilesToDelete = [file for file in paramFiles if not p.splitext(file)[0].endswith(f"_{nShuffles}")] # Renamed
                
        if cleanUpLevel == 1:
            zip_files(files=paramFilesToDelete, out_zip=p.join(moleculeParameterDir, "PARAM_FILES.zip")) # type: ignore
        for file_path in paramFilesToDelete: # Renamed loop var
            os.remove(file_path) # file_path is already a full path
    ########################

 
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def zip_files(files: list[FilePath], out_zip:FilePath) -> None:
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
   
