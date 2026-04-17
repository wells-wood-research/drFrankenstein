import os
from os import path as p

from parmed.charmm import CharmmParameterSet

from .Agnostic_Protocol import Agnostic_Assembly_Monster
from .Agnostic_Protocol import Agnostic_Assembly_Assistant
from .Antechamber_Protocol import Antechamber_Assembly_Monster
from .CGenFF_Protocol import CGenFF_Assembly_Monster
from .CGenFF_Protocol import CGenFF_Assembly_Assistant
from .Antechamber_Protocol import Antechamber_Assembly_Assistant
from OperatingTools import drLogger


class FilePath:
    pass

class DirectoryPath:
    pass



@drLogger.experiment_logger("Parameter Assembly")
def parameter_assembly_protocol(config:dict) -> dict:
    assemblyProtocol = config["miscInfo"]["assemblyProtocol"]

    if assemblyProtocol.upper() == "AGNOSTIC":
        config = agnostic_assembly_protocol(config)
    elif assemblyProtocol.upper() == "ANTECHAMBER":
        config = amber_assembly_protocol(config)
    elif assemblyProtocol.upper() == "CGENFF":
        config = CGenFF_assembly_protocol(config)

    config.setdefault("checkpointInfo", {})["assemblyComplete"] = True
    return config
    

def amber_assembly_protocol(config:dict) -> dict:
    """
    Main protocol for assembly of starting parameters for AMBER parameters
    This protocol:
    1. Make a MOl2 with antechamber
    2. Make a FRCMOD file with tleap
    3. Optionally maps backbone atoms from gaff2 to AMBER atoms
    """

    ## unpack config
    moleculeName = config["moleculeInfo"]["moleculeName"]
    amberHome = config["pathInfo"]["amberHome"]

    ## find amber parameter sets
    gaff2Dat, parm19Dat = Antechamber_Assembly_Assistant.find_default_amber_parameters(amberHome)

    ## create an empty dict in config
    config["runtimeInfo"]["madeByAssembly"] = {}

    ## make a dir for assembly of initial parameters
    assemblyDir = p.join(config["pathInfo"]["outputDir"], "04_parameter_assembly")
    os.makedirs(assemblyDir, exist_ok=True)
    config["runtimeInfo"]["madeByAssembly"]["assemblyDir"] = assemblyDir

    ## create a mol2 file with antechamber
    moleculeMol2 = p.join(assemblyDir, f"{moleculeName}_capped.mol2")
    Antechamber_Assembly_Assistant.make_mol2_with_antechamber(config["runtimeInfo"]["madeByCapping"]["cappedPdb"], moleculeMol2, assemblyDir, config)

    if config["torsionScanInfo"]["runScansOn"]["phiPsi"]:
        ## reset capping atom types to gaff2 defaults
        Antechamber_Assembly_Monster.change_capping_types_amber(moleculeMol2, config)
    else:
        ## reset backbone atom types to AMBER defaults
        Antechamber_Assembly_Monster.change_backbone_types_amber(moleculeMol2, config)

    ## create frcmod file
    moleculeFrcmod = p.join(assemblyDir, f"{moleculeName}_capped.frcmod")
    Antechamber_Assembly_Assistant.create_frcmod_file(moleculeMol2, moleculeFrcmod, gaff2Dat)
    ## replace gaff2 params with amber19sb parameters
    Antechamber_Assembly_Monster.replace_parameters(moleculeFrcmod, parm19Dat)
    ## create a PRMTOP file
    moleculePrmtop = Antechamber_Assembly_Assistant.run_tleap_to_make_params(moleculeMol2, moleculeFrcmod, assemblyDir, "amber")

    ## update config
    config["runtimeInfo"]["madeByAssembly"]["cappedMol2"] = moleculeMol2
    config["runtimeInfo"]["madeByAssembly"]["assembledFrcmod"] = moleculeFrcmod
    config["runtimeInfo"]["madeByAssembly"]["assembledPrmtop"] = moleculePrmtop

    return config

def agnostic_assembly_protocol(config:dict) -> dict:
    """
    Main protocol for assembly of starting parameters for agnostic parameters
    This protocol:
    [] Find .hess files from charges calculation (might depend on protocol used for charges)
    [] loop through .hess files
        [] construct hessian matrix from .hess file
        [] invert to get compliance matrix
        [] construct adjacancy and angle matrix frpm cappedPdb
        [] convert adjacency and angle matrix to bond and angle lists
        [] calculate distance and angle strengths (k values)
        [] extract optimised bond and angle values from cappedPdb

    [] average over r0 and k0 values for each bond and angle
    [] create a dict with atomIds as keys and:
        [] Element
        [] nBonds
        [] Charge
        [] r0 values
        [] k0 values
    
    """

    ## create an empty dict in config
    config["runtimeInfo"]["madeByAssembly"] = {}

    ## make a dir for assembly of initial parameters
    assemblyDir = p.join(config["pathInfo"]["outputDir"], "04_parameter_assembly")
    os.makedirs(assemblyDir, exist_ok=True)
    config["runtimeInfo"]["madeByAssembly"]["assemblyDir"] = assemblyDir

    ## find hessian files
    hessianFiles, xyzFiles = Agnostic_Assembly_Assistant.find_hess_xyz_files(config)

    bonds, angles, dihedrals, impropers = Agnostic_Assembly_Monster.get_molecule_connectivity(config)

    allKbValues = []
    allKThetaValues = []
    allr0values = []
    allTheta0values = []
    for hessianFile, xyzFile in zip(hessianFiles, xyzFiles):
        kbValues, r0values, kThetaValues,  theta0values = Agnostic_Assembly_Monster.process_hessian(hessianFile=hessianFile,
                                                                                                bonds=bonds,
                                                                                                angles=angles,
                                                                                                xyzFile=xyzFile,
                                                                                                config=config)
        allKbValues.append(kbValues)
        allr0values.append(r0values)
        allKThetaValues.append(kThetaValues)
        allTheta0values.append(theta0values)

    kbValues = Agnostic_Assembly_Monster.average_dict_values(allKbValues)
    r0values = Agnostic_Assembly_Monster.average_dict_values(allr0values)
    kThetaValues = Agnostic_Assembly_Monster.average_dict_values(allKThetaValues)
    theta0values = Agnostic_Assembly_Monster.average_dict_values(allTheta0values)

    atomFeaturesDf, atomTypesDict = Agnostic_Assembly_Monster.construct_atom_features(kbValues, r0values, kThetaValues, theta0values, config=config)


    atomTypesDf = Agnostic_Assembly_Monster.assign_atom_types_by_clustering(atomFeaturesDf)

    atomTypesDf = Agnostic_Assembly_Monster.reassign_capping_types(atomTypesDf, config=config)

    if not config["torsionScanInfo"]["runScansOn"]["phiPsi"]:
        atomTypesDf = Agnostic_Assembly_Monster.reassign_backbone_types(atomTypesDf, config=config)

    atomTypeMasses = Agnostic_Assembly_Monster.get_atom_type_masses(atomTypesDf)

    (bondParamsByType,
      angleParamsByType,
        groupedDihedrals,
          groupedImpropers) = Agnostic_Assembly_Monster.group_params_by_atom_type(atomTypesDf,
                                                                          kbValues, r0values,
                                                                            kThetaValues, theta0values, 
                                                                                dihedrals, impropers)
    
    nonBondedParams = Agnostic_Assembly_Monster.assign_non_bonded_by_analogy(atomTypesDf)

    config = Agnostic_Assembly_Monster.construct_frcmod(atomTypeMasses,
                                        bondParamsByType, 
                                        angleParamsByType,
                                           groupedDihedrals, groupedImpropers,
                                           nonBondedParams,
                                                 config=config)

    config = Agnostic_Assembly_Monster.construct_mol2(atomTypesDf, config=config)

    config = Agnostic_Assembly_Assistant.run_tleap_to_make_params(config)

    return config

def CGenFF_assembly_protocol(config:dict) -> dict:
    """
    Main protocol for assembly of starting parameters for CHARMM parameters
    Before running this, a stream (.STR) file has to have been created 
    using CGenFF.
    This protocol:
    1. Split STR -> PRM + RTF
    2. Create PSF file from PRM and RTF file
    3. Optionally maps backbone atoms to CHARMM36m atoms
    
    Args: 
        config (dict): contains all run info
    Returns:
        config (dict): updated config   
    """
    ## create an empty dict in config
    config["runtimeInfo"]["madeByAssembly"] = {}
    runScansOn = config["torsionScanInfo"]["runScansOn"]

    ## make a dir for assembly of initial parameters
    assemblyDir = p.join(config["pathInfo"]["outputDir"], "04_parameter_assembly")
    os.makedirs(assemblyDir, exist_ok=True)
    config["runtimeInfo"]["madeByAssembly"]["assemblyDir"] = assemblyDir
    ## Split STR file into RTF and PRM
    cappedRtf, cappedPrm = CGenFF_Assembly_Monster.split_charmm_str(config)
    ## Find default CHARMM parameters in drFrankenstein SRC folder
    charmmDefaultParams = CGenFF_Assembly_Assistant.find_default_charmm_parameters()
    ## create PSF file using RTF from capped RTF and CGenFF RTF
    cappedPsf = CGenFF_Assembly_Monster.make_charmm_psf(rtfFiles=[cappedRtf,
                                                  charmmDefaultParams["cgenffRtf"]], config = config, 
                                                  suffix = "capped"
                                                    )

    parmedPsf = CGenFF_Assembly_Assistant.load_psf_with_params(psfFile = cappedPsf,
                                                        params = (charmmDefaultParams["cgenffRtf"],
                                                                charmmDefaultParams["cgenffPrm"],
                                                                cappedRtf, cappedPrm))
    

    parmedPsf, nameToCgenffType, nameToDesiredType = CGenFF_Assembly_Monster.set_backbone_types_psf(parmedPsf,
                                                                        config)


    missingPrm = CGenFF_Assembly_Monster.create_missing_prm(parmedPsf = parmedPsf,
                                                cappedRtf = cappedRtf,
                                                cappedPrm = cappedPrm,
                                                charmmDefaultParams = charmmDefaultParams,
                                                nameToCgenffType = nameToCgenffType,
                                                config = config)
    
    completeParameterSet = CharmmParameterSet(cappedRtf, cappedPrm, missingPrm, *charmmDefaultParams.values())

    ## load parameters back in to parmed PSF  
    parmedPsf.load_parameters(completeParameterSet)

    assembledPrm, assembledPsf = CGenFF_Assembly_Assistant.save_modified_prm_file(parmedPsf, config)
    assembledRtf = CGenFF_Assembly_Assistant.update_rtf_types(cappedRtf, nameToDesiredType, config)
    assembledPsf = CGenFF_Assembly_Monster.make_charmm_psf(rtfFiles = [assembledRtf,], config = config, suffix = "assembled")

    config["runtimeInfo"]["madeByAssembly"]["assembledPrm"] = assembledPrm
    config["runtimeInfo"]["madeByAssembly"]["assembledRtf"] = assembledRtf
    config["runtimeInfo"]["madeByAssembly"]["assembledPsf"] = assembledPsf


    return config
