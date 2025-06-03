import os
from os import path as p
import warnings
from parmed.charmm import CharmmParameterSet
from . import Assembly_Monster
from . import Assembly_Assistant

class FilePath:
    pass

class DirectoryPath:
    pass

def amber_assembly_protocol(config: dict) -> dict:
    """
    Main protocol for assembly of starting parameters for AMBER parameters
    This protocol:
    1. Make a MOl2 with antechamber
    2. Make a FRCMOD file with tleap
    3. Optionally maps backbone atoms from gaff2 to AMBER atoms
    """
    moleculeName = config['moleculeInfo']['moleculeName']
    amberHome = config['pathInfo']['amberHome']
    (gaff2Dat, parm19Dat) = Assembly_Assistant.find_default_amber_parameters(amberHome)
    config['runtimeInfo']['madeByAssembly'] = {}
    assemblyDir = p.join(config['pathInfo']['outputDir'], '02_parameter_assembly')
    os.makedirs(assemblyDir, exist_ok=True)
    config['runtimeInfo']['madeByAssembly']['assemblyDir'] = assemblyDir
    moleculeMol2 = p.join(assemblyDir, f'{moleculeName}_capped.mol2')
    Assembly_Assistant.pdb2mol2(config['runtimeInfo']['madeByCapping']['cappedPdb'], moleculeMol2, assemblyDir, config)
    if config['torsionScanInfo']['preserveBackboneTorsions']:
        Assembly_Monster.change_backbone_types_amber(moleculeMol2, config)
    moleculeFrcmod = p.join(assemblyDir, f'{moleculeName}_capped.frcmod')
    Assembly_Assistant.create_frcmod_file(moleculeMol2, moleculeFrcmod, gaff2Dat)
    Assembly_Monster.replace_parameters(moleculeFrcmod, parm19Dat)
    moleculePrmtop = Assembly_Assistant.run_tleap_to_make_params(moleculeMol2, moleculeFrcmod, assemblyDir, 'amber')
    config['runtimeInfo']['madeByAssembly']['cappedMol2'] = moleculeMol2
    config['runtimeInfo']['madeByAssembly']['assembledFrcmod'] = moleculeFrcmod
    config['runtimeInfo']['madeByAssembly']['assembledPrmtop'] = moleculePrmtop
    return config

def charmm_assembly_protocol(config: dict) -> dict:
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
    config['runtimeInfo']['madeByAssembly'] = {}
    assemblyDir = p.join(config['pathInfo']['outputDir'], '02_parameter_assembly')
    os.makedirs(assemblyDir, exist_ok=True)
    config['runtimeInfo']['madeByAssembly']['assemblyDir'] = assemblyDir
    (cappedRtf, cappedPrm) = Assembly_Monster.split_charmm_str(config)
    charmmDefaultParams = Assembly_Assistant.find_default_charmm_parameters()
    cappedPsf = Assembly_Monster.make_charmm_psf(cappedRtf, charmmDefaultParams['cgenffRtf'], config)
    parmedPsf = Assembly_Assistant.load_psf_with_params(psfFile=cappedPsf, params=(charmmDefaultParams['cgenffRtf'], charmmDefaultParams['cgenffPrm'], cappedRtf, cappedPrm))
    (parmedPsf, nameToCgenffType, nameToDesiredType) = Assembly_Monster.set_backbone_types_psf(parmedPsf, config)
    missingPrm = Assembly_Monster.create_missing_prm(parmedPsf=parmedPsf, cappedRtf=cappedRtf, cappedPrm=cappedPrm, charmmDefaultParams=charmmDefaultParams, nameToCgenffType=nameToCgenffType, config=config)
    completeParameterSet = CharmmParameterSet(cappedRtf, cappedPrm, missingPrm, *charmmDefaultParams.values())
    parmedPsf.load_parameters(completeParameterSet)
    config = Assembly_Assistant.save_modified_parameter_files(parmedPsf, config)
    config = Assembly_Assistant.update_rtf_types(cappedRtf, nameToDesiredType, config)
    return config