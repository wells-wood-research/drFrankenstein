import os
from os import path as p
import warnings
import parmed as pmd
from parmed.charmm import CharmmParameterSet, CharmmPsfFile
from parmed.exceptions import ParameterWarning
from psfgen import PsfGen

# Suppress ParameterWarning
warnings.filterwarnings('ignore', category=ParameterWarning)


from . import Assembly_Monster
from . import Assembly_Assistant
# Placeholder classes (extend if needed)
class FilePath:
    pass

class DirectoryPath:
    pass



def charmm_assembly_protocol(config:dict) -> dict:
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
    config["runtimeInfo"]["madeByAssembly"] = {}

    ## make a dir for assembly of initial parameters
    assemblyDir = p.join(config["pathInfo"]["outputDir"], "02_assembly")
    os.makedirs(assemblyDir, exist_ok=True)
    config["runtimeInfo"]["madeByAssembly"]["assemblyDir"] = assemblyDir
    ## Split STR file into RTF and PRM
    cappedRtf, cappedPrm = Assembly_Monster.split_charmm_str(config)
    ## Find default CHARMM parameters in drFrankenstein SRC folder
    charmmDefaultParams = Assembly_Assistant.find_default_charmm_parameters()
    ## create PSF file using RTF from capped RTF and CGenFF RTF
    cappedPsf = Assembly_Monster.make_charmm_psf(cappedRtf,
                                                  charmmDefaultParams["cgenffRtf"],
                                                    config)

    parmedPsf = Assembly_Assistant.load_psf_with_params(psfFile = cappedPsf,
                                                        params = (charmmDefaultParams["cgenffRtf"],
                                                                charmmDefaultParams["cgenffPrm"],
                                                                cappedRtf, cappedPrm))
    
    parmedPsf, nameToCgenffType = Assembly_Monster.set_backbone_types_psf(parmedPsf,
                                                                           config)

    missingPrm = Assembly_Monster.create_missing_prm(parmedPsf = parmedPsf,
                                                cappedRtf = cappedRtf,
                                                cappedPrm = cappedPrm,
                                                charmmDefaultParams = charmmDefaultParams,
                                                nameToCgenffType = nameToCgenffType,
                                                config = config)
    ##NOTE try to skip a couple of steps here - can't remember why I did this!
    completeParameterSet = CharmmParameterSet(cappedRtf, cappedPrm, missingPrm, *charmmDefaultParams.values())
    parmedPsf.load_parameters(completeParameterSet)


    config = Assembly_Assistant.save_modified_parameter_files(parmedPsf, config)

    return config