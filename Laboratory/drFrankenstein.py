import os
from os import path as p
import traceback
import ruamel.yaml as ruamel
## drFRANKENSTEIN LIBRARIES ##
from Experiments.Protocol_1_Capping import Capping_Doctor
from Experiments.Protocol_2_Assembly import Assembly_Doctor
from Experiments.Protocol_3_Wriggling import Wriggling_Doctor
from Experiments.Protocol_4_Twisting import Twisted_Doctor
from Experiments.Protocol_5_Charging import Charged_Doctor
from Experiments.Protocol_6_Stitching import Stitching_Doctor
from Experiments.Protocol_7_Creation import drCreator
from Experiments.Protocol_8_Reporter import Reporting_Doctor

from OperatingTools import drYaml
from OperatingTools import drSplash
# from OperatingTools import Timer # Timer is not directly used in this file
from OperatingTools import validate_config
from OperatingTools import handle_CGenFF_dependancy
from typing import Dict, Any, Optional, List # Added for type hints
from ruamel.yaml.comments import CommentedMap # For config object type

## CLEAN CODE ##
class FilePath: # Assuming FilePath might be a str or an object with a path attribute
    pass
class DirectoryPath: # Assuming DirectoryPath might be a str or an object
    pass

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def main() -> None:
    """
    Main protocol for drFrankenstein

    Args:
        None
    Returns:
        None
    
    """
    ## get config yaml file from argpass command-line-argument
    configYaml: FilePath = drYaml.get_config_input_arg()
    ## load into dict, check for bad formatting
    config: CommentedMap = drYaml.read_input_yaml(config_file=configYaml) # type: ignore
    ## check config for errors
    ##TODO: re-do once config is settled
    validated_config = validate_config.validate_config(config=config) # type: ignore
    if validated_config is None:
        # Validation failed, errors already printed by validate_config
        return # Or handle error appropriately
    config = validated_config


    ## unpack config to find outputDir, make directory
    outputDir: DirectoryPath = config["pathInfo"]["outputDir"] # type: ignore
    os.makedirs(outputDir,exist_ok=True)
    ## deal with checkpointing, lets us skip steps if already done
    config = drYaml.read_config_with_checkpoints(config=config, out_dir=outputDir)
    config = drYaml.init_config_checkpoints(config=config, out_dir=outputDir)

    ## initialise runtimeInfo
    config = drYaml.initialize_runtime_info(config=config) # Corrected call
    ## save config back to yaml
    drYaml.write_config_to_yaml(config=config, out_dir=outputDir)

    ## add capping groups to input molecule
    ## TODO: cope with staples etc. that need 2 or more capping groups
    checkpointInfo: Dict[str, bool] = config["checkpointInfo"] # type: ignore
    if not checkpointInfo["cappingComplete"]:
        print("Running Capping Protocol")
        config = Capping_Doctor.capping_protocol(config=config)
        drYaml.write_config_to_yaml(config=config, out_dir=outputDir)

    ## run assembly protocol
    if not checkpointInfo["assemblyComplete"]:
        if config["parameterFittingInfo"]["forceField"] == "CHARMM": # type: ignore
            config = handle_CGenFF_dependancy.handle_cgenff_dependency(config=config) # Corrected call
            drYaml.write_config_to_yaml(config=config, out_dir=outputDir)
            config = Assembly_Doctor.charmm_assembly_protocol(config=config)
            drYaml.write_config_to_yaml(config=config, out_dir=outputDir)  
        elif config["parameterFittingInfo"]["forceField"] == "AMBER": # type: ignore
            config = Assembly_Doctor.amber_assembly_protocol(config=config)
            drYaml.write_config_to_yaml(config=config, out_dir=outputDir)

    ## run conformer generation protocol
    if not checkpointInfo["conformersComplete"]:
        drSplash.show_wriggle_splash()
        config = Wriggling_Doctor.conformer_generation_protocol(config=config)
        drYaml.write_config_to_yaml(config=config, out_dir=outputDir)

    ## run torsion scanning
    if not checkpointInfo["scanningComplete"]:
        drSplash.show_twist_splash()
        config = Twisted_Doctor.twist_protocol(config=config)
        drYaml.write_config_to_yaml(config=config, out_dir=outputDir)

    ## run charge calculations
    if not checkpointInfo["chargesComplete"]:
        drSplash.show_charge_splash()
        config = Charged_Doctor.charge_protocol(config=config)
        drYaml.write_config_to_yaml(config=config, out_dir=outputDir)

    ## run torsion parameter fitting
    if not checkpointInfo["torsionFittingComplete"]:
        drSplash.show_stitch_splash()
        if config["parameterFittingInfo"]["forceField"] == "AMBER": # type: ignore
            config = Stitching_Doctor.torsion_fitting_protocol_amber(config=config) # Corrected call
        elif config["parameterFittingInfo"]["forceField"] == "CHARMM": # type: ignore
            config = Stitching_Doctor.torsion_fitting_protocol_charmm(config=config) # Corrected call
        drYaml.write_config_to_yaml(config=config, out_dir=outputDir)

    ## run final creation
    if not checkpointInfo["finalCreationComplete"]:
        drSplash.show_creation_splash()
        config = drCreator.create_the_monster(config=config)
        drYaml.write_config_to_yaml(config=config, out_dir=outputDir)

    ## run reporting to make html
    if not checkpointInfo["reportingComplete"]:
        config = Reporting_Doctor.reporter_protocol(config=config) # type: ignore
        drYaml.write_config_to_yaml(config=config, out_dir=outputDir)

    ## show what we have created SPLASH
    drSplash.show_what_have_we_created(config=config)

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def handle_exceptions(e: Exception, pdb_name: str) -> Dict[str, Any]:
    tb_list = traceback.extract_tb(e.__traceback__) # Renamed tb to tb_list to avoid confusion with module
    fullTraceBack: List[str] = []
    functionName = 'Unknown'
    lineNumber = 'Unknown'
    lineOfCode = 'Unknown' # This might be hard to get accurately or might be empty
    scriptName = 'Unknown'

    if tb_list:
        # tb_list.reverse() # Reversing is not standard for traceback presentation. Usually newest call last.
        fullTraceBack = [f"{frame.filename}:{frame.lineno} in {frame.name}" for frame in tb_list]
        last_frame = tb_list[-1] # Get the most recent frame for specific details
        functionName = last_frame.name
        lineNumber = str(last_frame.lineno) # Ensure string
        lineOfCode = last_frame.line if last_frame.line else 'Not available' # Ensure not None
        scriptName = last_frame.filename
    
    errorType = type(e).__name__
    errorData: Dict[str, Any] = {
        "pdbName": pdb_name,
        "errorType": errorType,
        "errorMessage": str(e),
        "functionName": functionName,
        "lineNumber": lineNumber,
        "lineOfCode": lineOfCode,
        "scriptName": scriptName,
        "fullTraceBack": fullTraceBack
    }
    return errorData

# 🗲𗲲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        errorReport = handle_exceptions(e=e, pdb_name="drFrankenstein")
        drSplash.print_botched(error_report=errorReport)
