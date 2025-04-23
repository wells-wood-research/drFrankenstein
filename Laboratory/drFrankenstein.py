import os
from os import path as p
import traceback
import ruamel.yaml as ruamel
## drFRANKENSTEIN LIBRARIES ##
from Experiments.Protocol_1_Capping import Capping_Doctor
from Experiments.Protocol_2_Wriggling import Wriggling_Doctor
from Experiments.Protocol_3_Twisting import Twisted_Doctor
from Experiments.Protocol_4_Charging import Charged_Doctor
from Experiments.Protocol_5_Stitching import Stitching_Doctor
from Experiments.Protocol_6_Creation import drCreator

from OperatingTools import drYaml
from OperatingTools import drSplash
from OperatingTools import Timer 
from OperatingTools import validate_config
## CLEAN CODE ##
class FilePath:
    pass
class DirectoryPath:
    pass

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def main():
    ## get config yaml file from argpass command-line-argument
    configYaml = drYaml.get_config_input_arg()
    ## load into dict, check for bad formatting
    config = drYaml.read_input_yaml(configYaml)
    ## check config for errors
    config = validate_config.validate_config(config)
    ## unpack config to find outputDir, make directory
    outputDir = config["pathInfo"]["outputDir"]
    os.makedirs(outputDir,exist_ok=True)
    ## deal with checkpointing, lets us skip steps if already done
    config = drYaml.read_config_with_checkpoints(config, outputDir)
    config = drYaml.init_config_checkpoints(config, outputDir)

    ## initialise runtimeInfo
    config = drYaml.initialise_runtime_info(config)
    ## set up function timers
    config = Timer.sort_output_directory(config)
    ## save config back to yaml
    drYaml.write_config_to_yaml(config, outputDir)

    ## add capping groups to input molecule
    ## TODO: cope with staples etc. that need 2 or more capping groups
    checkpointInfo = config["checkpointInfo"]
    if not checkpointInfo["cappingComplete"]:
        print("Running Capping Protocol")
        config = Capping_Doctor.capping_protocol(config=config)
        drYaml.write_config_to_yaml(config, outputDir)
    ## run conformer generation protocol
    if not checkpointInfo["conformersComplete"]:
        drSplash.show_wriggle_splash()
        config = Wriggling_Doctor.conformer_generation_protocol(config=config)
        drYaml.write_config_to_yaml(config, outputDir)
    ## run torsion scanning
    if not checkpointInfo["scanningComplete"]:
        drSplash.show_twist_splash()
        config = Twisted_Doctor.twist_protocol(config=config)
        drYaml.write_config_to_yaml(config, outputDir)

    ## run charge calculations
    if not checkpointInfo["chargesComplete"]:
        drSplash.show_charge_splash()
        config = Charged_Doctor.charge_protocol(config=config)
        drYaml.write_config_to_yaml(config, outputDir)

    ## run torsion parameter fitting
    if not checkpointInfo["torsionFittingComplete"]:
        drSplash.show_stitch_splash()
        if config["parameterFittingInfo"]["forceField"] == "AMBER":
            config = Stitching_Doctor.torsion_fitting_protocol_AMBER(config=config)
        elif config["parameterFittingInfo"]["forceField"] == "CHARMM":
            config = Stitching_Doctor.torsion_fitting_protocol_CHARMM(config=config)

        drYaml.write_config_to_yaml(config, outputDir)

    if not checkpointInfo["finalCreationComplete"]:
        drSplash.show_creation_splash()
        drCreator.create_the_monster(config=config)
    
    
    drSplash.show_what_have_we_created(config["moleculeInfo"]["moleculeName"])

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def handle_exceptions(e, pdbName):
    tb = traceback.extract_tb(e.__traceback__)
    if tb:
        tb.reverse()
        fullTraceBack = [f"{frame.filename}:{frame.lineno} in {frame.name}" for frame in tb]
        last_frame = tb[-1]
        functionName = last_frame.name
        lineNumber = last_frame.lineno
        lineOfCode = last_frame.line
        scriptName = last_frame.filename
    else:
        functionName = 'Unknown'
        lineNumber = 'Unknown'
        lineOfCode = 'Unknown'
    
    errorType = type(e).__name__
    errorData: dict = {
        "pdbName": pdbName,
        "errorType": errorType,
        "errorMessage": str(e),
        "functionName": functionName,
        "lineNumber": lineNumber,
        "lineOfCode": lineOfCode,
        "scriptName": scriptName,
        "fullTraceBack": fullTraceBack
    }
    return errorData

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        errorReport = handle_exceptions(e, "drFrankenstein")
        drSplash.print_botched(errorReport)
