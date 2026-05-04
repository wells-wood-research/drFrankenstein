import os
from os import path as p
import traceback
import ruamel.yaml as ruamel
## drFRANKENSTEIN LIBRARIES ##
from Experiments.Protocol_1_Capping import Capping_Doctor
from Experiments.Protocol_4_Assembly import Assembly_Doctor
from Experiments.Protocol_2_Wriggling import Wriggling_Doctor
from Experiments.Protocol_5_Twisting import Twisted_Doctor
from Experiments.Protocol_3_Charging import Charged_Doctor
from Experiments.Protocol_6_Stitching import Stitching_Doctor
from Experiments.Protocol_7_Creation import drCreator
from Experiments.Protocol_8_Reporter import Reporting_Doctor

from OperatingTools import drYaml
from OperatingTools import drSplash
from OperatingTools import Timer
from OperatingTools import drLogger
from OperatingTools import set_config_defaults
from OperatingTools import validate_config
from OperatingTools import handle_CGenFF_dependancy
from OperatingTools import pdb_checker
## CLEAN CODE ##
class FilePath:
    pass
class DirectoryPath:
    pass

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def main():
    """
    Main protocol for drFrankenstein

    Args:
        None
    Returns:
        None
    
    """
    logger = None
    try:
        drSplash.show_mad_man()

        ## get config yaml file from argpass command-line-argument
        configYaml = drYaml.get_config_input_arg()
        ## load into dict, check for bad formatting
        config = drYaml.read_input_yaml(configYaml)

        config = set_config_defaults.apply_defaults_and_validate(config)

        ## check config for errors#
        config = validate_config.validate_config(config)

        ## unpack config to find outputDir, make directory
        outputDir = config["pathInfo"]["outputDir"]
        os.makedirs(outputDir,exist_ok=True)
        
        ## Initialize logging system
        log_dir = p.join(outputDir, "00_logs")
        logger = drLogger.ExperimentLogger(log_dir)
        drLogger.set_logger(logger)  # Set global logger instance
        logger.logger.info(f"drFrankenstein Parameterisation started")
        
        ## deal with checkpointing, lets us skip steps if already done
        config = drYaml.read_config_with_checkpoints(config, outputDir)
        config = drYaml.init_config_checkpoints(config, outputDir)
        ## check input PDB
        pdb_checker.check_pdb(config)

        ## initialise runtimeInfo
        config = drYaml.initialise_runtime_info(config)
        ## save config back to yaml
        drYaml.write_config_to_yaml(config, outputDir)

        ## add capping groups to input molecule
        ## TODO: cope with staples etc. that need 2 or more capping groups
        checkpointInfo = config["checkpointInfo"]
        if not checkpointInfo["cappingComplete"]:
            drSplash.show_capping_splash()
            config = Capping_Doctor.capping_protocol(config=config)
            drYaml.write_config_to_yaml(config, outputDir)
        
        ## run conformer generation protocol
        if not checkpointInfo["conformersComplete"]:
            drSplash.show_wriggle_splash()
            config = Wriggling_Doctor.conformer_generation_protocol(config=config)
            drYaml.write_config_to_yaml(config, outputDir)

        ## run charge calculations
        if not checkpointInfo["chargesComplete"]:
            drSplash.show_charge_splash()
            config = Charged_Doctor.charge_protocol(config=config)
            drYaml.write_config_to_yaml(config, outputDir)

        if not checkpointInfo["assemblyComplete"]:
            ## TODO: drSplash?
            config = Assembly_Doctor.parameter_assembly_protocol(config=config)
            drYaml.write_config_to_yaml(config, outputDir)


        ## run torsion scanning
        if not checkpointInfo["scanningComplete"]:
            drSplash.show_twist_splash()
            config = Twisted_Doctor.twist_protocol(config=config)
            drYaml.write_config_to_yaml(config, outputDir)


        ## run torsion parameter fitting
        if not checkpointInfo["torsionFittingComplete"]:
            drSplash.show_stitch_splash()
            config = Stitching_Doctor.torsion_fitting_protocol(config=config)
            drYaml.write_config_to_yaml(config, outputDir)

        ## run final creation
        if not checkpointInfo["finalCreationComplete"]:
            drSplash.show_creation_splash()
            config = drCreator.create_the_monster(config=config)
            drYaml.write_config_to_yaml(config, outputDir)

        ## run reporting to make html
        if not checkpointInfo["reportingComplete"]:
            config = Reporting_Doctor.reporter_protocol(config=config)
            drYaml.write_config_to_yaml(config, outputDir)

        ## show what we have created SPLASH
        drSplash.show_what_have_we_created(config)
        
        ## Log completion
        logger.log_completion_tag()
    
    except Exception as e:
        errorReport = handle_exceptions(e, "drFrankenstein")
        # Log the error if logger was initialized
        if logger:
            logger.log_error_report(errorReport)
        raise

# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def handle_exceptions(e, pdbName):
    tb = traceback.extract_tb(e.__traceback__)
    if tb:
        tb.reverse()
        fullTraceBack = [f"{frame.filename}:{frame.lineno} in {frame.name}" for frame in tb]
        lastFrame = tb[-1]
        functionName = lastFrame.name
        lineNumber = lastFrame.lineno
        lineOfCode = lastFrame.line
        scriptName = lastFrame.filename
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
        
        # Log the error to file if logger was initialized
        try:
            # Try to get logger from config (if it was created)
            # This might fail if error occurred before logger initialization
            import sys
            logger = None
            # Check if we can access the logger through the exception context
            # If main() was called, we would have created a logger
            # For now, create a new logger to ensure error is logged
            import tempfile
            temp_log_dir = tempfile.gettempdir()
            from OperatingTools import drLogger
            error_logger = drLogger.ExperimentLogger(temp_log_dir)
            error_logger.log_error_report(errorReport)
        except Exception as log_err:
            # If logging fails, that's okay - still print to console
            pass
        
        drSplash.print_botched(errorReport)
