import os
from os import path as p
import traceback
import ruamel.yaml as ruamel
## drFRANKENSTEIN LIBRARIES ##
import drInputs
import drCapper.capping_protocol as capping_protocol
import drTwist.Doctor
import drCharge
from drHybrid import parameter_fitting_protocol
import drSplash
import drCreator
import drConformers
## CLEAN CODE ##
class FilePath:
    pass
class DirectoryPath:
    pass


def main():
    ## get config yaml file from argpass command-line-argument
    configYaml = drInputs.get_config_input_arg()
    ## load into dict, check for bad formatting
    ## TODO: write config_checker for bad args
    config = drInputs.read_input_yaml(configYaml)
    ## unpack config to find outputDir, make directory
    outputDir = config["pathInfo"]["outputDir"]
    os.makedirs(outputDir,exist_ok=True)

    config = read_config_with_checkpoints(config, outputDir)
    config = init_config_checkpoints(config, outputDir)

    write_config_to_yaml(config, outputDir)

    ## add capping groups to input molecule
    ## TODO: cope with staples etc. that need 2 or more capping groups
    checkpointInfo = config["checkpointInfo"]
    if not checkpointInfo["cappingComplete"]:
        print("Running Capping Protocol")
        config = capping_protocol.capping_protocol(config)
        write_config_to_yaml(config, outputDir)

    ## run conformer generation protocol
    if not checkpointInfo["conformersComplete"]:
        print("Running Conformer Generation Protocol")
        config = drConformers.conformer_generation_protocol(config)
        write_config_to_yaml(config, outputDir)

    if not checkpointInfo["scanningComplete"]:
        print("Running Twist Protocol")
        config = drTwist.Doctor.twist_protocol(config)
        write_config_to_yaml(config, outputDir)

    if not checkpointInfo["chargesComplete"]:
        print("Running Charge Protocol")
        config = drCharge.charge_protocol(config)
        write_config_to_yaml(config, outputDir)

    if not checkpointInfo["torsionFittingComplete"]:
        print("Running Parameter Fitting Protocol")
        config = parameter_fitting_protocol.torsion_fitting_protocol(config)
        write_config_to_yaml(config, outputDir)

    if not checkpointInfo["finalCreationComplete"]:
        drCreator.create_the_monster(config)
    
    print("WHAT HAVE WE DONE??")
    print("WHAT HAVE WE CREATED???")
    print("")
    print("Parameters for your non-natural amino acid")



def read_config_with_checkpoints(config, outDir):
    drFrankensteinYaml = p.join(outDir, "drFrankenstein.yaml")
    ruamelParser = ruamel.YAML()

    if p.isfile(drFrankensteinYaml):
        with open(drFrankensteinYaml, "r") as f:
            config = ruamelParser.load(f)
    return config

def init_config_checkpoints(config, outDir):
    drFrankensteinYaml = p.join(outDir, "drFrankenstein.yaml")

    if not p.exists(drFrankensteinYaml):
        config["checkpointInfo"] = {
            "cappingComplete": False,
            "conformersComplete": False,
            "scanningComplete": False,
            "chargesComplete": False,
            "torsionFittingComplete": False,
            "finalCreationComplete": False
        }
        return config
    else:
        return config


def write_config_to_yaml(config, outDir):
    drFrankensteinYaml = p.join(outDir, "drFrankenstein.yaml")
    ruamelParser = ruamel.YAML()
    with open(drFrankensteinYaml, "w") as f:
        ruamelParser.dump(config, f)

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


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        errorReport = handle_exceptions(e, "drFrankestein")
        drSplash.print_botched(errorReport)
