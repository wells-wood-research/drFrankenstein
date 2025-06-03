import os
from os import path as p
import traceback
import ruamel.yaml as ruamel
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
from OperatingTools import Timer
from OperatingTools import validate_config
from OperatingTools import handle_CGenFF_dependancy

class FilePath:
    pass

class DirectoryPath:
    pass

def main():
    """
    Main protocol for drFrankenstein

    Args:
        None
    Returns:
        None
    
    """
    configYaml = drYaml.get_config_input_arg()
    config = drYaml.read_input_yaml(configYaml)
    outputDir = config['pathInfo']['outputDir']
    os.makedirs(outputDir, exist_ok=True)
    config = drYaml.read_config_with_checkpoints(config, outputDir)
    config = drYaml.init_config_checkpoints(config, outputDir)
    config = drYaml.initialise_runtime_info(config)
    drYaml.write_config_to_yaml(config, outputDir)
    checkpointInfo = config['checkpointInfo']
    if not checkpointInfo['cappingComplete']:
        print('Running Capping Protocol')
        config = Capping_Doctor.capping_protocol(config=config)
        drYaml.write_config_to_yaml(config, outputDir)
    if not checkpointInfo['assemblyComplete']:
        if config['parameterFittingInfo']['forceField'] == 'CHARMM':
            config = handle_CGenFF_dependancy.handle_CGenFF_dependancy(config)
            drYaml.write_config_to_yaml(config, outputDir)
            config = Assembly_Doctor.charmm_assembly_protocol(config=config)
            drYaml.write_config_to_yaml(config, outputDir)
        elif config['parameterFittingInfo']['forceField'] == 'AMBER':
            config = Assembly_Doctor.amber_assembly_protocol(config=config)
            drYaml.write_config_to_yaml(config, outputDir)
    if not checkpointInfo['conformersComplete']:
        drSplash.show_wriggle_splash()
        config = Wriggling_Doctor.conformer_generation_protocol(config=config)
        drYaml.write_config_to_yaml(config, outputDir)
    if not checkpointInfo['scanningComplete']:
        drSplash.show_twist_splash()
        config = Twisted_Doctor.twist_protocol(config=config)
        drYaml.write_config_to_yaml(config, outputDir)
    if not checkpointInfo['chargesComplete']:
        drSplash.show_charge_splash()
        config = Charged_Doctor.charge_protocol(config=config)
        drYaml.write_config_to_yaml(config, outputDir)
    if not checkpointInfo['torsionFittingComplete']:
        drSplash.show_stitch_splash()
        if config['parameterFittingInfo']['forceField'] == 'AMBER':
            config = Stitching_Doctor.torsion_fitting_protocol_AMBER(config=config)
        elif config['parameterFittingInfo']['forceField'] == 'CHARMM':
            config = Stitching_Doctor.torsion_fitting_protocol_CHARMM(config=config)
        drYaml.write_config_to_yaml(config, outputDir)
    if not checkpointInfo['finalCreationComplete']:
        drSplash.show_creation_splash()
        config = drCreator.create_the_monster(config=config)
        drYaml.write_config_to_yaml(config, outputDir)
    if not checkpointInfo['reportingComplete']:
        config = Reporting_Doctor.reporter_protocol(config=config)
        drYaml.write_config_to_yaml(config, outputDir)
    drSplash.show_what_have_we_created(config)

def handle_exceptions(e, pdbName):
    tb = traceback.extract_tb(e.__traceback__)
    if tb:
        tb.reverse()
        fullTraceBack = [f'{frame.filename}:{frame.lineno} in {frame.name}' for frame in tb]
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
    errorData: dict = {'pdbName': pdbName, 'errorType': errorType, 'errorMessage': str(e), 'functionName': functionName, 'lineNumber': lineNumber, 'lineOfCode': lineOfCode, 'scriptName': scriptName, 'fullTraceBack': fullTraceBack}
    return errorData
if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        errorReport = handle_exceptions(e, 'drFrankenstein')
        drSplash.print_botched(errorReport)