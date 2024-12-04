import os
from os import path as p

## drFRANKENSTEIN LIBRARIES ##
import drInputs
import drCapper
import drTwist

## CLEAN CODE ##
class FilePath:
    pass
class DirectoryPath:
    pass


def main():
    configYaml = drInputs.get_config_input_arg()
    config = drInputs.read_input_yaml(configYaml)
    outputDir = config["pathInfo"]["outputDir"]
    os.makedirs(outputDir,exist_ok=True)

    cappedMolPdb = drCapper.capping_protocol(config)
    config["moleculeInfo"]["moleculePdb"] = cappedMolPdb

    drTwist.twist_protocol(config)
    


if __name__ == "__main__":
    main()