import os
from os import path as p

## drFRANKENSTEIN LIBRARIES ##
import drInputs
import drCapper
import drTwist
import drCharge

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
    ## add capping groups to input molecule
    ## TODO: cope with staples etc. that need 2 or more capping groups
    cappedMolPdb = drCapper.capping_protocol(config)
    ## update config
    config["moleculeInfo"]["cappedPdb"] = cappedMolPdb
    ## run torsion scanning protocol
    drTwist.twist_protocol(config)
    ## run charge fitting
    drCharge.charge_protocol(config)

    


if __name__ == "__main__":
    main()