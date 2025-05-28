import os
from os import path as p
import argpass
import ruamel.yaml as ruamel
from ruamel.yaml.comments import CommentedMap

## CLEAN CODE ##
class FilePath:
    pass
class DirectoryPath:
    pass

from typing import Union, Dict, Any # For type hints
from ruamel.yaml.comments import CommentedMap # For specific dict type

def initialize_runtime_info(config: CommentedMap) -> CommentedMap: # Renamed function
    """
    Initializes the 'runtimeInfo' section in the configuration if it doesn't exist,
    and adds a descriptive banner comment before it.

    Args:
        config: The configuration object (ruamel.yaml.comments.CommentedMap).

    Returns:
        The modified configuration object.
    """
    runtimeBanner = '\n'.join([
        "##########################################################################################################################",                                     
        "                                    _     _                      ___            __                                       #",         
        "             _ __   _   _   _ __   | |_  (_)  _ __ ___     ___  |_ _|  _ __    / _|   ___                                #",     
        "            | '__| | | | | | '_ \  | __| | | | '_ ` _ \   / _ \  | |  | '_ \  | |_   / _ \                               #",     
        "            | |    | |_| | | | | | | |_  | | | | | | | | |  __/  | |  | | | | |  _| | (_) |                              #",     
        "            |_|     \__,_| |_| |_|  \__| |_| |_| |_| |_|  \___| |___| |_| |_| |_|    \___/                               #",     
        "                                                                                                                         #",     
        "##########################################################################################################################"                                     
    ])


    runtimeInfo = config.get("runtimeInfo", None)
    if runtimeInfo is None:
        config["runtimeInfo"] = {}
        # Use the banner variable for the comment
        config.yaml_set_comment_before_after_key("runtimeInfo", before=runtimeBanner)

    return config



# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def read_config_with_checkpoints(config: CommentedMap, out_dir: DirectoryPath) -> CommentedMap:
    """
    Reads the drFrankenstein.yaml file from the output directory if it exists,
    effectively loading checkpointed configuration.

    Args:
        config: The current configuration object.
        out_dir: The output directory where drFrankenstein.yaml might be located.

    Returns:
        The loaded configuration if the file exists, otherwise the original config.
    """
    drFrankensteinYaml: FilePath = p.join(out_dir, "drFrankenstein.yaml") # type: ignore
    ruamelParser = ruamel.YAML()
    ruamelParser.preserve_quotes = True  # Ensure quotes are preserved

    if p.isfile(drFrankensteinYaml):
        with open(drFrankensteinYaml, "r") as f:
            loaded_config = ruamelParser.load(f)
            if isinstance(loaded_config, dict): # Ensure loaded content is a dict/CommentedMap
                 return loaded_config # type: ignore
            # Optionally handle cases where the file is empty or not valid YAML
            # For now, returning original config if load fails to produce a dict
            return config # type: ignore
    return config
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def init_config_checkpoints(config: CommentedMap, out_dir: DirectoryPath) -> CommentedMap:
    """
    Initializes the 'checkpointInfo' section in the configuration if a
    drFrankenstein.yaml file doesn't already exist in the output directory.
    Adds a descriptive banner comment before it.

    Args:
        config: The configuration object.
        out_dir: The output directory.

    Returns:
        The modified configuration object.
    """
    drFrankensteinYaml: FilePath = p.join(out_dir, "drFrankenstein.yaml") # type: ignore
    checkpointBanner = '\n'.join([
        "##########################################################################################################################",                                     
        "                                                                                                                         #",
        "               _                     _                      _           _     ___            __                          #",
        "         ___  | |__     ___    ___  | | __  _ __     ___   (_)  _ __   | |_  |_ _|  _ __    / _|   ___                   #",
        "        / __| | '_ \   / _ \  / __| | |/ / | '_ \   / _ \  | | | '_ \  | __|  | |  | '_ \  | |_   / _ \                  #",
        "       | (__  | | | | |  __/ | (__  |   <  | |_) | | (_) | | | | | | | | |_   | |  | | | | |  _| | (_) |                 #",
        "        \___| |_| |_|  \___|  \___| |_|\_\ | .__/   \___/  |_| |_| |_|  \__| |___| |_| |_| |_|    \___/                  #",
        "                                           |_|                                                                           #",
        "                                                                                                                         #",
        "##########################################################################################################################"                                      
    ])                             


    if not p.exists(drFrankensteinYaml):
        config["checkpointInfo"] = {
            "cappingComplete": False,
            "assemblyComplete": False,
            "conformersComplete": False,
            "scanningComplete": False,
            "chargesComplete": False,
            "torsionFittingComplete": False,
            "finalCreationComplete": False,
            "reportingComplete": False
        }
        config.yaml_set_comment_before_after_key("checkpointInfo", before=checkpointBanner)

        return config # type: ignore
    else:
        return config # type: ignore
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def write_config_to_yaml(config: CommentedMap, out_dir: DirectoryPath) -> None:
    """
    Writes the configuration object to drFrankenstein.yaml in the output directory.

    Args:
        config: The configuration object to write.
        out_dir: The output directory.
    """
    drFrankensteinYaml: FilePath = p.join(out_dir, "drFrankenstein.yaml") # type: ignore
    ruamelParser = ruamel.YAML()
    ruamelParser.indent(mapping=2, sequence=4, offset=2)
    ruamelParser.default_flow_style = None
    ruamelParser.allow_unicode = True
    ruamelParser.preserve_quotes = True

    # monkey patch:
    ruamelParser.representer.ignore_aliases = lambda x: True

    with open(drFrankensteinYaml, "w") as f:
        ruamelParser.dump(config, f)
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

def get_config_input_arg() -> FilePath:
    """
    Sets up argpass to read the config.yaml file from command line
    Reads a YAML file using the "--config" flag with argpass

    Returns:
    - configFile (FilePath)
    """
    # create an argpass parser, read config file,
    parser = argpass.ArgumentParser()
    parser.add_argument(f"--config")
    args = parser.parse_args()

    configFile: FilePath = args.config

    return configFile
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

def read_input_yaml(config_file: FilePath) -> CommentedMap:
    """
    Reads YAML file into a dict

    Args:
    - configFile (str): Path to the YAML configuration file.

    Returns:
    # config (CommentedMap): Parsed YAML content as a CommentedMap.
    """
    yellow = "\033[33m"
    reset = "\033[0m"
    teal = "\033[38;5;37m"
    try:
        ruamelParser = ruamel.YAML()
        ruamelParser.preserve_quotes = True  # Ensure quotes are preserved
        with open(config_file, "r") as yamlFile: # type: ignore
            config: CommentedMap = ruamelParser.load(yamlFile)
            if not isinstance(config, dict): # Check if loading returned a dict-like object
                raise ruamel.YAMLError(f"YAML content in {config_file} did not parse as a dictionary.")
            return config
        

    except FileNotFoundError:
        print(f"-->{' '*4}Config file {config_file} not found.")
        exit(1) # Consider raising FileNotFoundError for library usage
    except ruamel.YAMLError as exc:
        print(f"-->{' '*4}{yellow}Error while parsing YAML file: {config_file}{reset}")
        if hasattr(exc, 'problem_mark'):
            mark = exc.problem_mark
            print(f"{' '*7}Problem found at line {mark.line + 1}, column {mark.column + 1}:")
            if exc.context:
                print(f"{' '*7}{exc.problem} {exc.context}")
            else:
                print(f"{' '*7}{exc.problem}")
            print(f"{' '*7}Please correct the data and retry.")
        else:
            print(f"{' '*7}Something went wrong while parsing the YAML file.")
            print(f"{' '*7}Please check the file for syntax errors.")
        print(f"\n{teal}TIP:{reset} Large language models (LLMs) like GPT-4 can be helpful for debugging YAML files.")
        print(f"{' '*5}If you get stuck with the formatting, ask a LLM for help!")
        exit(1)