import os
from os import path as p
import argpass
import ruamel.yaml as ruamel
from ruamel.yaml.comments import CommentedMap
import numpy as np

## CLEAN CODE ##
class FilePath:
    pass
class DirectoryPath:
    pass


def initialise_runtime_info(config: dict) -> dict:
    """Ensure runtimeInfo exists and attach its banner comment on first creation."""
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
def read_config_with_checkpoints(config: dict, outDir: DirectoryPath) -> dict:
    """Reload the output YAML when it already exists so checkpoints are preserved."""
    drFrankensteinYaml = p.join(outDir, "drFrankenstein.yaml")
    ruamelParser = ruamel.YAML()
    ruamelParser.preserve_quotes = True  # Ensure quotes are preserved

    if p.isfile(drFrankensteinYaml):
        with open(drFrankensteinYaml, "r") as f:
            config = ruamelParser.load(f)
    return config
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def init_config_checkpoints(config: dict, outDir: DirectoryPath) -> dict:
    """Create the checkpoint structure when no output YAML has been written yet."""
    drFrankensteinYaml = p.join(outDir, "drFrankenstein.yaml")
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
            "conformersComplete": False,
            "chargesComplete": False,
            "assemblyComplete": False,
            "scanningComplete": False,
            "torsionFittingComplete": False,
            "finalCreationComplete": False,
            "reportingComplete": False
        }
        config.yaml_set_comment_before_after_key("checkpointInfo", before=checkpointBanner)

        return config
    else:
        return config
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def write_config_to_yaml(config: dict, outDir: DirectoryPath) -> None:
    """Write the current config to YAML after sanitising non-serialisable values."""
    drFrankensteinYaml = p.join(outDir, "drFrankenstein.yaml")
    ruamelParser = ruamel.YAML()
    ruamelParser.indent(mapping=2, sequence=4, offset=2)
    ruamelParser.default_flow_style = None
    ruamelParser.allow_unicode = True
    ruamelParser.preserve_quotes = True

    # monkey patch:
    ruamelParser.representer.ignore_aliases = lambda x: True

    # Helper to sanitize values (convert numpy types, arrays, bytes, and other non-serializable objects)
    def _sanitize_value(v: object) -> object:
        """Recursively convert values into YAML-safe primitives."""
        # Avoid importing numpy repeatedly
        try:
            import numpy as _np
        except Exception:
            _np = None

        # Handle CommentedMap/dict in-place
        if isinstance(v, (dict, CommentedMap)):
            for k in list(v.keys()):
                v[k] = _sanitize_value(v[k])
            return v

        # Lists/tuples/sets -> lists
        if isinstance(v, (list, tuple, set)):
            return [_sanitize_value(i) for i in v]

        # Numpy arrays -> lists
        if _np is not None and isinstance(v, _np.ndarray):
            return _sanitize_value(v.tolist())

        # Numpy scalar types
        if _np is not None and isinstance(v, (_np.integer,)):
            return int(v)
        if _np is not None and isinstance(v, (_np.floating,)):
            return float(v)
        if _np is not None and isinstance(v, (_np.bool_ ,)):
            return bool(v)

        # Bytes -> decode
        if isinstance(v, (bytes, bytearray)):
            try:
                return v.decode("utf-8")
            except Exception:
                return v.decode("utf-8", errors="replace")

        # Objects exposing tolist/item
        try:
            if hasattr(v, "tolist") and not isinstance(v, (str, bytes, bytearray)):
                return _sanitize_value(v.tolist())
            if hasattr(v, "item") and callable(getattr(v, "item")):
                return _sanitize_value(v.item())
        except Exception:
            pass

        # Functions, classes and other callables -> repr
        if callable(v) or isinstance(v, type):
            return repr(v)

        # Default: return as-is
        return v

    # Create a sanitized copy in-place where possible to preserve CommentedMap/comment objects
    sanitized_config = _sanitize_value(config)

    with open(drFrankensteinYaml, "w") as f:
        try:
            ruamelParser.dump(sanitized_config, f)
        except Exception as e:
            # Provide a helpful error message including a sample repr of the problematic object
            raise RuntimeError(f"Failed to write YAML config to {drFrankensteinYaml}: {e}")
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

def read_input_yaml(configFile: FilePath) -> dict:
    """
    Reads YAML file into a dict

    Args:
    - configFile (str): Path to the YAML configuration file.

    Returns:
    - config (dict): Parsed YAML content as a dictionary.
    """
    yellow = "\033[33m"
    reset = "\033[0m"
    teal = "\033[38;5;37m"
    try:
        ruamelParser = ruamel.YAML()
        ruamelParser.preserve_quotes = True  # Ensure quotes are preserved
        with open(configFile, "r") as yamlFile:
            config: dict = ruamelParser.load(yamlFile)
            return config
        

    except FileNotFoundError:
        print(f"-->{' '*4}Config file {configFile} not found.")
        exit(1)
    except ruamel.YAMLError as exc:
        print(f"-->{' '*4}{yellow}Error while parsing YAML file:{reset}")
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
