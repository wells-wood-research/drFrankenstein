
import argpass
import ruamel.yaml as ruamel
## CLEAN CODE ##
class FilePath:
    pass
class DirectoryPath:
    pass


# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def read_config_with_checkpoints(config, outDir):
    drFrankensteinYaml = p.join(outDir, "drFrankenstein.yaml")
    ruamelParser = ruamel.YAML()

    if p.isfile(drFrankensteinYaml):
        with open(drFrankensteinYaml, "r") as f:
            config = ruamelParser.load(f)
    return config
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
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
# 🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def write_config_to_yaml(config, outDir):
    drFrankensteinYaml = p.join(outDir, "drFrankenstein.yaml")
    ruamelParser = ruamel.YAML()
    ruamelParser.indent(mapping=2, sequence=4, offset=2)
    ruamelParser.top_comment = True
    with open(drFrankensteinYaml, "w") as f:
        ruamelParser.dump(config, f)

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