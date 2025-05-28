import textwrap
from textwrap import dedent
import ruamel.yaml as ruamel
from ruamel.yaml.comments import CommentedMap
class FilePath(str):
    pass


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
############################################################################
def methods_writer_protocol(config: dict) -> dict:
    """
    Runs the methods writer protocol

    Args:
        config (dict): Dictionary containing all information needed

    Returns:
        None
    """
    conformerMethods : str = write_conformer_method(config)


    scanningMethods: str = write_torsion_scanning_method(config)


    chargeMethods: str = write_charge_calculation_method(config)

    fittingMethods = write_fitting_methods(config)

    methods = {
        "conformerMethods": conformerMethods,
        "scanningMethods": scanningMethods,
        "chargeMethods": chargeMethods,
        "fittingMethods": fittingMethods
    }

    return methods
############################################################################


def write_conformer_method(config: dict) -> str:
    """
    Writes the conformer method to be used in the methods writer protocol

    Args:
        config (dict): Dictionary containing all information needed

    Returns:
        conformerMethod (str): The conformer method to be used in the methods writer protocol
    """
    ## unpack config ##
    moleculeName: str = config["moleculeInfo"]["moleculeName"]
    nTermini: list = config["moleculeInfo"]["nTermini"]
    nConformersGenerated = config["runtimeInfo"]["madeByConformers"]["nConformersGenerated"]
    isCapped = False
    if len(nTermini) != 0:
        isCapped = True

    
    conformerMethod = dedent(f"""
A library of {nConformersGenerated} low-energy conformers of {"capped" if isCapped else ""} {moleculeName} was generated using ORCA's GOAT protocol. \
These calculations were performed using the GFN2-XTB semi-empirical method. \
The conformers generated in this step were used as inputs for the subsequent torsion-scanning and charge calculation protocols.
                            """)
    
    return conformerMethod

############################################################################
def write_torsion_scanning_method(config: dict) -> str:
    """
    Writes the method section for the torsion-scanning protocol to be used in the methods writer protocol

    Args:
        config (dict): Dictionary containing all information needed

    Returns:
        torsionMethod (str): The torsion-scanning method to be used in the methods writer protocol
    """
    ## unpack config ##
    moleculeName: str = config["moleculeInfo"]["moleculeName"]
    nConformersGenerated = config["runtimeInfo"]["madeByConformers"]["nConformersGenerated"]
    nConformers = config["torsionScanInfo"]["nConformers"]
    preserveBackboneTorsions = config["torsionScanInfo"]["preserveBackboneTorsions"]
    scanMethod = config["torsionScanInfo"]["scanMethod"]
    scanSolvationMethod = config["torsionScanInfo"]["scanSolvationMethod"]
    singlePointMethod = config["torsionScanInfo"]["singlePointMethod"]
    singlePointSolvationMethod = config["torsionScanInfo"]["singlePointSolvationMethod"]

    if nConformers == -1:
        nConformersUsed = nConformersGenerated
    elif nConformers > nConformersGenerated:
        nConformersUsed = nConformersGenerated
    else:
        nConformersUsed = nConformers

    torsionMethod = dedent(f"""
Torsion scans were performed on all rotatable bonds in {moleculeName}{", except those represented by backbone Phi and Psi angles" if preserveBackboneTorsions else ""}. \
For each rotatable dihedral full 360 degree relaxed scan was performed in the forwards and backwards directions. Each scan was conducted using 10 degree increments. \
Scans were performed at the {scanMethod} level{" with the " + scanSolvationMethod + " solvation method" if scanSolvationMethod is not None else ""}. \
                            """)
    
    if singlePointMethod:
        singlePointText = dedent(f"""
The energy profiles of each torsion scan were refined using single-point calculations at the {singlePointMethod} level \
{"with the " + singlePointSolvationMethod + " solvation method" if singlePointSolvationMethod is not None else ""}. \
                                 """)
        torsionMethod = torsionMethod + singlePointText
    
    averagingText = dedent(f"""
The torsion scanning process was performed {nConformersUsed} times, each with a different conformer geometry as a starting point. \
After discarding scans with unphysically high energy barriers, the final energy profile of each torsion was calculated by averaging each scan. \
                            """)

    torsionMethod = torsionMethod + averagingText
    
    return torsionMethod
############################################################################
def write_charge_calculation_method(config: dict) -> str:
    """
    Writes the method section for the charge calculation protocol 
    Args:
        config (dict): Dictionary containing all information needed

    Returns:
        conformerMethod (str): The charge calculation  method to be used in the methods writer protocol
    """
    ## unpack config ##
    moleculeName: str = config["moleculeInfo"]["moleculeName"]
    nConformersGenerated = config["runtimeInfo"]["madeByConformers"]["nConformersGenerated"]
    nConformers = config["chargeFittingInfo"]["nConformers"]
    nWaters = config["runtimeInfo"]["madeByCharges"]["nWaters"]
    chargeFittingProtocol = config["chargeFittingInfo"]["chargeFittingProtocol"]

    optMethod  = config["chargeFittingInfo"]["optMethod"]
    optSolvationMethod = config["chargeFittingInfo"]["optSolvationMethod"]
    singlePointMethod = config["chargeFittingInfo"]["singlePointMethod"]
    singlePointSolvationMethod = config["chargeFittingInfo"]["singlePointSolvationMethod"]

    if nConformers == -1:
        nConformersUsed = nConformersGenerated
    elif nConformers > nConformersGenerated:
        nConformersUsed = nConformersGenerated
    else:
        nConformersUsed = nConformers

    ## write charge method
    chargesMethod = dedent(f"""
The partial charges of {moleculeName} were calculated using the {chargeFittingProtocol} charge fitting protocol. \
                             """)
    
    if chargeFittingProtocol == "SOLVATOR":
        chargesMethod = chargesMethod + dedent(f"""
This protocol used ORCA's SOLVATOR method to place {nWaters} water molecules around the {moleculeName} molecule. \
The SOLVATOR protocol was performed at the GFN2-XTB level and using the stochastic solvation procedure. \
The stochastic procedure is extremely fast, but produces sub-optimal water interaction networks. \
To obtain better placements of the solvating water molecules, a QM/MM geometry optimisation was performed, \
where the {moleculeName} molecule was treated at the {optMethod} level and the solvating water molecules were treated at the \
MM level as TIP3P waters. \
Single-point energy calculations were then performed using QM/MM calculations were the {moleculeName} molecule was treated at the \
{singlePointMethod} level and the solvating water molecules were treated at the MM level as TIP3P waters. \
                                                   """)

    elif chargeFittingProtocol == "RESP":
        chargesMethod = chargesMethod + dedent(f"""
Initially a geometry optimisation was performed at the {optMethod} level \
{"with the " + optSolvationMethod + " solvation method" if optSolvationMethod is not None else ""}. \
Single-point energy calculations were then performed using single-point calculations at the {singlePointMethod} level \
{"with the " + singlePointSolvationMethod + " solvation method" if singlePointSolvationMethod is not None else ""}. \
                                                   """)

    elif chargeFittingProtocol == "RESP2":
        chargesMethod = chargesMethod + dedent(f"""
Initially a geometry optimisation was performed at the {optMethod} level \
{"with the " + optSolvationMethod + " solvation method" if optSolvationMethod is not None else ""}. 
A pair of single-point energy calculations were then performed using single-point calculations at the {singlePointMethod} level. \
One single-point was performed using the {singlePointSolvationMethod} solvation method and the other in the gas-phase. \
                                                   """)

    chargesMethod = chargesMethod + dedent(f"""
The above process was performed on {nConformersUsed} conformers of {moleculeName}, generated previously using ORCA's GOAT protocol. \
                                               """)

    chargesMethod = chargesMethod + dedent(f"""
MultiWFN was then used to assign partial charges to each atom in {moleculeName} \
using the restrained electrostatic potential (RESP) method. \
The relative energies of each conformer were used to weight their contributions to the final partial charges. \
                                                   """)
        
    if chargeFittingProtocol == "RESP2":
        chargesMethod = chargesMethod + dedent(f"""
This charge fitting process was performed once for the solvated systems and once for the gas-phase systems.
The final partial charges were then calculated using a 60:40weighted average between the solvated and gas-phase partial charges. \
                                               """)


    return chargesMethod


def write_fitting_methods(config: dict) -> str:
    ## unpack config
    maxCosineFunctions = config["parameterFittingInfo"]["maxCosineFunctions"]
    nShuffles = config["parameterFittingInfo"]["nShuffles"]
    fittingMethods = dedent(f"""
To obtain MM force field parameters that accurately represent the energies obtained in the previous torsion-scanning step, \
an iterative fitting procedure was performed. \
For each rotatable dihedral the QM(torsion) energy is calculated as:\n\
QM(torsion) = QM(total) - MM(total) + MM(torsion)\n\
QM(total) is the scan energies obtained in the previous torsion-scanning step. \
MM(total) is obtained by performing single-point energy calculations at the MM level on the geometries from the torsion-scanning step, \
these calculations were performed in OpenMM. \
MM(torsion) is obtained directly as a sum of cosine components described in the current parameter file for the torsion of interest. \
Once the QM(torsion) energies have been calculated, a maximum of {maxCosineFunctions} cosine functions were fit to the energy profile using the Fourier transform method. 
The amplitudes, phases and periodicities of the cosine functions were then used to update the parameters for the torsion of interest. \
As this iterative fitting process proceeds, the evaluation of QM(torsion) for each torsion is calculated using the updated parameters. \
This process was repeated {nShuffles} times, each time the order of the torsions was shuffled to remove any artifacts from the arbitrary order of the torsions. \
                            """)

    return fittingMethods


if __name__ == "__main__":
    configFile = "/home/esp/scriptDevelopment/drFrankenstein/04_GLY_CHARMM/drFrankenstein.yaml"

    config = read_input_yaml(configFile)

    methods_writer_protocol(config)