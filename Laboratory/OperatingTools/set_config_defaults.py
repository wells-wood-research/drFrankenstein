import os
import io
import multiprocessing as mp
from ruamel.yaml import YAML

from OperatingTools import drSplash

def _has_fatal_errors(d):
    """Recursively checks a dictionary for fatal error messages."""
    if isinstance(d, dict):
        return any(_has_fatal_errors(v) for v in d.values())
    elif isinstance(d, str) and "default" not in d.lower():
        return True
    return False


def apply_defaults_and_validate(config):
    """Reads the yaml file, applies defaults, and validates all parameters."""
    errors = {}

    # --- Section: moleculeInfo ---
    errors['moleculeInfo'] = {}
    moleculeInfo = config.get('moleculeInfo', {})
    for key, expected_type in [('charge', int), ('multiplicity', int), ('moleculeName', str)]:
        if key not in moleculeInfo:
            errors['moleculeInfo'][key] = "Missing required key."
        elif not isinstance(moleculeInfo[key], expected_type):
            errors['moleculeInfo'][key] = f"Must be type {expected_type}, not {type(moleculeInfo[key])}."
        else:
            errors['moleculeInfo'][key] = None
    for key, expected_type in [('chargeGroups', dict), ('backboneAliases', dict)]:
        if key not in moleculeInfo:
             errors['moleculeInfo'][key] = f"Default Used: {None}"
        elif not isinstance(moleculeInfo[key], (expected_type, type(None))):
             errors['moleculeInfo'][key] = f"Must be type {expected_type} or None."
        else:
            errors['moleculeInfo'][key] = None
    moleculeInfo.setdefault('chargeGroups', None)
    moleculeInfo.setdefault('backboneAliases', None)

    # --- Section: pathInfo ---
    errors['pathInfo'] = {}
    pathInfo = config.setdefault('pathInfo', {})
    if 'inputDir' not in pathInfo:
        errors['pathInfo']['inputDir'] = f"Default Used: {os.getcwd()}"
    pathInfo.setdefault('inputDir', os.getcwd())
    
    if 'outputDir' not in pathInfo:
        if moleculeInfo.get('moleculeName') and errors['moleculeInfo']['moleculeName'] is None:
            default_path = os.path.join(os.getcwd(), f"{moleculeInfo['moleculeName']}_FrankenParams")
            pathInfo['outputDir'] = default_path
            errors['pathInfo']['outputDir'] = f"Default Used: {default_path}"
        else:
            errors['pathInfo']['outputDir'] = "Cannot set default, 'moleculeName' is missing or invalid."

    if 'amberHome' not in pathInfo:
        errors['pathInfo']['amberHome'] = f"Default Used: {os.environ.get('AMBERHOME')}"
    pathInfo.setdefault('amberHome', os.environ.get('AMBERHOME'))
    
    if 'cgenffExe' not in pathInfo:
        errors['pathInfo']['cgenffExe'] = f"Default Used: {None}"
    pathInfo.setdefault('cgenffExe', None)

    for key, expected_type in [('inputDir', str), ('outputDir', str), ('amberHome', str), ('cgenffExe', str)]:
        if key not in errors['pathInfo']:
             errors['pathInfo'][key] = None
        if key in pathInfo and not isinstance(pathInfo[key], (expected_type, type(None))):
             errors['pathInfo'][key] = f"Must be type {expected_type} or None."

    for key, is_dir in [('multiWfnDir', True), ('orcaExe', False)]:
        path_val = pathInfo.get(key)
        if not path_val:
            errors['pathInfo'][key] = "Missing required key."
        elif not isinstance(path_val, str):
            errors['pathInfo'][key] = f"Must be type str, not {type(path_val)}."
        elif is_dir and not os.path.isdir(path_val):
            errors['pathInfo'][key] = f"Directory not found: {path_val}"
        elif not is_dir and not os.path.isfile(path_val):
            errors['pathInfo'][key] = f"File not found: {path_val}"
        else:
            if key not in errors['pathInfo']:
                 errors['pathInfo'][key] = None

    # --- Section: torsionScanInfo ---
    errors['torsionScanInfo'] = {}
    torsionScanInfo = config.get('torsionScanInfo', {})
    if not torsionScanInfo:
        errors['torsionScanInfo'] = "Missing required key."
    runScansOn = torsionScanInfo.setdefault('runScansOn', {})
    errors['torsionScanInfo']['runScansOn'] = {}
    for key in ['phiPsi', 'polarProtons', 'nonPolarProtons', 'amides', 'nonAromaticRings']:
        if key not in runScansOn:
            default_val = key in ['phiPsi', 'polarProtons']
            runScansOn[key] = default_val
            errors['torsionScanInfo']['runScansOn'][key] = f"Default Used: {default_val}"
        elif not isinstance(runScansOn[key], bool):
            errors['torsionScanInfo']['runScansOn'][key] = f"Must be type bool."
        else:
            errors['torsionScanInfo']['runScansOn'][key] = None

    if 'scanMethod' not in torsionScanInfo:
        errors['torsionScanInfo']['scanMethod'] = "Required: use a SIMPLE INPUT line from the ORCA input library"
    # --- Simplified setdefault and error logging for remaining sections ---
    def set_default(section, key, default_value, error_dict):
        if key not in section:
            error_dict[key] = f"Default Used: {default_value}"
        section.setdefault(key, default_value)
        
    set_default(torsionScanInfo, 'nConformers', -1, errors['torsionScanInfo'])
    set_default(torsionScanInfo, 'scanSolvationMethod', None, errors['torsionScanInfo'])
    set_default(torsionScanInfo, 'singlePointMethod', None, errors['torsionScanInfo'])
    set_default(torsionScanInfo, 'singlePointSolvationMethod', None, errors['torsionScanInfo'])

    # --- Section: chargeFittingInfo ---
    errors['chargeFittingInfo'] = {}
    chargeFittingInfo = config.setdefault('chargeFittingInfo', {})
    for key in ['chargeFittingProtocol', 'optMethod', 'singlePointMethod']:
        if key not in chargeFittingInfo:
            errors['chargeFittingInfo'][key] = "Missing required key."
        else:
            errors['chargeFittingInfo'][key] = None
    
    set_default(chargeFittingInfo, 'nConformers', -1, errors['chargeFittingInfo'])
    set_default(chargeFittingInfo, 'nCoresPerCalculation', 1, errors['chargeFittingInfo'])
    if chargeFittingInfo.get('chargeFittingProtocol') == 'SOLVATOR':
        set_default(chargeFittingInfo, 'waterDensity', 10, errors['chargeFittingInfo'])

    # --- Section: parameterFittingInfo ---
    errors['parameterFittingInfo'] = {}
    parameterFittingInfo = config.setdefault('parameterFittingInfo', {})
    if 'forceField' not in parameterFittingInfo:
        errors['parameterFittingInfo']['forceField'] = "Missing required key."
    else:
        errors['parameterFittingInfo']['forceField'] = None
    
    set_default(parameterFittingInfo, 'maxCosineFunctions', 3, errors['parameterFittingInfo'])
    set_default(parameterFittingInfo, 'nShuffles', 50, errors['parameterFittingInfo'])
    set_default(parameterFittingInfo, 'l2DampingFactor', 0.1, errors['parameterFittingInfo'])
    set_default(parameterFittingInfo, 'sagvolSmoothing', True, errors['parameterFittingInfo'])

    # --- Section: miscInfo ---
    errors['miscInfo'] = {}
    miscInfo = config.setdefault('miscInfo', {})
    set_default(miscInfo, 'availableCpus', mp.cpu_count(), errors['miscInfo'])
    set_default(miscInfo, 'cleanUpLevel', 1, errors['miscInfo'])

    # --- Final Error Check ---
    if _has_fatal_errors(errors):
        drSplash.print_config_error(errors)

    return config

# --- Main Execution ---
if __name__ == "__main__":
    raise NotImplementedError