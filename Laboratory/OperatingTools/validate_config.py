import os
from typing import Dict, Any, List, Union, Optional, Tuple
from OperatingTools import drSplash


def _add_error(errors: Dict[str, str], keyPath: str, message: str):
    """Adds an error message to the errors dictionary."""
    errors[keyPath] = message

def _check_key_exists(data: Dict[str, Any], key: str, sectionPath: str, errors: Dict[str, str]) -> bool:
    """Checks if a key exists in the dictionary. Adds error if missing."""
    if key not in data:
        _add_error(errors, f"{sectionPath}.{key}", f"Missing required key '{key}'")
        return False
    return True

def _validate_type(value: Any, expectedType: Union[type, Tuple[type, ...]], keyPath: str, errors: Dict[str, str]) -> bool:
    """Checks if the value's type matches the expected type(s). Adds error if mismatched."""
    if not isinstance(value, expectedType):
        # Improved error message formatting for None/NoneType
        if isinstance(expectedType, type):
            expectedTypeStr = 'None' if expectedType is type(None) else expectedType.__name__
        else: # It's a tuple
            typeNames = []
            for t in expectedType:
                typeNames.append('None' if t is type(None) else t.__name__)
            expectedTypeStr = ' or '.join(typeNames)

        _add_error(errors, keyPath, f"Invalid type. Expected {expectedTypeStr}, but got {type(value).__name__}")
        return False
    return True

def _validate_allowed_values(value: Any, allowedValues: List[Any], keyPath: str, errors: Dict[str, str]) -> bool:
    """Checks if the value is within the list of allowed values. Adds error if not."""
    # Skip check if value is None, as None might be allowed by type but not explicitly in allowedValues list
    if value is None:
        return True # Assume type check already validated None if it was allowed
    if value not in allowedValues:
        # Filter out None from displayed allowed values if it's present, as it's handled by type check
        displayAllowed = [v for v in allowedValues if v is not None]
        _add_error(errors, keyPath, f"Invalid value '{value}'. Allowed values are: {displayAllowed}")
        return False
    return True

# --- Section Validators (snake_case names) ---

def _validate_path_info(sectionData: Optional[Dict[str, Any]], sectionName: str, errors: Dict[str, str]):
    """Validates the pathInfo section."""
    if sectionData is None:
        _add_error(errors, sectionName, "Missing required section 'pathInfo'")
        return

    if not isinstance(sectionData, dict):
         _add_error(errors, sectionName, f"Section '{sectionName}' should be a dictionary, but got {type(sectionData).__name__}")
         return

    # Keys here match the expected config structure, not changing them to camelCase
    expectedKeys = {
        "inputDir": str,
        "outputDir": str,
        "multiWfnDir": str,
        "orcaExe": str,
        "gaff2Dat": str,
    }

    for key, expectedType in expectedKeys.items():
        keyPath = f"{sectionName}.{key}"
        if _check_key_exists(sectionData, key, sectionName, errors):
            value = sectionData[key]
            _validate_type(value, expectedType, keyPath, errors)
            # Optional: Add os.path.exists checks here if needed
            # if key in ["inputDir", "multiWfnDir", "orcaExe", "gaff2Dat"]:
            #     if isinstance(value, str) and not os.path.exists(value):
            #          _add_error(errors, keyPath, f"Path does not exist: '{value}'")
            # elif key == "outputDir":
            #     parentDir = os.path.dirname(os.path.abspath(value))
            #     if isinstance(value, str) and not os.access(parentDir, os.W_OK):
            #          _add_error(errors, keyPath, f"Output directory's parent is not writable: '{parentDir}'")


def _validate_molecule_info(sectionData: Optional[Dict[str, Any]], sectionName: str, errors: Dict[str, str]):
    """Validates the moleculeInfo section."""
    if sectionData is None:
        _add_error(errors, sectionName, "Missing required section 'moleculeInfo'")
        return

    if not isinstance(sectionData, dict):
         _add_error(errors, sectionName, f"Section '{sectionName}' should be a dictionary, but got {type(sectionData).__name__}")
         return

    # Keys here match the expected config structure
    expectedKeys = {
        "charge": int,
        "multiplicity": int,
        "moleculePdb": str,
        "moleculeName": str,
        "chargeGroups": dict,
        "backboneAliases": dict,
    }

    for key, expectedType in expectedKeys.items():
        keyPath = f"{sectionName}.{key}"
        if _check_key_exists(sectionData, key, sectionName, errors):
            value = sectionData[key]
            if _validate_type(value, expectedType, keyPath, errors):
                # Specific validations for list/dict contents
                if key == "chargeGroups":
                    # value is the chargeGroups dictionary
                    for groupName, groupData in value.items():
                        groupKeyPath = f"{keyPath}.{groupName}"
                        if not isinstance(groupData, dict):
                            _add_error(errors, groupKeyPath, f"Charge group '{groupName}' must be a dictionary.")
                            continue # Skip further checks for this malformed group

                        # Check 'atoms' key within the group
                        if _check_key_exists(groupData, "atoms", groupKeyPath, errors):
                            atomsValue = groupData["atoms"]
                            atomsKeyPath = f"{groupKeyPath}.atoms"
                            if _validate_type(atomsValue, list, atomsKeyPath, errors):
                                if not all(isinstance(item, str) for item in atomsValue):
                                    _add_error(errors, atomsKeyPath, "All items in 'atoms' list must be strings.")

                        # Check 'charge' key within the group
                        if _check_key_exists(groupData, "charge", groupKeyPath, errors):
                            chargeValue = groupData["charge"]
                            chargeKeyPath = f"{groupKeyPath}.charge"
                            _validate_type(chargeValue, int, chargeKeyPath, errors)


def _validate_torsion_scan_info(sectionData: Optional[Dict[str, Any]], sectionName: str, errors: Dict[str, str]):
    """Validates the torsionScanInfo section."""
    if sectionData is None:
        _add_error(errors, sectionName, "Missing required section 'torsionScanInfo'")
        return

    if not isinstance(sectionData, dict):
         _add_error(errors, sectionName, f"Section '{sectionName}' should be a dictionary, but got {type(sectionData).__name__}")
         return

    # Keys here match the expected config structure
    # Use tuple (str, type(None)) for solvation methods to allow None
    expectedKeys = {
        "nConformers": int,
        "scanMethod": str,
        "scanSolvationMethod": (str, type(None)),  # Allow str or None
        "singlePointMethod": str,
        "singlePointSolvationMethod": (str, type(None)), # Allow str or None
        "skipPhiPSi": bool,
    }

    for key, expectedType in expectedKeys.items():
        keyPath = f"{sectionName}.{key}"
        if _check_key_exists(sectionData, key, sectionName, errors):
            value = sectionData[key]
            if _validate_type(value, expectedType, keyPath, errors):
                # Value range/specific value checks (only if type validation passed)
                if key == "nConformers" and isinstance(value, int) and value < -1:
                    _add_error(errors, keyPath, f"Value for '{key}' must be -1 or greater, but got {value}.")



def _validate_charge_fitting_info(sectionData: Optional[Dict[str, Any]], sectionName: str, errors: Dict[str, str]):
    """Validates the chargeFittingInfo section."""
    if sectionData is None:
        _add_error(errors, sectionName, "Missing required section 'chargeFittingInfo'")
        return

    if not isinstance(sectionData, dict):
         _add_error(errors, sectionName, f"Section '{sectionName}' should be a dictionary, but got {type(sectionData).__name__}")
         return

    # Keys here match the expected config structure
    # Use tuple (str, type(None)) for solvation methods to allow None
    expectedKeys = {
        "chargeFittingProtocol": str,
        "nConformers": int,
        "nCoresPerCalculation": int,
        "optMethod": str,
        "optSolvationMethod": (str, type(None)), # Allow str or None
        "singlePointMethod": str,
        "singlePointSolvationMethod": (str, type(None)), # Allow str or None
    }

    allowedProtocols = ["RESP", "RESP2", "SOLVATOR"]
    # allowedQmmmSolvation = ["TIP3P"] # Context specific check removed for generality

    for key, expectedType in expectedKeys.items():
        keyPath = f"{sectionName}.{key}"
        if _check_key_exists(sectionData, key, sectionName, errors):
            value = sectionData[key]
            if _validate_type(value, expectedType, keyPath, errors):
                 # Value range/specific value checks (only if type validation passed)
                if key == "chargeFittingProtocol" and isinstance(value, str):
                    _validate_allowed_values(value, allowedProtocols, keyPath, errors)
                elif key == "nConformers" and isinstance(value, int) and value < -1:
                    _add_error(errors, keyPath, f"Value for '{key}' must be -1 or greater, but got {value}.")
                elif key == "nCoresPerCalculation" and isinstance(value, int) and value <= 0:
                     _add_error(errors, keyPath, f"Value for '{key}' must be a positive integer, but got {value}.")


def _validate_parameter_fitting_info(sectionData: Optional[Dict[str, Any]], sectionName: str, errors: Dict[str, str]):
    """Validates the parameterFittingInfo section."""
    if sectionData is None:
        _add_error(errors, sectionName, "Missing required section 'parameterFittingInfo'")
        return

    if not isinstance(sectionData, dict):
         _add_error(errors, sectionName, f"Section '{sectionName}' should be a dictionary, but got {type(sectionData).__name__}")
         return

    # Keys here match the expected config structure
    expectedKeys = {
        "forceField": str,
        "maxCosineFunctions": int,
        "maxShuffles": int,
        "minShuffles": int,
        "l2DampingFactor": float | None,
        "sagvolSmoothing": bool
    }
    allowedForceFields = ["CHARMM", "AMBER"]

    for key, expectedType in expectedKeys.items():
        keyPath = f"{sectionName}.{key}"
        if _check_key_exists(sectionData, key, sectionName, errors):
            value = sectionData[key]
            if _validate_type(value, expectedType, keyPath, errors):
                # Value range/specific value checks (only if type validation passed)
                if key == "forceField" and isinstance(value, str):
                    _validate_allowed_values(value, allowedForceFields, keyPath, errors)
                elif key == "maxCosineFunctions" and isinstance(value, int) and value <= 0:
                    _add_error(errors, keyPath, f"Value for '{key}' must be a positive integer, but got {value}.")
                elif key == "maxShuffles" and isinstance(value, int) and value <= 0:
                    _add_error(errors, keyPath, f"Value for '{key}' must be a positive integer, but got {value}.")
                elif key == "minShuffles" and isinstance(value, int) and value <= 0:
                    _add_error(errors, keyPath, f"Value for '{key}' must be a positive integer, but got {value}.")
                elif key == "l2DampingFactor" and isinstance(value, float) and value <= 0:
                    _add_error(errors, keyPath, f"Value for '{key}' must be a positive float, but got {value}.")
                elif key == "sagvolSmoothing" and not isinstance(value, bool):
                    _add_error(errors, keyPath, f"Value for '{key}' must be a boolean, but got {value}.")
    if sectionData["minShuffles"] > sectionData["maxShuffles"]:
        _add_error(errors, f"{sectionName}.minShuffles", f"Value for 'minShuffles' must be less than or equal to 'maxShuffles'.")

def _validate_misc_info(sectionData: Optional[Dict[str, Any]], sectionName: str, errors: Dict[str, str]):
    """Validates the miscInfo section."""
    if sectionData is None:
        _add_error(errors, sectionName, "Missing required section 'miscInfo'")
        return

    if not isinstance(sectionData, dict):
         _add_error(errors, sectionName, f"Section '{sectionName}' should be a dictionary, but got {type(sectionData).__name__}")
         return

    # Keys here match the expected config structure
    expectedKeys = {
        "availableCpus": int,
        "cleanUpLevel": str,
    }
    allowedCleanupLevels = ["None", "basic", "full", "brutal"]
    # currentlyWorkingCleanupLevels = ["None", "basic"] # Based on comment

    for key, expectedType in expectedKeys.items():
        keyPath = f"{sectionName}.{key}"
        if _check_key_exists(sectionData, key, sectionName, errors):
            value = sectionData[key]
            if _validate_type(value, expectedType, keyPath, errors):
                # Value range/specific value checks (only if type validation passed)
                if key == "availableCpus" and isinstance(value, int) and value <= 0:
                     _add_error(errors, keyPath, f"Value for '{key}' must be a positive integer, but got {value}.")
                elif key == "cleanUpLevel" and isinstance(value, str):
                    if not _validate_allowed_values(value, allowedCleanupLevels, keyPath, errors):
                        pass # Error already added by _validate_allowed_values
                    # Optional: Add a warning for levels mentioned as not working yet
                    # elif value not in currentlyWorkingCleanupLevels:
                    #     print(f"Warning: Config specifies cleanUpLevel='{value}', but comment indicates only {currentlyWorkingCleanupLevels} are implemented.")


# --- Main Validation Function (snake_case name) ---

def validate_config(config: Dict[str, Any]) -> Union[Dict[str, Any], None]:
    """
    Validates the structure, types, and allowed values of the configuration dictionary.

    Args:
        config: The configuration dictionary to validate.

    Returns:
        The original config dictionary if validation passes.
        Returns None if validation fails (after calling drSplash.show_config_error).
    """
    errors: Dict[str, str] = {}

    if not isinstance(config, dict):
        # Use a top-level key for this fundamental error
        _add_error(errors, "config", f"Configuration must be a dictionary, but got {type(config).__name__}")
        # Cannot proceed if the top level isn't a dict, show error and exit
        drSplash.show_config_error(errors)
        return None

    # Define expected top-level sections
    expectedSections = [
        "pathInfo",
        "moleculeInfo",
        "torsionScanInfo",
        "chargeFittingInfo",
        "parameterFittingInfo",
        "miscInfo",
    ]

    # Check if all expected sections are present at the top level
    for sectionName in expectedSections:
        if sectionName not in config:
            _add_error(errors, sectionName, f"Missing required section '{sectionName}'")
        # We'll still pass config.get(sectionName) to validators,
        # they handle None case gracefully.

    # Validate each section using helper functions
    _validate_path_info(config.get("pathInfo"), "pathInfo", errors)
    _validate_molecule_info(config.get("moleculeInfo"), "moleculeInfo", errors)
    _validate_torsion_scan_info(config.get("torsionScanInfo"), "torsionScanInfo", errors)
    _validate_charge_fitting_info(config.get("chargeFittingInfo"), "chargeFittingInfo", errors)
    _validate_parameter_fitting_info(config.get("parameterFittingInfo"), "parameterFittingInfo", errors)
    _validate_misc_info(config.get("miscInfo"), "miscInfo", errors)

    if len(errors) == 0:
        return config # Validation successful
    else:
        drSplash.show_config_error(errors)


# --- Example Usage (variables also changed to camelCase) ---