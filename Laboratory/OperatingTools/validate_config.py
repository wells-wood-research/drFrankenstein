import os
from typing import Dict, Any, List, Union, Optional, Tuple
from OperatingTools import drSplash
from typing import Dict, Any, List, Union, Optional, Tuple # Ensure all necessary types are imported

# Placeholder classes (extend if needed)
class FilePath:
    pass
class DirectoryPath:
    pass

def _add_error(errors: Dict[str, str], key_path: str, message: str) -> None: # Added return type hint
    """Adds an error message to the errors dictionary."""
    errors[key_path] = message

def _check_key_exists(data: Dict[str, Any], key: str, section_path: str, errors: Dict[str, str]) -> bool:
    """Checks if a key exists in the dictionary. Adds error if missing."""
    if key not in data:
        _add_error(errors, f"{section_path}.{key}", f"Missing required key '{key}'")
        return False
    return True

def _validate_type(value: Any, expected_type: Union[type, Tuple[type, ...]], key_path: str, errors: Dict[str, str]) -> bool:
    """Checks if the value's type matches the expected type(s). Adds error if mismatched."""
    if not isinstance(value, expected_type):
        # Improved error message formatting for None/NoneType
        if isinstance(expected_type, type):
            expectedTypeStr = 'None' if expected_type is type(None) else expected_type.__name__
        else: # It's a tuple
            typeNames = []
            for t in expected_type:
                typeNames.append('None' if t is type(None) else t.__name__)
            expectedTypeStr = ' or '.join(typeNames)

        _add_error(errors, key_path, f"Invalid type. Expected {expectedTypeStr}, but got {type(value).__name__}")
        return False
    return True

def _validate_allowed_values(value: Any, allowed_values: List[Any], key_path: str, errors: Dict[str, str]) -> bool:
    """Checks if the value is within the list of allowed values. Adds error if not."""
    # Skip check if value is None, as None might be allowed by type but not explicitly in allowedValues list
    if value is None:
        return True # Assume type check already validated None if it was allowed
    if value not in allowed_values:
        # Filter out None from displayed allowed values if it's present, as it's handled by type check
        displayAllowed = [v for v in allowed_values if v is not None]
        _add_error(errors, key_path, f"Invalid value '{value}'. Allowed values are: {displayAllowed}")
        return False
    return True

# --- Section Validators (snake_case names) ---

def _validate_path_info(section_data: Optional[Dict[str, Any]], section_name: str, errors: Dict[str, str]) -> None:
    """Validates the pathInfo section."""
    if section_data is None:
        _add_error(errors, section_name, "Missing required section 'pathInfo'")
        return

    if not isinstance(section_data, dict):
         _add_error(errors, section_name, f"Section '{section_name}' should be a dictionary, but got {type(section_data).__name__}")
         return

    # Keys here match the expected config structure, not changing them to camelCase
    expectedKeys = {
        "inputDir": str,
        "outputDir": str,
        "multiWfnDir": str,
        "orcaExe": str,
        "gaff2Dat": str,
    }

    for key, expectedType_val in expectedKeys.items(): # Renamed expectedType to avoid conflict
        keyPath = f"{section_name}.{key}"
        if _check_key_exists(section_data, key, section_name, errors):
            value = section_data[key]
            _validate_type(value, expectedType_val, keyPath, errors)
            # Optional: Add os.path.exists checks here if needed
            # if key in ["inputDir", "multiWfnDir", "orcaExe", "gaff2Dat"]:
            #     if isinstance(value, str) and not os.path.exists(value):
            #          _add_error(errors, keyPath, f"Path does not exist: '{value}'")
            # elif key == "outputDir":
            #     parentDir = os.path.dirname(os.path.abspath(value))
            #     if isinstance(value, str) and not os.access(parentDir, os.W_OK):
            #          _add_error(errors, keyPath, f"Output directory's parent is not writable: '{parentDir}'")


def _validate_molecule_info(section_data: Optional[Dict[str, Any]], section_name: str, errors: Dict[str, str]) -> None:
    """Validates the moleculeInfo section."""
    if section_data is None:
        _add_error(errors, section_name, "Missing required section 'moleculeInfo'")
        return

    if not isinstance(section_data, dict):
         _add_error(errors, section_name, f"Section '{section_name}' should be a dictionary, but got {type(section_data).__name__}")
         return

    # Keys here match the expected config structure
    expectedKeys = {
        "charge": int,
        "multiplicity": int,
        "moleculePdb": str,
        "moleculeName": str,
        "nTermini": list,
        "cTermini": list,
        "chargeGroups": dict,
    }

    for key, expectedType_val in expectedKeys.items(): # Renamed expectedType
        keyPath = f"{section_name}.{key}"
        if _check_key_exists(section_data, key, section_name, errors):
            value = section_data[key]
            if _validate_type(value, expectedType_val, keyPath, errors):
                # Specific validations for list/dict contents
                if key in ["nTermini", "cTermini"]:
                    if not all(isinstance(item, str) for item in value): # type: ignore
                         _add_error(errors, keyPath, f"All items in list '{key}' must be strings.")
                elif key == "chargeGroups":
                    # value is the chargeGroups dictionary
                    for groupName, groupData in value.items(): # type: ignore
                        groupKeyPath = f"{keyPath}.{groupName}"
                        if not isinstance(groupData, dict):
                            _add_error(errors, groupKeyPath, f"Charge group '{groupName}' must be a dictionary.")
                            continue # Skip further checks for this malformed group

                        # Check 'atoms' key within the group
                        if _check_key_exists(groupData, "atoms", groupKeyPath, errors):
                            atomsValue = groupData["atoms"]
                            atomsKeyPath = f"{groupKeyPath}.atoms"
                            if _validate_type(atomsValue, list, atomsKeyPath, errors):
                                if not all(isinstance(item, str) for item in atomsValue): # type: ignore
                                    _add_error(errors, atomsKeyPath, "All items in 'atoms' list must be strings.")

                        # Check 'charge' key within the group
                        if _check_key_exists(groupData, "charge", groupKeyPath, errors):
                            chargeValue = groupData["charge"]
                            chargeKeyPath = f"{groupKeyPath}.charge"
                            _validate_type(chargeValue, int, chargeKeyPath, errors)


def _validate_torsion_scan_info(section_data: Optional[Dict[str, Any]], section_name: str, errors: Dict[str, str]) -> None:
    """Validates the torsionScanInfo section."""
    if section_data is None:
        _add_error(errors, section_name, "Missing required section 'torsionScanInfo'")
        return

    if not isinstance(section_data, dict):
         _add_error(errors, section_name, f"Section '{section_name}' should be a dictionary, but got {type(section_data).__name__}")
         return

    # Keys here match the expected config structure
    # Use tuple (str, type(None)) for solvation methods to allow None
    expectedKeys = {
        "nConformers": int,
        "nScanSteps": int,
        "scanMethod": str,
        "scanSolvationMethod": (str, type(None)),  # Allow str or None
        "singlePointMethod": str,
        "singlePointSolvationMethod": (str, type(None)), # Allow str or None
        "scanSinglePointsOn": str,
        "skipPhiPSi": bool,
    }

    for key, expectedType_val in expectedKeys.items(): # Renamed expectedType
        keyPath = f"{section_name}.{key}"
        if _check_key_exists(section_data, key, section_name, errors):
            value = section_data[key]
            if _validate_type(value, expectedType_val, keyPath, errors):
                # Value range/specific value checks (only if type validation passed)
                if key == "nConformers" and isinstance(value, int) and value < -1:
                    _add_error(errors, keyPath, f"Value for '{key}' must be -1 or greater, but got {value}.")
                elif key == "nScanSteps" and isinstance(value, int) and value != 37:
                    _add_error(errors, keyPath, f"Value for '{key}' must be 37 (based on config comment), but got {value}.")
                elif key == "scanSinglePointsOn" and isinstance(value, str):
                     _validate_allowed_values(value, ["all"], keyPath, errors)


def _validate_charge_fitting_info(section_data: Optional[Dict[str, Any]], section_name: str, errors: Dict[str, str]) -> None:
    """Validates the chargeFittingInfo section."""
    if section_data is None:
        _add_error(errors, section_name, "Missing required section 'chargeFittingInfo'")
        return

    if not isinstance(section_data, dict):
         _add_error(errors, section_name, f"Section '{section_name}' should be a dictionary, but got {type(section_data).__name__}")
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

    for key, expectedType_val in expectedKeys.items(): # Renamed expectedType
        keyPath = f"{section_name}.{key}"
        if _check_key_exists(section_data, key, section_name, errors):
            value = section_data[key]
            if _validate_type(value, expectedType_val, keyPath, errors):
                 # Value range/specific value checks (only if type validation passed)
                if key == "chargeFittingProtocol" and isinstance(value, str):
                    _validate_allowed_values(value, allowedProtocols, keyPath, errors)
                elif key == "nConformers" and isinstance(value, int) and value < -1:
                    _add_error(errors, keyPath, f"Value for '{key}' must be -1 or greater, but got {value}.")
                elif key == "nCoresPerCalculation" and isinstance(value, int) and value <= 0:
                     _add_error(errors, keyPath, f"Value for '{key}' must be a positive integer, but got {value}.")


def _validate_parameter_fitting_info(section_data: Optional[Dict[str, Any]], section_name: str, errors: Dict[str, str]) -> None:
    """Validates the parameterFittingInfo section."""
    if section_data is None:
        _add_error(errors, section_name, "Missing required section 'parameterFittingInfo'")
        return

    if not isinstance(section_data, dict):
         _add_error(errors, section_name, f"Section '{section_name}' should be a dictionary, but got {type(section_data).__name__}")
         return

    # Keys here match the expected config structure
    expectedKeys = {
        "forceField": str,
        "maxCosineFunctions": int,
        "nShuffles": int,
    }
    allowedForceFields = ["CHARMM", "AMBER"]

    for key, expectedType_val in expectedKeys.items(): # Renamed expectedType
        keyPath = f"{section_name}.{key}"
        if _check_key_exists(section_data, key, section_name, errors):
            value = section_data[key]
            if _validate_type(value, expectedType_val, keyPath, errors):
                # Value range/specific value checks (only if type validation passed)
                if key == "forceField" and isinstance(value, str):
                    _validate_allowed_values(value, allowedForceFields, keyPath, errors)
                elif key == "maxCosineFunctions" and isinstance(value, int) and value <= 0:
                    _add_error(errors, keyPath, f"Value for '{key}' must be a positive integer, but got {value}.")
                elif key == "nShuffles" and isinstance(value, int) and value <= 0:
                    _add_error(errors, keyPath, f"Value for '{key}' must be a positive integer, but got {value}.")


def _validate_misc_info(section_data: Optional[Dict[str, Any]], section_name: str, errors: Dict[str, str]) -> None:
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