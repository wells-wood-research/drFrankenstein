import os
import yaml # Assuming the config is loaded from YAML, otherwise use dict directly
from typing import Dict, Any, List, Union, Optional, Tuple

def _add_error(errors: Dict[str, str], key_path: str, message: str):
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
        expected_str = (
            expected_type.__name__ 
            if isinstance(expected_type, type) 
            else ' or '.join(t.__name__ for t in expected_type)
        )
        _add_error(errors, key_path, f"Invalid type. Expected {expected_str}, but got {type(value).__name__}")
        return False
    return True

def _validate_allowed_values(value: Any, allowed_values: List[Any], key_path: str, errors: Dict[str, str]) -> bool:
    """Checks if the value is within the list of allowed values. Adds error if not."""
    if value not in allowed_values:
        _add_error(errors, key_path, f"Invalid value '{value}'. Allowed values are: {allowed_values}")
        return False
    return True

def _validate_path_info(section_data: Optional[Dict[str, Any]], section_name: str, errors: Dict[str, str]):
    """Validates the pathInfo section."""
    if section_data is None:
        _add_error(errors, section_name, "Missing required section 'pathInfo'")
        return

    if not isinstance(section_data, dict):
         _add_error(errors, section_name, f"Section '{section_name}' should be a dictionary, but got {type(section_data).__name__}")
         return

    expected_keys = {
        "inputDir": str,
        "outputDir": str,
        "multiWfnDir": str,
        "orcaExe": str,
        "gaff2Dat": str,
    }

    for key, expected_type in expected_keys.items():
        key_path = f"{section_name}.{key}"
        if _check_key_exists(section_data, key, section_name, errors):
            value = section_data[key]
            _validate_type(value, expected_type, key_path, errors)
            # Optional: Add os.path.exists checks here if needed, but prompt didn't require it
            # Example:
            # if key in ["inputDir", "multiWfnDir", "orcaExe", "gaff2Dat"]:
            #     if isinstance(value, str) and not os.path.exists(value):
            #          _add_error(errors, key_path, f"Path does not exist: '{value}'")
            # elif key == "outputDir":
            #     # outputDir might not exist yet, maybe check parent dir writability?
            #     pass


def _validate_molecule_info(section_data: Optional[Dict[str, Any]], section_name: str, errors: Dict[str, str]):
    """Validates the moleculeInfo section."""
    if section_data is None:
        _add_error(errors, section_name, "Missing required section 'moleculeInfo'")
        return

    if not isinstance(section_data, dict):
         _add_error(errors, section_name, f"Section '{section_name}' should be a dictionary, but got {type(section_data).__name__}")
         return

    # Check top-level keys and types
    expected_keys = {
        "charge": int,
        "multiplicity": int,
        "moleculePdb": str,
        "moleculeName": str,
        "nTermini": list,
        "cTermini": list,
        "chargeGroups": dict,
    }

    for key, expected_type in expected_keys.items():
        key_path = f"{section_name}.{key}"
        if _check_key_exists(section_data, key, section_name, errors):
            value = section_data[key]
            if _validate_type(value, expected_type, key_path, errors):
                # Specific validations for list/dict contents
                if key in ["nTermini", "cTermini"]:
                    if not all(isinstance(item, str) for item in value):
                         _add_error(errors, key_path, f"All items in list '{key}' must be strings.")
                elif key == "chargeGroups":
                    for group_name, group_data in value.items():
                        group_key_path = f"{key_path}.{group_name}"
                        if not isinstance(group_data, dict):
                            _add_error(errors, group_key_path, f"Charge group '{group_name}' must be a dictionary.")
                            continue # Skip further checks for this malformed group

                        # Check 'atoms' key within the group
                        if _check_key_exists(group_data, "atoms", group_key_path, errors):
                            atoms_value = group_data["atoms"]
                            atoms_key_path = f"{group_key_path}.atoms"
                            if _validate_type(atoms_value, list, atoms_key_path, errors):
                                if not all(isinstance(item, str) for item in atoms_value):
                                    _add_error(errors, atoms_key_path, "All items in 'atoms' list must be strings.")

                        # Check 'charge' key within the group
                        if _check_key_exists(group_data, "charge", group_key_path, errors):
                            charge_value = group_data["charge"]
                            charge_key_path = f"{group_key_path}.charge"
                            _validate_type(charge_value, int, charge_key_path, errors)


def _validate_torsion_scan_info(section_data: Optional[Dict[str, Any]], section_name: str, errors: Dict[str, str]):
    """Validates the torsionScanInfo section."""
    if section_data is None:
        _add_error(errors, section_name, "Missing required section 'torsionScanInfo'")
        return

    if not isinstance(section_data, dict):
         _add_error(errors, section_name, f"Section '{section_name}' should be a dictionary, but got {type(section_data).__name__}")
         return

    expected_keys = {
        "nConformers": int,
        "nScanSteps": int,
        "scanMethod": str,
        "scanSolvationMethod": str, # Expecting string, e.g., "ALPB(water)" or "Null"
        "singlePointMethod": str,   # Expecting string, e.g., "revPBE..." or "Null"
        "singlePointSolvationMethod": str, # Expecting string, e.g., "CPCM(water)" or "Null"
        "scanSinglePointsOn": str,
        "skipPhiPSi": bool,
    }

    for key, expected_type in expected_keys.items():
        key_path = f"{section_name}.{key}"
        if _check_key_exists(section_data, key, section_name, errors):
            value = section_data[key]
            if _validate_type(value, expected_type, key_path, errors):
                # Value range/specific value checks
                if key == "nConformers" and value < -1:
                    _add_error(errors, key_path, f"Value for '{key}' must be -1 or greater, but got {value}.")
                elif key == "nScanSteps" and value != 37:
                     # Warning in config implies only 37 works, strict check here
                    _add_error(errors, key_path, f"Value for '{key}' must be 37 (based on config comment), but got {value}.")
                elif key == "scanSinglePointsOn":
                     # Warning implies only 'all' works, strict check
                     _validate_allowed_values(value, ["all"], key_path, errors)


def _validate_charge_fitting_info(section_data: Optional[Dict[str, Any]], section_name: str, errors: Dict[str, str]):
    """Validates the chargeFittingInfo section."""
    if section_data is None:
        _add_error(errors, section_name, "Missing required section 'chargeFittingInfo'")
        return

    if not isinstance(section_data, dict):
         _add_error(errors, section_name, f"Section '{section_name}' should be a dictionary, but got {type(section_data).__name__}")
         return

    expected_keys = {
        "chargeFittingProtocol": str,
        "nConformers": int,
        "nCoresPerCalculation": int,
        "optMethod": str,
        "optSolvationMethod": str, # Expecting string, e.g., "ALPB(water)" or "Null"
        "singlePointMethod": str,
        "singlePointSolvationMethod": str, # Expecting string, e.g., "CPCM(water)" or "TIP3P" or "Null"
    }

    allowed_protocols = ["RESP", "RESP2", "SOLVATOR"]
    # Allowed solvation for QM/MM protocol
    allowed_qmmm_solvation = ["TIP3P"]

    for key, expected_type in expected_keys.items():
        key_path = f"{section_name}.{key}"
        if _check_key_exists(section_data, key, section_name, errors):
            value = section_data[key]
            if _validate_type(value, expected_type, key_path, errors):
                 # Value range/specific value checks
                if key == "chargeFittingProtocol":
                    _validate_allowed_values(value, allowed_protocols, key_path, errors)
                elif key == "nConformers" and value < -1:
                    _add_error(errors, key_path, f"Value for '{key}' must be -1 or greater, but got {value}.")
                elif key == "nCoresPerCalculation" and value <= 0:
                     _add_error(errors, key_path, f"Value for '{key}' must be a positive integer, but got {value}.")
                # Check QM/MM solvation compatibility
                elif key == "singlePointSolvationMethod":
                     protocol = section_data.get("chargeFittingProtocol")
                     if protocol == "SOLVATOR" and value not in allowed_qmmm_solvation:
                         _add_error(errors, key_path, f"For 'SOLVATOR' protocol, '{key}' must be one of {allowed_qmmm_solvation}, but got '{value}'.")
                     # Note: Could add checks for other protocols needing specific solvation, if known.


def _validate_parameter_fitting_info(section_data: Optional[Dict[str, Any]], section_name: str, errors: Dict[str, str]):
    """Validates the parameterFittingInfo section."""
    if section_data is None:
        _add_error(errors, section_name, "Missing required section 'parameterFittingInfo'")
        return

    if not isinstance(section_data, dict):
         _add_error(errors, section_name, f"Section '{section_name}' should be a dictionary, but got {type(section_data).__name__}")
         return

    expected_keys = {
        "forceField": str,
        "maxCosineFunctions": int,
        "nShuffles": int,
    }
    allowed_ff = ["CHARMM", "AMBER"]

    for key, expected_type in expected_keys.items():
        key_path = f"{section_name}.{key}"
        if _check_key_exists(section_data, key, section_name, errors):
            value = section_data[key]
            if _validate_type(value, expected_type, key_path, errors):
                # Value range/specific value checks
                if key == "forceField":
                    _validate_allowed_values(value, allowed_ff, key_path, errors)
                elif key == "maxCosineFunctions" and value <= 0:
                    _add_error(errors, key_path, f"Value for '{key}' must be a positive integer, but got {value}.")
                elif key == "nShuffles" and value <= 0:
                    _add_error(errors, key_path, f"Value for '{key}' must be a positive integer, but got {value}.")


def _validate_misc_info(section_data: Optional[Dict[str, Any]], section_name: str, errors: Dict[str, str]):
    """Validates the miscInfo section."""
    if section_data is None:
        _add_error(errors, section_name, "Missing required section 'miscInfo'")
        return

    if not isinstance(section_data, dict):
         _add_error(errors, section_name, f"Section '{section_name}' should be a dictionary, but got {type(section_data).__name__}")
         return

    expected_keys = {
        "availableCpus": int,
        "cleanUpLevel": str,
    }
    allowed_cleanup = ["None", "basic", "full", "brutal"]
    # currently_working_cleanup = ["None", "basic"] # Based on comment

    for key, expected_type in expected_keys.items():
        key_path = f"{section_name}.{key}"
        if _check_key_exists(section_data, key, section_name, errors):
            value = section_data[key]
            if _validate_type(value, expected_type, key_path, errors):
                # Value range/specific value checks
                if key == "availableCpus" and value <= 0:
                     _add_error(errors, key_path, f"Value for '{key}' must be a positive integer, but got {value}.")
                elif key == "cleanUpLevel":
                    if not _validate_allowed_values(value, allowed_cleanup, key_path, errors):
                        pass # Error already added by _validate_allowed_values
                    # Optional: Add a warning for levels mentioned as not working yet
                    # elif value not in currently_working_cleanup:
                    #     print(f"Warning: Config specifies cleanUpLevel='{value}', but comment indicates only {currently_working_cleanup} are implemented.")


def validate_config(config: Dict[str, Any]) -> Dict[str, str]:
    """
    Validates the structure, types, and allowed values of the configuration dictionary.

    Args:
        config: The configuration dictionary to validate.

    Returns:
        A dictionary containing validation errors. Keys are path-like strings
        indicating the location of the error (e.g., 'pathInfo.inputDir'),
        and values are strings describing the problem. An empty dictionary
        is returned if the configuration is valid according to the checks.
    """
    errors: Dict[str, str] = {}

    if not isinstance(config, dict):
        errors["config"] = f"Configuration must be a dictionary, but got {type(config).__name__}"
        return errors # Cannot proceed if the top level isn't a dict

    # Define expected top-level sections
    expected_sections = [
        "pathInfo",
        "moleculeInfo",
        "torsionScanInfo",
        "chargeFittingInfo",
        "parameterFittingInfo",
        "miscInfo",
    ]

    # Check if all expected sections are present
    for section_name in expected_sections:
        if section_name not in config:
            _add_error(errors, section_name, f"Missing required section '{section_name}'")

    # Validate each section using helper functions
    _validate_path_info(config.get("pathInfo"), "pathInfo", errors)
    _validate_molecule_info(config.get("moleculeInfo"), "moleculeInfo", errors)
    _validate_torsion_scan_info(config.get("torsionScanInfo"), "torsionScanInfo", errors)
    _validate_charge_fitting_info(config.get("chargeFittingInfo"), "chargeFittingInfo", errors)
    _validate_parameter_fitting_info(config.get("parameterFittingInfo"), "parameterFittingInfo", errors)
    _validate_misc_info(config.get("miscInfo"), "miscInfo", errors)

    return errors

