import os
from typing import Dict, Any, List, Union, Optional, Tuple
from OperatingTools import drSplash

def _add_error(errors: Dict[str, str], keyPath: str, message: str):
    """Adds an error message to the errors dictionary."""
    errors[keyPath] = message

def _check_key_exists(data: Dict[str, Any], key: str, sectionPath: str, errors: Dict[str, str]) -> bool:
    """Checks if a key exists in the dictionary. Adds error if missing."""
    if key not in data:
        _add_error(errors, f'{sectionPath}.{key}', f"Missing required key '{key}'")
        return False
    return True

def _validate_type(value: Any, expectedType: Union[type, Tuple[type, ...]], keyPath: str, errors: Dict[str, str]) -> bool:
    """Checks if the value's type matches the expected type(s). Adds error if mismatched."""
    if not isinstance(value, expectedType):
        if isinstance(expectedType, type):
            expectedTypeStr = 'None' if expectedType is type(None) else expectedType.__name__
        else:
            typeNames = []
            for t in expectedType:
                typeNames.append('None' if t is type(None) else t.__name__)
            expectedTypeStr = ' or '.join(typeNames)
        _add_error(errors, keyPath, f'Invalid type. Expected {expectedTypeStr}, but got {type(value).__name__}')
        return False
    return True

def _validate_allowed_values(value: Any, allowedValues: List[Any], keyPath: str, errors: Dict[str, str]) -> bool:
    """Checks if the value is within the list of allowed values. Adds error if not."""
    if value is None:
        return True
    if value not in allowedValues:
        displayAllowed = [v for v in allowedValues if v is not None]
        _add_error(errors, keyPath, f"Invalid value '{value}'. Allowed values are: {displayAllowed}")
        return False
    return True

def _validate_path_info(sectionData: Optional[Dict[str, Any]], sectionName: str, errors: Dict[str, str]):
    """Validates the pathInfo section."""
    if sectionData is None:
        _add_error(errors, sectionName, "Missing required section 'pathInfo'")
        return
    if not isinstance(sectionData, dict):
        _add_error(errors, sectionName, f"Section '{sectionName}' should be a dictionary, but got {type(sectionData).__name__}")
        return
    expectedKeys = {'inputDir': str, 'outputDir': str, 'multiWfnDir': str, 'orcaExe': str, 'gaff2Dat': str}
    for (key, expectedType) in expectedKeys.items():
        keyPath = f'{sectionName}.{key}'
        if _check_key_exists(sectionData, key, sectionName, errors):
            value = sectionData[key]
            _validate_type(value, expectedType, keyPath, errors)

def _validate_molecule_info(sectionData: Optional[Dict[str, Any]], sectionName: str, errors: Dict[str, str]):
    """Validates the moleculeInfo section."""
    if sectionData is None:
        _add_error(errors, sectionName, "Missing required section 'moleculeInfo'")
        return
    if not isinstance(sectionData, dict):
        _add_error(errors, sectionName, f"Section '{sectionName}' should be a dictionary, but got {type(sectionData).__name__}")
        return
    expectedKeys = {'charge': int, 'multiplicity': int, 'moleculePdb': str, 'moleculeName': str, 'nTermini': list, 'cTermini': list, 'chargeGroups': dict}
    for (key, expectedType) in expectedKeys.items():
        keyPath = f'{sectionName}.{key}'
        if _check_key_exists(sectionData, key, sectionName, errors):
            value = sectionData[key]
            if _validate_type(value, expectedType, keyPath, errors):
                if key in ['nTermini', 'cTermini']:
                    if not all((isinstance(item, str) for item in value)):
                        _add_error(errors, keyPath, f"All items in list '{key}' must be strings.")
                elif key == 'chargeGroups':
                    for (groupName, groupData) in value.items():
                        groupKeyPath = f'{keyPath}.{groupName}'
                        if not isinstance(groupData, dict):
                            _add_error(errors, groupKeyPath, f"Charge group '{groupName}' must be a dictionary.")
                            continue
                        if _check_key_exists(groupData, 'atoms', groupKeyPath, errors):
                            atomsValue = groupData['atoms']
                            atomsKeyPath = f'{groupKeyPath}.atoms'
                            if _validate_type(atomsValue, list, atomsKeyPath, errors):
                                if not all((isinstance(item, str) for item in atomsValue)):
                                    _add_error(errors, atomsKeyPath, "All items in 'atoms' list must be strings.")
                        if _check_key_exists(groupData, 'charge', groupKeyPath, errors):
                            chargeValue = groupData['charge']
                            chargeKeyPath = f'{groupKeyPath}.charge'
                            _validate_type(chargeValue, int, chargeKeyPath, errors)

def _validate_torsion_scan_info(sectionData: Optional[Dict[str, Any]], sectionName: str, errors: Dict[str, str]):
    """Validates the torsionScanInfo section."""
    if sectionData is None:
        _add_error(errors, sectionName, "Missing required section 'torsionScanInfo'")
        return
    if not isinstance(sectionData, dict):
        _add_error(errors, sectionName, f"Section '{sectionName}' should be a dictionary, but got {type(sectionData).__name__}")
        return
    expectedKeys = {'nConformers': int, 'nScanSteps': int, 'scanMethod': str, 'scanSolvationMethod': (str, type(None)), 'singlePointMethod': str, 'singlePointSolvationMethod': (str, type(None)), 'scanSinglePointsOn': str, 'skipPhiPSi': bool}
    for (key, expectedType) in expectedKeys.items():
        keyPath = f'{sectionName}.{key}'
        if _check_key_exists(sectionData, key, sectionName, errors):
            value = sectionData[key]
            if _validate_type(value, expectedType, keyPath, errors):
                if key == 'nConformers' and isinstance(value, int) and (value < -1):
                    _add_error(errors, keyPath, f"Value for '{key}' must be -1 or greater, but got {value}.")
                elif key == 'nScanSteps' and isinstance(value, int) and (value != 37):
                    _add_error(errors, keyPath, f"Value for '{key}' must be 37 (based on config comment), but got {value}.")
                elif key == 'scanSinglePointsOn' and isinstance(value, str):
                    _validate_allowed_values(value, ['all'], keyPath, errors)

def _validate_charge_fitting_info(sectionData: Optional[Dict[str, Any]], sectionName: str, errors: Dict[str, str]):
    """Validates the chargeFittingInfo section."""
    if sectionData is None:
        _add_error(errors, sectionName, "Missing required section 'chargeFittingInfo'")
        return
    if not isinstance(sectionData, dict):
        _add_error(errors, sectionName, f"Section '{sectionName}' should be a dictionary, but got {type(sectionData).__name__}")
        return
    expectedKeys = {'chargeFittingProtocol': str, 'nConformers': int, 'nCoresPerCalculation': int, 'optMethod': str, 'optSolvationMethod': (str, type(None)), 'singlePointMethod': str, 'singlePointSolvationMethod': (str, type(None))}
    allowedProtocols = ['RESP', 'RESP2', 'SOLVATOR']
    for (key, expectedType) in expectedKeys.items():
        keyPath = f'{sectionName}.{key}'
        if _check_key_exists(sectionData, key, sectionName, errors):
            value = sectionData[key]
            if _validate_type(value, expectedType, keyPath, errors):
                if key == 'chargeFittingProtocol' and isinstance(value, str):
                    _validate_allowed_values(value, allowedProtocols, keyPath, errors)
                elif key == 'nConformers' and isinstance(value, int) and (value < -1):
                    _add_error(errors, keyPath, f"Value for '{key}' must be -1 or greater, but got {value}.")
                elif key == 'nCoresPerCalculation' and isinstance(value, int) and (value <= 0):
                    _add_error(errors, keyPath, f"Value for '{key}' must be a positive integer, but got {value}.")

def _validate_parameter_fitting_info(sectionData: Optional[Dict[str, Any]], sectionName: str, errors: Dict[str, str]):
    """Validates the parameterFittingInfo section."""
    if sectionData is None:
        _add_error(errors, sectionName, "Missing required section 'parameterFittingInfo'")
        return
    if not isinstance(sectionData, dict):
        _add_error(errors, sectionName, f"Section '{sectionName}' should be a dictionary, but got {type(sectionData).__name__}")
        return
    expectedKeys = {'forceField': str, 'maxCosineFunctions': int, 'nShuffles': int}
    allowedForceFields = ['CHARMM', 'AMBER']
    for (key, expectedType) in expectedKeys.items():
        keyPath = f'{sectionName}.{key}'
        if _check_key_exists(sectionData, key, sectionName, errors):
            value = sectionData[key]
            if _validate_type(value, expectedType, keyPath, errors):
                if key == 'forceField' and isinstance(value, str):
                    _validate_allowed_values(value, allowedForceFields, keyPath, errors)
                elif key == 'maxCosineFunctions' and isinstance(value, int) and (value <= 0):
                    _add_error(errors, keyPath, f"Value for '{key}' must be a positive integer, but got {value}.")
                elif key == 'nShuffles' and isinstance(value, int) and (value <= 0):
                    _add_error(errors, keyPath, f"Value for '{key}' must be a positive integer, but got {value}.")

def _validate_misc_info(sectionData: Optional[Dict[str, Any]], sectionName: str, errors: Dict[str, str]):
    """Validates the miscInfo section."""
    if sectionData is None:
        _add_error(errors, sectionName, "Missing required section 'miscInfo'")
        return
    if not isinstance(sectionData, dict):
        _add_error(errors, sectionName, f"Section '{sectionName}' should be a dictionary, but got {type(sectionData).__name__}")
        return
    expectedKeys = {'availableCpus': int, 'cleanUpLevel': str}
    allowedCleanupLevels = ['None', 'basic', 'full', 'brutal']
    for (key, expectedType) in expectedKeys.items():
        keyPath = f'{sectionName}.{key}'
        if _check_key_exists(sectionData, key, sectionName, errors):
            value = sectionData[key]
            if _validate_type(value, expectedType, keyPath, errors):
                if key == 'availableCpus' and isinstance(value, int) and (value <= 0):
                    _add_error(errors, keyPath, f"Value for '{key}' must be a positive integer, but got {value}.")
                elif key == 'cleanUpLevel' and isinstance(value, str):
                    if not _validate_allowed_values(value, allowedCleanupLevels, keyPath, errors):
                        pass

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
        _add_error(errors, 'config', f'Configuration must be a dictionary, but got {type(config).__name__}')
        drSplash.show_config_error(errors)
        return None
    expectedSections = ['pathInfo', 'moleculeInfo', 'torsionScanInfo', 'chargeFittingInfo', 'parameterFittingInfo', 'miscInfo']
    for sectionName in expectedSections:
        if sectionName not in config:
            _add_error(errors, sectionName, f"Missing required section '{sectionName}'")
    _validate_path_info(config.get('pathInfo'), 'pathInfo', errors)
    _validate_molecule_info(config.get('moleculeInfo'), 'moleculeInfo', errors)
    _validate_torsion_scan_info(config.get('torsionScanInfo'), 'torsionScanInfo', errors)
    _validate_charge_fitting_info(config.get('chargeFittingInfo'), 'chargeFittingInfo', errors)
    _validate_parameter_fitting_info(config.get('parameterFittingInfo'), 'parameterFittingInfo', errors)
    _validate_misc_info(config.get('miscInfo'), 'miscInfo', errors)
    if len(errors) == 0:
        return config
    else:
        drSplash.show_config_error(errors)