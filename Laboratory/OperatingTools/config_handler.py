import os
import io
import multiprocessing as mp
from typing import Dict, Any, List, Union, Optional, Tuple

from ruamel.yaml import YAML
import openmm as _openmm

try:
    from OperatingTools import drSplash
except ImportError:  # pragma: no cover - fallback for package-style imports
    from Laboratory.OperatingTools import drSplash


def _has_fatal_errors(d: object) -> bool:
    """Recursively checks a dictionary for fatal error messages."""
    if isinstance(d, dict):
        return any(_has_fatal_errors(v) for v in d.values())
    elif isinstance(d, str) and "default" not in d.lower():
        return True
    return False


def _detect_gpu_platform() -> str:
    """Detects available GPU platforms in order of preference: CUDA > HIP > OpenCL."""
    for candidate in ["CUDA", "HIP", "OpenCL"]:
        try:
            _openmm.Platform.getPlatformByName(candidate)
            return candidate
        except Exception:
            continue
    return "CPU"


def set_default(section: dict, key: str, default_value: object, error_dict: dict) -> None:
    """Apply a default value and record it in the error dictionary when used."""
    if key not in section:
        error_dict[key] = f"Default Used: {default_value}"
    section.setdefault(key, default_value)


def apply_defaults_and_validate(config: dict) -> dict:
    """Apply config defaults and run lightweight validation plus default reporting."""
    errors = {}

    # --- Section: moleculeInfo ---
    errors["moleculeInfo"] = {}
    moleculeInfo = config.get("moleculeInfo", {})
    for key, expected_type in [("charge", int), ("multiplicity", int), ("moleculeName", str)]:
        if key not in moleculeInfo:
            errors["moleculeInfo"][key] = "Missing required key."
        elif not isinstance(moleculeInfo[key], expected_type):
            errors["moleculeInfo"][key] = f"Must be type {expected_type}, not {type(moleculeInfo[key])}."
        else:
            errors["moleculeInfo"][key] = None
    for key, expected_type in [("chargeGroups", dict), ("backboneAliases", dict)]:
        if key not in moleculeInfo:
            errors["moleculeInfo"][key] = f"Default Used: {None}"
        elif not isinstance(moleculeInfo[key], (expected_type, type(None))):
            errors["moleculeInfo"][key] = f"Must be type {expected_type} or None."
        else:
            errors["moleculeInfo"][key] = None
    moleculeInfo.setdefault("chargeGroups", None)
    moleculeInfo.setdefault("backboneAliases", None)

    # --- Section: conformerGenerationInfo ---
    errors["conformerGenerationInfo"] = {}
    conformerGenerationInfo = config.setdefault("conformerGenerationInfo", {})
    if not isinstance(conformerGenerationInfo, dict):
        errors["conformerGenerationInfo"]["section"] = "Must be a dictionary."
        conformerGenerationInfo = {}
        config["conformerGenerationInfo"] = conformerGenerationInfo
    set_default(conformerGenerationInfo, "goatMode", "GOAT", errors["conformerGenerationInfo"])
    set_default(conformerGenerationInfo, "energyCutoff", 6.0, errors["conformerGenerationInfo"])
    set_default(conformerGenerationInfo, "goatMethod", "XTB2", errors["conformerGenerationInfo"])
    set_default(conformerGenerationInfo, "conformerSelction", "ENERGY", errors["conformerGenerationInfo"])

    # --- Section: pathInfo ---
    errors["pathInfo"] = {}
    pathInfo = config.setdefault("pathInfo", {})
    if "inputDir" not in pathInfo:
        errors["pathInfo"]["inputDir"] = f"Default Used: {os.getcwd()}"
    pathInfo.setdefault("inputDir", os.getcwd())

    if "outputDir" not in pathInfo:
        if moleculeInfo.get("moleculeName") and errors["moleculeInfo"]["moleculeName"] is None:
            default_path = os.path.join(os.getcwd(), f"{moleculeInfo['moleculeName']}_FrankenParams")
            pathInfo["outputDir"] = default_path
            errors["pathInfo"]["outputDir"] = f"Default Used: {default_path}"
        else:
            errors["pathInfo"]["outputDir"] = "Cannot set default, 'moleculeName' is missing or invalid."

    if "amberHome" not in pathInfo:
        amber_home = os.environ.get("AMBERHOME")
        errors["pathInfo"]["amberHome"] = f"Default Used: {amber_home}"
        if amber_home is not None:
            pathInfo["amberHome"] = amber_home

    if "cgenffExe" not in pathInfo:
        errors["pathInfo"]["cgenffExe"] = f"Default Used: {None}"
    pathInfo.setdefault("cgenffExe", None)

    for key, expected_type in [("inputDir", str), ("outputDir", str), ("amberHome", str), ("cgenffExe", str)]:
        if key not in errors["pathInfo"]:
            errors["pathInfo"][key] = None
        if key in pathInfo and not isinstance(pathInfo[key], (expected_type, type(None))):
            errors["pathInfo"][key] = f"Must be type {expected_type} or None."

    for key, is_dir in [("multiWfnDir", True), ("orcaExe", False)]:
        path_val = pathInfo.get(key)
        if not path_val:
            errors["pathInfo"][key] = "Missing required key."
        elif not isinstance(path_val, str):
            errors["pathInfo"][key] = f"Must be type str, not {type(path_val)}."
        elif is_dir and not os.path.isdir(path_val):
            errors["pathInfo"][key] = f"Directory not found: {path_val}"
        elif not is_dir and not os.path.isfile(path_val):
            errors["pathInfo"][key] = f"File not found: {path_val}"
        else:
            if key not in errors["pathInfo"]:
                errors["pathInfo"][key] = None

    # --- Section: torsionScanInfo ---
    errors["torsionScanInfo"] = {}
    torsionScanInfo = config.setdefault("torsionScanInfo", {})
    if not isinstance(torsionScanInfo, dict):
        errors["torsionScanInfo"]["section"] = "Must be a dictionary."
        torsionScanInfo = {}
        config["torsionScanInfo"] = torsionScanInfo
    runScansOn = torsionScanInfo.setdefault("runScansOn", {})
    if not isinstance(runScansOn, dict):
        errors["torsionScanInfo"]["runScansOn"] = "Must be type dict."
        runScansOn = {}
        torsionScanInfo["runScansOn"] = runScansOn
    errors["torsionScanInfo"]["runScansOn"] = {}
    for key in ["phiPsi", "polarProtons", "nonPolarProtons", "amides", "nonAromaticRings"]:
        if key not in runScansOn:
            default_val = key in ["phiPsi", "polarProtons"]
            runScansOn[key] = default_val
            errors["torsionScanInfo"]["runScansOn"][key] = f"Default Used: {default_val}"
        elif not isinstance(runScansOn[key], bool):
            errors["torsionScanInfo"]["runScansOn"][key] = "Must be type bool."
        else:
            errors["torsionScanInfo"]["runScansOn"][key] = None

    if "scanMethod" not in torsionScanInfo:
        errors["torsionScanInfo"]["scanMethod"] = "Required: use a SIMPLE INPUT line from the ORCA input library"
    elif not isinstance(torsionScanInfo["scanMethod"], str):
        errors["torsionScanInfo"]["scanMethod"] = "Must be type str."
    else:
        errors["torsionScanInfo"]["scanMethod"] = None

    set_default(torsionScanInfo, "nConformers", -1, errors["torsionScanInfo"])
    set_default(torsionScanInfo, "nCoresPerCalculation", 1, errors["torsionScanInfo"])
    set_default(torsionScanInfo, "scanSolvationMethod", None, errors["torsionScanInfo"])
    set_default(torsionScanInfo, "singlePointMethod", None, errors["torsionScanInfo"])
    set_default(torsionScanInfo, "singlePointSolvationMethod", None, errors["torsionScanInfo"])

    # --- Section: chargeFittingInfo ---
    errors["chargeFittingInfo"] = {}
    chargeFittingInfo = config.setdefault("chargeFittingInfo", {})
    for key in ["chargeFittingProtocol", "optMethod", "singlePointMethod"]:
        if key not in chargeFittingInfo:
            errors["chargeFittingInfo"][key] = "Missing required key."
        else:
            errors["chargeFittingInfo"][key] = None

    set_default(chargeFittingInfo, "nConformers", -1, errors["chargeFittingInfo"])
    set_default(chargeFittingInfo, "nCoresPerCalculation", 1, errors["chargeFittingInfo"])
    set_default(chargeFittingInfo, "enforceDefaultBackboneCharges", False, errors["chargeFittingInfo"])
    if chargeFittingInfo.get("chargeFittingProtocol") == "SOLVATOR":
        set_default(chargeFittingInfo, "waterDensity", 10, errors["chargeFittingInfo"])
    else:
        set_default(chargeFittingInfo, "waterDensity", None, errors["chargeFittingInfo"])

    if not isinstance(chargeFittingInfo.get("enforceDefaultBackboneCharges"), bool):
        errors["chargeFittingInfo"]["enforceDefaultBackboneCharges"] = "Must be type bool."

    # --- Section: parameterFittingInfo ---
    errors["parameterFittingInfo"] = {}
    parameterFittingInfo = config.setdefault("parameterFittingInfo", {})
    if "forceField" not in parameterFittingInfo:
        errors["parameterFittingInfo"]["forceField"] = "Missing required key."
    else:
        errors["parameterFittingInfo"]["forceField"] = None

    set_default(parameterFittingInfo, "maxCosineFunctions", 4, errors["parameterFittingInfo"])
    set_default(parameterFittingInfo, "maxShuffles", 50, errors["parameterFittingInfo"])
    set_default(parameterFittingInfo, "minShuffles", 10, errors["parameterFittingInfo"])
    set_default(parameterFittingInfo, "maeTolTotal", None, errors["parameterFittingInfo"])
    set_default(parameterFittingInfo, "maeTolTorsion", None, errors["parameterFittingInfo"])
    set_default(parameterFittingInfo, "l2DampingFactor", 0.1, errors["parameterFittingInfo"])
    set_default(parameterFittingInfo, "sagvolSmoothing", True, errors["parameterFittingInfo"])
    set_default(parameterFittingInfo, "cosineMinScoreImprovement", 1e-4, errors["parameterFittingInfo"])

    # --- Section: miscInfo ---
    errors["miscInfo"] = {}
    miscInfo = config.setdefault("miscInfo", {})
    set_default(miscInfo, "availableCpus", mp.cpu_count(), errors["miscInfo"])
    set_default(miscInfo, "cleanUpLevel", 1, errors["miscInfo"])
    set_default(miscInfo, "seed", 1818, errors["miscInfo"])
    set_default(miscInfo, "debug", False, errors["miscInfo"])

    if "gpuPlatform" not in miscInfo or miscInfo.get("gpuPlatform") is None:
        detected = _detect_gpu_platform()
        miscInfo["gpuPlatform"] = detected
        errors["miscInfo"]["gpuPlatform"] = f"Default Used: {detected}"
    else:
        errors["miscInfo"]["gpuPlatform"] = None

    if _has_fatal_errors(errors):
        drSplash.print_config_error(errors)

    return config


def _add_error(errors: Dict[str, str], keyPath: str, message: str) -> None:
    """Record a validation error for a specific config path."""
    errors[keyPath] = message


def _check_key_exists(data: Dict[str, Any], key: str, sectionPath: str, errors: Dict[str, str]) -> bool:
    """Check whether a required key exists and record an error if it does not."""
    if key not in data:
        _add_error(errors, f"{sectionPath}.{key}", f"Missing required key '{key}'")
        return False
    return True


def _validate_type(value: Any, expectedType: Union[type, Tuple[type, ...]], keyPath: str, errors: Dict[str, str]) -> bool:
    """Check a value's type against the expected type or tuple of types."""
    if not isinstance(value, expectedType):
        if isinstance(expectedType, type):
            expectedTypeStr = "None" if expectedType is type(None) else expectedType.__name__
        else:
            typeNames = []
            for t in expectedType:
                typeNames.append("None" if t is type(None) else t.__name__)
            expectedTypeStr = " or ".join(typeNames)

        _add_error(errors, keyPath, f"Invalid type. Expected {expectedTypeStr}, but got {type(value).__name__}")
        return False
    return True


def _validate_allowed_values(value: Any, allowedValues: List[Any], keyPath: str, errors: Dict[str, str]) -> bool:
    """Check whether a value is one of the allowed configuration values."""
    if value is None:
        return True
    if value not in allowedValues:
        displayAllowed = [v for v in allowedValues if v is not None]
        _add_error(errors, keyPath, f"Invalid value '{value}'. Allowed values are: {displayAllowed}")
        return False
    return True


def _validate_path_info(sectionData: Optional[Dict[str, Any]], sectionName: str, errors: Dict[str, str]) -> None:
    """Validate the `pathInfo` section."""
    if sectionData is None:
        _add_error(errors, sectionName, "Missing required section 'pathInfo'")
        return
    if not isinstance(sectionData, dict):
        _add_error(errors, sectionName, f"Section '{sectionName}' should be a dictionary, but got {type(sectionData).__name__}")
        return
    requiredKeys = {"inputDir": str, "outputDir": str, "multiWfnDir": str, "orcaExe": str}
    optionalKeys = {"amberHome": (str, type(None)), "cgenffExe": (str, type(None))}
    for key, expectedType in requiredKeys.items():
        keyPath = f"{sectionName}.{key}"
        if _check_key_exists(sectionData, key, sectionName, errors):
            _validate_type(sectionData[key], expectedType, keyPath, errors)
    for key, expectedType in optionalKeys.items():
        if key in sectionData:
            keyPath = f"{sectionName}.{key}"
            _validate_type(sectionData[key], expectedType, keyPath, errors)


def _validate_molecule_info(sectionData: Optional[Dict[str, Any]], sectionName: str, errors: Dict[str, str]) -> None:
    """Validate the `moleculeInfo` section."""
    if sectionData is None:
        _add_error(errors, sectionName, "Missing required section 'moleculeInfo'")
        return
    if not isinstance(sectionData, dict):
        _add_error(errors, sectionName, f"Section '{sectionName}' should be a dictionary, but got {type(sectionData).__name__}")
        return
    requiredKeys = {"charge": int, "multiplicity": int, "moleculeName": str}
    optionalKeys = {"chargeGroups": (dict, type(None)), "backboneAliases": (dict, type(None))}
    for key, expectedType in requiredKeys.items():
        keyPath = f"{sectionName}.{key}"
        if _check_key_exists(sectionData, key, sectionName, errors):
            _validate_type(sectionData[key], expectedType, keyPath, errors)
    for key, expectedType in optionalKeys.items():
        if key in sectionData:
            keyPath = f"{sectionName}.{key}"
            value = sectionData[key]
            if _validate_type(value, expectedType, keyPath, errors):
                if key == "chargeGroups" and isinstance(value, dict):
                    for groupName, groupData in value.items():
                        groupKeyPath = f"{keyPath}.{groupName}"
                        if not isinstance(groupData, dict):
                            _add_error(errors, groupKeyPath, f"Charge group '{groupName}' must be a dictionary.")
                            continue
                        if _check_key_exists(groupData, "atoms", groupKeyPath, errors):
                            atomsValue = groupData["atoms"]
                            atomsKeyPath = f"{groupKeyPath}.atoms"
                            if _validate_type(atomsValue, list, atomsKeyPath, errors):
                                if not all(isinstance(item, str) for item in atomsValue):
                                    _add_error(errors, atomsKeyPath, "All items in 'atoms' list must be strings.")
                        if _check_key_exists(groupData, "charge", groupKeyPath, errors):
                            _validate_type(groupData["charge"], int, f"{groupKeyPath}.charge", errors)


def _validate_torsion_scan_info(sectionData: Optional[Dict[str, Any]], sectionName: str, errors: Dict[str, str]) -> None:
    """Validate the `torsionScanInfo` section."""
    if sectionData is None:
        _add_error(errors, sectionName, "Missing required section 'torsionScanInfo'")
        return
    if not isinstance(sectionData, dict):
        _add_error(errors, sectionName, f"Section '{sectionName}' should be a dictionary, but got {type(sectionData).__name__}")
        return
    requiredKeys = {"runScansOn": dict, "nConformers": int, "scanMethod": str}
    optionalKeys = {
        "scanSolvationMethod": (str, type(None)),
        "singlePointMethod": (str, type(None)),
        "singlePointSolvationMethod": (str, type(None)),
        "nCoresPerCalculation": int,
    }
    for key, expectedType in requiredKeys.items():
        keyPath = f"{sectionName}.{key}"
        if _check_key_exists(sectionData, key, sectionName, errors):
            value = sectionData[key]
            if _validate_type(value, expectedType, keyPath, errors):
                if key == "nConformers" and isinstance(value, int) and value < -1:
                    _add_error(errors, keyPath, f"Value for '{key}' must be -1 or greater, but got {value}.")
    for key, expectedType in optionalKeys.items():
        if key in sectionData:
            keyPath = f"{sectionName}.{key}"
            value = sectionData[key]
            if _validate_type(value, expectedType, keyPath, errors):
                if key == "nCoresPerCalculation" and isinstance(value, int) and value <= 0:
                    _add_error(errors, keyPath, f"Value for '{key}' must be a positive integer, but got {value}.")
    if "runScansOn" in sectionData and isinstance(sectionData["runScansOn"], dict):
        runScansOn = sectionData["runScansOn"]
        runScansKeyPath = f"{sectionName}.runScansOn"
        for scanType, enabled in runScansOn.items():
            keyPath = f"{runScansKeyPath}.{scanType}"
            if not isinstance(enabled, bool):
                _add_error(errors, keyPath, f"Scan type '{scanType}' must be a boolean, but got {type(enabled).__name__}")


def _validate_conformer_generation_info(sectionData: Optional[Dict[str, Any]], sectionName: str, errors: Dict[str, str]) -> None:
    """Validate the `conformerGenerationInfo` section."""
    if sectionData is None:
        _add_error(errors, sectionName, "Missing required section 'conformerGenerationInfo'")
        return
    if not isinstance(sectionData, dict):
        _add_error(errors, sectionName, f"Section '{sectionName}' should be a dictionary, but got {type(sectionData).__name__}")
        return
    requiredKeys = {"goatMode": str, "energyCutoff": (int, float), "goatMethod": str, "conformerSelction": str}
    allowedGoatModes = ["GOAT", "GOAT-ENTROPY"]
    allowedSelectionMethods = ["ENERGY", "DIVERSE"]
    allowedGoatMethods = ["XTB2", "GFN-FF"]
    for key, expectedType in requiredKeys.items():
        keyPath = f"{sectionName}.{key}"
        if _check_key_exists(sectionData, key, sectionName, errors):
            value = sectionData[key]
            if _validate_type(value, expectedType, keyPath, errors):
                if key == "goatMode":
                    _validate_allowed_values(value, allowedGoatModes, keyPath, errors)
                elif key == "goatMethod":
                    _validate_allowed_values(value, allowedGoatMethods, keyPath, errors)
                elif key == "conformerSelction":
                    _validate_allowed_values(value, allowedSelectionMethods, keyPath, errors)
                elif key == "energyCutoff" and isinstance(value, (int, float)) and value < 0:
                    _add_error(errors, keyPath, f"Value for '{key}' must be 0 or greater, but got {value}.")


def _validate_charge_fitting_info(sectionData: Optional[Dict[str, Any]], sectionName: str, errors: Dict[str, str]) -> None:
    """Validate the `chargeFittingInfo` section."""
    if sectionData is None:
        _add_error(errors, sectionName, "Missing required section 'chargeFittingInfo'")
        return
    if not isinstance(sectionData, dict):
        _add_error(errors, sectionName, f"Section '{sectionName}' should be a dictionary, but got {type(sectionData).__name__}")
        return
    requiredKeys = {
        "chargeFittingProtocol": str,
        "nConformers": int,
        "nCoresPerCalculation": int,
        "optMethod": str,
        "singlePointMethod": str,
    }
    optionalKeys = {
        "optSolvationMethod": (str, type(None)),
        "singlePointSolvationMethod": (str, type(None)),
        "waterDensity": (int, float, type(None)),
        "enforceDefaultBackboneCharges": bool,
    }
    allowedProtocols = ["RESP", "RESP2", "SOLVATOR"]
    for key, expectedType in requiredKeys.items():
        keyPath = f"{sectionName}.{key}"
        if _check_key_exists(sectionData, key, sectionName, errors):
            value = sectionData[key]
            if _validate_type(value, expectedType, keyPath, errors):
                if key == "chargeFittingProtocol" and isinstance(value, str):
                    _validate_allowed_values(value, allowedProtocols, keyPath, errors)
                elif key == "nConformers" and isinstance(value, int) and value < -1:
                    _add_error(errors, keyPath, f"Value for '{key}' must be -1 or greater, but got {value}.")
                elif key == "nCoresPerCalculation" and isinstance(value, int) and value <= 0:
                    _add_error(errors, keyPath, f"Value for '{key}' must be a positive integer, but got {value}.")
    for key, expectedType in optionalKeys.items():
        if key in sectionData:
            keyPath = f"{sectionName}.{key}"
            _validate_type(sectionData[key], expectedType, keyPath, errors)


def _validate_parameter_fitting_info(sectionData: Optional[Dict[str, Any]], sectionName: str, errors: Dict[str, str]) -> None:
    """Validates the parameterFittingInfo section."""
    if sectionData is None:
        _add_error(errors, sectionName, "Missing required section 'parameterFittingInfo'")
        return
    if not isinstance(sectionData, dict):
        _add_error(errors, sectionName, f"Section '{sectionName}' should be a dictionary, but got {type(sectionData).__name__}")
        return
    expectedKeys = {
        "forceField": str,
        "maxCosineFunctions": int,
        "maxShuffles": int,
        "minShuffles": int,
        "l2DampingFactor": float | None,
        "sagvolSmoothing": (bool, dict, type(None)),
    }
    allowedForceFields = ["CHARMM", "AMBER"]
    for key, expectedType in expectedKeys.items():
        keyPath = f"{sectionName}.{key}"
        if _check_key_exists(sectionData, key, sectionName, errors):
            value = sectionData[key]
            if _validate_type(value, expectedType, keyPath, errors):
                if key == "forceField" and isinstance(value, str):
                    _validate_allowed_values(value, allowedForceFields, keyPath, errors)
                elif key in ("maxCosineFunctions", "maxShuffles", "minShuffles") and isinstance(value, int) and value <= 0:
                    _add_error(errors, keyPath, f"Value for '{key}' must be a positive integer, but got {value}.")
                elif key == "l2DampingFactor" and isinstance(value, float) and value <= 0:
                    _add_error(errors, keyPath, f"Value for '{key}' must be a positive float, but got {value}.")
    if (
        "minShuffles" in sectionData
        and "maxShuffles" in sectionData
        and isinstance(sectionData["minShuffles"], int)
        and isinstance(sectionData["maxShuffles"], int)
        and sectionData["minShuffles"] > sectionData["maxShuffles"]
    ):
        _add_error(errors, f"{sectionName}.minShuffles", "Value for 'minShuffles' must be less than or equal to 'maxShuffles'.")


def _validate_misc_info(sectionData: Optional[Dict[str, Any]], sectionName: str, errors: Dict[str, str]) -> None:
    """Validates the miscInfo section."""
    if sectionData is None:
        _add_error(errors, sectionName, "Missing required section 'miscInfo'")
        return
    if not isinstance(sectionData, dict):
        _add_error(errors, sectionName, f"Section '{sectionName}' should be a dictionary, but got {type(sectionData).__name__}")
        return
    requiredKeys = {"availableCpus": int}
    optionalKeys = {
        "cleanUpLevel": int,
        "assemblyProtocol": str,
        "seed": int,
        "conformerSelectionMethods": str,
        "debug": bool,
        "gpuPlatform": (str, type(None)),
    }
    allowedAssemblyProtocols = ["ANTECHAMBER", "CGENFF", "AGNOSTIC"]
    allowedSelectionMethods = ["ENERGY", "DIVERSE"]
    for key, expectedType in requiredKeys.items():
        keyPath = f"{sectionName}.{key}"
        if _check_key_exists(sectionData, key, sectionName, errors):
            value = sectionData[key]
            if _validate_type(value, expectedType, keyPath, errors):
                if key == "availableCpus" and isinstance(value, int) and value <= 0:
                    _add_error(errors, keyPath, f"Value for '{key}' must be a positive integer, but got {value}.")
    for key, expectedType in optionalKeys.items():
        if key in sectionData:
            keyPath = f"{sectionName}.{key}"
            value = sectionData[key]
            if _validate_type(value, expectedType, keyPath, errors):
                if key == "cleanUpLevel" and isinstance(value, int) and (value < 0 or value > 3):
                    _add_error(errors, keyPath, f"Value for '{key}' must be between 0 and 3, but got {value}.")
                elif key == "assemblyProtocol" and isinstance(value, str):
                    _validate_allowed_values(value, allowedAssemblyProtocols, keyPath, errors)
                elif key == "seed" and isinstance(value, int) and value < 0:
                    _add_error(errors, keyPath, f"Value for '{key}' must be 0 or greater, but got {value}.")
                elif key == "conformerSelectionMethods" and isinstance(value, str):
                    if value.upper() in allowedSelectionMethods:
                        continue
                    _validate_allowed_values(value, allowedSelectionMethods, keyPath, errors)
                elif key == "gpuPlatform" and isinstance(value, str):
                    allowedGpu = ["CUDA", "HIP", "OpenCL", "CPU"]
                    if value.upper() not in allowedGpu:
                        _add_error(errors, keyPath, f"Invalid value '{value}'. Allowed values are: {allowedGpu}")


def _validate_backbone_charge_enforcement(config: Dict[str, Any], errors: Dict[str, str]) -> None:
    """Cross-section validation for `enforceDefaultBackboneCharges`."""
    chargeFittingInfo = config.get("chargeFittingInfo")
    moleculeInfo = config.get("moleculeInfo")
    if not isinstance(chargeFittingInfo, dict) or not isinstance(moleculeInfo, dict):
        return
    if not chargeFittingInfo.get("enforceDefaultBackboneCharges", False):
        return
    backboneAliases = moleculeInfo.get("backboneAliases")
    keyPathPrefix = "moleculeInfo.backboneAliases"
    if not isinstance(backboneAliases, dict):
        _add_error(
            errors,
            keyPathPrefix,
            "When chargeFittingInfo.enforceDefaultBackboneCharges is True, moleculeInfo.backboneAliases must be a dictionary containing N, CA, HA, C, O (and H for non-proline residues).",
        )
        return
    ## Proline-type residues have a ring nitrogen with no backbone amide H, so an
    ## "H" alias may be absent or empty; only require it for non-proline residues.
    isProline = not backboneAliases.get("H")
    requiredBackboneKeys = ["N", "CA", "HA", "C", "O"] if isProline else ["N", "H", "CA", "HA", "C", "O"]
    for backboneKey in requiredBackboneKeys:
        aliasKeyPath = f"{keyPathPrefix}.{backboneKey}"
        if backboneKey not in backboneAliases:
            _add_error(errors, aliasKeyPath, f"Missing required key '{backboneKey}' when enforceDefaultBackboneCharges is True.")
            continue
        aliases = backboneAliases[backboneKey]
        if not isinstance(aliases, list):
            _add_error(errors, aliasKeyPath, f"Invalid type. Expected list, but got {type(aliases).__name__}")
            continue
        if len(aliases) == 0:
            _add_error(errors, aliasKeyPath, "Alias list cannot be empty when enforceDefaultBackboneCharges is True.")
            continue
        if not all(isinstance(alias, str) for alias in aliases):
            _add_error(errors, aliasKeyPath, "All aliases must be strings.")


def _validate_torsion_scan_cores(config: Dict[str, Any], errors: Dict[str, str]) -> None:
    """Cross-section validation for `torsionScanInfo.nCoresPerCalculation`."""
    torsionScanInfo = config.get("torsionScanInfo")
    miscInfo = config.get("miscInfo")
    if not isinstance(torsionScanInfo, dict) or not isinstance(miscInfo, dict):
        return
    nCoresPerCalculation = torsionScanInfo.get("nCoresPerCalculation")
    availableCpus = miscInfo.get("availableCpus")
    if isinstance(nCoresPerCalculation, int) and isinstance(availableCpus, int) and nCoresPerCalculation > availableCpus:
        _add_error(errors, "torsionScanInfo.nCoresPerCalculation", "nCoresPerCalculation cannot exceed miscInfo.availableCpus.")


def validate_config(config: Dict[str, Any]) -> Union[Dict[str, Any], None]:
    """Validate the structure, types, and allowed values of the configuration dictionary."""
    errors: Dict[str, str] = {}

    if not isinstance(config, dict):
        _add_error(errors, "config", f"Configuration must be a dictionary, but got {type(config).__name__}")
        drSplash.show_config_error(errors)
        return None

    expectedSections = [
        "pathInfo",
        "moleculeInfo",
        "conformerGenerationInfo",
        "torsionScanInfo",
        "chargeFittingInfo",
        "parameterFittingInfo",
        "miscInfo",
    ]
    for sectionName in expectedSections:
        if sectionName not in config:
            _add_error(errors, sectionName, f"Missing required section '{sectionName}'")

    _validate_path_info(config.get("pathInfo"), "pathInfo", errors)
    _validate_molecule_info(config.get("moleculeInfo"), "moleculeInfo", errors)
    _validate_conformer_generation_info(config.get("conformerGenerationInfo"), "conformerGenerationInfo", errors)
    _validate_torsion_scan_info(config.get("torsionScanInfo"), "torsionScanInfo", errors)
    _validate_charge_fitting_info(config.get("chargeFittingInfo"), "chargeFittingInfo", errors)
    _validate_parameter_fitting_info(config.get("parameterFittingInfo"), "parameterFittingInfo", errors)
    _validate_misc_info(config.get("miscInfo"), "miscInfo", errors)
    _validate_backbone_charge_enforcement(config, errors)
    _validate_torsion_scan_cores(config, errors)

    if len(errors) == 0:
        return config
    drSplash.show_config_error(errors)
    return None
