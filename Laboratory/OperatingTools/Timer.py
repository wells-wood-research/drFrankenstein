import time
import csv
import functools
import os
from os import path as p
from datetime import datetime

from . import drYaml
from typing import Callable, Any # Added for more specific type hints

def time_function(function_alias: str, function_group: str = None) -> Callable:
    """
    Simple decorator to time the execution of a function.
    Writes results to the config dictionary.
    If functionName exists in config's timeInfo, adds to existing execution time.

    Args:
        function_alias: An alias for the function, used for reporting.
        function_group: An optional group name for the function, used for reporting.

    Returns:
        The decorator function.
    """
    def decorator(func: Callable) -> Callable:
        @functools.wraps(func)
        def wrapper(*args: Any, **kwargs: Any) -> Any:
            # Assuming 'config' is always passed as a keyword argument to the wrapped function
            if "config" not in kwargs:
                # Or handle error: raise ValueError("Config not provided to timed function")
                # Alternatively, if 'config' can be an arg:
                # config = args[0] if args and isinstance(args[0], dict) else kwargs.get("config")
                # This part depends on how 'config' is expected to be passed.
                # Forcing it as kwarg as per original comment.
                raise ValueError("The 'config' dictionary must be passed as a keyword argument to the timed function.")

            config = kwargs["config"]
            startTime = time.time()
            result = func(*args, **kwargs)  # Execute the function
            executionTime = time.time() - startTime
            
            functionName = func.__name__

            timeInfo = config["runtimeInfo"].get("timeInfo", None)
            if timeInfo is None: # Simplified check
                timeInfo = {}
                config["runtimeInfo"]["timeInfo"] = timeInfo

            if functionName not in timeInfo: # Simplified check
                config["runtimeInfo"]["timeInfo"][functionName] = {}
                config["runtimeInfo"]["timeInfo"][functionName]["functionAlias"] = function_alias
                config["runtimeInfo"]["timeInfo"][functionName]["executionTime"] = executionTime
                config["runtimeInfo"]["timeInfo"][functionName]["functionGroup"] = function_group
            else:
                # Ensure "executionTime" exists before trying to add to it
                currentExecutionTime = config["runtimeInfo"]["timeInfo"][functionName].get("executionTime", 0.0)
                totalExecutionTime = currentExecutionTime + executionTime
                config["runtimeInfo"]["timeInfo"][functionName]["executionTime"] = totalExecutionTime
                # Update alias and group if they are not set, or if you want to allow overriding
                if "functionAlias" not in config["runtimeInfo"]["timeInfo"][functionName]:
                     config["runtimeInfo"]["timeInfo"][functionName]["functionAlias"] = function_alias
                if "functionGroup" not in config["runtimeInfo"]["timeInfo"][functionName] and function_group is not None:
                     config["runtimeInfo"]["timeInfo"][functionName]["functionGroup"] = function_group


            return result  # Return the function's result unchanged
        return wrapper
    return decorator