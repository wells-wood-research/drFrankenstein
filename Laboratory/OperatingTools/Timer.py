import time
import csv
import functools
import os
from os import path as p
from datetime import datetime

from . import drYaml

def time_function(functionAlias: str, functionGroup: str = None):
    """
    Simple decorator to time the execution of a function
    Writes results to a csv, location is found in config
    If functionName exists, adds to existing execution time
    """
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            config = kwargs["config"]  # Only works when args are explicitly defined as a kwarg
            startTime = time.time()
            result = func(*args, **kwargs)  # Execute the function
            executionTime = time.time() - startTime
            
            functionName = func.__name__

            timeInfo = config["runtimeInfo"].get("timeInfo", None)
            if not config["runtimeInfo"].get("timeInfo", None):
                timeInfo = {}
                config["runtimeInfo"]["timeInfo"] = timeInfo

            if not config["runtimeInfo"]["timeInfo"].get(functionName, None):
                config["runtimeInfo"]["timeInfo"][functionName] = {}
                config["runtimeInfo"]["timeInfo"][functionName]["functionAlias"] = functionAlias
                config["runtimeInfo"]["timeInfo"][functionName]["executionTime"] = executionTime
                config["runtimeInfo"]["timeInfo"][functionName]["functionGroup"] = functionGroup
            else:
                totalExecutionTime =    config["runtimeInfo"]["timeInfo"][functionName]["executionTime"] + executionTime
                config["runtimeInfo"]["timeInfo"][functionName]["executionTime"] += totalExecutionTime

            return result  # Return the function's result unchanged
        return wrapper
    return decorator