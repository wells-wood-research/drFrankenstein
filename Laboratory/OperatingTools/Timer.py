import time
import csv
import functools
import os
from os import path as p
from datetime import datetime

from . import drYaml

def time_function(functionAlias: str, functionGroup: str | None = None):
    """Time a function call and store the results in config runtime info."""
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            # Support extracting the config either as a kwarg "config" or as the first positional arg
            if "config" in kwargs:
                config = kwargs["config"]
            elif len(args) > 0:
                config = args[0]
            else:
                raise KeyError("config")

            startTime = time.time()
            startCpuTimes = os.times()
            result = func(*args, **kwargs)  # Execute the function
            executionTime = time.time() - startTime
            endCpuTimes = os.times()
            cpuTime = (
                (endCpuTimes.user + endCpuTimes.system + endCpuTimes.children_user + endCpuTimes.children_system)
                - (startCpuTimes.user + startCpuTimes.system + startCpuTimes.children_user + startCpuTimes.children_system)
            )
            
            functionName = func.__name__

            timeInfo = config["runtimeInfo"].get("timeInfo", None)
            if not config["runtimeInfo"].get("timeInfo", None):
                timeInfo = {}
                config["runtimeInfo"]["timeInfo"] = timeInfo

            if not config["runtimeInfo"]["timeInfo"].get(functionName, None):
                config["runtimeInfo"]["timeInfo"][functionName] = {}
                config["runtimeInfo"]["timeInfo"][functionName]["functionAlias"] = functionAlias
                config["runtimeInfo"]["timeInfo"][functionName]["executionTime"] = executionTime
                config["runtimeInfo"]["timeInfo"][functionName]["cpuTime"] = cpuTime
                config["runtimeInfo"]["timeInfo"][functionName]["functionGroup"] = functionGroup
            else:
                config["runtimeInfo"]["timeInfo"][functionName]["executionTime"] += executionTime
                existingCpuTime = config["runtimeInfo"]["timeInfo"][functionName].get("cpuTime", 0.0)
                config["runtimeInfo"]["timeInfo"][functionName]["cpuTime"] = existingCpuTime + cpuTime

            return result  # Return the function's result unchanged
        return wrapper
    return decorator
