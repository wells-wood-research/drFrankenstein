import time
import csv
import functools
import os
from os import path as p
from datetime import datetime

def sort_output_directory(config):
    outputDir = config["pathInfo"]["outputDir"]
    timeDir = p.join(outputDir, "00_Timers")
    os.makedirs(timeDir, exist_ok=True)
    config["pathInfo"]["timeDir"] = timeDir
    return config



def time_function():
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            config = args[0]       ### THIS WILL ONLY WORK FOR ONE ARGUMENT THAT IS CONFIG
            timeDir = config["pathInfo"]["timeDir"]
            timeCsv = p.join(timeDir, "timeLog.csv")

            print(timeCsv)
            startTime = time.time()
            result = func(*args, **kwargs)  # Execute the function and preserve its return value
            executionTime = time.time() - startTime
            
            timeStamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            functionName = func.__name__
            fileExists = p.isfile(timeCsv)
                
            with open(timeCsv, 'a', newline='') as csvFile:
                fieldNames = ['timeStamp', 'functionName', 'executionTimeSeconds']
                writer = csv.DictWriter(csvFile, fieldnames=fieldNames)
                
                if not fileExists:
                    writer.writeheader()
                
                writer.writerow({
                    'timeStamp': timeStamp,
                    'functionName': functionName,
                    'executionTimeSeconds': executionTime
                })
            return result  # Return the result unchanged (e.g., a ruamel.yaml dict)
        return wrapper
    return decorator
