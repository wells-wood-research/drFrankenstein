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
    config["runtimeInfo"]["timeDir"] = timeDir
    return config



def time_function():
    """
    Simple decorator to time the execution of a function
    Writes results to a csv, location is found in config
    If functionName exists, adds to existing execution time
    """
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            config = kwargs["config"]  # Only works when args are explicitly defined
            timeDir = config["runtimeInfo"]["timeDir"]
            timeCsv = p.join(timeDir, "timeLog.csv")

            startTime = time.time()
            result = func(*args, **kwargs)  # Execute the function
            executionTime = time.time() - startTime
            
            timeStamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            functionName = func.__name__
            fileExists = p.isfile(timeCsv)
            
            # Read existing data if file exists
            existing_data = []
            if fileExists:
                with open(timeCsv, 'r', newline='') as csvFile:
                    reader = csv.DictReader(csvFile)
                    existing_data = list(reader)
            
            # Check if functionName exists and update execution time
            function_found = False
            for row in existing_data:
                if row['functionName'] == functionName:
                    function_found = True
                    # Add new execution time to existing
                    row['executionTimeSeconds'] = str(float(row['executionTimeSeconds']) + executionTime)
                    row['timeStamp'] = timeStamp  # Update timestamp
                    break
            
            # Write back all data
            with open(timeCsv, 'w', newline='') as csvFile:
                fieldNames = ['timeStamp', 'functionName', 'executionTimeSeconds']
                writer = csv.DictWriter(csvFile, fieldnames=fieldNames)
                
                writer.writeheader()
                
                if function_found:
                    # Write updated existing data
                    writer.writerows(existing_data)
                else:
                    # Write existing data plus new entry
                    writer.writerows(existing_data)
                    writer.writerow({
                        'timeStamp': timeStamp,
                        'functionName': functionName,
                        'executionTimeSeconds': executionTime
                    })
            
            return result  # Return the function's result unchanged
        return wrapper
    return decorator