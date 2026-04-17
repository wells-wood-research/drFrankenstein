# drFrankenstein Logging System

## Overview

A robust logging system has been implemented for drFrankenstein to track experiment execution, subprocess calls, and runtime metrics. All logs are written to a timestamped log file in the output directory.

## Features

### 1. **Experiment Tracking**
- **Start/End Logging**: Each protocol experiment is logged with start and end timestamps
- **Duration Tracking**: Automatic calculation and logging of execution time in `DD:HH:MM:SS` format
- **Success/Failure Status**: Track whether each experiment completed successfully or failed

### 2. **Subprocess Call Logging**
- **Command Tracking**: Every subprocess call is logged with the full command string
- **Success Indicator**: Return code is logged for each subprocess invocation
- **Error Handling**: Failed subprocess calls are logged with error details

### 3. **Completion Tag**
- **Easy Search**: Upon successful completion, a searchable tag is written to the log:
  ```
  ***drFrankenstein_Parameterisation_Complete***
  ```

## Log File Location

Logs are stored in:
```
<outputDir>/logs/drFrankenstein.log
```

where `<outputDir>` is specified in your drFrankenstein configuration file.

## Log File Format

Each log entry follows this format:
```
YYYY-MM-DD HH:MM:SS - [LEVEL] - message
```

### Example Log Output - Success

```
2026-04-17 13:07:56 - [INFO] - drFrankenstein Parameterisation started
2026-04-17 13:07:58 - [INFO] - START EXPERIMENT: Termini Capping
2026-04-17 13:08:15 - [INFO] - SUBPROCESS CALL: /path/to/orca input.inp | Status: SUCCESS (code: 0)
2026-04-17 13:08:45 - [INFO] - END EXPERIMENT: Termini Capping | Status: SUCCESS | Duration: 00:00:47:12
2026-04-17 13:08:47 - [INFO] - START EXPERIMENT: Conformer Generation
...
2026-04-17 15:30:22 - [INFO] - END EXPERIMENT: Report Generation | Status: SUCCESS | Duration: 00:02:05:33
2026-04-17 15:30:23 - [INFO] - ***drFrankenstein_Parameterisation_Complete***
```

### Example Log Output - Error

When an error occurs, a detailed traceback is logged in **text-friendly format** (no ANSI colors):

```
2026-04-17 14:21:24 - [ERROR] - ================================================================================
2026-04-17 14:21:24 - [ERROR] - DRFRANKENSTEIN ERROR REPORT
2026-04-17 14:21:24 - [ERROR] - ================================================================================
2026-04-17 14:21:24 - [ERROR] - 
For System: ALA_molecule
2026-04-17 14:21:24 - [ERROR] - 
In Script: .../Protocol_5_Twisting/Twisted_Monster.py
2026-04-17 14:21:24 - [ERROR] - In Function: calculate_torsion_angles
2026-04-17 14:21:24 - [ERROR] - At Line Number: 187
2026-04-17 14:21:24 - [ERROR] - 
Error Type: ValueError
2026-04-17 14:21:24 - [ERROR] - Error Message: Torsion angle out of valid range: 450 degrees
2026-04-17 14:21:24 - [ERROR] - Line of Code: if angle > 360: raise ValueError(...)
2026-04-17 14:21:24 - [ERROR] - 
2026-04-17 14:21:24 - [ERROR] - FULL DEBUG TRACEBACK
2026-04-17 14:21:24 - [ERROR] - LINE NUMBER    FUNCTION                           FILE PATH
2026-04-17 14:21:24 - [ERROR] - 105            main                               .../drFrankenstein.py
2026-04-17 14:21:24 - [ERROR] - 42             twist_protocol                     .../Twisted_Doctor.py
2026-04-17 14:21:24 - [ERROR] - 187            calculate_torsion_angles           .../Twisted_Monster.py
2026-04-17 14:21:24 - [ERROR] - 150            validate_angle                     .../Twisted_Monster.py
2026-04-17 14:21:24 - [ERROR] - ================================================================================
```

## Tracked Experiments

The following protocol experiments are automatically logged:

1. **Termini Capping** - Protocol 1
2. **Conformer Generation** - Protocol 2
3. **Charge Calculation** - Protocol 3
4. **Parameter Assembly** - Protocol 4
5. **Torsion Scanning** - Protocol 5
6. **Torsion Parameter Fitting** - Protocol 6
7. **Final Creation** - Protocol 7
8. **Report Generation** - Protocol 8

## Using the Logging System

### Important: Global Logger (No Config Pollution)

The logger is now a **global singleton** to avoid YAML serialization issues. Do NOT store it in the config dict.

### In Experiment Code

The logging system is automatically integrated into all protocol functions via decorators:

```python
from OperatingTools import drLogger

@drLogger.experiment_logger("My Experiment Name")
def my_protocol(config: dict) -> dict:
    # Your protocol code here
    # No need to access logger - it's automatic!
    return config
```

The decorator automatically uses the global logger. No config modification needed.

### In Subprocess Calls

Use the provided logging-aware subprocess wrappers:

```python
from OperatingTools import drSubprocess

# Instead of: subprocess.call(cmd)
drSubprocess.logged_call(cmd)

# Instead of: subprocess.run(cmd)
drSubprocess.logged_run(cmd)
```

Or log manually:

```python
from OperatingTools import drLogger

result = subprocess.call(cmd)
logger = drLogger.get_logger()
if logger:
    logger.log_subprocess_call(cmd, result)
```

### Accessing the Logger Globally

The logger is accessible from anywhere without polluting the config:

```python
from OperatingTools import drLogger

logger = drLogger.get_logger()
if logger:
    exp_id = logger.start_experiment("Custom Experiment")
    # ... do work ...
    logger.end_experiment(exp_id, success=True)
```

## API Reference

### ExperimentLogger Class

#### `__init__(log_dir: str)`
Initialize the logger with a log directory.

#### `start_experiment(experiment_name: str) -> str`
Log the start of an experiment. Returns a unique experiment ID for tracking.

#### `end_experiment(experiment_id: str, success: bool = True) -> str`
Log the end of an experiment. Returns formatted duration in `DD:HH:MM:SS` format.

#### `log_subprocess_call(command: List[str], return_code: int, output: str = None) -> None`
Log a subprocess call with its return code.

#### `log_subprocess_error(command: List[str], exception: Exception, output: str = None) -> None`
Log a subprocess error with exception details.

#### `log_completion_tag() -> None`
Write the completion tag to the log file.

#### `log_error_report(error_report: dict) -> None`
Log a detailed error report in text-friendly format (no ANSI colors). Takes an error report dict with keys: pdbName, errorType, errorMessage, functionName, lineNumber, lineOfCode, scriptName, fullTraceBack.

#### `get_log_file_path() -> str`
Return the path to the log file.

#### `_format_duration(seconds: int) -> str`
Format seconds into DD:HH:MM:SS format.

### Module-Level Functions

#### `get_logger() -> ExperimentLogger | None`
Get the global logger instance. Returns None if not initialized.

#### `set_logger(logger: ExperimentLogger) -> None`
Set the global logger instance. Called automatically by main() to initialize logging.

### Decorators

#### `@experiment_logger(experiment_name: str)`
Decorator for automatic experiment start/end logging with timing. Uses the global logger instance automatically.

### Utilities

#### `drSubprocess.logged_call(command, **kwargs) -> int`
Logging-aware wrapper for `subprocess.call()`. Uses global logger instance.

#### `drSubprocess.logged_run(command, **kwargs) -> CompletedProcess`
Logging-aware wrapper for `subprocess.run()`. Uses global logger instance.

## Duration Format

Durations are reported in `DD:HH:MM:SS` format:
- **DD** = Days (00-99+)
- **HH** = Hours (00-23)
- **MM** = Minutes (00-59)
- **SS** = Seconds (00-59)

Examples:
- `00:00:00:30` = 30 seconds
- `00:00:05:45` = 5 minutes 45 seconds
- `00:02:30:15` = 2 hours 30 minutes 15 seconds
- `01:06:15:30` = 1 day 6 hours 15 minutes 30 seconds

## Checking Logs

### View the entire log:
```bash
cat <outputDir>/logs/drFrankenstein.log
```

### Search for a specific experiment:
```bash
grep "Termini Capping" <outputDir>/logs/drFrankenstein.log
```

### Find the completion tag:
```bash
grep "drFrankenstein_Parameterisation_Complete" <outputDir>/logs/drFrankenstein.log
```

### Get timing information:
```bash
grep "END EXPERIMENT" <outputDir>/logs/drFrankenstein.log
```

### Check for errors:
```bash
grep "ERROR\|FAILED" <outputDir>/logs/drFrankenstein.log
```

## Integration with Existing Code

### Already Integrated

All main protocol functions in the following files have been instrumented with logging:

- `Laboratory/drFrankenstein.py` (main entry point)
- `Laboratory/Experiments/Protocol_1_Capping/Capping_Doctor.py`
- `Laboratory/Experiments/Protocol_2_Wriggling/Wriggling_Doctor.py`
- `Laboratory/Experiments/Protocol_3_Charging/Charged_Doctor.py`
- `Laboratory/Experiments/Protocol_4_Assembly/Assembly_Doctor.py`
- `Laboratory/Experiments/Protocol_5_Twisting/Twisted_Doctor.py`
- `Laboratory/Experiments/Protocol_6_Stitching/Stitching_Doctor.py`
- `Laboratory/Experiments/Protocol_7_Creation/drCreator.py`
- `Laboratory/Experiments/Protocol_8_Reporter/Reporting_Doctor.py`
- `Laboratory/OperatingTools/drOrca.py` (subprocess calls)

### Extending Logging

To add logging to additional functions:

1. Import the logging decorator:
   ```python
   from OperatingTools import drLogger
   ```

2. Add the decorator to your function:
   ```python
   @drLogger.experiment_logger("Your Experiment Name")
   def your_protocol(config):
       # code here - no need to touch config
       pass
   ```

3. For subprocess calls:
   ```python
   from OperatingTools import drSubprocess
   result = drSubprocess.logged_call(cmd)  # No config needed!
   ```

4. To access logger directly (optional):
   ```python
   from OperatingTools import drLogger
   logger = drLogger.get_logger()
   if logger:
       # logger available
       pass
   ```

## Testing

A comprehensive test suite is available at:
```
tests/test_logging_system.py
```

Run tests with:
```bash
cd /home/esp/scriptDevelopment/drFrankenstein
python -m pytest tests/test_logging_system.py -v
```

## Files Modified

### New Files
- `Laboratory/OperatingTools/drLogger.py` - Main logging module
- `Laboratory/OperatingTools/drSubprocess.py` - Subprocess logging utilities
- `tests/test_logging_system.py` - Logging system tests
- `Documentation/LOGGING_SYSTEM.md` - This documentation

### Modified Files
- `Laboratory/drFrankenstein.py` - Initialize logger, add completion tag
- `Laboratory/Experiments/Protocol_1_Capping/Capping_Doctor.py` - Added logging decorator
- `Laboratory/Experiments/Protocol_2_Wriggling/Wriggling_Doctor.py` - Added logging decorator
- `Laboratory/Experiments/Protocol_3_Charging/Charged_Doctor.py` - Added logging decorator
- `Laboratory/Experiments/Protocol_4_Assembly/Assembly_Doctor.py` - Added logging decorator
- `Laboratory/Experiments/Protocol_5_Twisting/Twisted_Doctor.py` - Added logging decorator
- `Laboratory/Experiments/Protocol_6_Stitching/Stitching_Doctor.py` - Added logging decorator
- `Laboratory/Experiments/Protocol_7_Creation/drCreator.py` - Added logging decorator
- `Laboratory/Experiments/Protocol_8_Reporter/Reporting_Doctor.py` - Added logging decorator
- `Laboratory/OperatingTools/drOrca.py` - Added subprocess call logging
