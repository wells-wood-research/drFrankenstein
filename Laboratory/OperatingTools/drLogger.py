import os
import logging
import functools
import time
import subprocess
from datetime import datetime, timedelta
from typing import Optional, List, Any
from os import path as p


class ExperimentLogger:
    """
    Robust logging system for drFrankenstein experiments.
    Tracks experiment lifecycle, subprocess calls, and provides formatted output.
    """

    def __init__(self, log_dir: str):
        """
        Initialize the experiment logger.
        
        Args:
            log_dir: Directory where log files will be stored
        """
        self.log_dir = log_dir
        os.makedirs(log_dir, exist_ok=True)
        
        self.log_file = p.join(log_dir, "drFrankenstein.log")
        self.experiments = {}  # Track active experiments
        self._setup_logger()

    def _setup_logger(self):
        """Configure the logger with both file and console handlers."""
        self.logger = logging.getLogger("drFrankenstein")
        self.logger.setLevel(logging.DEBUG)
        
        # Remove existing handlers to avoid duplicates
        self.logger.handlers = []
        
        # File handler
        file_handler = logging.FileHandler(self.log_file)
        file_handler.setLevel(logging.DEBUG)
        
        # Console handler
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.INFO)
        
        # Formatter
        formatter = logging.Formatter(
            '%(asctime)s - [%(levelname)s] - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        
        file_handler.setFormatter(formatter)
        console_handler.setFormatter(formatter)
        
        self.logger.addHandler(file_handler)
        self.logger.addHandler(console_handler)

    def start_experiment(self, experiment_name: str) -> str:
        """
        Log the start of an experiment.
        
        Args:
            experiment_name: Name of the experiment (e.g., "Termini Capping", "Charge Calculation")
            
        Returns:
            experiment_id: Unique identifier for tracking this experiment instance
        """
        experiment_id = f"{experiment_name}_{datetime.now().strftime('%Y%m%d_%H%M%S_%f')}"
        start_time = time.time()
        
        self.experiments[experiment_id] = {
            'name': experiment_name,
            'start_time': start_time,
            'start_datetime': datetime.now(),
            'end_time': None,
            'status': 'running'
        }
        
        self.logger.info(f"START EXPERIMENT: {experiment_name}")
        self.logger.debug(f"Experiment ID: {experiment_id}")
        return experiment_id

    def end_experiment(self, experiment_id: str, success: bool = True) -> Optional[str]:
        """
        Log the end of an experiment with duration.
        
        Args:
            experiment_id: The experiment_id returned from start_experiment()
            success: Whether the experiment completed successfully
            
        Returns:
            formatted_duration: Duration in DD:HH:MM:SS format
        """
        if experiment_id not in self.experiments:
            self.logger.warning(f"Unknown experiment ID: {experiment_id}")
            return None
        
        exp = self.experiments[experiment_id]
        end_time = time.time()
        exp['end_time'] = end_time
        
        duration_seconds = int(end_time - exp['start_time'])
        formatted_duration = self._format_duration(duration_seconds)
        
        status_str = "SUCCESS" if success else "FAILED"
        exp['status'] = 'success' if success else 'failed'
        
        self.logger.info(
            f"END EXPERIMENT: {exp['name']} | "
            f"Status: {status_str} | Duration: {formatted_duration}"
        )
        
        return formatted_duration

    def log_subprocess_call(self, 
                          command: List[str], 
                          return_code: int,
                          output: Optional[str] = None) -> None:
        """
        Log a subprocess call with command and result.
        
        Args:
            command: The command list that was executed
            return_code: The return code from the subprocess
            output: Optional stdout/stderr output from the command
        """
        command_str = ' '.join(str(c) for c in command)
        status_str = "SUCCESS" if return_code == 0 else "FAILED"
        
        self.logger.info(f"SUBPROCESS CALL: {command_str} | Status: {status_str} (code: {return_code})")
        
        if return_code != 0 and output:
            self.logger.debug(f"Subprocess output:\n{output}")

    def log_subprocess_error(self,
                           command: List[str],
                           exception: Exception,
                           output: Optional[str] = None) -> None:
        """
        Log a subprocess error.
        
        Args:
            command: The command list that was attempted
            exception: The exception that was raised
            output: Optional stdout/stderr output
        """
        command_str = ' '.join(str(c) for c in command)
        self.logger.error(f"SUBPROCESS ERROR: {command_str} | Exception: {str(exception)}")
        
        if output:
            self.logger.debug(f"Subprocess output:\n{output}")

    def log_completion_tag(self) -> None:
        """Write the completion tag to the log file."""
        self.logger.info("***drFrankenstein_Parameterisation_Complete***")

    def log_error_report(self, error_report: dict) -> None:
        """
        Log a detailed error report in text-friendly format (no ANSI colors).
        
        Args:
            error_report: Dictionary with error details (as created by drFrankenstein.handle_exceptions)
                         Should contain: pdbName, errorType, errorMessage, functionName, 
                         lineNumber, lineOfCode, scriptName, fullTraceBack
        """
        if not error_report:
            return
        
        # Header
        self.logger.error("=" * 80)
        self.logger.error("DRFRANKENSTEIN ERROR REPORT")
        self.logger.error("=" * 80)
        
        # Summary info
        self.logger.error(f"\nFor System: {error_report.get('pdbName', 'UNKNOWN')}")
        self.logger.error(f"\nIn Script: {error_report.get('scriptName', 'UNKNOWN')}")
        self.logger.error(f"In Function: {error_report.get('functionName', 'UNKNOWN')}")
        self.logger.error(f"At Line Number: {error_report.get('lineNumber', 'UNKNOWN')}")
        
        # Error details
        self.logger.error(f"\nError Type: {error_report.get('errorType', 'UNKNOWN')}")
        self.logger.error(f"Error Message: {error_report.get('errorMessage', 'UNKNOWN')}")
        self.logger.error(f"Line of Code: {error_report.get('lineOfCode', 'UNKNOWN')}")
        
        # Full traceback
        traceback_list = error_report.get('fullTraceBack', [])
        if traceback_list:
            self.logger.error("\n" + "-" * 80)
            self.logger.error("FULL DEBUG TRACEBACK")
            self.logger.error("-" * 80)
            self.logger.error(f"{'LINE NUMBER':<15}{'FUNCTION':<35}FILE PATH")
            self.logger.error(f"{'---':<15}{'---':<35}---")
            
            for traceback_line in traceback_list:
                # Parse the traceback line format: "/path/to/file.py:LINE in FUNCTION"
                try:
                    parts = traceback_line.split(":")
                    script_path = parts[0]
                    line_and_func = ":".join(parts[1:])  # In case path has colons
                    
                    # Extract line number and function name
                    if " in " in line_and_func:
                        line_num = line_and_func.split(" in ")[0].strip()
                        func_name = line_and_func.split(" in ")[1].strip()
                    else:
                        line_num = line_and_func.strip()
                        func_name = "UNKNOWN"
                    
                    self.logger.error(f"{line_num:<15}{func_name:<35}{script_path}")
                except Exception as e:
                    # If parsing fails, just log the raw line
                    self.logger.error(f"  {traceback_line}")
        
        self.logger.error("=" * 80)

    def _format_duration(self, seconds: int) -> str:
        """
        Format duration in DD:HH:MM:SS format.
        
        Args:
            seconds: Total seconds
            
        Returns:
            Formatted duration string
        """
        td = timedelta(seconds=seconds)
        days = td.days
        hours, remainder = divmod(td.seconds, 3600)
        minutes, secs = divmod(remainder, 60)
        return f"{days:02d}:{hours:02d}:{minutes:02d}:{secs:02d}"

    def get_log_file_path(self) -> str:
        """Return the path to the log file."""
        return self.log_file


def experiment_logger(experiment_name: str):
    """
    Decorator to automatically log experiment start/end with timing.
    
    Usage:
        @experiment_logger("My Experiment Name")
        def my_protocol(config):
            ...
    
    Note: Requires 'config' to be passed as a keyword argument.
          Logger is accessed via config['logger'] (or can be created if needed).
    """
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            config = kwargs.get("config")
            if not config:
                # If config not found, just run the function normally
                return func(*args, **kwargs)
            
            logger = config.get("logger")
            if not logger:
                # Logger not initialized yet
                return func(*args, **kwargs)
            
            exp_id = logger.start_experiment(experiment_name)
            try:
                result = func(*args, **kwargs)
                logger.end_experiment(exp_id, success=True)
                return result
            except Exception as e:
                logger.end_experiment(exp_id, success=False)
                raise
        
        return wrapper
    return decorator


def subprocess_wrapper(logger: 'ExperimentLogger'):
    """
    Wrapper factory for subprocess.call and subprocess.run that logs all invocations.
    
    Usage:
        from OperatingTools import drLogger
        logger = drLogger.ExperimentLogger(log_dir)
        safe_call = drLogger.subprocess_wrapper(logger)(subprocess.call)
        safe_call(['ls', '-la'])
    
    Returns:
        A wrapped subprocess function that logs calls and returns
    """
    def wrapper(original_func):
        @functools.wraps(original_func)
        def wrapped_subprocess(*args, **kwargs):
            # Extract command from args
            if args:
                command = args[0]
            else:
                command = kwargs.get('args', kwargs.get('cmd', 'unknown'))
            
            # Convert command to list if it's a string
            if isinstance(command, str):
                command = [command]
            
            try:
                result = original_func(*args, **kwargs)
                
                # Handle different return types
                if hasattr(result, 'returncode'):
                    return_code = result.returncode
                else:
                    return_code = result
                
                logger.log_subprocess_call(command, return_code)
                return result
            
            except Exception as e:
                logger.log_subprocess_error(command, e)
                raise
        
        return wrapped_subprocess
    
    return wrapper
