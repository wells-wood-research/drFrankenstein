"""
Subprocess wrappers that integrate with drLogger for command tracking.
These functions wrap subprocess.call and subprocess.run to automatically log
command execution and results. Uses global logger instance to avoid
polluting the config dict (which gets written to YAML).
"""

import subprocess
from typing import List, Optional


def logged_call(command: List[str], **kwargs) -> int:
    """
    Wrapper for subprocess.call that logs the command and result.
    
    Args:
        command: List of command and arguments
        **kwargs: Additional arguments to pass to subprocess.call
        
    Returns:
        The return code from subprocess.call
    """
    from . import drLogger
    logger = drLogger.get_logger()
    
    try:
        result = subprocess.call(command, **kwargs)
        if logger:
            logger.log_subprocess_call(command, result)
        return result
    except Exception as e:
        if logger:
            logger.log_subprocess_error(command, e)
        raise


def logged_run(command: List[str], **kwargs) -> subprocess.CompletedProcess:
    """
    Wrapper for subprocess.run that logs the command and result.
    
    Args:
        command: List of command and arguments
        **kwargs: Additional arguments to pass to subprocess.run
        
    Returns:
        The CompletedProcess object from subprocess.run
    """
    from . import drLogger
    logger = drLogger.get_logger()
    
    try:
        result = subprocess.run(command, **kwargs)
        if logger:
            logger.log_subprocess_call(command, result.returncode)
        return result
    except Exception as e:
        if logger:
            logger.log_subprocess_error(command, e)
        raise
