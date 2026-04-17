"""
Subprocess wrappers that integrate with drLogger for command tracking.
These functions wrap subprocess.call and subprocess.run to automatically log
command execution and results.
"""

import subprocess
from typing import List, Optional, Any


def logged_call(command: List[str], 
                config: Optional[dict] = None,
                **kwargs) -> int:
    """
    Wrapper for subprocess.call that logs the command and result.
    
    Args:
        command: List of command and arguments
        config: Optional config dict containing logger
        **kwargs: Additional arguments to pass to subprocess.call
        
    Returns:
        The return code from subprocess.call
    """
    logger = None
    if config:
        logger = config.get("logger")
    
    try:
        result = subprocess.call(command, **kwargs)
        if logger:
            logger.log_subprocess_call(command, result)
        return result
    except Exception as e:
        if logger:
            logger.log_subprocess_error(command, e)
        raise


def logged_run(command: List[str],
               config: Optional[dict] = None,
               **kwargs) -> subprocess.CompletedProcess:
    """
    Wrapper for subprocess.run that logs the command and result.
    
    Args:
        command: List of command and arguments
        config: Optional config dict containing logger
        **kwargs: Additional arguments to pass to subprocess.run
        
    Returns:
        The CompletedProcess object from subprocess.run
    """
    logger = None
    if config:
        logger = config.get("logger")
    
    try:
        result = subprocess.run(command, **kwargs)
        if logger:
            logger.log_subprocess_call(command, result.returncode)
        return result
    except Exception as e:
        if logger:
            logger.log_subprocess_error(command, e)
        raise
