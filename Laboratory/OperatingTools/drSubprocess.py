"""
Subprocess wrappers that integrate with drLogger for command tracking.
These functions wrap subprocess.call and subprocess.run to automatically log
command execution and results. Uses global logger instance to avoid
polluting the config dict (which gets written to YAML).
"""

import subprocess
from typing import Any, List


def logged_call(command: List[str], **kwargs: Any) -> int:
    """Call a subprocess command and log the result if logging is enabled."""
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


def logged_run(command: List[str], **kwargs: Any) -> subprocess.CompletedProcess:
    """Run a subprocess command and log the result if logging is enabled."""
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
