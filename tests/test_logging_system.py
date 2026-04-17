"""
Tests for the drLogger logging system
"""
import os
import tempfile
import subprocess
from pathlib import Path
import sys

# Add Laboratory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'Laboratory'))

from OperatingTools import drLogger
from OperatingTools import drSubprocess


def test_experiment_logger_initialization():
    """Test that the logger initializes correctly"""
    with tempfile.TemporaryDirectory() as tmpdir:
        logger = drLogger.ExperimentLogger(tmpdir)
        
        # Check log directory was created
        assert os.path.isdir(tmpdir)
        assert os.path.isfile(logger.get_log_file_path())
        print("✓ Logger initialization test passed")


def test_experiment_start_end():
    """Test experiment start and end logging with timing"""
    with tempfile.TemporaryDirectory() as tmpdir:
        logger = drLogger.ExperimentLogger(tmpdir)
        
        exp_id = logger.start_experiment("Test Experiment")
        assert exp_id is not None
        assert "Test Experiment" in exp_id
        
        # Simulate some work
        import time
        time.sleep(0.5)
        
        duration = logger.end_experiment(exp_id, success=True)
        assert duration is not None
        assert "00:00:00" in duration  # At least seconds format
        
        # Check log file contains experiment info
        log_file = logger.get_log_file_path()
        with open(log_file, 'r') as f:
            content = f.read()
            assert "Test Experiment" in content
            assert "SUCCESS" in content
        
        print(f"✓ Experiment logging test passed (duration: {duration})")


def test_completion_tag():
    """Test that completion tag is written to log"""
    with tempfile.TemporaryDirectory() as tmpdir:
        logger = drLogger.ExperimentLogger(tmpdir)
        logger.log_completion_tag()
        
        log_file = logger.get_log_file_path()
        with open(log_file, 'r') as f:
            content = f.read()
            assert "***drFrankenstein_Parameterisation_Complete***" in content
        
        print("✓ Completion tag test passed")


def test_subprocess_logging():
    """Test subprocess call logging"""
    with tempfile.TemporaryDirectory() as tmpdir:
        logger = drLogger.ExperimentLogger(tmpdir)
        drLogger.set_logger(logger)  # Set global logger
        
        # Test successful command
        result = drSubprocess.logged_call(['echo', 'test'])
        assert result == 0
        
        log_file = logger.get_log_file_path()
        with open(log_file, 'r') as f:
            content = f.read()
            assert "SUBPROCESS CALL" in content
            assert "echo" in content
            assert "SUCCESS" in content
        
        print("✓ Subprocess logging test passed")


def test_subprocess_run_logging():
    """Test subprocess.run logging"""
    with tempfile.TemporaryDirectory() as tmpdir:
        logger = drLogger.ExperimentLogger(tmpdir)
        drLogger.set_logger(logger)  # Set global logger
        
        # Test successful command
        result = drSubprocess.logged_run(['ls', '-l'], capture_output=True)
        assert result.returncode == 0
        
        log_file = logger.get_log_file_path()
        with open(log_file, 'r') as f:
            content = f.read()
            assert "SUBPROCESS CALL" in content
            assert "ls" in content
        
        print("✓ Subprocess.run logging test passed")


def test_duration_formatting():
    """Test duration formatting in DD:HH:MM:SS format"""
    with tempfile.TemporaryDirectory() as tmpdir:
        logger = drLogger.ExperimentLogger(tmpdir)
        
        # Test various durations
        test_cases = [
            (0, "00:00:00:00"),
            (60, "00:00:01:00"),
            (3661, "00:01:01:01"),
            (86400, "01:00:00:00"),
            (90061, "01:01:01:01"),
        ]
        
        for seconds, expected in test_cases:
            formatted = logger._format_duration(seconds)
            assert formatted == expected, f"Expected {expected}, got {formatted}"
        
        print("✓ Duration formatting test passed")


def test_experiment_without_logger():
    """Test that decorator handles missing logger gracefully"""
    from OperatingTools.drLogger import experiment_logger
    
    @experiment_logger("Test Protocol")
    def test_protocol(config):
        return config
    
    # Should work even without logger in config
    config = {}
    result = test_protocol(config=config)
    assert result == config
    
    print("✓ Missing logger gracefully handled test passed")


def test_error_report_logging():
    """Test error report logging in text-friendly format"""
    with tempfile.TemporaryDirectory() as tmpdir:
        logger = drLogger.ExperimentLogger(tmpdir)
        
        # Create a mock error report (similar to drFrankenstein.handle_exceptions output)
        error_report = {
            "pdbName": "test_molecule",
            "errorType": "ValueError",
            "errorMessage": "Invalid parameter value provided",
            "functionName": "calculate_charges",
            "lineNumber": "125",
            "lineOfCode": "result = invalid_function()",
            "scriptName": "/home/user/drFrankenstein/charges.py",
            "fullTraceBack": [
                "/home/user/drFrankenstein/main.py:42 in main",
                "/home/user/drFrankenstein/charges.py:125 in calculate_charges",
                "/home/user/drFrankenstein/utils.py:50 in helper_function",
            ]
        }
        
        logger.log_error_report(error_report)
        
        log_file = logger.get_log_file_path()
        with open(log_file, 'r') as f:
            content = f.read()
            # Verify key error information is in the log
            assert "DRFRANKENSTEIN ERROR REPORT" in content
            assert "test_molecule" in content
            assert "ValueError" in content
            assert "Invalid parameter value provided" in content
            assert "calculate_charges" in content
            assert "main.py" in content
            # Verify no ANSI color codes
            assert "\033[" not in content
        
        print("✓ Error report logging test passed")


if __name__ == '__main__':
    print("Running drLogger tests...")
    print()
    
    test_experiment_logger_initialization()
    test_experiment_start_end()
    test_completion_tag()
    test_subprocess_logging()
    test_subprocess_run_logging()
    test_duration_formatting()
    test_experiment_without_logger()
    test_error_report_logging()
    
    print()
    print("✓ All logging tests passed!")
