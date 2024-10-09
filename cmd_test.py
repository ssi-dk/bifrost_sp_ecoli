import pytest
import os
import subprocess
from bifrostlib import datahandling
from bifrost_sp_ecoli import launcher

@pytest.fixture
def test_db_connection():
    """Ensure there is a database connection and it's a test database."""
    # see https://github.com/ssi-dk/bifrost_whats_my_species/blob/master/tests/test_simple.py
    assert datahandling.has_a_database_connection(), "No database connection!"
    assert "TEST" in os.environ["BIFROST_DB_KEY"].upper(), "BIFROST_DB_KEY does not contain 'TEST'"

def test_launcher_info(test_db_connection):
    """
    Test using the launcher to run --info.
    """
    launcher.run_pipeline(["--info"])

def test_launcher_help(test_db_connection):
    """
    Test using the launcher to run --help.
    """
    launcher.run_pipeline(["--help"])

def test_bifrost_help_command(test_db_connection):
    """
    Test that the bifrost_sp_ecoli command runs and returns help info without errors.
    """
    bifrost_install_dir = os.environ.get("BIFROST_INSTALL_DIR")
    bifrost_db_key = os.environ.get("BIFROST_DB_KEY")
    
    # Ensure environment variables are set correctly
    assert bifrost_install_dir is not None, "BIFROST_INSTALL_DIR is not set!"
    assert bifrost_db_key is not None, "BIFROST_DB_KEY is not set!"

    # Run the bifrost_sp_ecoli help command
    cmd = f"BIFROST_DB_KEY={bifrost_db_key} python -m bifrost_sp_ecoli -h"
    result = subprocess.run(cmd, shell=True, capture_output=True)

    # Check if the command ran successfully
    assert result.returncode == 0, f"Command failed with error: {result.stderr.decode()}"
    assert "usage:" in result.stdout.decode(), "Help output not returned!"

# This allows pytest to run all test functions in this file
if __name__ == "__main__":
    pytest.main()