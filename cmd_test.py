import pytest
import os
import subprocess
from bifrostlib import datahandling
from bifrostlib import database_interface
from bifrost_sp_ecoli import launcher
import pymongo
import shutil

def connect_to_database():
    """Connect to the MongoDB database using the BIFROST_DB_KEY environment variable."""
    db_key = os.getenv("BIFROST_DB_KEY")

    if not db_key:
        print("BIFROST_DB_KEY is not set!")  # the DB key is missing
        raise ValueError("BIFROST_DB_KEY is not set!")  # stop execution


    try:
        # create a database instance
        client = pymongo.MongoClient(db_key)

        # Attempt to get server info to verify the connection
        server_info = client.server_info()  # Raises an exception if the connection fails
        print("MongoDB client created successfully!")
        print("Server Info:", server_info)  # Print the server information
    except Exception as e:
        print(f"Failed to create MongoDB client: {e}")  # Print a detailed error message
        raise  # Re-raise the exception to indicate failure to the caller
    
    return client

@pytest.fixture
def test_db_connection():
    """Ensure there is a database connection and it's a test database."""
    client = None
    try:
        client = connect_to_database()
        # Check if we can fetch server info to confirm the connection
        client.server_info()  
        print("Database connection successful!")
        assert "TEST" in os.environ["BIFROST_DB_KEY"].upper(), "BIFROST_DB_KEY does not contain 'TEST'"

    except Exception as e:
        assert False, f"Failed to connect to the database: {e}"

    finally:
        if client:
            client.close()  

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

# This allows pytest to run all test functions in this file
if __name__ == "__main__":
    pytest.main()