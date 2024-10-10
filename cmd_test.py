import pytest
import os
import subprocess
from bifrostlib import datahandling
from bifrostlib import database_interface
from bifrost_sp_ecoli import launcher
import pymongo
import shutil

#@pytest.fixture

def test_mongo_connection():
    # Get the MongoDB connection string from the environment variable
    
    client = pymongo.MongoClient(os.environ["BIFROST_DB_KEY"])

    print("Attempt to connect to the database")
    ping = client.admin.command('ping')
    print("Successfully connected to the MongoDB database.")

    # Clean up: close the client
    client.close()

if __name__ == "__main__":
    # Run the test
    test_mongo_connection()