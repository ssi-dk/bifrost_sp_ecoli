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
    mongo_connection_string = os.environ["BIFROST_DB_KEY"]

    if not mongo_connection_string:
        print("No MongoDB connection string found in environment variables.")
        return False

    try:
        # Create a MongoDB client
        client = pymongo.MongoClient(mongo_connection_string)

        print("Attempt to connect to the database")
        ping = client.admin.command('ping')
        print("Successfully connected to the MongoDB database.")
        
        except ConnectionFailure:
        print("Server not available")
        except pymongo.errors.PyMongoError as e:
        print("Pymongo error: " + str(e))
        except Exception as exc:
        print("Exception: " + str(exc))

        # Clean up: close the client
        client.close()
        return True

    except pymongo.errors.ServerSelectionTimeoutError as e:
        print("Failed to connect to the MongoDB database:", e)
        return False

if __name__ == "__main__":
    # Run the test
    test_mongo_connection()