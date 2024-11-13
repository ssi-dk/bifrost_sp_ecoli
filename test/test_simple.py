import pytest
import os
import pymongo

# Test to check database connection using the BIFROST_DB_KEY
@pytest.fixture
def test_connection():
    # Check if BIFROST_DB_KEY is set in environment
    db_url = os.getenv("BIFROST_DB_KEY")
    assert db_url is not None, "Database connection string (BIFROST_DB_KEY) is not set."
    
    # Try to connect to the database
    client = pymongo.MongoClient(db_url)
    try:
        client.server_info()  # This triggers a connection and checks server status
    except pymongo.errors.PyMongoError as e:
        pytest.fail(f"Database connection failed: {e}")
    finally:
        client.close()

# Test to check if the BIFROST_INSTALL_DIR is set correctly
@pytest.fixture
def test_install_dir():
    # Ensure that the install directory is set
    install_dir = os.getenv("BIFROST_INSTALL_DIR")
    assert install_dir is not None, "BIFROST_INSTALL_DIR is not set."
    assert os.path.isdir(install_dir), f"The install directory '{install_dir}' does not exist."

# A simple test class to run these tests
class TestBifrostSetup:
    #print("test connection")
    #def test_connection(self, test_connection):
    #    pass

    print("test install dir")
    def test_install_dir(self, test_install_dir):
        # This is handled by the fixture, so just a placeholder test
        pass
