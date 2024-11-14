import pytest
from bifrostlib import datahandling
from bifrostlib import database_interface
from bifrost_sp_ecoli import launcher
import pymongo
import os
import shutil

@pytest.fixture
def test_cwd_in_class():
    """Test that checks the bifrost install directory."""
    print(f"bifrost cwd: {self.bifrost_install_dir}")
    assert self.bifrost_install_dir != ""
    
def test_connection():
    """Test that checks the bifrost install directory."""
    assert datahandling.has_a_database_connection()

# Test class to check the Bifrost setup
class TestBifrostSetup:
    component_name = "bifrost_sp_ecoli_v.0.0.1"
    # component_name = component_name + "__171019"

    @pytest.fixture(autouse=True)
    def setup(self):
        """Fixture to set up necessary environment variables for testing."""
        self.bifrost_install_dir = os.environ["BIFROST_INSTALL_DIR"]
        self.test_dir = f"{self.bifrost_install_dir}/bifrost/test_data/output/test__whats_my_species/"
        self.r1 = f"{self.bifrost_install_dir}/bifrost/test_data/samples/S1_R1.fastq.gz"
        self.r2 = f"{self.bifrost_install_dir}/bifrost/test_data/samples/S1_R2.fastq.gz"
        self.json_entries = [
            {
                "_id": {"$oid": "000000000000000000000001"},
                "name": "S1",
                "components": [],
                "categories": {"paired_reads": {"summary": {"data": [self.r1, self.r2]}}},
            }
        ]
        self.bson_entries = [database_interface.json_to_bson(i) for i in self.json_entries]
