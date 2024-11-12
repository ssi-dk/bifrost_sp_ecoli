import os
from pathlib import Path

#import pymongo
import pytest
from bifrost_sp_ecoli import launcher
from bifrostlib import common, database_interface, datahandling
from bifrostlib.datahandling import (Component, ComponentReference, Run,
                                     RunReference, Sample, SampleReference)


bifrost_install_dir = os.environ['BIFROST_INSTALL_DIR']
bifrost_config_and_data_path = Path(f"{bifrost_install_dir}/bifrost/test_data")

@pytest.fixture
def sample(use_collection):
    json_entries = [
        {
            "_id": {"$oid": "000000000000000000000001"},
            "display_name": "S1",
            "name": "S1",
            "components": [],
            "categories": {
                "paired_reads": {
                    "summary": {
                        "data": [f"{bifrost_config_and_data_path}/samples/S1_R1.fastq.gz",
                                 f"{bifrost_config_and_data_path}/samples/S1_R2.fastq.gz"]
                    }
                }
            }
        }
    ]
    bson_entries = [database_interface.json_to_bson(i) for i in json_entries]
    use_collection("samples")
    sample = Sample(value=json_entries[0])
    sample.save()
    return sample


class TestBifrostSpEcoli:
    component_name = "bifrost_sp_ecoli_v0.0.1"
    current_dir = os.getcwd()
    json_entries = [
        {
            "_id": {"$oid": "000000000000000000000001"},
            "display_name": "S1",
            "name": "S1",
            "components": [],
            "categories": {
                "paired_reads": {
                    "summary": {
                        "data": [f"{bifrost_config_and_data_path}/samples/S1_R1.fastq.gz",
                                 f"{bifrost_config_and_data_path}/samples/S1_R2.fastq.gz"]
                    }
                }
            }
        }
    ]
    bson_entries = [database_interface.json_to_bson(i) for i in json_entries]

    def test_info(self):
        launcher.main(["launcher.py", "--info"])

    def test_help(self):
        launcher.main(["launcher.py", "--help"])
