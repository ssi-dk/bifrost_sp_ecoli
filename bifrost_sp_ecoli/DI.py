import json
import os
from bifrostlib import common
from bifrostlib.datahandling import SampleReference, Sample, ComponentReference, Component, SampleComponentReference, SampleComponent


class Dependencies:
    def __init__(self, config):
        self.config = config
        self.sample_ref = None
        self.sample = None
        self.component_ref = None
        self.component = None
        self.samplecomponent_ref = None
        self.samplecomponent = None

    def initialize_references(self):
        self.sample_ref = Sample.Reference(_id=self.config.get('sample_id', None), name=self.config.get('sample_name', None))
        self.sample = Sample.load(self.sample_ref)
        if not self.sample:
            raise ValueError("Invalid sample provided")

        self.component_ref = Component.Reference(name=self.config['component_name'])
        self.component = Component.load(self.component_ref)
        if not self.component:
            raise ValueError("Invalid component provided")

        self.samplecomponent_ref = SampleComponent.Reference(
            name=SampleComponent.Reference.name_generator(self.sample.to_reference(), self.component.to_reference())
        )
        self.samplecomponent = SampleComponent.load(self.samplecomponent_ref)
        if not self.samplecomponent:
            self.samplecomponent = SampleComponent(
                sample_reference=self.sample.to_reference(),
                component_reference=self.component.to_reference()
            )

    def get_database_path(self, resources_dir=None, db_name=None):

        if resources_dir is None:
            resources_dir = f"{os.environ['BIFROST_INSTALL_DIR']}/bifrost/components/bifrost_{self.component['display_name']}"

        database_name = db_name if db_name else self.component['resources']['db']

        if db_name is None:
            print(f"Default db_name used: {database_name}")

        db_path = f"{resources_dir}/bifrost_sp_ecoli/{database_name}"
        print(f"Constructed database path: {db_path}")

        return db_path

    def get_sample_id(self):
        return self.sample['name']

    def get_paired_reads(self):
        try:
            return self.sample['categories']['paired_reads']['summary']['data']
        except KeyError:
            raise ValueError("Paired reads data is missing or improperly formatted in the sample")

    def set_status_running(self):
        common.set_status_and_save(self.sample, self.samplecomponent, "Running")

    def set_status_failure(self):
        common.set_status_and_save(self.sample, self.samplecomponent, "Failure")


# If the script is executed directly, test its behavior
if __name__ == "__main__":

    # Testing
    test_config = {
        "sample_id": "sample123",
        "sample_name": "Sample 123",
        "component_name": "ecolityping_component",
    }

    os.environ["BIFROST_INSTALL_DIR"] = "/path/to/bifrost"
    os.environ["CONDA_PREFIX"] = "/path/to/conda"

    # Create Dependencies instance
    deps = Dependencies(test_config)

    # Initialize
    try:
        deps.initialize_references()
        print("Initialized references successfully!")
    except Exception as e:
        print(f"Failed to initialize references: {e}")

    # Test database path retrieval
    try:
        default_db_path = deps.get_database_path()
        custom_db_path = deps.get_database_path(resources_dir="/custom/resources", db_name="custom_db.fasta")
        print(f"Default database path: {default_db_path}")
        print(f"Custom database path: {custom_db_path}")
    except Exception as e:
        print(f"Failed to retrieve database path: {e}")

