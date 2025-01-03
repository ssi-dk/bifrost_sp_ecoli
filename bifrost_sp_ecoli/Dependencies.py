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
        self.sample_ref = SampleReference(_id=self.config.get('sample_id', None), name=self.config.get('sample_name', None))
        self.sample:Sample = Sample.load(self.sample_ref)
        if not self.sample:
            raise ValueError("Invalid sample provided")

        self.component_ref = ComponentReference(name=self.config['component_name'])
        self.component:Component = Component.load(self.component_ref) #it was component:component in snakemake - consider this and check if it influence anything
        if not self.component:
            raise ValueError("Invalid component provided")

        self.samplecomponent_ref = SampleComponentReference(
            name=SampleComponentReference.name_generator(self.sample.to_reference(), self.component.to_reference())
        )
        self.samplecomponent = SampleComponent.load(self.samplecomponent_ref)
        if not self.samplecomponent:
            self.samplecomponent = SampleComponent(
                sample_reference=self.sample.to_reference(),
                component_reference=self.component.to_reference()
            )

     def get_bifrost_install_dir(self,bifrost_dir:str=None) -> str:

        if bifrost_dir is None:
            bifrost_dir = os.environ.get("BIFROST_INSTALL_DIR", "")
            if not bifrost_dir:
                raise EnvironmentError("BIFROST_INSTALL_DIR environment variable is not set")
        
        return bifrost_dir

    def get_resources_dir(self,resources_dir:str=None,bifrost_dir:str=None)->str:

        if resources_dir is None:
            bifrost_install_dir = self.get_bifrost_install_dir(bifrost_dir)
            resources_dir = f"{bifrost_install_dir}/bifrost/components/bifrost_{self.component['display_name']}"
        
        return resources_dir
                                                   
    def get_database_path(self,db_name:str=None, db_path:str=None,resources_dir:str=None, bifrost_dir:str=None)->str:

        resources = self.get_resources_dir(resources_dir,bifrost_dir)

        database_name = db_name if db_name else self.component['resources']['db']

        if db_path is None:
            db_path = f"{resources}/bifrost_{self.component['display_name']}/{database_name}"
            print(f"Constructed database path: {db_path}")

        return db_path

    def get_kma_path(self,conda_path:str=None)->str:
        
        if conda_path is None:
            conda_prefix = os.environ.get("CONDA_PREFIX", "")
            
            if not conda_prefix:
                raise EnvironmentError("CONDA_PREFIX environment variable is not set")
            
            kma_path = f"{conda_prefix}/bin"
        else:
            kma_path = f"{conda_path}/bin"
        
        return kma_path

    def get_sample_id(self)->str:
        return self.sample['name']

    def get_paired_reads(self)->list:
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

