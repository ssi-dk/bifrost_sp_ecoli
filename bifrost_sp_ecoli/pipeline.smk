#- Templated section: start ------------------------------------------------------------------------
import os
import sys
import traceback

from bifrostlib import common
from bifrostlib.datahandling import SampleReference
from bifrostlib.datahandling import Sample
from bifrostlib.datahandling import ComponentReference
from bifrostlib.datahandling import Component
from bifrostlib.datahandling import SampleComponentReference
from bifrostlib.datahandling import SampleComponent
from Dependencies import Dependencies #importing the class from our dependencies python file

os.umask(0o2)

try:
    #print(config)
    
    # Dependency injection
    DI = Dependencies(config)
    
    # sample_ref, sample, component_ref, component, samplecomponent_ref, samplecomponent
    DI.initialize_references()
    
    # status 
    DI.set_status("Running")
    
except Exception as error:
    print(traceback.format_exc(), file=sys.stderr)
    raise Exception("failed to set sample, component and/or samplecomponent")

onerror:
    if not DI.samplecomponent.has_requirements():
        DI.set_status("Requirements not met")
    if DI.samplecomponent['status'] == "Running":
        DI.set_status("Failure")

# Define dynamic paths and parameters using Dependencies
resources_dir = DI.get_resources_dir() # Utilize the BIFROST_INSTALL_DIR
database_path = DI.get_database_path() 
kma_path = DI.get_kma_path() # Utilize the CONDA_PREFIX
sample_id = DI.get_sample_id()
paired_reads = DI.get_paired_reads()

rule all:
    input:
        f"{DI.component['name']}/datadump_complete"
    run:
        DI.set_status("Success")

#touch(temp(f"{component['name']}/initialized")),
rule setup:
    output:
        init_file = touch(temp(f"{DI.component['name']}/initialized")),
    params:
        folder = DI.component['name']
    run:
        DI.samplecomponent['path'] = os.path.join(os.getcwd(), DI.component['name'])
        DI.samplecomponent.save()


rule_name = "check_requirements"
rule check_requirements:
    message:
        f"Running step:{rule_name}"
    log:
        out_file = f"{DI.component['name']}/log/{rule_name}.out.log",
        err_file = f"{DI.component['name']}/log/{rule_name}.err.log",
    benchmark:
        f"{DI.component['name']}/benchmarks/{rule_name}.benchmark"
    input:
        folder = rules.setup.output.init_file,
    output:
        check_file = f"{DI.component['name']}/requirements_met",
    run:
        if DI.samplecomponent.has_requirements():
            with open(output.check_file, "w") as fh:
                fh.write("")

#- Templated section: end --------------------------------------------------------------------------

#* Dynamic section: start **************************************************************************

rule_name = "run_ecolityping"
rule run_ecolityping:
    message:
        f"Running step:{rule_name}"
    log:
        out_file = f"{DI.component['name']}/log/{rule_name}.out.log",
        err_file = f"{DI.component['name']}/log/{rule_name}.err.log",
    benchmark:
        f"{DI.component['name']}/benchmarks/{rule_name}.benchmark",
    input:  # files
        rules.check_requirements.output.check_file,
        reads = paired_reads,
        db = database_path,
    params:  # values
        id = sample_id,
        update = "no",
        kma_dir = kma_path
    output:
        _aln = f"{rules.setup.params.folder}/ecoli_analysis/{id}/sp_ecoli_fbi/colipost.aln",
        _frag = f"{rules.setup.params.folder}/ecoli_analysis/{id}/sp_ecoli_fbi/colipost.frag.gz",
        _fsa = f"{rules.setup.params.folder}/ecoli_analysis/{id}/sp_ecoli_fbi/colipost.fsa",
        _mat = f"{rules.setup.params.folder}/ecoli_analysis/{id}/sp_ecoli_fbi/colipost.mat.gz",
        _res = f"{rules.setup.params.folder}/ecoli_analysis/{id}/sp_ecoli_fbi/colipost.res",
    shell:
        """
        # Type
        python3 {resources_dir}/bifrost_sp_ecoli/ecoli_fbi/ecolityping.py -i {params.id} -R1 {input.reads[0]} -R2 {input.reads[1]} -db {input.db} -k {params.kma_dir} --update \
{params.update} -o {rules.setup.params.folder}/ecoli_analysis 1> {log.out_file} 2> {log.err_file}
	"""

rule_name = "run_postecolityping"
rule run_postecolityping:
    message:
        f"Running step:{rule_name}"
    log:
	out_file = f"{DI.component['name']}/log/{rule_name}.out.log",
        err_file = f"{DI.component['name']}/log/{rule_name}.err.log",
    benchmark:
        f"{DI.component['name']}/benchmarks/{rule_name}.benchmark",
    input:
        rules.check_requirements.output.check_file,
        #folder = rules.run_ecolityping.output.folder,
        _aln = rules.run_ecolityping.output._aln,
        _frag = rules.run_ecolityping.output._frag,
        _fsa = rules.run_ecolityping.output._fsa,
        _mat = rules.run_ecolityping.output._mat,
        _res = rules.run_ecolityping.output._res,
    output:
        _file = f"{rules.setup.params.folder}/ecoli_analysis/{rules.run_ecolityping.params.id}/{rules.run_ecolityping.params.id}.json",
        _tsv = f"{rules.setup.params.folder}/ecoli_analysis/{rules.run_ecolityping.params.id}/{rules.run_ecolityping.params.id}.tsv",
    params:  # values
        sample_id = rules.run_ecolityping.params.sample_id,
    shell:
        """
        # Process
        python3 {resources_dir}/bifrost_sp_ecoli/ecoli_fbi/postecolityping.py -i {params.id} -d {rules.setup.params.folder}/ecoli_analysis 1> {log.out_file} 2> {log.err_file}
        """

#* Dynamic section: end ****************************************************************************

#- Templated section: start ------------------------------------------------------------------------
rule_name = "datadump"
rule datadump:
    message:
        f"Running step:{rule_name}"
    log:
        out_file = f"{DI.component['name']}/log/{rule_name}.out.log",
        err_file = f"{DI.component['name']}/log/{rule_name}.err.log",
    benchmark:
        f"{DI.component['name']}/benchmarks/{rule_name}.benchmark",
    input:
        ecoli_analysis_output_file = rules.run_postecolityping.output._file,
        ecoli_analysis_output_tsv = rules.run_postecolityping.output._tsv,
    output:
        complete = rules.all.input
    params:
        samplecomponent_ref_json = DI.samplecomponent.to_reference().json
    script:
        f"{resources_dir}/bifrost_sp_ecoli/datadump.py"
#- Templated section: end --------------------------------------------------------------------------
