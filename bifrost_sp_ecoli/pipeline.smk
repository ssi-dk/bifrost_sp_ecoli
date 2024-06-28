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
os.umask(0o2)

try:
    #print(config)
    sample_ref = SampleReference(_id=config.get('sample_id', None), name=config.get('sample_name', None))
    sample:Sample = Sample.load(sample_ref) # schema 2.1
    if sample is None:
        raise Exception("invalid sample passed")
    component_ref = ComponentReference(name=config['component_name'])
    component:Component = Component.load(reference=component_ref) # schema 2.1
    if component is None:
        raise Exception("invalid component passed")
    samplecomponent_ref = SampleComponentReference(name=SampleComponentReference.name_generator(sample.to_reference(), component.to_reference()))
    samplecomponent = SampleComponent.load(samplecomponent_ref)
    if samplecomponent is None:
        samplecomponent:SampleComponent = SampleComponent(sample_reference=sample.to_reference(), component_reference=component.to_reference()) # schema 2.1
    common.set_status_and_save(sample, samplecomponent, "Running")
except Exception as error:
    print(traceback.format_exc(), file=sys.stderr)
    raise Exception("failed to set sample, component and/or samplecomponent")

onerror:
    if not samplecomponent.has_requirements():
        common.set_status_and_save(sample, samplecomponent, "Requirements not met")
    if samplecomponent['status'] == "Running":
        common.set_status_and_save(sample, samplecomponent, "Failure")

envvars:
    "BIFROST_INSTALL_DIR",
    "CONDA_PREFIX"

resources_dir=f"{os.environ['BIFROST_INSTALL_DIR']}/bifrost/components/bifrost_{component['display_name']}"

rule all:
    input:
        f"{component['name']}/datadump_complete"
    run:
        common.set_status_and_save(sample, samplecomponent, "Success")

rule setup:
    output:
        init_file = touch(temp(f"{component['name']}/initialized")),
    params:
        folder = component['name']
    run:
        samplecomponent['path'] = os.path.join(os.getcwd(), component['name'])
        samplecomponent.save()


rule_name = "check_requirements"
rule check_requirements:
    message:
        f"Running step:{rule_name}"
    log:
        out_file = f"{component['name']}/log/{rule_name}.out.log",
        err_file = f"{component['name']}/log/{rule_name}.err.log",
    benchmark:
        f"{component['name']}/benchmarks/{rule_name}.benchmark"
    input:
        folder = rules.setup.output.init_file,
    output:
        check_file = f"{component['name']}/requirements_met",
    params:
        samplecomponent
    run:
        if samplecomponent.has_requirements():
            with open(output.check_file, "w") as fh:
                fh.write("")

#- Templated section: end --------------------------------------------------------------------------

#* Dynamic section: start **************************************************************************
rule_name = "run_ecolityping"
rule run_ecolityping:
    message:
        f"Running step:{rule_name}"
    log:
        out_file = f"{component['name']}/log/{rule_name}.out.log",
        err_file = f"{component['name']}/log/{rule_name}.err.log",
    benchmark:
        f"{component['name']}/benchmarks/{rule_name}.benchmark",
    input:  # files
        rules.check_requirements.output.check_file,
        reads = reads,
        db = f"{resources_dir}/{component['resources']['db']}",
    params:  # values
        sample_id = sample.name,
        update = "no",
        kma_path = #TODO
    output:
        folder = directory(rules.setup.params.folder + "/ecoli_analysis"),
    shell:
        """
        # Type
        python3 ecoli_fbi/ecolityping.py -i {params.sample_id} -R1 {input.reads[0]} -R2 {input.reads[1]} -o {output.folder} -db {input.db} -k {params.kma_path} --update {params.update} 1> {log.out_file} 2> {log.err_file}
        """


rule_name = "run_postecolityping"
rule run_postecolityping:
    message:
        f"Running step:{rule_name}"
    log:
        out_file = f"{component['name']}/log/{rule_name}.out.log",
        err_file = f"{component['name']}/log/{rule_name}.err.log",
    benchmark:
        f"{component['name']}/benchmarks/{rule_name}.benchmark",
    input:  # files
        rules.check_requirements.output.check_file,
        reads = sample['categories']['paired_reads']['summary']['data'],
    params:  # values
        sample_id = sample.name,
    output:
        folder = directory(rules.setup.params.folder + "/ecoli_analysis"),
        _file = f"{folder}/{sample_id}.json"
    shell:
        """
        # Process
        python3 ecoli_fbi/postecolityping.py -i {params.sample_id} -d {output.folder} 1> {log.out_file} 2> {log.err_file}"
        """


# rule_name = "run_qc_ecoli_summary"
# rule run_qc_ecoli_summary:
#     message:
#         f"Running step:{rule_name}"
#     log:
#         out_file = f"{component['name']}/log/{rule_name}.out.log",
#         err_file = f"{component['name']}/log/{rule_name}.err.log",
#     benchmark:
#         f"{component['name']}/benchmarks/{rule_name}.benchmark",
#     input:  # files
#         rules.check_requirements.output.check_file,
#     output:
#         folder = directory(rules.setup.params.folder + "/ecoli_analysis"),
#     shell:
#         """
#         # Summarize
#         python3 ecoli_fbi/qc_ecoli_summary.py -i {output.folder} -o {output.folder} 1> {log.out_file} 2> {log.err_file}
#         """


#* Dynamic section: end ****************************************************************************

#- Templated section: start ------------------------------------------------------------------------
rule_name = "datadump"
rule datadump:
    message:
        f"Running step:{rule_name}"
    log:
        out_file = f"{component['name']}/log/{rule_name}.out.log",
        err_file = f"{component['name']}/log/{rule_name}.err.log",
    benchmark:
        f"{component['name']}/benchmarks/{rule_name}.benchmark"
    input:
        #* Dynamic section: start ******************************************************************
        # TODO
        #dtartrate = rules.run_dtartrate.output._file,  # Needs to be output of final rule
        ecoli_analysis_output_file = run_postecolityping.output._file
        #subspecies = rules.run_subspecies.output._file  # Needs to be output of final rule
        #* Dynamic section: end ********************************************************************
    output:
        complete = rules.all.input
    params:
        samplecomponent_ref_json = samplecomponent.to_reference().json
    script:
        os.path.join(os.path.dirname(workflow.snakefile), "datadump.py")
#- Templated section: end --------------------------------------------------------------------------
