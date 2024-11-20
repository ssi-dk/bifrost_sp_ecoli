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
    sample_id=sample['name']
    print(f"sample_ref {sample_ref} and sample id {sample_id}")
    if sample is None:
        raise Exception("invalid sample passed")
    component_ref = ComponentReference(name=config['component_name'])
    component:Component = Component.load(reference=component_ref) # schema 2.1
    print(f"Component ref: {component_ref}")
    print(f"Component ref: {component}")
    if component is None:
        raise Exception("invalid component passed")
    samplecomponent_ref = SampleComponentReference(name=SampleComponentReference.name_generator(sample.to_reference(), component.to_reference()))
    samplecomponent = SampleComponent.load(samplecomponent_ref)
    print(f"sample component_ref {samplecomponent_ref}")
    print(f"sample component {samplecomponent}")
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
print(f"resource dir {resources_dir}")

print(f"component db {component['resources']['db']}")
print(f"component {component['resources']}")

rule all:
    input:
        f"{component['name']}/datadump_complete"
    run:
        common.set_status_and_save(sample, samplecomponent, "Success")

#touch(temp(f"{component['name']}/initialized")),
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

print(f"kma database {resources_dir}/bifrost_sp_ecoli/{component['resources']['db']}")
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
        reads = sample['categories']['paired_reads']['summary']['data'],
        db = f"{resources_dir}/bifrost_sp_ecoli/{component['resources']['db']}",
    params:  # values
        sample_id = sample_id,
        update = "no",
        kma_path = f"{os.environ['CONDA_PREFIX']}/bin"
    output:
        folder = directory(rules.setup.params.folder + "/ecoli_analysis"),
        _aln = f"{rules.setup.params.folder}/ecoli_analysis/{sample_id}/sp_ecoli_fbi/colipost.aln",
        _frag = f"{rules.setup.params.folder}/ecoli_analysis/{sample_id}/sp_ecoli_fbi/colipost.frag.gz",
        _fsa = f"{rules.setup.params.folder}/ecoli_analysis/{sample_id}/sp_ecoli_fbi/colipost.fsa",
        _mat = f"{rules.setup.params.folder}/ecoli_analysis/{sample_id}/sp_ecoli_fbi/colipost.mat.gz",
        _res = f"{rules.setup.params.folder}/ecoli_analysis/{sample_id}/sp_ecoli_fbi/colipost.res",
    shell:
        """
        echo {rules.run_ecolityping.output.folder}/{rules.run_ecolityping.params.sample_id}/sp_ecoli_fbi/colipost.aln
        echo {output._aln}
        # Type
        python3 {resources_dir}/bifrost_sp_ecoli/ecoli_fbi/ecolityping.py -i {params.sample_id} -R1 {input.reads[0]} -R2 {input.reads[1]} -db {input.db} -k {params.kma_path} --update \
{params.update} -o {output.folder} 1> {log.out_file} 2> {log.err_file}
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
    input:
        rules.check_requirements.output.check_file,
        folder = rules.run_ecolityping.output.folder,
        _aln = rules.run_ecolityping.output._aln,
        _frag = rules.run_ecolityping.output._frag,
        _fsa = rules.run_ecolityping.output._fsa,
        _mat = rules.run_ecolityping.output._mat,
        _res = rules.run_ecolityping.output._res,
    output:
        _file = f"{rules.run_ecolityping.output.folder}/{rules.run_ecolityping.params.sample_id}/{rules.run_ecolityping.params.sample_id}.json",
        _tsv = f"{rules.run_ecolityping.output.folder}/{rules.run_ecolityping.params.sample_id}/{rules.run_ecolityping.params.sample_id}.tsv",
    params:  # values
        sample_id = rules.run_ecolityping.params.sample_id,
    shell:
        """
        # Process
        python3 {resources_dir}/bifrost_sp_ecoli/ecoli_fbi/postecolityping.py -i {params.sample_id} -d {input.folder} 1> {log.out_file} 2> {log.err_file}
        """

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
        ecoli_analysis_output_file = rules.run_postecolityping.output._file,
        ecoli_analysis_output_tsv = rules.run_postecolityping.output._tsv,
    output:
        complete = rules.all.input
    params:
        samplecomponent_ref_json = samplecomponent.to_reference().json
    script:
        f"{resources_dir}/bifrost_sp_ecoli/datadump.py"
#- Templated section: end --------------------------------------------------------------------------
