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
from pathlib import Path

os.umask(0o2)

try:
    #print(config)
    sample_ref = SampleReference(_id=config.get('sample_id', None), name=config.get('sample_name', None))
    sample:Sample = Sample.load(sample_ref) # schema 2.1

    if sample is None:
        raise Exception("invalid sample passed")

    sample_id=sample['name']
    print(f"sample_ref {sample_ref} and sample id {sample_id}")

    component_ref = ComponentReference(name=config['component_name'])
    component:Component = Component.load(reference=component_ref) # schema 2.1

    if component is None:
        raise Exception("invalid component passed")

    print(f"Component ref: {component_ref}")
    print(f"Component ref: {component}")

    samplecomponent_ref = SampleComponentReference(
        name=SampleComponentReference.name_generator(sample.to_reference(), component.to_reference())
    )
    
    samplecomponent = SampleComponent.load(samplecomponent_ref)
    print(f"sample component_ref {samplecomponent_ref}")
    print(f"sample component {samplecomponent}")
    
    if samplecomponent is None:
        samplecomponent:SampleComponent = SampleComponent(
            sample_reference=sample.to_reference(), 
            component_reference=component.to_reference()
        ) # schema 2.1
    
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


#NB! -> this is simply because i need somewhere else to alter the db directory stored within the mongoDB
component['resources']['db'] = "resources/ecoligenes"
print(f"component db {component['resources']['db']}")
print(f"component {component['resources']}")

resources_dir=f"{os.environ['BIFROST_INSTALL_DIR']}/bifrost/components/bifrost_{component['display_name']}"
print(f"resource dir {resources_dir}")

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
        req_species = component['requirements']['sample']['categories']['species_detection']['summary'].get('species')
        sample_species = sample["categories"]['species_detection']['summary'].get('species')

        print("Requirement species:", req_species)
        print("Sample species:", sample_species)
        print("has_requirements():", samplecomponent.has_requirements())

        if samplecomponent.has_requirements():
            with open(output.check_file, "w") as fh:
                fh.write("")


# Paths for KMA
BASE = Path(workflow.basedir)
DB_FASTA = f"{BASE}/resources/ecoligenes.fasta"

# Build index in run dir (writable + avoids races on shared install resources)
DB_RUN_DIR = f"{component['name']}/db"
DB_PREFIX  = f"{DB_RUN_DIR}/ecoligenes"  # produces ecoligenes.name, ecoligenes.seq.b, etc.

ANALYSIS_DIR = f"{component['name']}/ecoli_analysis/{sample_id}/sp_ecoli_fbi"
KMA_OUT_PREFIX = f"{ANALYSIS_DIR}/colipost"
TYPING_DIR = f"{ANALYSIS_DIR}/ecolitype"

#- Templated section: end --------------------------------------------------------------------------
rule_name = "run_kma_index"
rule run_kma_index:
    message:
        f"Running step:{rule_name}"
    conda:
        "../envs/kma.yaml"
    log:
        out_file = f"{component['name']}/log/{rule_name}.out.log",
        err_file = f"{component['name']}/log/{rule_name}.err.log",
    input:
        req = rules.check_requirements.output.check_file,
        fasta = DB_FASTA,
    params:
        prefix = DB_PREFIX
    output:
        name = f"{DB_PREFIX}.name",
    benchmark:
        f"{component['name']}/benchmarks/{rule_name}.benchmark"
    shell:
        """
        set -euo pipefail

        kma_index -i {input.fasta} -o {params.prefix} 1> {log.out_file} 2> {log.err_file}
        """

rule_name = "run_kma"
rule run_kma:
    message:
        f"Running step:{rule_name}"
    conda:
        "../envs/kma.yaml"
    log:
        out_file = f"{component['name']}/log/{rule_name}.out.log",
        err_file = f"{component['name']}/log/{rule_name}.err.log",
    benchmark:
        f"{component['name']}/benchmarks/{rule_name}.benchmark",
    input:
        req = rules.check_requirements.output.check_file,
        db_idx = rules.run_kma_index.output.name,  # forces indexing first
        reads = sample['categories']['paired_reads']['summary']['data'],
    output:
        aln  = f"{KMA_OUT_PREFIX}.aln",
        frag = f"{KMA_OUT_PREFIX}.frag.gz",
        fsa  = f"{KMA_OUT_PREFIX}.fsa",
        mat  = f"{KMA_OUT_PREFIX}.mat.gz",
        res  = f"{KMA_OUT_PREFIX}.res",
    params:
        db_prefix = DB_PREFIX,
        output_prefix = KMA_OUT_PREFIX
    shell:
        """
        set -euo pipefail
        mkdir -p "$(dirname "{params.output_prefix}")"
        
        kma -ipe {input.reads[0]} {input.reads[1]} -matrix -t_db {params.prefix} -o {params.output_prefix} 1> {log.out_file} 2> {log.err_file}            
        """

rule_name = "run_ecolityping"
rule run_ecolityping:
    message:
        f"Running step:{rule_name}"
    conda:
        "../envs/ecolityping.yaml"
    log:
        out_file = f"{component['name']}/log/{rule_name}.out.log",
        err_file = f"{component['name']}/log/{rule_name}.err.log",
    benchmark:
        f"{component['name']}/benchmarks/{rule_name}.benchmark",
    input:
        req = rules.check_requirements.output.check_file,
        res = rules.run_kma.output.res,
        reads = sample['categories']['paired_reads']['summary']['data'],
    output:
        output_prefix = f"{TYPING_DIR}"
        output_tsv = f"{TYPING_DIR}_final.tsv",
    params:
        genefilter = "GeneFilter.yaml",
        species = sample['categories']['species_detection']['summary']['species'],
        sample_name = sample_id,
        options="--store_tmp --verbose"
    shell:
        """
        set -euo pipefail
        mkdir -p "{output.output_dir}"

        python ecolityping_filter.py --KMA_res {input.res} --config {params.genefilter} --organism {params.species} --sample_id {params.sample_name} {params.options} --output {output.output_prefix}
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
        final_tsv = rules.run_ecolityping.output.output_tsv,
    output:
        complete = rules.all.input
    params:
        samplecomponent_ref_json = samplecomponent.to_reference().json
    script:
        f"{resources_dir}/bifrost_sp_ecoli/datadump.py"
#- Templated section: end --------------------------------------------------------------------------
