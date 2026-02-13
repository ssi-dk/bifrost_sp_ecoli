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
    sample_ref = SampleReference(
        _id=config.get('sample_id', None),
        name=config.get('sample_name', None)
    )
    sample: Sample = Sample.load(sample_ref)
    if sample is None:
        raise Exception("invalid sample passed")

    component_ref = ComponentReference(name=config['component_name'])
    component: Component = Component.load(reference=component_ref)
    if component is None:
        raise Exception("invalid component passed")

    samplecomponent_ref = SampleComponentReference(
        name=SampleComponentReference.name_generator(
            sample.to_reference(),
            component.to_reference()
        )
    )
    samplecomponent = SampleComponent.load(samplecomponent_ref)
    if samplecomponent is None:
        samplecomponent = SampleComponent(
            sample_reference=sample.to_reference(),
            component_reference=component.to_reference()
        )

    common.set_status_and_save(sample, samplecomponent, "Running")

except Exception:
    print(traceback.format_exc(), file=sys.stderr)
    raise Exception("failed to set sample, component and/or samplecomponent")

onerror:
    if not samplecomponent.has_requirements():
        common.set_status_and_save(sample, samplecomponent, "Requirements not met")
    if samplecomponent["status"] == "Running":
        common.set_status_and_save(sample, samplecomponent, "Failure")

envvars:
    "BIFROST_INSTALL_DIR",
    "CONDA_PREFIX"

rule all:
    input:
        f"{component['name']}/datadump_complete"
    run:
        common.set_status_and_save(sample, samplecomponent, "Success")

rule setup:
    output:
        init_file = touch(temp(f"{component['name']}/initialized")),
    params:
        folder = component["name"]
    run:
        samplecomponent["path"] = os.path.join(os.getcwd(), component["name"])
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

# ------------------------------------------------------------------
# KMA database shipped with the component
# ------------------------------------------------------------------

COMPONENT_DIR = os.path.dirname(workflow.snakefile)
DB_PREFIX = os.path.join(COMPONENT_DIR, "resources", "ecoligenes")

ANALYSIS_DIR = f"{component['name']}/ecoli_analysis"
KMA_OUT_PREFIX = f"{ANALYSIS_DIR}/kma"
TYPING_DIR = f"{ANALYSIS_DIR}/ecolitype"

# ------------------------------------------------------------------
# run_kma — use prebuilt DB, no indexing
# ------------------------------------------------------------------
rule_name = "run_kma"
rule run_kma:
    message:
        f"Running step:{rule_name}"
    log:
        out_file = f"{component['name']}/log/{rule_name}.out.log",
        err_file = f"{component['name']}/log/{rule_name}.err.log",
    benchmark:
        f"{component['name']}/benchmarks/{rule_name}.benchmark"
    input:
        req   = rules.check_requirements.output.check_file,
        db_idx = [
            f"{DB_PREFIX}.name",
            f"{DB_PREFIX}.seq.b",
            f"{DB_PREFIX}.length.b",
            f"{DB_PREFIX}.comp.b",
        ],
        reads = sample['categories']['trimmed_reads']['summary']['data'],
    output:
        aln  = f"{KMA_OUT_PREFIX}.aln",
        frag = f"{KMA_OUT_PREFIX}.frag.gz",
        fsa  = f"{KMA_OUT_PREFIX}.fsa",
        mat  = f"{KMA_OUT_PREFIX}.mat.gz",
        res  = f"{KMA_OUT_PREFIX}.res",
    params:
        db_prefix     = DB_PREFIX,
        output_prefix = KMA_OUT_PREFIX
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{params.output_prefix}")"

        kma -ipe {input.reads[0]} {input.reads[1]} \
            -matrix \
            -t_db {params.db_prefix} \
            -o {params.output_prefix} \
            1> {log.out_file} 2> {log.err_file}
        """

# ------------------------------------------------------------------
# ecolityping
# ------------------------------------------------------------------
rule_name = "run_ecolityping"
rule run_ecolityping:
    message:
        f"Running step:{rule_name}"
    log:
        out_file = f"{component['name']}/log/{rule_name}.out.log",
        err_file = f"{component['name']}/log/{rule_name}.err.log",
    benchmark:
        f"{component['name']}/benchmarks/{rule_name}.benchmark"
    input:
        req   = rules.check_requirements.output.check_file,
        res   = rules.run_kma.output.res,
        reads = sample['categories']['paired_reads']['summary']['data'],
    output:
        output_prefix = directory(f"{TYPING_DIR}"),
        output_tsv    = f"{TYPING_DIR}_final.tsv",
    params:
        genefilter    = os.path.join(os.path.dirname(workflow.snakefile), "GeneFilter.yaml"),
        species       = sample['categories']['species_detection']['summary']['species'],
        sample_name   = sample['name'],
        options       = "--store_tmp --verbose",
        typing_script = os.path.join(os.path.dirname(workflow.snakefile), "rule__ecolityping.py")
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{output.output_prefix}"

        python {params.typing_script} \
            --KMA_res {input.res} \
            --config {params.genefilter} \
            --organism "{params.species}" \
            --sample_id "{params.sample_name}" \
            {params.options} \
            --output "{output.output_prefix}"
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
        os.path.join(os.path.dirname(workflow.snakefile), "datadump.py")

#- Templated section: end --------------------------------------------------------------------------

