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

DB_PREFIX = os.path.join(os.path.dirname(workflow.snakefile), "resources/ecoligenes")

# ------------------------------------------------------------------
# run_kma use prebuilt DB, no indexing
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
        aln  = f"{component['name']}/kma.aln",
        frag = f"{component['name']}/kma.frag.gz",
        fsa  = f"{component['name']}/kma.fsa",
        mat  = f"{component['name']}/kma.mat.gz",
        res  = f"{component['name']}/kma.res",
    params:
        db_prefix     = DB_PREFIX,
        output_prefix = f"{component['name']}/kma"
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
        req = rules.check_requirements.output.check_file,
        res = rules.run_kma.output.res,
    output:
        output_tsv = f"{component['name']}/{sample['name']}_typing_final.tsv"
    params:
        component_prefix = f"{component['name']}/{sample['name']}_typing",
        genefilter    = os.path.join(os.path.dirname(workflow.snakefile), "GeneFilter.yaml"),
        species       = sample['categories']['species_detection']['summary']['species'],
        sample_name   = sample['name'],
        options       = "--store_tmp --verbose",
        typing_script = os.path.join(os.path.dirname(workflow.snakefile), "rule__ecolityping.py")
    shell:
        r"""
        set -euo pipefail

        python {params.typing_script} \
            --KMA_res {input.res} \
            --config {params.genefilter} \
            --organism "{params.species}" \
            --sample_id "{params.sample_name}" \
            {params.options} \
            --output "{params.component_prefix}"
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

