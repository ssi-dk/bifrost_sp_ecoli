#!/usr/bin/env python3
import os
from typing import Dict
from pprint import pprint
import json


from bifrostlib import common
from bifrostlib.datahandling import (Category, Component, Sample,
                                     SampleComponent, SampleComponentReference)


def extract_results_from_json(serotype: Category, results: Dict, component_name: str, file_name: str) -> None:
    with open(file_name) as fd:
        results = json.load(fd)
    serotype["summary"] = results


def datadump(input: object, output: object, samplecomponent_ref_json: Dict):
    samplecomponent_ref = SampleComponentReference(value=samplecomponent_ref_json)
    samplecomponent = SampleComponent.load(samplecomponent_ref)
    sample = Sample.load(samplecomponent.sample)
    component = Component.load(samplecomponent.component)
    
    serotype = sample.get_category("serotype")
    if serotype is None:
        serotype = Category(value={
            "name": "bifrost_sp_ecoli",
            "component": {"id": samplecomponent["component"]["_id"], "name": samplecomponent["component"]["name"]},
            "summary": {
                "ecoli_fbi": "",
            },
            "report": {}
        })
    extract_results_from_json(
        serotype,
        samplecomponent["results"],
        samplecomponent["component"]["name"],
        input.ecoli_analysis_output_file)
    
    samplecomponent.set_category(serotype)
    sample.set_category(serotype)
    samplecomponent.save_files()
    common.set_status_and_save(sample, samplecomponent, "Success")
    pprint(output.complete)
    with open(output.complete[0], "w+") as fh:
        fh.write("done")


datadump(
    snakemake.input,
    snakemake.output,
    snakemake.params.samplecomponent_ref_json,
)
