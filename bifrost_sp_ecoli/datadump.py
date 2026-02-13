#!/usr/bin/env python3
import pandas as pd
from typing import Dict, Any

from bifrostlib import common
from bifrostlib.datahandling import (
    Category,
    Component,
    Sample,
    SampleComponent,
    SampleComponentReference,
)


# ----------------------------------------------------------------------
# Convert new ecolitype TSV row → old Bifrost serotype summary structure
# ----------------------------------------------------------------------
def convert_ecolitype_row_to_summary(row: Dict[str, Any]) -> Dict[str, Any]:
    """
    Convert the ecolityping TSV row into the classic Bifrost serotype summary.
    """

    # sero_serotype_finder looks like: "O157;H7"
    o_type, h_type = row["sero_serotype_finder"].split(";")

    summary: Dict[str, Any] = {
        "isolate": row["sample_id"],

        # Old structure expected separate wzx / wzy even if identical
        "wzx": o_type if o_type != "-" else "-",
        "wzy": o_type if o_type != "-" else "-",

        # H-type
        "fliC": h_type if h_type != "-" else "-",

        # Combined OH field
        "OH": f"{o_type}:{h_type}" if o_type != "-" and h_type != "-" else "-",

        # stx subtype(s)
        "stx": row["Toxin"] if row["Toxin"] != "-" else "-",

        # Adhesin / virulence
        "eae": row["Adhesin"],
        "ehx": row["Virulence"],

        # Other genes (your TSV encodes this inside toxin_details)
        "other": "-",

        # Static fields (same as old component)
        "wgsrun": "-",
        "wgsdate": "-",
        "ST": "NA",
        "STgenes": "NA",

        # The new verbose field is your toxin_details
        "verbose": row["toxin_details"],
    }

    return summary


# ----------------------------------------------------------------------
# Main datadump entry point
# ----------------------------------------------------------------------
def datadump(
    input: Any,
    output: Any,
    samplecomponent_ref_json: Dict[str, Any]
) -> None:
    """
    Load ecolityping TSV, convert to old-style summary, store in Bifrost.
    """

    # Load Bifrost objects
    samplecomponent_ref = SampleComponentReference(value=samplecomponent_ref_json)
    samplecomponent: SampleComponent = SampleComponent.load(samplecomponent_ref)
    sample: Sample = Sample.load(samplecomponent.sample)
    component: Component = Component.load(samplecomponent.component)

    # Get or create serotype category
    serotype: Category | None = sample.get_category("serotype")
    if serotype is None:
        serotype = Category(value={
            "name": "bifrost_sp_ecoli",
            "component": {
                "id": samplecomponent["component"]["_id"],
                "name": samplecomponent["component"]["name"],
            },
            "summary": {},
            "report": {},
        })

    # Read ecolitype_final.tsv
    final_tsv_path: str = input.final_tsv
    df = pd.read_csv(final_tsv_path, sep="\t")

    # Convert first (and only) row to old-style summary
    row_dict: Dict[str, Any] = df.to_dict(orient="records")[0]
    summary_dict: Dict[str, Any] = convert_ecolitype_row_to_summary(row_dict)

    # Store summary in category
    serotype["summary"] = summary_dict

    # Save back to Bifrost
    samplecomponent.set_category(serotype)
    sample.set_category(serotype)
    samplecomponent.save_files()
    common.set_status_and_save(sample, samplecomponent, "Success")

    # Write completion flag
    with open(output.complete[0], "w") as fh:
        fh.write("done")

datadump(
    snakemake.input,
    snakemake.output,
    snakemake.params.samplecomponent_ref_json,
)
