# bifrost_sp_ecoli

This component is used to analyze specific samples belonging to the *Escherichia coli* or *Shigella* species and determine shiga toxin subtype(s), O type, H type, Adhesin markers, Virulence markers and other detected genes not classified within these main groups.

## Requirements
- The component is alignment-based using the tool [KMA](https://bitbucket.org/genomicepidemiology/kma).
- The versions are described in the [environment.yaml](https://github.com/ssi-dk/bifrost_sp_ecoli/blob/main/environment.yml)
- The alignment uses the trimmed reads from another [component](https://github.com/ssi-dk/bifrost_min_read_check/tree/master).
- The alignment uses a curated reference [database](https://github.com/ssi-dk/bifrost_sp_ecoli/tree/main/bifrost_sp_ecoli/resources).
- This species-specific analysis can only run on samples that belong to the two (currently) previously mentioned species, which have been determined from another [component](https://github.com/ssi-dk/bifrost_whats_my_species).

## Download
```bash
git clone https://github.com/ssi-dk/bifrost_sp_ecoli.git
cd bifrost_sp_ecoli
git submodule init
git submodule update
bash install.sh -i LOCAL
conda activate bifrost_sp_ecoli_vx.x.x
export BIFROST_INSTALL_DIR='/your/path/'
BIFROST_DB_KEY="/your/key/here/" python -m bifrost_sp_ecoli -h
```

## Run the snakemake analysis
Each component can be run on each sample individually using one snakemake command, replacing the config sample_name with the appropriate dataset name. The component name needs to align with the tag used for the conda environment from the most updated GitHub tag
```bash
snakemake -p --nolock --cores 5 -s <path>/pipeline.smk --config sample_name="sample_name" component_name=sp_ecoli__v0.0.1 --rerun-incomplete
```

## Add a species
If any additional species need this typing component, the files [config.yaml](https://github.com/ssi-dk/bifrost_sp_ecoli/blob/main/bifrost_sp_ecoli/config.yaml) (describing the species name and abbreviation) and [GeneFilter.yaml](https://github.com/ssi-dk/bifrost_sp_ecoli/blob/main/bifrost_sp_ecoli/GeneFilter.yaml) (containing the allele-specific threshold) needs to be altered. 

## Analysis
To determine the various groups described above, the [rule_ecolityping](https://github.com/ssi-dk/bifrost_sp_ecoli/blob/main/bifrost_sp_ecoli/rule__ecolityping.py) defined in the [pipeline](https://github.com/ssi-dk/bifrost_sp_ecoli/blob/main/bifrost_sp_ecoli/pipeline.smk) wrangles the procuded kma results using the following.

### KMA result
One example of the aligned kma results *kma.res*
```bash
#Template       Score   Expected        Template_length Template_Identity       Template_Coverage       Query_Identity  Query_Coverage  Depth   q_value p_value
1__wzx__O157__JH959508     78656            1730            1392          100.00          100.00          100.00          100.00           56.60        73613.56        1.0e-26
2__wzy__O157__JH953200     24223            1532            1185          100.00          100.00          100.00          100.00           20.27        19989.30        1.0e-26
5__fliC__H7__AF228487     359574            1726            1758          100.00          100.00          100.00          100.00          204.77        354427.40       1.0e-26
12__eae__eae-42__AF071034         440253            2546            2805          100.00          100.00          100.00          100.00          157.42        432671.71       1.0e-26
32__ehxA__ehxA-3__AB011549        134647            3574            2997           99.97          100.00           99.97          100.00           45.08        124291.98       1.0e-26
6__stx2__stx2-a__X07865   305905            1280            1241           99.92          100.00           99.92          100.00          246.78        302085.84       1.0e-26
7__stx2__stx2-c__AB015057          70033            1552            1241          100.00          100.00          100.00          100.00           56.33        65510.17        1.0e-26
```
Each template describes the gene belonging to one of the categories (shiga toxin subtype(s), O type, H type, Adhesin markers, Virulence markers, other genes). 

### Gene filtering
The following columns *template_cov*, *query_ident*, *Depth* from the KMA result are used for downstream filtering. These values are compared against gene-specific thresholds defined in the [configuration](https://github.com/ssi-dk/bifrost_sp_ecoli/blob/main/bifrost_sp_ecoli/GeneFilter.yaml) file with a list of thresholds for the three columns [#template_cov, query_ident, Depth] as shown below. 

e.g.
```bash
Escherichia:
  stx1: [98,98,10]
  stx2: [98,98,10]
.
.
.
  ehxA: [98,98,10]
  other: [90,90,10]
```
To determine any of the described categories, in this example, the Template needs template coverage and query identity above 98%  with a depth of coverage above 10. Only hits that pass the relevant thresholds are retained for downstream typing, helping ensure that the final calls (e.g. *1__wzx__O157__JH959508* → gene = *wzx*, allele = *O157*) accurately reflect the allele pathogenic profile of the analyzed isolate. 

Once the templates have been filtered according to the described thresholds, the various categories are set according to the following criteria:

- **Shiga Toxin:** Only PASSED hits for _stx*_ genes are used to make the toxin call, and alleles are collected separately for *stx1* and *stx2*. Toxin calling is not based on a single best hit — it reports all passing toxin alleles, and only the unique alleles are reported in the final result. Thus, the *Toxin* field (see below) is the union of all unique *stx1* and *stx2* alleles that passed.

- **O-type:** Only PASSED hits for _wz*_ genes (`wzx`, `wzy`, `wzt`, `wzm`) are used to build O-type candidates. The candidates are separated into pairs `wzx/wzy` and `wzt/wzm`, and the final column `O_type` is only assigned if all valid candidate pairs agree on the same allele. Thus, non-concordant PASS hits among _wz*_ genes do not produce an O-type call, but are retained in the `toxin_details` column.

    Examples:
    1. `wzx/wzy = O157` and `wzt/wzm = O157` → `O_type = O157`
    2. `wzx/wzy = O157` and `wzt/wzm` is invalid or ambiguous → `O_type = O157`
    3. `wzx/wzy` is invalid or ambiguous and `wzt/wzm = O157` → `O_type = O157`
    4. `wzx/wzy = O157` and `wzt/wzm = O26` → `O_type = -`
    5. Both pairs are ambiguous or mismatching → `O_type = -`

- **H-type:** Only PASSED hits for *fliC* or any _fl*_ alleles used to build the H-type candidates. If *fliC* is present, the H-type is collected, and any _fl*_ allele is disregarded. If *fliC* is absent, the H-type is determined from any  _fl*_ alleles. Thus, the H-typing includes a preference rule using *fliC* with _fl*_ as a fallback and only makes the final call when the results are unambiguous. 

 - **Adhesin-type:** Only PASSED hits for *eae* genes are used to determine Adhesin. The *Adhesin* column is treated as a presence/absence marker, not subtype calls, thus if one or several PASS eae alleles are present, the column is set to *positive* 

 - **Virulence-type:** Only PASSED hits for *ehx* genes are used to determine Adhesin. The *Virulence* column is treated as a presence/absence marker, not subtype calls, thus if one or several PASS eae alleles are present, the column is set to *positive*

### Final results
The pipeline generates numerous output files with filtered genes for each category (e.g. *<sample_name>__typing_H.tsv*) , while also creating a final concatenated output file (*<sample_name>_typing_final.tsv*), which is stored within the database (see [datadump](https://github.com/ssi-dk/bifrost_sp_ecoli/blob/main/bifrost_sp_ecoli/datadump.py)).

e.g.
```bash
sample_id       Toxin   sero_serotype_finder    Adhesin Virulence       toxin_details
"<sample_name>"  stx2-a;stx2-c   O157;H7 positive        positive        toxin:pass:stx2-a_100.00_99.92_246.78;stx2-c_100.00_100.00_56.33|fail:-||Otype:pass:wzx_O157_100.00_100.00_56.60;wzy_O157_100.00_100.00_20.27|fail:-||Htype:pass:fliC_H7_100.00_100.00_204.77|fail:-||adhesin_virulence:pass:eae_eae-42_100.00_100.00_157.42;ehxA_ehxA-3_100.00_99.97_45.08|fail:-||other:pass:-|fail:-
```

### Conclusion
From the final results above for this test example *"<sample_name>"*, it is possible to determine this strain *Escherichia coli* O157:H7 is a highly virulent Shiga toxin-producing E. coli (STEC).
