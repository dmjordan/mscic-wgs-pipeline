# MSCIC COVID-19 Biobank WGS Pipeline

This pipeline is written and maintained by D M Jordan, <daniel.jordan1@mssm.edu>.
Contact me with any questions.

## Files

The scripts in this directory are mostly called by ../../files/WGS/dodo.py to
process data in the ../../files/WGS directory; see the README file in that
directory for more information. (../../files/WGS/dodo.py is symlinked to
dodo.py in this directory.)

The following files in this directory are actually called by doit in the
current pipeline:

* **dodo.py**: doit task definition file. Symlinked as ../../files/WGS/dodo.py, and intended
to be called by running doit from that directory.
* **hail_wgs.py**: contains functions to perform analyses with hail. These need to be run in a
    spark session, and the `start_hail()` function in that file needs to be run
    before any of them. (doit enforces this automatically.) Most of the
    functions in this file perform discrete tasks defined in dodo.py.
* **seqarray_genesis.R**: contains functions to perform analysis with SeqArray and GENESIS in R.
    These functions are designed to be run from within bsub, and will run in
    parallel up to 64 processes if `start_cluster()` is called first. (again,
    doit enforces this.) Tasks that require more parallelism than this
    are broken out into separate files and run with mpirun. Again, most of
    these functions perform discrete tasks defined in dodo.py, which calls out
    to R using rpy2.
* **mpi_vcf2gds.R**: converts a VCF "shards" directory containing many small VCF files into a
    directory containing many small SeqArray GDS files. Run using `mpirun Rscript mpi_vcf2gds.R`, with
    the path to the shards directory as a command line argument.
* **mpi_genesis_gwas.R**: runs GWAS with GENESIS, loading an already-computed null model. Run using
    `mpirun Rscript mpi_genesis_gwas.R`, with the name of the trait as a command line argument.
* **race_prediction.py**: defines a single function, called from dodo.py, to infer races from PCA.
* **build_design_matrix.py**: defines a single function, called from dodo.py, to construct a design
    matrix for genetic association studies from the clinical data table and
    other data sources.

## The Pipeline

The general structure of the pipeline is that there is a series of "endpoints," each of which represents a specific
filtering or transformation of the data, and analysis steps that can be done on each endpoint.
Each endpoint can exist in a number of different data formats, and typically each analysis step requires a specific
format. The primary reason for using doit is to seamlessly handle which formats and which endpoints are required for 
each analysis. 

### Endpoints

| Endpoint Name | Description                                                          | Produced by Task |
|---------------|----------------------------------------------------------------------|------------------|
| `full`        | Complete data set, after QC and matching to clinical data            | `match_samples`  |
| `lof`         | Filtered to only predicted loss of function variants                 | `lof_filter`     |
| `gwas`        | Filtered to variants with MAF > 1% and in Hardy-Weinberg equilibrium | `gwas_filter`    |
| `ld`          | Same as `gwas`, but pruned so no variants are in LD                  | `ld_prune`       |

Each endpoint can also be prefixed with a racial group (`white_`, `black_`, `hispanic_`, or `asian_`) to subset 
samples to those predicted to be that race. (Currently this doesn't work for `lof`, not for any technical reason but
just because it isn't required in the pipeline.) They can also be postfixed with a chromosome (e.g. `_chr3` or `_chrX`)
to subset sites to one chromosome. (This currently doesn't work for every analysis task, again because it isn't required
for most of them.)


### Format Conversions

| Task Name           | From Format        | To Format                     | Implementation                            |
|---------------------|--------------------|-------------------------------|-------------------------------------------|
| `vcf2mt`            | VCF                | Hail (.mt)                    | `convert_vcf_to_mt` in hail_wgs.py        |
| `mt2plink`          | Hail               | PLINK (.bed/.bim/.fam)        | `convert_mt_to_plink` in hail_wgs.py      |
| `build_snp_gds`     | PLINK              | SNPRelate GDS (.snp.gds)      | `build_snp_gds` in seqarray_genesis.R     |
| `mt2vcfshards`      | Hail               | Split VCFs                    | `convert_mt_to_vcf_shards` in hail_wgs.py |
| `build_vcf`         | Split VCFs         | Combined VCF                  | Runs bcftools                             |
| `vcf2gds_shards`    | Split VCFs         | Split SeqArray GDS (.seq.gds) | mpi_vcf2gds.R                             |
| `build_seq_gds`     | Split SeqArray GDS | Combined SeqArray GDS         | `build_seq_gds` in seqarray_genesis.R     |
 
Note that there are two different GDS formats: one used by the SNPRelate library, and one used by the SeqArray library. 
They are not compatible and are treated as separate formats. (Ask me how I figured that out.)

### Preprocessing, Filtering, and Annotation

| Task Name           | Description                                                        | Implementation                       |
|---------------------|--------------------------------------------------------------------|--------------------------------------|
| `qc`                | Perform sample and variant QC                                      | `run_hail_qc` in  hail_wgs.py        |
| `match_samples`     | Match samples to the clinical data table and perform sex inference | `match_samples` in hail_wgs.py       |
| `king`              | Calculate an uncorrected kinship matrix                            | Runs king                            |
| `gwas_filter`       | Filter by MAF and HWE p-value for GWAS                             | `gwas_filter` in hail_wgs.py         |
| `ld_prune`          | Prune SNPs in LD                                                   | `ld_prune` in hail_wgs.py            |
| `pcair`             | Calculate kinship-corrected PCA                                    | `run_pcair` in seqarray_genesis.R    |
| `pcrelate`          | Calculate ancestry-corrected kinship matrix                        | `run_pcrelate` in seqarray_genesis.R |
| `race_prediction`   | Calculate inferred race categories from PCA                        | race_prediction.py                   |
| `split_races`       | Split data files by inferred ancestry                              | `subset_mt_samples` in hail_wgs.py   |
| `split_chromosomes` | Split data files by chromosome                                     | `split_chromosomes` in hail_wgs.py   |
| `vep`               | Run VEP to annotate variants                                       | `run_vep` in hail_wgs.py             |
| `lof_filter`        | Filter to only predicted loss of function sites                    | `filter_lof_hc` in hail_wgs.py       |

### Association Analyses
Most of the actual association analyses have subtasks for individual phenotypes, plus special subtasks
`all` and `phenotypes_of_interest`.

| Task Name       | Description                                                | Implementation                                        |
|-----------------|------------------------------------------------------------|-------------------------------------------------------|
| `design_matrix` | Generate phenotypes and covariates for association studies | build_design_matrix.py                                |
| `null_model`    | Calculate GENESIS null models                              | `genesis_null_model_exhaustive` in seqarray_genesis.R |
| `gwas_to_run`   | Look up which null models were successfully generated      | `gwas_to_run` in dodo.py                              |
| `run_gwas`      | Run GENESIS GWAS analysis                                  | mpi_genesis_gwas.R                                    |
| `run_smmat`     | Run GENESIS rare variant analysis on pLOF variants         | `run_smmat` in seqarray_genesis.R                     |

### Postprocessing

| Task Name             | Description                                                                    | Implementation                          |
|-----------------------|--------------------------------------------------------------------------------|-----------------------------------------|
| `gwas_plots`          | Produces Manhattan and Q-Q plots for all GWAS traits                           | `make_gwas_plots` in seqarray_genesis.R |
| `ld_scores`           | Calculates LD scores for the biobank cohort                                    | Runs ldsc.py                            |
| `munge_sumstats`      | Convert GENESIS's output tables into the .sumstats.gz format required by LDSC  | Runs munge_sumstats.py from LDSC        |
| `ld_score_regression` | Runs LD score regression between all pairs of GWAS traits                      | Runs ldsc.py                            | 

## Hail Scripts

There are also a series of utility scripts for working with hail and spark.
While these may be useful for diagnostics or external analyses, they are not
actually necessary for the doit version of the pipeline: dodo.py currently
launches a spark cluster for running hail all on its own. This may change soon, though, as I am
in the process of converting the doit pipeline to submit tasks through bsub rather than
expecting to be submitted by bsub.

* **lsf_submit_hail_jupyter.sh**: When submitted with bsub, starts a jupyter notebook server running hail. 
    You must then open an SSH tunnel to this notebook in order to interact with
    it.
* **lsf_submit_hail.sh**:
    Call with the command to run as a command-line argument. Submits that 
    command to bsub, along with the setup for a spark cluster suitable for 
    running hail.
* **lsf_submit_hail_schade01.sh**:
    Same as above, but specifically requests the schade01-1 private node.
* **launch_local_hail_notebook.sh**:
    Starts a jupyter notebook running hail on the local machine. You should
    probably only do this on a private node that has many cores free.
* **launch_local_spark_cluster.sh**:
    Call with the command to run as a command-line argument. Starts a spark 
    cluster suitable for running hail on the local machine, and then runs
    that command on it. Again, only do this if you're on a node that has
    many cores free.

