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
| `split_chromosomes` | Hail               | Single-chromosome Hail        | `split_chromosomes` in hail_wgs.py        | 

These tasks have subtasks representing different pipeline steps.
Note that there are two different GDS formats: one used by the SNPRelate library, and one used by the SeqArray library. 
They are not compatible and are treated as separate formats. (Ask me how I feel about that.)
 

### QC Steps

2. Use Hail to perform sample and variant QC. -- hail_wgs.py:run_hail_qc()
3. Use Hail to match samples to the clinical data table and perform sex
    inference. -- hail_wgs.py:match_samples()
Ancestry and Kinship:
4. Dump Hail data to PLINK format. -- hail_wgs.py:convert_mt_to_plink()
5. Use KING to calculate a kinship matrix. -- dodo.py
6. Use Hail to filter by MAF and HWE p-value for GWAS. -- hail_wgs.py:gwas_filter()
7. Use Hail to prune SNPs in LD. -- hail_wgs.py:ld_prune()
8. Repeat step 4 to dump LD-pruned data to PLINK format.
9. Convert PLINK format to SNPArray GDS format. -- seqarray_genesis.R:build_snp_gds()
10. Run PCAir to get kinship-corrected PCs. -- seqarray_genesis.R:run_pcair()
11. Run PCRelate to get ancestry-corrected kinship. -- seqarray_genesis.R:run_pcrelate()
12. Calculate inferred race categories from PCs. -- race_prediction.py
13. Split LD-pruned data files by inferred ancestry. -- hail_wgs.py:subset_mt_samples()
14. Repeat steps 8-10 to run PCAir on each single ancestry.
15. Repeat step 13 to split unpruned GWAS-filtered data by ancestry.
Association Analyses:
16. Generate design matrix for association studies. -- build_design_matrix.py
17. Calculate GENESIS null model. -- seqarray_genesis.R:genesis_null_model()
18. Dump unpruned GWAS-filtered data to VCF shards. -- hail_wgs.py:convert_mt_to_vcf_shards()
19. Convert VCF shards to SeqArray GDS shards. -- mpi_vcf2gds.R
20. Run GENESIS GWAS analysis. -- mpi_genesis_gwas.R
21. Repeat steps 18-19 to create GDS shards for unfiltered rare variant data.
22. Merge GDS shards to single SeqArray GDS file -- seqarray_genesis.R:build_seq_gds()
23. Run GENESIS rare variant analysis. -- <not yet finalized>
VCF file output:
24. Merge VCF shards to single VCF files using bcftools. -- dodo.py
25. Repeat steps 19 and 24 to create VCF files for each individual ancestry.


There are also a series of utility scripts for working with hail and spark.
While these may be useful for diagnostics or external analyses, they are not
actually necessary for the doit version of the pipeline: dodo.py currently
launches a spark cluster for running hail all on its own.

lsf_submit_hail_jupyter.sh
    When submitted with bsub, starts a jupyter notebook server running hail. 
    You must then open an SSH tunnel to this notebook in order to interact with
    it.
lsf_submit_hail.sh
    Call with the command to run as a command-line argument. Submits that 
    command to bsub, along with the setup for a spark cluster suitable for 
    running hail.
lsf_submit_hail_schade01.sh
    Same as above, but specifically requests the schade01-1 private node.
launch_local_hail_notebook.sh
    Starts a jupyter notebook running hail on the local machine. You should
    probably only do this on a private node that has many cores free.
launch_local_spark_cluster.sh
    Call with the command to run as a command-line argument. Starts a spark 
    cluster suitable for running hail on the local machine, and then runs
    that command on it. Again, only do this if you're on a node that has
    many cores free.

All other files in this directory are preliminary or experimental versions of 
the files named above. They may be informative to look at, but aren't in
current use.
