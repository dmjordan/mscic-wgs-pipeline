from pathlib import Path

ORIGINAL_VCF = Path('/sc/arion/projects/mscic1/techdev.incoming/DNA/all_625_samples_cohort_vcf/625_Samples.cohort.vcf.gz')
COVARIATES_FILE = '/sc/arion/projects/mscic1/data/covariates/clinical_data_deidentified_allsamples/Biobank_clinical_data_table_by_blood_sample_deidentified_UNCONSENTED.csv.gz'
EXOME_BED = "padded_Twist_ComprehensiveExome_targets_hg38.bed"
DESIGN_MATRIX = "/sc/arion/projects/mscic1/data/covariates/clinical_data_deidentified_allsamples/jordad05/625_Samples.cohort.QC_filtered.sample_matched.age_flowcell_PCAir_dmatrix.csv"

COHORT_STEM = ORIGINAL_VCF.stem
QC_STEM = COHORT_STEM + ".QC_filtered"
SAMPLE_MATCHED_STEM = QC_STEM + ".sample_matched"
GWAS_STEM = SAMPLE_MATCHED_STEM + ".GWAS_filtered"
LD_STEM = GWAS_STEM + ".LD_pruned"

SCRIPTSDIR = Path('../../scripts/WGS/').resolve()
HAIL_WGS = SCRIPTSDIR / "hail_wgs.py"
SEQARRAY_GENESIS = SCRIPTSDIR / "seqarray_genesis.R"

TRAITS_OF_INTEREST = ["max_severity_moderate", "severity_ever_severe", "severity_ever_eod", "max_who",
        "severity_ever_increased", "who_ever_increased", "who_ever_decreased",
        "recovered", "highest_titer_irnt", "days_onset_to_encounter_log", "covid_encounter_days_log"]

rule gwas_traits_of_interest:
    input: expand("{phenotype}.GENESIS.assoc.txt", phenotype=TRAITS_OF_INTEREST)

# format conversions

rule vcf2mt:
    input: ORIGINAL_VCF
    output: directory(f"{COHORT_STEM}.mt")
    resources: hail=1
    shell: f"{HAIL_WGS} convert-vcf-to-mt {ORIGINAL_VCF} {COHORT_STEM}.mt"

rule mt2plink:
    input: "{prefix}.mt"
    output: multiext("{prefix}", ".bed", ".bim", ".fam")
    resources: hail=1
    shell: "{HAIL_WGS} convert-mt-to-plink {prefix}.mt" # TODO: fix chr prefix

rule plink2snpgds:
    input: multiext("{prefix}", ".bed", ".bim", ".fam")
    output: "{prefix}.snp.gds"
    shell: "Rscript {SEQARRAY_GENESIS} plink2snpgds {prefix}"

rule mt2vcfshards:
    input: "{prefix}.mt"
    output: directory("{prefix}.shards.vcf.bgz")
    resources: hail=1
    shell: "{HAIL_WGS} convert-mt-to-vcf-shards {prefix}.mt {prefix}.shards.vcf.bgz"

rule build_vcf:
    input: "{prefix}.shards.vcf.bgz"
    output: "{prefix}.vcf.bgz"
    shell:
        """
        ml bcftools
        ml htslib
        bcftools concat --naive -Oz -o {prefix}.vcf.bgz {prefix}.shards.vcf.bgz/part-*.bgz
        tabix {prefix}.vcf.bgz
        """

rule vcf2seqgds_shards:
    input: "{prefix}.shards.vcf.bgz"
    output: directory("{prefix}.shards.seq.gds")
    resources:
        threads=128, mem_mb=16000
    shell:
        """
        ml openmpi
        mpirun --mca mpi_warn_on_fork 0 Rscript {SCRIPTSDIR}/mpi_vcf2gds.R {prefix}.shards.vcf.bgz
        """

rule vcf2seqgds_single:
    input: "{prefix}.shards.vcf.bgz"
    output: "{prefix}.seq.gds"
    resources: threads=64
    shell: "Rscript {SCRIPTSDIR}/seqvcf2gds.R {prefix}.shards.vcf.bgz {prefix}.seq.gds"

# qc steps

rule qc:
    input: f"{COHORT_STEM}.mt"
    output: directory(f"{QC_STEM}.mt")
    resources: hail=1
    shell: "{HAIL_WGS} run-hail-qc 625_Samples.cohort.mt"

rule match_samples:
    input: f"{QC_STEM}.mt", COVARIATES_FILE
    output: directory(f"{SAMPLE_MATCHED_STEM}.mt")
    resources: hail=1
    shell: "{HAIL_WGS} match-samples {COVARIATES_FILE} {QC_STEM}.mt"

# race and ancestry steps

rule king:
    input: multiext("{prefix}", ".bed", ".bim", ".fam")
    output: multiext("{prefix}", ".kin0", "X.kin", "X.kin0")
    resources: threads=16
    shell: "ml king && king -b {prefix}.bed --kinship --cpus {resources.threads} --prefix {prefix}"

rule pcair:
    input:
        gds=f"{LD_STEM}.snp.gds",
        king=f"{SAMPLE_MATCHED_STEM}.kin0"
    output: multiext(f"{SAMPLE_MATCHED_STEM}.PCAir", ".RDS", ".txt")
    params:
        output_stem=lambda wildcards, output: Path(output[0]).with_suffix("")
    shell: "Rscript {SEQARRAY_GENESIS} pcair {input.gds} {params.output_stem}"


use rule pcair as pcair_race with:
    input:
        gds=f"{LD_STEM}.{{race}}_only.snp.gds",
        king=f"{SAMPLE_MATCHED_STEM}.kin0"
    output:
        multiext(f"{SAMPLE_MATCHED_STEM}.{{race}}_only.PCAir", ".RDS", ".txt")

rule race_prediction:
    input:
        COVARIATES_FILE,
        f"{SAMPLE_MATCHED_STEM}.PCAir"
    output:
        expand("{race}.indiv_list.txt", race=["WHITE", "BLACK", "HISPANIC", "ASIAN"]),
        f"{SAMPLE_MATCHED_STEM}.race_and_PCA.txt"
    script: SCRIPTSDIR / "race_prediction.py"


rule split_races:
    input:
        indiv_list="{race}.indiv_list.txt",
        mt=f"{SAMPLE_MATCHED_STEM}.mt"
    output:
        directory(f"{SAMPLE_MATCHED_STEM}.{{race}}_only.mt")
    resources: hail=1
    shell: "{HAIL_WGS} subset-mt-samples {input.mt} {input.indiv_list} {output[0]}"

rule pcrelate:
    input:
        f"{SAMPLE_MATCHED_STEM}.PCAir.RDS",
        f"{LD_STEM}.snp.gds"
    output:
        f"{SAMPLE_MATCHED_STEM}.PCRelate.RDS",
    shell: "Rscript {SEQARRAY_GENESIS} pcrelate {prefix}"

# variant subsets

rule gwas_filter:
    input: "{prefix}.mt"
    output: directory("{prefix}.GWAS_filtered.mt")
    resources: hail=1
    shell: "{HAIL_WGS} gwas-filter {prefix}.mt"

rule rare_filter:
    input: "{prefix}.mt"
    output: directory("{prefix}.rare_filtered.mt")
    resources: hail=1
    shell: "{HAIL_WGS} rare-filter {prefix}.mt"

rule exome_filter:
    input: "{prefix}.mt"
    output: directory("{prefix}.exome_filtered.mt")
    resources: hail=1
    shell: "{HAIL_WGS} restrict-to-bed {prefix}.mt {EXOME_BED}"

rule prune_ld:
    input: "{prefix}.mt"
    output: directory("{prefix}.LD_pruned.mt")
    resources: hail=1
    shell: "{HAIL_WGS} ld-prune {prefix}.mt"

# association tests

rule design_matrix:
    input:
        COVARIATES_FILE,
        "flowcells.csv",
        f"{SAMPLE_MATCHED_STEM}.PCAir.txt",
        "MSCIC_blood_viral_load_predictions.csv"
    output:
        DESIGN_MATRIX
    script:
        SCRIPTSDIR / "build_design_matrix.py"

rule null_model:
    input:
        DESIGN_MATRIX,
        f"{SAMPLE_MATCHED_STEM}.PCRelate.RDS"
    output:
        f"{SAMPLE_MATCHED_STEM}.{{phenotype}}.null.RDS"
    resources:
        threads=128, mem_mb=16000
    shell:
        """
        ml openmpi
        mpirun --mca mpi_warn_on_fork 0 Rscript {SCRIPTSDIR}/mpi_null_model_exhaustive.R {SAMPLE_MATCHED_STEM}.{phenotype}.null.RDS {phenotype}"
        """

rule run_gwas:
    input:
        f"{GWAS_STEM}.shards.seq.gds",
        f"{SAMPLE_MATCHED_STEM}.{{phenotype}}.null.RDS"
    output:
        "{phenotype}.GENESIS.assoc.txt"
    resources:
        threads=128, mem_mb=16000
    shell:
        """
        ml openmpi 
        mpirun --mca mpi_warn_on_fork 0 Rscript {SCRIPTSDIR}/mpi_genesis_gwas.R {SAMPLE_MATCHED_STEM} {phenotype}"
        """