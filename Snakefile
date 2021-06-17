from pathlib import Path

ORIGINAL_VCF = Path('/sc/arion/projects/mscic1/techdev.incoming/DNA/all_625_samples_cohort_vcf/625_Samples.cohort.vcf.gz')
COVARIATES_FILE = '/sc/arion/projects/mscic1/data/covariates/clinical_data_deidentified_allsamples/Biobank_clinical_data_table_by_blood_sample_deidentified_UNCONSENTED.csv.gz'
EXOME_BED = "padded_Twist_ComprehensiveExome_targets_hg38.bed"
DESIGN_MATRIX = "/sc/arion/projects/mscic1/data/covariates/clinical_data_deidentified_allsamples/jordad05/625_Samples.cohort.QC_filtered.sample_matched.age_flowcell_PCAir_dmatrix.csv"

COHORT_STEM = ORIGINAL_VCF.with_suffix('').stem  # because original VCF is .vcf.bgz
QC_STEM = COHORT_STEM + ".QC_filtered"
SAMPLE_MATCHED_STEM = QC_STEM + ".sample_matched"
GWAS_STEM = SAMPLE_MATCHED_STEM + ".GWAS_filtered"
LD_STEM = GWAS_STEM + ".LD_pruned"

SCRIPTSDIR = Path('../../scripts/WGS/').resolve()
HAIL_WGS = str(SCRIPTSDIR / "hail_wgs.py")
HAIL_WRAPPER = str(SCRIPTSDIR / "lsf_hail_wrapper.py")
SEQARRAY_GENESIS = str(SCRIPTSDIR / "seqarray_genesis.R")

TRAITS_OF_INTEREST = ["max_severity_moderate", "severity_ever_severe", "severity_ever_eod", "max_who",
        "severity_ever_increased", "who_ever_increased", "who_ever_decreased",
        "recovered", "highest_titer_irnt", "days_onset_to_encounter_log", "covid_encounter_days_log"]

rule gwas_traits_of_interest:
    input: expand("{phenotype}.GENESIS.assoc.txt", phenotype=TRAITS_OF_INTEREST)

# format conversions

rule vcf2mt:
    input:
        vcf=ORIGINAL_VCF
    output:
        mt=directory(f"{COHORT_STEM}.mt")
    resources: 
        cpus=128
    params:
        hail_script=HAIL_WGS,
        pass_output=True,
        hail_cmd="convert-vcf-to-mt"
    script: HAIL_WRAPPER

rule mt2plink:
    input:
        mt="{prefix}.mt"
    output:
        multiext("{prefix}", ".bed", ".bim", ".fam")
    resources:
        cpus=128
    params:
        hail_script=HAIL_WGS,
        hail_cmd="convert-mt-to-plink"
    script: HAIL_WRAPPER
    # TODO: fix chr prefix

rule plink2snpgds:
    input:
        multiext("{prefix}", ".bed", ".bim", ".fam")
    output:
        gds="{prefix}.snp.gds"
    shell: "Rscript {SEQARRAY_GENESIS} plink2snpgds {wildcards.prefix}"

rule mt2vcfshards:
    input:
        mt="{prefix}.mt",
        vcf=ORIGINAL_VCF
    output:
        shards_dir=directory("{prefix}.shards.vcf.bgz")
    resources:
        cpus=128
    params:
        hail_script=HAIL_WGS,
        hail_cmd="convert-mt-to-vcf-shards"
    script: HAIL_WRAPPER

rule build_vcf:
    input:
        shards_dir="{prefix}.shards.vcf.bgz"
    output:
        vcf="{prefix}.vcf.bgz"
    shell:
        """
        ml bcftools
        ml htslib
        bcftools concat --naive -Oz -o {output.vcf} {input.shards_dir}/part-*.bgz
        tabix {output.vcf}
        """

rule vcf2seqgds_shards:
    input:
        shards_dir="{prefix}.shards.vcf.bgz"
    output:
        shards_dir=directory("{prefix}.shards.seq.gds")
    resources:
        mem_mb=16000,
        cpus=128
    shell:
        """
        ml openmpi
        mpirun --mca mpi_warn_on_fork 0 Rscript {SCRIPTSDIR}/mpi_vcf2gds.R {input.shards_dir}
        """

rule vcf2seqgds_single:
    input:
        shards_dir="{prefix}.shards.vcf.bgz"
    output:
        gds="{prefix}.seq.gds"
    resources:
        cpus=48,
        single_host=1
    shell: "Rscript {SCRIPTSDIR}/seqvcf2gds.R {input.shards_dir} {output.gds}"

# qc steps

rule qc:
    input:
        mt=f"{COHORT_STEM}.mt"
    output:
        mt=directory(f"{QC_STEM}.mt")
    resources:
        cpus=128
    params:
        hail_script=HAIL_WGS,
        hail_cmd="run-hail-qc"
    script: HAIL_WRAPPER

rule match_samples:
    input:
        covariates = COVARIATES_FILE,
        mt = f"{QC_STEM}.mt"
    output:
        mt=directory(f"{SAMPLE_MATCHED_STEM}.mt")
    resources:
        cpus=128
    params:
        hail_script=HAIL_WGS,
        hail_cmd="match-samples"
    script: HAIL_WRAPPER

# race and ancestry steps

rule king:
    input:
        multiext("{prefix}", ".bed", ".bim", ".fam")
    output:
        multiext("{prefix}", ".kin0", "X.kin", "X.kin0")
    resources:
        cpus=16,
        single_host=1
    shell: "ml king && king -b {input[0]} --kinship --cpus {threads} --prefix {wildcards.prefix}"

rule pcair:
    input:
        gds=f"{LD_STEM}.snp.gds",
        king=f"{SAMPLE_MATCHED_STEM}.kin0"
    output:
        multiext(f"{SAMPLE_MATCHED_STEM}.PCAir", ".RDS", ".txt")
    params: output_stem=lambda wildcards, output: Path(output[0]).with_suffix('').stem
    shell: "Rscript {SEQARRAY_GENESIS} pcair {input.gds} {params.output_stem}"


use rule pcair as pcair_race with:
   input:
       gds=f"{LD_STEM}.{{race}}_only.snp.gds",
       king=f"{SAMPLE_MATCHED_STEM}.kin0"
   output:
       multiext(f"{SAMPLE_MATCHED_STEM}.{{race}}_only.PCAir", ".RDS", ".txt")

rule race_prediction:
    input:
        covariates=COVARIATES_FILE,
        pcair=f"{SAMPLE_MATCHED_STEM}.PCAir.txt"
    output:
        expand("{race}.indiv_list.txt", race=["WHITE", "BLACK", "HISPANIC", "ASIAN"]),
        table=f"{SAMPLE_MATCHED_STEM}.race_and_PCA.txt"
    script: SCRIPTSDIR / "race_prediction.py"


rule split_races:
    input:
        indiv_list="{race}.indiv_list.txt",
        mt=f"{SAMPLE_MATCHED_STEM}.mt"
    output:
        mt=directory(f"{SAMPLE_MATCHED_STEM}.{{race}}_only.mt")
    resources:
        cpus=128
    params:
        hail_script=HAIL_WGS,
        hail_cmd="subset-mt-samples",
        pass_output=True
    script: HAIL_WRAPPER

rule pcrelate:
    input:
        pcair=f"{SAMPLE_MATCHED_STEM}.PCAir.RDS",
        gds=f"{LD_STEM}.snp.gds"
    output:
        rds=f"{SAMPLE_MATCHED_STEM}.PCRelate.RDS",
    shell: "Rscript {SEQARRAY_GENESIS} pcrelate {SAMPLE_MATCHED_STEM}"

# variant subsets

rule gwas_filter:
    input:
        mt="{prefix}.mt"
    output:
        mt=directory("{prefix}.GWAS_filtered.mt")
    resources:
        cpus=128
    params:
        hail_script=HAIL_WGS,
        hail_cmd="gwas-filter"
    script: HAIL_WRAPPER

rule rare_filter:
    input:
        mt="{prefix}.mt"
    output:
        mt=directory("{prefix}.rare_filtered.mt")
    resources:
        cpus=128
    params:
        hail_script=HAIL_WGS,
        hail_cmd="rare-filter"
    script: HAIL_WRAPPER

rule exome_filter:
    input:
        mt="{prefix}.mt",
        bed=EXOME_BED
    output:
        mt=directory("{prefix}.exome_filtered.mt")
    resources:
        cpus=128,
        mem_mb=12000
    params:
        hail_script=HAIL_WGS,
        hail_cmd="restrict-to-bed",
        pass_output=True
    script: HAIL_WRAPPER

rule prune_ld:
    input:
        mt="{prefix}.mt"
    output:
        mt=directory("{prefix}.LD_pruned.mt")
    resources:
        cpus=128
    params:
        hail_script=HAIL_WGS,
        hail_cmd="ld-prune"
    script: HAIL_WRAPPER

# association tests

rule design_matrix:
    input:
        COVARIATES_FILE,
        flowcells="flowcells.csv",
        pcair=f"{SAMPLE_MATCHED_STEM}.PCAir.txt",
        bvl="MSCIC_blood_viral_load_predictions.csv"
    output:
        DESIGN_MATRIX
    script:
        SCRIPTSDIR / "build_design_matrix.py"

rule null_model:
    input:
        DESIGN_MATRIX,
        rds=f"{SAMPLE_MATCHED_STEM}.PCRelate.RDS"
    output:
        rds=f"{SAMPLE_MATCHED_STEM}.{{phenotype}}.null.RDS"
    resources:
        cpus=128,
        mem_mb=16000
    shell:
        """
        ml openmpi
        mpirun --mca mpi_warn_on_fork 0 Rscript {SCRIPTSDIR}/mpi_null_model_exhaustive.R {output.rds} {wildcards.phenotype}"
        """

rule run_gwas:
    input:
        gds=f"{GWAS_STEM}.shards.seq.gds",
        null_nodel=f"{SAMPLE_MATCHED_STEM}.{{phenotype}}.null.RDS"
    output:
        txt="{phenotype}.GENESIS.assoc.txt"
    resources:
        mem_mb=16000,
        cpus=128
    shell:
        """
        ml openmpi 
        mpirun --mca mpi_warn_on_fork 0 Rscript {SCRIPTSDIR}/mpi_genesis_gwas.R {SAMPLE_MATCHED_STEM} {wildcards.phenotype}"
        """
