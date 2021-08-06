import itertools
import os
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

REGEN_EXOME_PATTERN = Path("/sc/private/regen/data/Regeneron/SINAI_Freeze_Two_pVCF/data/pVCF/QC_passed/freeze2-ontarget/biallelic/SINAI_Freeze_Two.GL.pVCF.PASS.onTarget.biallelic.{chrom}.vcf.gz")
REGEN_WORKING_DIR = Path("/sc/private/regen/IPM-general/jordad05/mscic/")

BIOME_SPLITCHR_STEM = str(REGEN_WORKING_DIR / REGEN_EXOME_PATTERN.with_suffix('').stem)  # because original VCF is .vcf.bgz
BIOME_SPLITCHR_SAMPLE_MATCHED_STEM = BIOME_SPLITCHR_STEM + ".sample_matched"
BIOME_SPLITCHR_GWAS_STEM = BIOME_SPLITCHR_SAMPLE_MATCHED_STEM + ".GWAS_filtered"
BIOME_SPLITCHR_LD_STEM = BIOME_SPLITCHR_GWAS_STEM + ".LD_pruned"

BIOME_CHRALL_SAMPLE_MATCHED_STEM = BIOME_SPLITCHR_SAMPLE_MATCHED_STEM.replace('.{chrom}', '.chrall')
BIOME_CHRALL_GWAS_STEM = BIOME_SPLITCHR_GWAS_STEM.replace('.{chrom}', '.chrall')
BIOME_CHRALL_LD_STEM = BIOME_SPLITCHR_LD_STEM.replace('.{chrom}', '.chrall')

BIOME_GSA_VCF = Path("/sc/private/regen/data/Regeneron/GSA/imputed_tgp_p3_vcf/GSA_chr_all.vcf.gz")
BIOME_GSA_BGEN_PATTERN = Path("/sc/private/regen/data/Regeneron/GSA/imputed_tgp_p3_bgen8bit/GSA_{chrom}.8bit.bgen")
BIOME_GSA_BGEN_FILES = [BIOME_GSA_BGEN_PATTERN.with_name(BIOME_GSA_BGEN_PATTERN.name.format(chrom=f"chr{chrom:02}")) for chrom in range(1,23)]
BIOME_GSA_BGEN_INDEX_FILES = [REGEN_WORKING_DIR / (file.name + ".idx2") for file in BIOME_GSA_BGEN_FILES]
BIOME_GSA_STEM = str(REGEN_WORKING_DIR / BIOME_GSA_VCF.with_suffix('').stem)  # because original VCF is .vcf.bgz
BIOME_GSA_SAMPLE_MATCHED_STEM = BIOME_GSA_STEM + ".sample_matched"
BIOME_GSA_GWAS_STEM = BIOME_GSA_SAMPLE_MATCHED_STEM + ".GWAS_filtered"
BIOME_GSA_LD_STEM = BIOME_GSA_GWAS_STEM + ".LD_pruned"

BIOME_CLINICAL_TABLE = REGEN_WORKING_DIR / "clinical_severity_table_regenid.csv"
EXCLUDED_REGENIDS_TABLE = REGEN_WORKING_DIR / "Regen_Biobank.csv"
BIOME_DESIGN_MATRIX = REGEN_WORKING_DIR / "regenid_dmatrix.csv"
BIOME_GSA_DESIGN_MATRIX = REGEN_WORKING_DIR / "regenid_gsa_dmatrix.csv"

TRAITS_OF_INTEREST = ["max_severity_moderate", "severity_ever_severe", 
        "severity_ever_eod", "max_who",
        "severity_ever_increased", "severity_ever_decreased", "who_ever_increased", "who_ever_decreased",
        "recovered", "highest_titer_irnt", "days_onset_to_encounter_log", "covid_encounter_days_log"]

BIOME_TRAITS = ["max_severity_moderate", "severity_ever_severe", "severity_ever_eod",
                "max_who", "severity_ever_increased", "severity_ever_decreased",
                "who_ever_increased", "who_ever_decreased", "covid_encounter_days_log"]

GTEX_MODELS_DIR = Path(config["resourcesdir"]) / "metaxcan_data" / "models"
MASHR_MODELS_DIR = GTEX_MODELS_DIR / "eqtl" / "mashr"
ALL_TISSUES = [model.with_suffix("").name[6:] for model in MASHR_MODELS_DIR.glob("*.db")]

MSCIC_EQTL_DIR = Path("/sc/arion/projects/mscic1/results/Noam/MainCovid/QTL/")

wildcard_constraints:
    chrom=r"chr([0-9]{1,2}|[XYM])",
    race=r"WHITE|BLACK|ASIAN|HISPANIC",
    phenotype_untagged=r"[a-z_]+",
    phenotype_suffix="(_[A-Z_]+)?",
    tissue="|".join(ALL_TISSUES),
    prefix=".*(?<!shards)"


# aggregation rules


rule gwas_traits_of_interest:
    input:
        expand("{phenotype}.GENESIS.{suffix}",
            phenotype=TRAITS_OF_INTEREST, suffix=["assoc.txt", "assoc.for_locuszoom.txt.bgz", "qq.png", "manhattan.png"])

rule biome_gwas_traits_of_interest:
    input:
        expand("BIOME_{phenotype}.GENESIS.{suffix}",
            phenotype=BIOME_TRAITS, suffix=["assoc.txt", "qq.png", "manhattan.png"])

rule biome_gsa_gwas_traits_of_interest:
    input:
        expand("BIOME_GSA_{phenotype}.GENESIS.{suffix}",
            phenotype=BIOME_TRAITS, suffix=["assoc.txt", "qq.png", "manhattan.png"])

rule smmat_lof_traits_of_interest:
    input:
        expand("{phenotype}.LOF.GENESIS.SMMAT.assoc.txt",
            phenotype=TRAITS_OF_INTEREST)

rule metaxcan_eqtl_mashr_traits_of_interest:
    input:
        expand("spredixcan_results/eqtl/mashr/{phenotype}.smultixcan.txt", phenotype=TRAITS_OF_INTEREST)


rule metaxcan_eqtl_elastic_net_traits_of_interest:
    input:
        expand("spredixcan_results/eqtl/elastic_net/{phenotype}.smultixcan.txt", phenotype=TRAITS_OF_INTEREST)

rule metaxcan_go_eqtl_elastic_net_traits_of_interest:
    input:
        expand("spredixcan_results/eqtl/elastic_net/{phenotype}.{tissue}.GOFigure", phenotype=TRAITS_OF_INTEREST,
            tissue=ALL_TISSUES + ["smultixcan"])


rule coloc2_gtex_traits_of_interest:
    input:
        expand("coloc2/{phenotype}.{tissue}.full_table.txt", phenotype=TRAITS_OF_INTEREST, tissue=ALL_TISSUES)


rule coloc2_mscic_traits_of_interest:
    input:
        expand("coloc2/{phenotype}.Whole_Blood_mscic{race}.full_table.txt", phenotype=TRAITS_OF_INTEREST, race=["WHITE", "BLACK", "HISPANIC", "ASIAN", "ALL"])


rule coloc2_gtex_go_traits_of_interest:
    input:
        expand("coloc2/{phenotype}.{tissue}.GOFigure", phenotype=TRAITS_OF_INTEREST, tissue=ALL_TISSUES)

rule metal_traits_of_interest:
    input:
        expand("{phenotype}.race_meta_1.tbl", phenotype=TRAITS_OF_INTEREST)


# summary / report rules

rule coloc2_report_traits_of_interest:
    input:
        expand("coloc2/{phenotype}.{tissue}.full_table.txt",phenotype=TRAITS_OF_INTEREST,tissue=ALL_TISSUES)
    output:
        "coloc2_report.traits_of_interest.txt"
    script: os.path.join(config["scriptsdir"], "coloc2_report.py")

rule metaxcan_report_elastic_net_traits_of_interest:
    input:
        expand("spredixcan_results/eqtl/elastic_net/{phenotype}.{tissue}.csv", phenotype=TRAITS_OF_INTEREST,tissue=ALL_TISSUES),
        expand("spredixcan_results/eqtl/elastic_net/{phenotype}.smultixcan.txt", phenotype=TRAITS_OF_INTEREST)
    output:
        "metaxcan_report.traits_of_interest.txt"
    script: os.path.join(config["scriptsdir"], "metaxcan_report.py")


# format conversions

rule vcf2mt:
    input:
        vcf=ORIGINAL_VCF
    output:
        mt=directory(f"{COHORT_STEM}.mt")
    params:
        pass_output=True,
        hail_cmd="convert-vcf-to-mt"
    resources:
        cpus = 128,
        mem_mb = 12000
    script: os.path.join(config["scriptsdir"],"lsf_hail_wrapper.py")

rule biome_vcf2mt:
    input:
        vcf=REGEN_EXOME_PATTERN
    output:
        mt=directory(BIOME_SPLITCHR_STEM + ".mt")
    params:
        pass_output=True,
        hail_cmd="convert-vcf-to-mt",
        hail_extra_args="--filter-multi"
    resources:
        cpus = 128,
        mem_mb = 11500
    script: os.path.join(config["scriptsdir"],"lsf_hail_wrapper.py")


rule index_bgen:
    input:
        BIOME_GSA_BGEN_FILES
    output:
        [directory(idx) for idx in BIOME_GSA_BGEN_INDEX_FILES]
    params:
        pass_output=True,
        hail_cmd="index-bgen",
    resources:
        cpus = 22,
        mem_mb = 11500
    script: os.path.join(config["scriptsdir"],"lsf_hail_wrapper.py")


rule biome_gsa_bgen2mt:
    input:
        BIOME_GSA_BGEN_PATTERN.with_name("GSA_Regeneron_ID.sample"),
        BIOME_GSA_BGEN_FILES, BIOME_GSA_BGEN_INDEX_FILES
    output:
        mt=directory(BIOME_GSA_STEM + ".mt")
    params:
        pass_output=True,
        hail_cmd="convert-bgen-to-mt",
    resources:
        cpus = 128,
        mem_mb = 11500
    script: os.path.join(config["scriptsdir"],"lsf_hail_wrapper.py")


rule mt2plink:
    input:
        mt="{prefix}.mt"
    output:
        multiext("{prefix}", ".bed", ".bim", ".fam")
    params:
        hail_cmd="convert-mt-to-plink"
    resources:
        cpus=128,
        mem_mb=11500
    script: os.path.join(config["scriptsdir"],"lsf_hail_wrapper.py")

rule plink2snpgds:
    input:
        multiext("{prefix}", ".bed", ".bim", ".fam")
    output:
        gds="{prefix}.snp.gds"
    params:
        genesis_cmd="plink2snpgds"
    script: os.path.join(config["scriptsdir"],"seqarray_genesis.R")


def original_vcf(wildcards):
    path_stem = Path(wildcards.prefix).stem
    if path_stem.startswith("SINAI"):
        return str(REGEN_EXOME_PATTERN).format("chr21")
    elif path_stem.startswith("GSA_chr_all"):
        return str(BIOME_GSA_VCF)
    elif path_stem.startswith("625_Samples.cohort"):
        return ORIGINAL_VCF
    else:
        raise ValueError(f"Don't know where to find the original VCF for prefix {wildcards.prefix}")


rule mt2vcfshards:
    input:
        mt="{prefix}.mt",
        vcf=original_vcf
    output:
        shards_dir=directory("{prefix}.shards.vcf.bgz")
    params:
        hail_cmd="convert-mt-to-vcf-shards"
    resources:
        cpus=128,
        mem_mb=11500
    script: os.path.join(config["scriptsdir"],"lsf_hail_wrapper.py")

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
    params:
        script_path=os.path.join(config["scriptsdir"], "mpi_vcf2gds.R")
    shell:
        """
        ml openmpi
        mpirun --mca mpi_warn_on_fork 0 Rscript {params.script_path} {input.shards_dir}
        """

rule vcf2seqgds_single:
    input:
        shards_dir="{prefix}.shards.vcf.bgz"
    output:
        gds="{prefix}.seq.gds"
    resources:
        cpus=64
    script: os.path.join(config["scriptsdir"], "seqvcf2gds.R")

# qc steps

rule qc:
    input:
        mt=f"{COHORT_STEM}.mt"
    output:
        mt=directory(f"{QC_STEM}.mt")
    params:
        hail_cmd="run-hail-qc"
    resources:
        cpus=128,
        mem_mb=12000
    script: os.path.join(config["scriptsdir"],"lsf_hail_wrapper.py")

rule match_samples:
    input:
        covariates = COVARIATES_FILE,
        mt = f"{QC_STEM}.mt"
    output:
        mt=directory(f"{SAMPLE_MATCHED_STEM}.mt")
    params:
        hail_cmd="match-samples"
    resources:
        cpus = 128,
        mem_mb = 12000
    script: os.path.join(config["scriptsdir"],"lsf_hail_wrapper.py")

rule biome_sample_list:
    input:
        BIOME_CLINICAL_TABLE, EXCLUDED_REGENIDS_TABLE
    output:
        os.path.join(REGEN_WORKING_DIR, "covid19_hospitalized_regenids.indiv_list.txt")
    shell: """comm -23 <(cut -d, -f1 {input[0]} | sort | uniq) \\
                    <(sort {input[1]}) \\
                    | awk 'NR == 1 {{print "Subject_ID"}} NR > 1 {{print}}' > {output[0]}"""

rule biome_match_samples:
    input:
        mt = "{prefix}.mt",
        indiv_list = os.path.join(REGEN_WORKING_DIR, "covid19_hospitalized_regenids.indiv_list.txt")
    output:
        mt=directory("{prefix}.sample_matched.mt")
    params:
        hail_cmd="subset-mt-samples",
        pass_output=True
    resources:
            cpus = 128,
            mem_mb = 16000
    script: os.path.join(config["scriptsdir"],"lsf_hail_wrapper.py")


# race and ancestry steps

rule king:
    input:
        multiext("{prefix}", ".bed", ".bim", ".fam")
    output:
        "{prefix}.kin0"
    resources:
        cpus=16,
        single_host=1
    shell: "ml king && king -b {input[0]} --kinship --cpus {resources.cpus} --prefix {wildcards.prefix}"

rule pcair:
    input:
        gds="{prefix}.GWAS_filtered.LD_pruned.snp.gds",
        king="{prefix}.kin0"
    output:
        multiext("{prefix}.PCAir", ".RDS", ".txt")
    params:
        output_stem=lambda wildcards, output: str(Path(output[0]).with_suffix('').with_suffix('')),
        genesis_cmd="pcair"
    script: os.path.join(config["scriptsdir"],"seqarray_genesis.R")

rule pcair_race:
   input:
       gds=f"{LD_STEM}.{{race}}_only.snp.gds",
       king=f"{SAMPLE_MATCHED_STEM}.kin0"
   output:
       multiext(f"{SAMPLE_MATCHED_STEM}.{{race}}_only.PCAir", ".RDS", ".txt")
   params:
       output_stem=lambda wildcards, output: str(Path(output[0]).with_suffix('').with_suffix('')),
       genesis_cmd="pcair"
   script: os.path.join(config["scriptsdir"],"seqarray_genesis.R")


rule race_prediction:
    input:
        covariates=COVARIATES_FILE,
        pcair=f"{SAMPLE_MATCHED_STEM}.PCAir.txt"
    output:
        expand("{race}.indiv_list.txt", race=["WHITE", "BLACK", "HISPANIC", "ASIAN"]),
        table=f"{SAMPLE_MATCHED_STEM}.race_and_PCA.csv"
    script: os.path.join(config["scriptsdir"], "race_prediction.py")


rule split_races:
    input:
        mt=f"{SAMPLE_MATCHED_STEM}.mt",
        indiv_list="{race}.indiv_list.txt"
    output:
        mt=directory(f"{SAMPLE_MATCHED_STEM}.{{race}}_only.mt")
    params:
        hail_cmd="subset-mt-samples",
        pass_output=True
    resources:
        cpus=128,
        mem_mb=12000
    script: os.path.join(config["scriptsdir"],"lsf_hail_wrapper.py")

rule pcrelate:
    input:
        pcair="{prefix}.PCAir.RDS",
        gds="{prefix}.GWAS_filtered.LD_pruned.snp.gds"
    output:
        rds="{prefix}.PCRelate.RDS"
    params:
        genesis_cmd="pcrelate"
    script: os.path.join(config["scriptsdir"],"seqarray_genesis.R")

# variant subsets

rule gwas_filter:
    input:
        mt="{prefix}.mt"
    output:
        mt=directory("{prefix}.GWAS_filtered.mt")
    params:
        hail_cmd="gwas-filter"
    resources:
        cpus=128, 
        mem_mb=11500
    script: os.path.join(config["scriptsdir"],"lsf_hail_wrapper.py")

rule rare_filter:
    input:
        mt="{prefix}.mt"
    output:
        mt=directory("{prefix}.rare_filtered.mt")
    params:
        hail_cmd="rare-filter"
    resources:
        cpus=128,
        mem_mb=12000
    script: os.path.join(config["scriptsdir"],"lsf_hail_wrapper.py")

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
        hail_cmd="restrict-to-bed",
        pass_output=True
    script: os.path.join(config["scriptsdir"],"lsf_hail_wrapper.py")

rule prune_ld:
    input:
        mt="{prefix}.mt"
    output:
        mt=directory("{prefix}.LD_pruned.mt")
    resources:
        cpus = 128,
        mem_mb = 16000
    params:
        hail_cmd="ld-prune"
    script: os.path.join(config["scriptsdir"],"lsf_hail_wrapper.py")

rule lof_filter:
    input:
        mt=f"{SAMPLE_MATCHED_STEM}.VEP_annotated.mt"
    output:
        mt=directory(f"{SAMPLE_MATCHED_STEM}.LOF_filtered.mt")
    params:
        hail_cmd="filter-lof-hc"
    resources:
        cpus=128,
        mem_mb=12000
    script: os.path.join(config["scriptsdir"],"lsf_hail_wrapper.py")

rule functional_filter:
    input:
        mt=f"{SAMPLE_MATCHED_STEM}.VEP_annotated.mt"
    output:
        mt=directory(f"{SAMPLE_MATCHED_STEM}.functional_filtered.mt")
    params:
        hail_cmd="filter-functional-variation"
    resources:
        cpus=128,
        mem_mb=12000
    script: os.path.join(config["scriptsdir"],"lsf_hail_wrapper.py")

rule pext_filter:
    input:
        mt=f"{SAMPLE_MATCHED_STEM}.pext_annotated.mt"
    output:
        mt=f"{SAMPLE_MATCHED_STEM}.pext_filtered.mt"
    params:
        hail_cmd="filter-hi-pext"
    resources:
        cpus = 128,
        mem_mb = 12000
    script: os.path.join(config["scriptsdir"],"lsf_hail_wrapper.py")

rule pext_lof_filter:
    input:
        mt=f"{SAMPLE_MATCHED_STEM}.pext_annotated.mt"
    output:
        mt=f"{SAMPLE_MATCHED_STEM}.pext_LOF_filtered.mt"
    params:
        hail_cmd="filter-hi-pext-lof"
    resources:
        cpus = 128,
        mem_mb = 12000
    script: os.path.join(config["scriptsdir"],"lsf_hail_wrapper.py")

rule chrom_split:
    input:
        mt="{prefix}.mt"
    output:
        mt=directory("{prefix}.{chrom}.mt")
    params:
        hail_cmd="split-chromosomes",
        hail_extra_args="{chrom}"
    resources:
        cpus = 128,
        mem_mb = 12000
    script: os.path.join(config["scriptsdir"],"lsf_hail_wrapper.py")

rule chrom_merge:
    input:
        expand("{prefix}.chr{chrom}.{suffix}.mt", chrom=list(range(1,23)) + ["X"], allow_missing=True)
    output:
        directory("{prefix}.chrall.{suffix}.mt")
    resources:
        cpus = 128,
        mem_mb = 11500
    params:
        hail_cmd="merge-chromosomes",
        pass_output=True
    script: os.path.join(config["scriptsdir"], "lsf_hail_wrapper.py")

ruleorder: chrom_merge > biome_match_samples
ruleorder: chrom_merge > gwas_filter
ruleorder: chrom_merge > prune_ld


# association tests

rule design_matrix:
    input:
        COVARIATES_FILE,
        expand(f"{SAMPLE_MATCHED_STEM}.{{race}}_only.PCAir.txt", race=["WHITE", "BLACK", "ASIAN", "HISPANIC"]),
        flowcells="flowcells.csv",
        pcair=f"{SAMPLE_MATCHED_STEM}.PCAir.txt",
        bvl="MSCIC_blood_viral_load_predictions.csv"
    output:
        DESIGN_MATRIX
    script:
        os.path.join(config["scriptsdir"], "build_design_matrix.py")

rule biome_design_matrix:
    input:
        table=BIOME_CLINICAL_TABLE,
        excluded_ids=EXCLUDED_REGENIDS_TABLE,
        pcair="{prefix}.PCAir.txt"
    output:
        "{prefix}.biome_dmatrix.csv"
    script:
        os.path.join(config["scriptsdir"], "build_biome_design_matrix.py")

rule null_model:
    input:
        DESIGN_MATRIX,
        rds=f"{SAMPLE_MATCHED_STEM}.PCRelate.RDS"
    output:
        rds=f"{SAMPLE_MATCHED_STEM}.{{phenotype_untagged}}.null.RDS"
    resources:
        cpus=128,
        mem_mb=20000
    params:
        script_path=os.path.join(config["scriptsdir"], "mpi_null_model_exhaustive.R")
    shell:
        """
        ml openmpi
        mpirun --mca mpi_warn_on_fork 0 Rscript {params.script_path} {SAMPLE_MATCHED_STEM} {wildcards.phenotype_untagged}
        """

rule biome_null_model:
    input:
        "{prefix}.biome_dmatrix.csv",
        rds="{prefix}.PCRelate.RDS"
    output:
        rds="{prefix}.BIOME_{phenotype_untagged}.null.RDS"
    resources:
        cpus=128,
        mem_mb=20000
    params:
        script_path=os.path.join(config["scriptsdir"], "mpi_null_model_exhaustive_biome.R")
    shell:
        """
        ml openmpi
        mpirun --mca mpi_warn_on_fork 0 Rscript {params.script_path} {wildcards.prefix} {wildcards.phenotype_untagged}
        """


rule null_model_race:
    input:
        DESIGN_MATRIX,
        rds=f"{SAMPLE_MATCHED_STEM}.PCRelate.RDS",
        indiv_list="{race}.indiv_list.txt"
    output:
        rds=f"{SAMPLE_MATCHED_STEM}.{{phenotype_untagged}}_{{race}}.null.RDS"
    resources:
        cpus=128,
        mem_mb=16000
    params:
        script_path=os.path.join(config["scriptsdir"],"mpi_null_model_exhaustive_single_race.R")
    shell:
        """
        ml openmpi
        mpirun --mca mpi_warn_on_fork 0 Rscript {params.script_path} \\
                                                {SAMPLE_MATCHED_STEM} \\
                                                {wildcards.phenotype_untagged} \\
                                                {wildcards.race}        """

rule null_model_loo:
    input:
        DESIGN_MATRIX,
        rds=f"{SAMPLE_MATCHED_STEM}.PCRelate.RDS"
    output:
        rds=f"{SAMPLE_MATCHED_STEM}.{{phenotype_untagged}}_leave_{{sample}}_out.null.RDS"
    resources:
        cpus=128,
        mem_mb=16000
    params:
        script_path=os.path.join(config["scriptsdir"],"mpi_null_model_exhaustive.R")
    shell:
        """
        ml openmpi
        mpirun --mca mpi_warn_on_fork 0 Rscript {params.script_path} \\
                                                {SAMPLE_MATCHED_STEM} \\
                                                {wildcards.phenotype_untagged} \\
                                                _leave_{wildcards.sample}_out {wildcards.sample} 
        """

rule run_gwas:
    input:
        gds=f"{GWAS_STEM}.shards.seq.gds",
        null_model=f"{SAMPLE_MATCHED_STEM}.{{phenotype_untagged}}{{phenotype_suffix}}.null.RDS"
    output:
        txt="{{phenotype_untagged}}{{phenotype_suffix}}.GENESIS.assoc.txt"
    resources:
        cpus=128,
        mem_mb=16000
    params:
        script_path=os.path.join(config["scriptsdir"], "mpi_genesis_gwas.R")
    shell:
        """
        ml openmpi
        mpirun --mca mpi_warn_on_fork 0 Rscript {params.script_path} {SAMPLE_MATCHED_STEM} {wildcards.phenotype_untagged}{wildcards.phenotype_suffix}
        """

rule run_biome_gwas:
    input:
        gds=f"{BIOME_CHRALL_GWAS_STEM}.shards.seq.gds",
        null_model=f"{BIOME_CHRALL_SAMPLE_MATCHED_STEM}.BIOME_{{phenotype_untagged}}.null.RDS"
    output:
        txt="BIOME_{phenotype_untagged}.GENESIS.assoc.txt"
    resources:
        cpus=128,
        mem_mb=11500
    params:
        script_path=os.path.join(config["scriptsdir"], "mpi_genesis_gwas.R")
    shell:
        """
        ml openmpi
        mpirun --mca mpi_warn_on_fork 0 Rscript {params.script_path} {BIOME_CHRALL_SAMPLE_MATCHED_STEM} BIOME_{wildcards.phenotype_untagged}
        """

rule run_biome_gsa_gwas:
    input:
        gds=f"{BIOME_GSA_GWAS_STEM}.shards.seq.gds",
        null_model=f"{BIOME_GSA_SAMPLE_MATCHED_STEM}.BIOME_{{phenotype_untagged}}.null.RDS"
    output:
        txt="BIOME_GSA_{phenotype_untagged}.GENESIS.assoc.txt"
    resources:
        cpus=128,
        mem_mb=11500
    params:
        script_path=os.path.join(config["scriptsdir"], "mpi_genesis_gwas.R")
    shell:
        """
        ml openmpi
        mpirun --mca mpi_warn_on_fork 0 Rscript {params.script_path} {BIOME_GSA_SAMPLE_MATCHED_STEM} BIOME_{wildcards.phenotype_untagged} {output.txt}
        """

ruleorder: run_biome_gsa_gwas > run_biome_gwas > run_gwas

rule gwas_plots:
    input:
        "{phenotype}.GENESIS.assoc.txt"
    output:
        "{phenotype}.GENESIS.qq.png",
        "{phenotype}.GENESIS.manhattan.png"
    resources:
        mem_mb="16000"
    params:
        genesis_cmd="gwas_plots"
    script: os.path.join(config["scriptsdir"],"seqarray_genesis.R")

rule run_smmat:
    input:
        f"{SAMPLE_MATCHED_STEM}.{{subset}}_filtered.seq.gds",
        f"{SAMPLE_MATCHED_STEM}.{{phenotype}}.null.RDS"
    output:
        "{phenotype}.{subset}.GENESIS.SMMAT.assoc.txt"
    resources:
        cpus=128,
        mem_mb=16000
    params:
        script_path=os.path.join(config["scriptsdir"], "mpi_genesis_smmat.R")
    shell:
        """
        ml openmpi
        mpirun --mca mpi_warn_on_fork 0 Rscript {params.script_path} {input[0]} {input[1]} {wildcards.phenotype}.{wildcards.subset}
        """

rule metal_three_races:
    input:
        expand("{phenotype_untagged}_{race}.GENESIS.assoc.txt",
            race=["WHITE", "BLACK", "HISPANIC"],
            allow_missing=True)
    output:
        "{phenotype_untagged}.3races_meta_1.tbl"
    script: os.path.join(config["scriptsdir"], "run_metal.py")


rule metal_four_races:
    input:
        expand("{phenotype_untagged}_{race}.GENESIS.assoc.txt",
            race=["WHITE", "BLACK", "HISPANIC", "ASIAN"],
            allow_missing=True)
    output:
        "{phenotype_untagged}.race_meta_1.tbl"
    script: os.path.join(config["scriptsdir"], "run_metal.py")


rule locuszoom_format:
    input:
        "{phenotype}.GENESIS.assoc.txt"
    output:
        multiext("{phenotype}.GENESIS.assoc.for_locuszoom", ".txt.bgz", ".txt.bgz.tbi")
    shell: 
        """ml htslib
        awk -v OFS='\t' '(NR == 1) {{print "#CHR", "POS", "FREQ", "REF", "ALT", "PVAL", "BETA", "SE"}} (NR > 1) {{print $2, $3, $6, $16, $15, $11, $12, $13}}' {input[0]} | \\
        bgzip -c > {output[0]}
        tabix -s 1 -b 2 -e 2 {output[0]}
        """

# variant annotation tasks

rule vep:
    input:
        f"{SAMPLE_MATCHED_STEM}.mt",
        "vep/vep_config_script.json"
    output:
        directory(f"{SAMPLE_MATCHED_STEM}.VEP_annotated.mt")
    params:
        hail_cmd="run-vep"
    resources:
        cpus = 128,
        mem_mb = 12000
    script: os.path.join(config["scriptsdir"],"lsf_hail_wrapper.py")

rule isoform_expression_tsv:
    input:
        "/sc/arion/projects/mscic1/data/covariates/clinical_data_deidentified_allsamples/RNASeq_SummarizedExperiment/SummarizedExperiment_Kallisto_Transcripts_RNA_Samples.RDS"
    output:
        "transcript_isoform_tpm_table.tsv"
    params:
        genesis_cmd="isoform_table"
    script: os.path.join(config["scriptsdir"],"seqarray_genesis.R")

rule isoform_expression_ht:
    input:
        "transcript_isoform_tpm_table.tsv"
    output:
        "tx_annot.tx_summary.ht",
        "tx_annot.gene_maximums.ht"
    params:
        hail_cmd="transform-tpm-table",
        pass_output=True
    resources:
        cpus=128,
        mem_mb=12000
    script: os.path.join(config["scriptsdir"],"lsf_hail_wrapper.py")

rule annotate_pext:
    input:
        f"{SAMPLE_MATCHED_STEM}.VEP_annotated.mt",
        "tx_annot.tx_summary.ht",
        "tx_annot.gene_maximums.ht"
    output:
        f"{SAMPLE_MATCHED_STEM}.pext_annotated.mt"
    params:
        hail_cmd="annotate-pext"
    resources:
        cpus = 128,
        mem_mb = 12000
    script: os.path.join(config["scriptsdir"],"lsf_hail_wrapper.py")

# downstream analyses

# LD Score Regression

rule ld_scores:
    input:
        multiext(f"{GWAS_STEM}.{{chrom}}", ".bed", ".bim", ".fam")
    output:
        multiext(f"{GWAS_STEM}.ld_chr/{{chrom}}.l2", ".M", ".M_5_50", ".ldscore.gz")
    shell:
        """ml ldsc
        mkdir -p {GWAS_STEM}.ld_chr
        ldsc.py --bfile {GWAS_STEM}.{wildcards.chrom} \
        -l2 --ld-wind-kb 1000 --out {GWAS_STEM}/{wildcards.chrom}"""

rule munge_sumstats:
    input:
        assoc="{phenotype}.GENESIS.assoc.txt"
    output:
        sumstats="{phenotype}.GENESIS.assoc.sumstats.gz"
    shell:
        """ml ldsc
        munge_sumstats.py --sumstats {input.assoc} \
                          --out {output.sumstats} \
                          --snp variant.id --N-col n.obs --frq freq \
                          --a1 effect.allele --a2 other.allele \
                          --p Score.pval --signed-sumstats Est,0
        """

rule ldsc:
    input:
        "{phenotype_1}.GENESIS.assoc.sumstats.gz",
        "{phenotype_2}.GENESIS.assoc.sumstats.gz",
        expand(f"{GWAS_STEM}.ld_chr/{{chrom}}.l2.ldscore.gz", chrom=[f"chr{i}" for i in range(1,23)] + ["chrX"])
    output:
        "{phenotype_1}.{phenotype_2}.assoc.rg.log"
    resources:
        mem_mb=16000
    shell:
        """
        ml ldsc
        ldsc.py --rg {input[1]},{input[2]} \ 
                --ref-ld-chr {GWAS_STEM}.ld_chr/ \
                --w-ld-chr {GWAS_STEM}.ld_chr/ \
                --out {wildcards.phenotype_1}.{wildcards.phenotype_2}.assoc.rg"
        """

# MetaXcan


rule metaxcan_harmonize:
    input:
        assoc="{phenotype}.GENESIS.assoc.txt"
    output:
        harmonized="{phenotype}.GENESIS.assoc.metaxcan_harmonized.txt"
    resources:
        mem_mb=20000
    params:
        script_path=os.path.join(config["scriptsdir"], "summary-gwas-imputation", "src", "gwas_parsing.py"),
        metadata_file=os.path.join(config["resourcesdir"], "metaxcan_data", "reference_panel_1000G", "variant_metadata.txt.gz")
    shell:
        """python {params.script_path} \
            -separator ' '  \
            -gwas_file {input.assoc}  \
            -snp_reference_metadata {params.metadata_file} METADATA \
            -output_column_map variant.id variant_id  \
            -output_column_map other.allele non_effect_allele  \
            -output_column_map effect.allele effect_allele  \
            -output_column_map Est effect_size  \
            -output_column_map Est.SE standard_error  \
            -output_column_map Score.pval pvalue  \
            -output_column_map chr chromosome \
            -output_column_map pos position  \
            -output_column_map n.obs sample_size  \
            -output_column_map freq frequency  \
            -output_order variant_id panel_variant_id chromosome position effect_allele non_effect_allele frequency pvalue zscore effect_size standard_error sample_size  \
            -output {output.harmonized}
        """

rule spredixcan:
    input:
        harmonized="{phenotype}.GENESIS.assoc.metaxcan_harmonized.txt",
        model= GTEX_MODELS_DIR / "eqtl" / "{model_type}" / "{model_type}_{tissue}.db",
        covar = GTEX_MODELS_DIR / "eqtl" / "{model_type}" / "{model_type}_{tissue}.txt.gz"
    output:
        predixcan="spredixcan_results/eqtl/{model_type}/{phenotype}.{tissue}.csv"
    params:
        script_path=os.path.join(config["scriptsdir"], "MetaXcan", "software", "SPrediXcan.py"),
        pval_script_path=os.path.join(config["scriptsdir"], "p_adjust.py")
    shell:
        """mkdir -p spredixcan_results/eqtl/{wildcards.model_type}
        python {params.script_path}  \
                --gwas_file {input.harmonized}  \
                --snp_column panel_variant_id  \
                --chromosome_column chromosome  \
                --position_column position  \
                --effect_allele_column effect_allele  \
                --non_effect_allele_column non_effect_allele  \
                --beta_column effect_size  \
                --se_column standard_error  \
                --pvalue_column pvalue  \
                --model_db_snp_key varID  \
                --keep_non_rsid  \
                --additional_output  \
                --overwrite  \
                --throw  \
                --model_db_path {input.model}  \
                --covariance {input.covar}  \
                --output_file {output.predixcan}
        python {params.pval_script_path} {output.predixcan}
        """

rule smultixcan:
    input:
        expand("spredixcan_results/eqtl/{model_type}/{phenotype}.{tissue}.csv", tissue=ALL_TISSUES, allow_missing=True),
        harmonized="{phenotype}.GENESIS.assoc.metaxcan_harmonized.txt",
        covar=GTEX_MODELS_DIR / "gtex_v8_expression_{model_type}_snp_smultixcan_covariance.txt.gz"
    output:
        multixcan="spredixcan_results/eqtl/{model_type}/{phenotype}.smultixcan.txt"
    params:
        script_path=os.path.join(config["scriptsdir"], "MetaXcan", "software", "SMulTiXcan.py"),
        pval_script_path=os.path.join(config["scriptsdir"],"p_adjust.py")

    shell:
        r"""python {params.script_path}  \
                --models_folder {GTEX_MODELS_DIR}/eqtl/{wildcards.model_type}  \
                --models_name_pattern "{wildcards.model_type}_(.*)\.db" \
                --snp_covariance {input.covar}  \
                --metaxcan_folder spredixcan_results/eqtl/{wildcards.model_type}  \
                --metaxcan_filter "{wildcards.phenotype}\.(.*)\.csv" \
                --metaxcan_file_name_parse_pattern "(.*)\.(.*)\.csv" \
                --gwas_file {input.harmonized}  \
                --snp_column panel_variant_id  \
                --chromosome_column chromosome  \
                --position_column position  \
                --effect_allele_column effect_allele  \
                --non_effect_allele_column non_effect_allele  \
                --beta_column effect_size  \
                --se_column standard_error  \
                --pvalue_column pvalue  \
                --model_db_snp_key varID  \
                --keep_non_rsid  \
                --cutoff_condition_number 30  \
                --throw  \
                --output {output.multixcan}
        python {params.pval_script_path} {output.multixcan}
        """

rule coloc2_gtex:
    input:
        eqtl="/sc/arion/projects/mscic1/resources/GTEx_Analysis_v8_eQTL/{tissue}.v8.signif_variant_gene_pairs.txt.gz",
        assoc="{phenotype}.GENESIS.assoc.txt"
    output:
        "coloc2/{phenotype}.{tissue}.full_table.txt",
        directory("coloc2/{phenotype}.{tissue}.coloc.output.perSNP")
    params:
        prefix=lambda wildcards: f"{wildcards.phenotype}.{wildcards.tissue}"
    resources:
        mem_mb=16000

    script: os.path.join(config["scriptsdir"],"do_coloc2_gtex.R")

rule coloc2_mscic:
    input:
        eqtl=MSCIC_EQTL_DIR / "{race}/results/GE_MAINCOVID_SV_30_eQTL_permutations.all.chunks.txt.gz",
        bim=MSCIC_EQTL_DIR / "{race}" / SAMPLE_MATCHED_STEM + ".{race}_only.GWAS_filtered.vcf.bgz_qc4_only_{race}_gene.bim",
        fam=MSCIC_EQTL_DIR / "{race}" / SAMPLE_MATCHED_STEM + ".{race}_only.GWAS_filtered.vcf.bgz_qc4_only_{race}_gene.fam",
        assoc="{phenotype}.GENESIS.assoc.txt"
    output:
        "coloc2/{phenotype}.Whole_Blood_mscic{race}.full_table.txt",
        directory("coloc2/{phenotype}.Whole_Blood_mscic{race}.output.perSNP")
    params:
        prefix=lambda wildcards: f"{wildcards.phenotype}.Whole_Blood_mscic{wildcards.race}"
    resources:
        mem_mb=16000

    script: os.path.join(config["scriptsdir"],"do_coloc2_mscic.R")


rule coloc2_go_enrichment:
    input:
        "coloc2/{phenotype}.{tissue}.full_table.txt"
    output:
        "coloc2/{phenotype}.{tissue}.go_results.txt"
    params:
        data_type = "coloc2"
    script:
        os.path.join(config["scriptsdir"],"run_topGO_enrichment.R")

rule spredixcan_go_enrichment:
    input:
        "spredixcan_results/eqtl/{model_type}/{phenotype}.{tissue}.csv"
    output:
        "spredixcan_results/eqtl/{model_type}/{phenotype}.{tissue}.go_results.txt"
    params:
        data_type = "spredixcan"
    script:
        os.path.join(config["scriptsdir"],"run_topGO_enrichment.R")

rule smultixcan_go_enrichment:
    input:
        "spredixcan_results/eqtl/{model_type}/{phenotype}.smultixcan.txt"
    output:
        "spredixcan_results/eqtl/{model_type}/{phenotype}.smultixcan.go_results.txt"
    params:
        data_type = "smultixcan"
    script:
        os.path.join(config["scriptsdir"],"run_topGO_enrichment.R")


rule gofigure:
    input:
        "{path_prefix}.go_results.txt"
    output:
        directory("{path_prefix}.GOFigure")
    shell:
        r"""python /sc/arion/projects/mscic1/data/GOFigure/GO-Figure/gofigure.py \
            --input {input[0]} \
            --input_type standard \
            --output {output[0]} \
            --ontology all \
            --file_type pdf \
            --outfile_appendix GOFigure_0.5 \
            --similarity_cutoff 0.5
        """
