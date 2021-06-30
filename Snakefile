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

TRAITS_OF_INTEREST = ["max_severity_moderate", "severity_ever_severe", "severity_ever_eod", "max_who",
        "severity_ever_increased", "who_ever_increased", "who_ever_decreased",
        "recovered", "highest_titer_irnt", "days_onset_to_encounter_log", "covid_encounter_days_log"]

GTEX_MODELS_DIR = Path(config["resourcesdir"]) / "metaxcan_data" / "models"
MASHR_MODELS_DIR = GTEX_MODELS_DIR / "eqtl" / "mashr"
ALL_TISSUES = [model.with_suffix("").name[6:] for model in MASHR_MODELS_DIR.glob("*.db")]

wildcard_constraints:
    chrom=r"chr([0-9]{1,2}|[XYM])",
    race=r"[A-Z]+",
    phenotype_untagged=r"[a-z_]+"


# aggregation rules


rule gwas_traits_of_interest:
    input:
        expand("{phenotype}.GENESIS.{suffix}",
            phenotype=TRAITS_OF_INTEREST, suffix=["assoc.txt", "qq.png", "manhattan.png"])

rule smmat_lof_traits_of_interest:
    input:
        expand("{phenotype}.LOF.GENESIS.SMMAT.{suffix}",
            phenotype=TRAITS_OF_INTEREST, suffix=["assoc.txt", "qq.png", "manhattan.png"])


rule metaxcan_eqtl_mashr_traits_of_interest:
    input:
        expand("spredixcan_results/eqtl/mashr/{phenotype}.smultixcan.txt", phenotype=TRAITS_OF_INTEREST)


rule metaxcan_eqtl_elastic_net_traits_of_interest:
    input:
        expand("spredixcan_results/eqtl/elastic_net/{phenotype}.smultixcan.txt", phenotype=TRAITS_OF_INTEREST)

rule coloc2_traits_of_interest:
    input:
        expand("coloc2/{phenotype}.{tissue}.full_table.txt", phenotype=TRAITS_OF_INTEREST, tissue=ALL_TISSUES)


rule metal_traits_of_interest:
    input:
        expand("{phenotype}.race_meta_1.tbl", phenotype=TRAITS_OF_INTEREST)


# utility base tasks

rule hail_base:
    resources:
        cpus=128
    script: os.path.join(config["scriptsdir"], "lsf_hail_wrapper.py") 

rule genesis_base:
    script: os.path.join(config["scriptsdir"], "seqarray_genesis.R")

# format conversions

use rule hail_base as vcf2mt with:
    input:
        vcf=ORIGINAL_VCF
    output:
        mt=directory(f"{COHORT_STEM}.mt")
    params:
        pass_output=True,
        hail_cmd="convert-vcf-to-mt"

use rule hail_base as mt2plink with:
    input:
        mt="{prefix}.mt"
    output:
        multiext("{prefix}", ".bed", ".bim", ".fam")
    params:
        hail_cmd="convert-mt-to-plink"

use rule genesis_base as plink2snpgds with:
    input:
        multiext("{prefix}", ".bed", ".bim", ".fam")
    output:
        gds="{prefix}.snp.gds"
    params:
        genesis_cmd="plink2snpgds"

use rule hail_base as mt2vcfshards with:
    input:
        mt="{prefix}.mt",
        vcf=ORIGINAL_VCF
    output:
        shards_dir=directory("{prefix}.shards.vcf.bgz")
    params:
        hail_cmd="convert-mt-to-vcf-shards"

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

use rule hail_base as qc with:
    input:
        mt=f"{COHORT_STEM}.mt"
    output:
        mt=directory(f"{QC_STEM}.mt")
    params:
        hail_cmd="run-hail-qc"

use rule hail_base as match_samples with:
    input:
        covariates = COVARIATES_FILE,
        mt = f"{QC_STEM}.mt"
    output:
        mt=directory(f"{SAMPLE_MATCHED_STEM}.mt")
    params:
        hail_cmd="match-samples"

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

use rule genesis_base as pcair with:
    input:
        gds=f"{LD_STEM}.snp.gds",
        king=f"{SAMPLE_MATCHED_STEM}.kin0"
    output:
        multiext(f"{SAMPLE_MATCHED_STEM}.PCAir", ".RDS", ".txt")
    params:
        output_stem=lambda wildcards, output: Path(output[0]).with_suffix('').stem,
        genesis_cmd="pcair"

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
        table=f"{SAMPLE_MATCHED_STEM}.race_and_PCA.csv"
    script: os.path.join(config["scriptsdir"], "race_prediction.py")


use rule hail_base as split_races with:
    input:
        mt=f"{SAMPLE_MATCHED_STEM}.mt",
        indiv_list="{race}.indiv_list.txt"
    output:
        mt=directory(f"{SAMPLE_MATCHED_STEM}.{{race}}_only.mt")
    params:
        hail_cmd="subset-mt-samples",
        pass_output=True

use rule genesis_base as pcrelate with:
    input:
        pcair=f"{SAMPLE_MATCHED_STEM}.PCAir.RDS",
        gds=f"{LD_STEM}.snp.gds"
    output:
        rds=f"{SAMPLE_MATCHED_STEM}.PCRelate.RDS"
    params:
        genesis_cmd="pcrelate",
        prefix=SAMPLE_MATCHED_STEM

# variant subsets

use rule hail_base as gwas_filter with:
    input:
        mt="{prefix}.mt"
    output:
        mt=directory("{prefix}.GWAS_filtered.mt")
    params:
        hail_cmd="gwas-filter"

use rule hail_base as rare_filter with:
    input:
        mt="{prefix}.mt"
    output:
        mt=directory("{prefix}.rare_filtered.mt")
    params:
        hail_cmd="rare-filter"

use rule hail_base as exome_filter with:
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

use rule hail_base as prune_ld with:
    input:
        mt="{prefix}.mt"
    output:
        mt=directory("{prefix}.LD_pruned.mt")
    resources:
        cpus=128,
        mem_mb=8000
    params:
        hail_cmd="ld-prune"

use rule hail_base as lof_filter with:
    input:
        mt=f"{SAMPLE_MATCHED_STEM}.VEP_annotated.mt"
    output:
        mt=directory(f"{SAMPLE_MATCHED_STEM}.LOF_filtered.mt")
    params:
        hail_cmd="filter-lof-hc"

use rule hail_base as functional_filter with:
    input:
        mt=f"{SAMPLE_MATCHED_STEM}.VEP_annotated.mt"
    output:
        mt=directory(f"{SAMPLE_MATCHED_STEM}.functional_filtered.mt")
    params:
        hail_cmd="filter-functional-variation"

use rule hail_base as pext_filter with:
    input:
        mt=f"{SAMPLE_MATCHED_STEM}.pext_annotated.mt"
    output:
        mt=f"{SAMPLE_MATCHED_STEM}.pext_filtered.mt"
    params:
        hail_cmd="filter-hi-pext"


use rule hail_base as pext_lof_filter with:
    input:
        mt=f"{SAMPLE_MATCHED_STEM}.pext_annotated.mt"
    output:
        mt=f"{SAMPLE_MATCHED_STEM}.pext_LOF_filtered.mt"
    params:
        hail_cmd="filter-hi-pext-lof"


use rule hail_base as chrom_split with:
    input:
        mt="{prefix}.mt"
    output:
        mt=directory("{prefix}.{chrom}.mt")
    params:
        hail_cmd="split-chromosomes",
        hail_extra_args="{wildcards.chrom}"

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

rule null_model:
    input:
        DESIGN_MATRIX,
        rds=f"{SAMPLE_MATCHED_STEM}.PCRelate.RDS"
    output:
        rds=f"{SAMPLE_MATCHED_STEM}.{{phenotype_untagged}}.null.RDS"
    resources:
        cpus=128,
        mem_mb=16000
    params:
        script_path=os.path.join(config["scriptsdir"], "mpi_null_model_exhaustive.R")
    shell:
        """
        ml openmpi
        mpirun --mca mpi_warn_on_fork 0 Rscript {params.script_path} {SAMPLE_MATCHED_STEM} {wildcards.phenotype_untagged}
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
        null_nodel=f"{SAMPLE_MATCHED_STEM}.{{phenotype}}.null.RDS"
    output:
        txt="{phenotype}.GENESIS.assoc.txt"
    resources:
        cpus=128,
        mem_mb=16000
    params:
        script_path=os.path.join(config["scriptsdir"], "mpi_genesis_gwas.R")
    shell:
        """
        ml openmpi
        mpirun --mca mpi_warn_on_fork 0 Rscript {params.script_path} {SAMPLE_MATCHED_STEM} {wildcards.phenotype}
        """

use rule genesis_base as gwas_plots with:
    input:
        "{phenotype}.GENESIS.assoc.txt"
    output:
        "{phenotype}.GENESIS.qq.png",
        "{phenotype}.GENESIS.manhattan.png"
    params:
        genesis_cmd="gwas_plots"

rule run_smmat:
    input:
        f"{SAMPLE_MATCHED_STEM}.{{subset}}_filtered.seq.gds",
        f"{SAMPLE_MATCHED_STEM}.{{phenotype}}.null.RDS"
    output:
        multiext("{phenotype}.{subset}.GENESIS.SMMAT", ".assoc.txt", ".manhattan.png", ".qq.png")
    resources:
        cpus=128,
        mem_mb=16000
    params:
        script_path=os.path.join(config["scriptsdir"], "mpi_genesis_smmat.R")
    shell:
        """
        ml openmpi
        mpirun --mca mpi_warn_on_fork 0 Rscript {input[0]} {input[1]} {wildcards.phenotype}.{wildcards.subset}
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
# variant annotation tasks

use rule hail_base as vep with:
    input:
        f"{SAMPLE_MATCHED_STEM}.mt",
        "vep/vep_config_script.json"
    output:
        directory(f"{SAMPLE_MATCHED_STEM}.VEP_annotated.mt")
    params:
        hail_cmd="run-vep"


use rule genesis_base as isoform_expression_tsv with:
    input:
        "/sc/arion/projects/mscic1/data/covariates/clinical_data_deidentified_allsamples/RNASeq_SummarizedExperiment/SummarizedExperiment_Kallisto_Transcripts_RNA_Samples.RDS"
    output:
        "transcript_isoform_tpm_table.tsv"
    params:
        genesis_cmd="isoform_table"

use rule hail_base as isoform_expression_ht with:
    input:
        "transcript_isoform_tpm_table.tsv"
    output:
        "tx_annot.tx_summary.ht",
        "tx_annot.gene_maximums.ht"
    params:
        hail_cmd="transform-tpm-table",
        pass_output=True

use rule hail_base as annotate_pext with:
    input:
        f"{SAMPLE_MATCHED_STEM}.VEP_annotated.mt",
        "tx_annot.tx_summary.ht",
        "tx_annot.gene_maximums.ht"
    output:
        f"{SAMPLE_MATCHED_STEM}.pext_annotated.mt"
    params:
        hail_cmd="annotate-pext"


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
        "{phenotype_1}.{phenotype_2}.assoc.rg.log",
    resources:
        mem_mb=16000
    shell:
        """
        ml ldsc
        ldsc.py --rg {output[1]},{output[2]} \ 
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
        script_path=os.path.join(config["scriptsdir"], "MetaXcan", "software", "SPrediXcan.py")
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
        """

rule smultixcan:
    input:
        expand("spredixcan_results/eqtl/{model_type}/{phenotype}.{tissue}.csv", tissue=ALL_TISSUES, allow_missing=True),
        harmonized="{phenotype}.GENESIS.assoc.metaxcan_harmonized.txt",
        covar=GTEX_MODELS_DIR / "gtex_v8_expression_{model_type}_snp_smultixcan_covariance.txt.gz"
    output:
        multixcan="spredixcan_results/eqtl/{model_type}/{phenotype}.smultixcan.txt"
    params:
        script_path=os.path.join(config["scriptsdir"], "MetaXcan", "software", "SMulTiXcan.py")
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
                --output {output.multixcan} \
        """

rule coloc2:
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

    script: os.path.join(config["scriptsdir"], "do_coloc2.R")
