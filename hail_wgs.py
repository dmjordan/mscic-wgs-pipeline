import os
import sys
import shutil
from pathlib import Path
import click

import hail as hl
import pandas as pd
from bokeh.io import export_png

import tx_annotation

@click.group("hail-wgs")
def cli():
    hl.init(tmp_dir=os.environ["TMPDIR"], min_block_size=128, local_tmpdir="/local/tmp/",
                log=hl.utils.timestamp_path(os.path.join(os.environ["SPARK_LOG_DIR"], 'hail')))

# creating commands manually so that the python functions are still callable

class ClickPathlibPath(click.Path):
    def coerce_path_result(self, rv):
        return Path(super().coerce_path_result(rv))


def split_vcf_G_field_min(expr, a_index):
    """For a VCF field with number=G (1 value per possible genotype), take the smallest value among
    genotypes that downcode to this genotype. e.g. for triallelic -> biallelic,
    allele 1 gets min(0/0,0/2,2/2),min(0/1,1/2),1/1 and
    allele 2 gets min(0/0,0/1,1/1),min(0/2,1/2),2/2
    Taking the min is appropriate for phred-scaled likelihood values like GP, PL, and PRI.
    This logic is from https://hail.is/docs/0.2/methods/genetics.html#hail.methods.split_multi"""
    return hl.or_missing(hl.is_defined(expr),
                         (hl.range(0, 3).map(lambda i: hl.min(hl.range(0, hl.len(expr))
                           .filter(lambda j: hl.downcode(hl.unphased_diploid_gt_index_call(j), a_index) 
                               == hl.unphased_diploid_gt_index_call(i))
                           .map(lambda j: expr[j])))))


def split_vcf_R_field_sum(expr, a_index):
    """For a VCF field with number=R (1 value per allele incl reference), take the sum across
    all alleles that downcode to this allele, e.g. for triallelic -> biallelic,
    allele 1 gets sum(0,2),1 and allele 2 gets sum(0,1),2.
    Sum is appropriate for count-based fields such as AD.
    This logic is from https://hail.is/docs/0.2/methods/genetics.html#hail.methods.split_multi"""
    return hl.or_missing(hl.is_defined(expr),
                         [hl.sum(expr) - expr[a_index], expr[a_index]])


@cli.command("convert-vcf-to-mt")
@click.argument("vcf_path", type=ClickPathlibPath())
@click.argument("mt_path", type=ClickPathlibPath())
@click.option("--filter-multi/--allow-multi", default=False)
@click.option("--hg19/--hg38", default=False)
def convert_vcf_to_mt(vcf_path, mt_path, filter_multi=False, hg19=False):
    vcf = hl.import_vcf(str(vcf_path.resolve()),
                        force_bgz=True, 
                        reference_genome="GRCh37" if hg19 else "GRCh38",
                        array_elements_required=False,
                        contig_recoding={chrom: f"chr{chrom}" for chrom in [str(x) for x in range(1,23)] + ["X", "Y"]} if not hg19 else {})
    if filter_multi:
        vcf = vcf.filter_rows(hl.len(vcf.alleles) == 2)
    else:
        vcf = hl.split_multi(vcf)
        vcf = vcf.annotate_entries(
                        GT=hl.downcode(vcf.GT, vcf.a_index),
                        AD=split_vcf_R_field_sum(vcf.AD, vcf.a_index),
                        AF=[vcf.AF[vcf.a_index - 1]],
                        F2R1=split_vcf_R_field_sum(vcf.F2R1, vcf.a_index),
                        F1R2=split_vcf_R_field_sum(vcf.F1R2, vcf.a_index),
                        GP=split_vcf_G_field_min(vcf.GP, vcf.a_index),
                        PL=split_vcf_G_field_min(vcf.PL, vcf.a_index),
                        PRI=split_vcf_G_field_min(vcf.PRI, vcf.a_index))
        new_info = vcf.info.annotate(
                        AC=[vcf.info.AC[vcf.a_index - 1]],
                        AF=[vcf.info.AF[vcf.a_index - 1]],
                        MLEAC=[vcf.info.MLEAC[vcf.a_index - 1]],
                        MLEAF=[vcf.info.MLEAF[vcf.a_index - 1]])
        vcf = vcf.annotate_rows(info=new_info)

    vcf.write(str(mt_path.resolve()), overwrite=True)


@cli.command("index-bgen")
@click.argument("bgen_paths", type=click.Path(exists=True, resolve_path=True), nargs=22)
@click.argument("idx_paths", type=click.Path(writable=True, resolve_path=True), nargs=22)
@click.option("--reference-genome", type=click.Choice(["GRCh37", "GRCh38"]), default="GRCh37")
def index_bgen(bgen_paths, idx_paths, reference_genome):
    index_file_map = {bgen_path: idx_path for bgen_path, idx_path in zip(bgen_paths, idx_paths)}
    hl.index_bgen(bgen_paths, index_file_map=index_file_map, reference_genome=reference_genome)

@cli.command("convert-bgen-to-mt")
@click.argument("sample_path", type=click.Path(exists=True, resolve_path=True))
@click.argument("bgen_paths", type=click.Path(exists=True, resolve_path=True), nargs=22)
@click.argument("idx_paths", type=click.Path(writable=True, resolve_path=True), nargs=22)
@click.argument("mt_path", type=click.Path(writable=True, resolve_path=True))
def convert_bgen_to_mt(sample_path, bgen_paths, idx_paths, mt_path):
    index_file_map = {bgen_path: idx_path for bgen_path, idx_path in zip(bgen_paths, idx_paths)}
    mt = hl.import_bgen(bgen_paths, entry_fields=["GT", "dosage"], sample_file=sample_path, index_file_map=index_file_map)
    mt.write(mt_path, overwrite=True)


@cli.command("run-hail-qc")
@click.argument("mt_path", type=ClickPathlibPath())
def run_hail_qc(mt_path):
    mt_path = mt_path.resolve()
    mt = hl.read_matrix_table(str(mt_path))
    mt = hl.sample_qc(mt)
    p = hl.plot.histogram(mt.sample_qc.call_rate, legend="Call Rate", range=(0.8,1))
    export_png(p, str(mt_path.with_suffix(".sample_qc.callrate_hist.png")))
    p = hl.plot.histogram(mt.sample_qc.gq_stats.mean, range=(10,70), legend='Mean Sample GQ')
    export_png(p, str(mt_path.with_suffix(".sample_qc.gq_hist.png")))
    p = hl.plot.histogram(mt.sample_qc.dp_stats.mean, range=(20,50), legend='Mean Sample DP')
    export_png(p, str(mt_path.with_suffix(".sample_qc.dp_hist.png")))
    p = hl.plot.scatter(mt.sample_qc.dp_stats.mean, mt.sample_qc.call_rate, xlabel='Mean DP', ylabel='Call Rate')
    export_png(p, str(mt_path.with_suffix(".sample_qc.dp_v_callrate_scatter.png")))

    mt_filtered = mt.filter_cols((mt.sample_qc.dp_stats.mean > 24) & (mt.sample_qc.call_rate > 0.85) & (mt.sample_qc.gq_stats.mean > 40))
    mt_filtered = mt_filtered.annotate_rows(gt_stats=hl.agg.call_stats(mt_filtered.GT, mt_filtered.alleles))
    mt_filtered = mt_filtered.filter_rows(mt_filtered.gt_stats.AC.any(lambda ac: (ac > 0) & (ac < mt_filtered.gt_stats.AN)))
    print('After applying sample QC, {0}/{1} samples and {2}/{3} variants remain.'.format(mt_filtered.count_cols(), mt.count_cols(), mt_filtered.count_rows(), mt.count_rows()), file=sys.stderr)
    
    mt_filtered = hl.variant_qc(mt_filtered)
    p = hl.plot.histogram(mt_filtered.variant_qc.call_rate, legend="Call Rate", range=(0.8,1))
    export_png(p, str(mt_path.with_suffix(".variant_qc.callrate_hist.png")))
    p = hl.plot.histogram(mt_filtered.variant_qc.gq_stats.mean, range=(10,70), legend='Mean Variant GQ')
    export_png(p, str(mt_path.with_suffix(".variant_qc.gq_hist.png")))
    p = hl.plot.histogram(mt_filtered.variant_qc.dp_stats.mean, range=(20,50), legend='Mean Variant DP')
    export_png(p, str(mt_path.with_suffix(".variant_qc.dp_hist.png")))
    p = hl.plot.scatter(mt_filtered.variant_qc.dp_stats.mean, mt_filtered.variant_qc.call_rate, xlabel='Mean DP', ylabel='Call Rate')
    export_png(p, str(mt_path.with_suffix(".variant_qc.dp_v_callrate_scatter.png")))
    mt_filtered = mt_filtered.filter_rows((mt_filtered.variant_qc.dp_stats.mean > 25) & (mt_filtered.variant_qc.call_rate > 0.9) & (mt_filtered.variant_qc.gq_stats.mean > 40))
    print('After applying sample and variant QC, {0}/{1} variants remain.'.format(mt_filtered.count_rows(), mt.count_rows()), file=sys.stderr)
    mt_filtered.write(str(mt_path.with_suffix(".QC_filtered.mt")), overwrite=True)


@cli.command("match-samples")
@click.argument("covariates_path", type=ClickPathlibPath())
@click.argument("mt_path", type=ClickPathlibPath())
@click.argument("external_x", type=ClickPathlibPath(), default=None, required=False)
def match_samples(covariates_path, mt_path, external_x=None):
    mt_path = mt_path.resolve()
    covariates_path = covariates_path.resolve()
    output_path = mt_path.with_suffix(".sample_matched.mt")
    
    mt = hl.read_matrix_table(str(mt_path))

    # reindex to subjects rather than samples
    # which requires removing duplicated subjects
    mt = mt.filter_cols(~mt.s.endswith("a"))
    mt = mt.annotate_cols(Subject_ID=mt.s.split("_")[-1].split("T")[0])
    mt = mt.rename({"Subject_ID": "s", "s": "sample_id"})
    mt = mt.key_cols_by("s")

    # load covariates table
    covariates_table = hl.import_table(str(covariates_path),
                                  force=True, impute=True, delimiter=",", quote='"')

    # remove samples without matching clinical data
    all_clinical_samples = covariates_table.aggregate(hl.agg.collect_as_set(covariates_table.Subject_ID))
    missing_samples = mt.filter_cols(~hl.set(all_clinical_samples).contains(mt.s))
    mt = mt.anti_join_cols(missing_samples.cols())

    # impute sex
    if external_x is not None:
        external_x = external_x.resolve()
        chrX_mt = hl.read_matrix_table(str(external_x))
        sex_table = hl.impute_sex(chrX_mt.GT, male_threshold=0.6, female_threshold=0.5)
        mt = mt.annotate_cols(is_female=sex_table[mt.sample_id].is_female)
    else:
        sex_table = hl.impute_sex(mt.GT, male_threshold=0.6, female_threshold=0.5)
        mt = mt.annotate_cols(is_female=sex_table[mt.s].is_female)

    # match sex to clinical table
    reported_sex=covariates_table.group_by("Subject_ID").aggregate(sex=hl.literal("/").join(
                                            hl.agg.collect_as_set(covariates_table.Sex)))
    sex_table = sex_table.annotate(reported=reported_sex[sex_table.s])
    mismatched_sex = sex_table.filter((sex_table.is_female &  
        (sex_table.reported.sex != "Female")) |
        (~sex_table.is_female & (sex_table.reported.sex != "Male")))
    mt = mt.anti_join_cols(mismatched_sex)
    mt = mt.distinct_by_col()
    mt.write(str(output_path), overwrite=True)


@cli.command("gwas-filter")
@click.argument("mt_path", type=ClickPathlibPath())
def gwas_filter(mt_path):
    mt_path = mt_path.resolve()
    mt = hl.read_matrix_table(str(mt_path))
    mt = hl.variant_qc(mt)
    mt = mt.filter_rows(mt.variant_qc.AF.all(lambda af: af > 0.01) & (mt.variant_qc.p_value_hwe > 1e-6))
    mt.write(str(mt_path.with_suffix(".GWAS_filtered.mt")), overwrite=True)


@cli.command("rare-filter")
@click.argument("mt_path", type=ClickPathlibPath())
def rare_filter(mt_path):
    mt_path = mt_path.resolve()
    mt = hl.read_matrix_table(str(mt_path))
    mt = hl.variant_qc(mt)
    mt = mt.filter_rows(mt.variant_qc.AC.any(lambda ac: ac == 1))
    mt.write(str(mt_path.with_suffix(".rare_filtered.mt")), overwrite=True)


@cli.command("ld-prune")
@click.argument("mt_path", type=ClickPathlibPath())
def ld_prune(mt_path):
    mt_path = mt_path.resolve()
    mt = hl.read_matrix_table(str(mt_path))
    pruned_variant_table = hl.ld_prune(mt.GT, r2=0.2, bp_window_size=500000)
    mt = mt.filter_rows(hl.is_defined(pruned_variant_table[mt.row_key]))
    mt.write(str(mt_path.with_suffix(".LD_pruned.mt")), overwrite=True)


@cli.command("convert-mt-to-vcf-shards")
@click.argument("mt_path", type=ClickPathlibPath())
@click.argument("original_vcf_path", type=ClickPathlibPath())
def convert_mt_to_vcf_shards(mt_path, original_vcf_path):
    original_vcf_path = original_vcf_path.resolve()
    mt_path = mt_path.resolve()
    output_vcf_dir = mt_path.with_suffix(".shards.vcf.bgz")
    shutil.rmtree(output_vcf_dir, ignore_errors=True)
    vcf_metadata = hl.get_vcf_metadata(str(original_vcf_path))
    mt = hl.read_matrix_table(str(mt_path))
    rows = mt.count_rows()
    mt = mt.repartition(1000) 
    hl.export_vcf(mt, str(output_vcf_dir), tabix=True, parallel="header_per_shard", metadata=vcf_metadata)


@cli.command("convert-mt-to-plink")
@click.argument("mt_path", type=ClickPathlibPath())
def convert_mt_to_plink(mt_path):
    mt_path = mt_path.resolve()
    mt = hl.read_matrix_table(str(mt_path))
    hl.export_plink(mt, str(mt_path.with_suffix("")), fam_id=mt.s, ind_id=mt.s)
    bimfile = pd.read_csv(mt_path.with_suffix(".bim"), delimiter="\t", header=None,
                          names=["chrom", "varid", "pos_cm", "pos_bp", "a1", "a2"],
                          dtype={"chrom": str,
                                 "varid": str,
                                 "pos_cm": float,
                                 "pos_bp": int,
                                 "a1": str,
                                 "a2": str})
    if bimfile["chrom"].str.startswith("chr").all():
        bimfile["chrom"] = bimfile["chrom"].str.slice(3)
    bimfile.to_csv(mt_path.with_suffix(".bim"), sep="\t", header=False, index=False)


@cli.command("convert-mt-to-bgen")
@click.argument("mt_path", type=ClickPathlibPath())
def convert_mt_to_plink(mt_path):
    mt_path = mt_path.resolve()
    mt = hl.read_matrix_table(str(mt_path))
    hl.export_bgen(mt, str(mt_path.with_suffix("")))


@cli.command("subset-mt-samples")
@click.argument("mt_path", type=ClickPathlibPath())
@click.argument("indiv_list", type=ClickPathlibPath())
@click.argument("out_path", type=ClickPathlibPath())
def subset_mt_samples(mt_path, indiv_list, out_path):
    mt_path = mt_path.resolve()
    mt = hl.read_matrix_table(str(mt_path))
    indiv_list_table = hl.import_table(str(indiv_list)).key_by("Subject_ID")
    mt = mt.semi_join_cols(indiv_list_table)
    mt = mt.annotate_rows(gt_stats=hl.agg.call_stats(mt.GT, mt.alleles))
    mt = mt.filter_rows(mt.gt_stats.AC.any(lambda ac: (ac > 0) & (ac < mt.gt_stats.AN)))
    mt.write(str(out_path.resolve()), overwrite=True)


@cli.command("filter-ac-by-phenotype")
@click.argument("mt_path", type=ClickPathlibPath())
@click.argument("phenotype_table_path", type=ClickPathlibPath())
@click.argument("phenotype")
def filter_ac_by_phenotype(mt_path, phenotype_table_path, phenotype):
    mt_path = mt_path.resolve()
    phenotype_table_path = phenotype_table_path.resolve()
    mt = hl.read_matrix_table(str(mt_path))
    pandas_table = pd.read_csv(phenotype_table_path, index_col=0)
    pandas_table.index.name = "Subject_ID"
    pandas_table = pandas_table[[phenotype]].dropna().reset_index()
    mt = mt.semi_join_cols(hl.Table.from_pandas(pandas_table, key="Subject_ID"))
    mt = mt.annotate_rows(gt_stats=hl.agg.call_stats(mt.GT, mt.alleles))
    mt = mt.filter_rows(mt.gt_stats.AC.any(lambda ac: (ac > 5) & (ac < mt.gt_stats.AN - 5)))
    mt.write(str(mt_path.with_suffix(f".{phenotype}_filtered.mt")), overwrite=True)


@cli.command("filter-ac")
@click.argument("mt_path", type=ClickPathlibPath())
@click.argument("ac_value", type=int)
def filter_ac(mt_path, ac_value):
    mt_path = mt_path.resolve()
    mt = hl.read_matrix_table(str(mt_path))
    mt = mt.annotate_rows(gt_stats=hl.agg.call_stats(mt.GT, mt.alleles))
    mt = mt.filter_rows(mt.gt_stats.AC.any(lambda ac: (ac > ac_value) & (ac < mt.gt_stats.AN - ac_value)))
    mt.write(str(mt_path.with_suffix(f".AC_filtered.mt")), overwrite=True)



@cli.command("run-vep")
@click.argument("mt_path", type=ClickPathlibPath())
@click.argument("extra_args", nargs=-1)
def run_vep(mt_path, extra_args):
    mt_path = mt_path.resolve()
    mt = hl.read_matrix_table(str(mt_path))
    mt = hl.vep(mt, config="/sc/arion/projects/mscic1/files/WGS/vep/vep_config_script.json")
    mt.write(str(mt_path.with_suffix(".VEP_annotated.mt")), overwrite=True)


def replace_annotation_in_filename(mt_path, suffix):
    if mt_path.stem.endswith("_annotated"):
        return mt_path.with_suffix("").with_suffix(f".{suffix}.mt")
    else:
        return mt_path.with_suffix(f".{suffix}.mt")


@cli.command("filter-lof-hc")
@click.argument("mt_path", type=ClickPathlibPath())
def filter_lof_hc(mt_path):
    mt_path = mt_path.resolve()
    mt = hl.read_matrix_table(str(mt_path))
    mt = mt.filter_rows(mt.vep.transcript_consequences.any(lambda x: x.lof == "HC"))
    outpath = replace_annotation_in_filename(mt_path, "LOF_filtered")
    mt.write(str(outpath), overwrite=True)


@cli.command("filter-functional-variation")
@click.argument("mt_path", type=ClickPathlibPath())
def filter_functional_variation(mt_path):
    mt_path = mt_path.resolve()
    mt = hl.read_matrix_table(str(mt_path))
    mt = mt.filter_rows(mt.vep.transcript_consequences.any(
        lambda x: ((
                    (x.lof == "HC") |
                        (
                            (x.sift_prediction == "deleterious") &
                            (x.polyphen_prediction == "probably_damaging")
                        )
                    ) & (x.tsl == 1))
    ))
    outpath = replace_annotation_in_filename(mt_path, "functional_filtered")
    mt.write(str(outpath), overwrite=True)


@cli.command("split-chromosomes")
@click.argument("mt_path", type=ClickPathlibPath())
@click.argument("chrom", type=str)
def split_chromosomes(mt_path, chrom):
    mt_path = mt_path.resolve()
    mt = hl.read_matrix_table(str(mt_path))
    mt_filtered = mt.filter_rows(mt.locus.contig == chrom)
    mt_filtered.write(str(mt_path.with_suffix(f".{chrom}.mt")), overwrite=True)


@cli.command("merge-chromosomes")
@click.argument("infiles", type=ClickPathlibPath(), nargs=-1)
@click.argument("outfile", type=ClickPathlibPath())
def merge_chromosomes(infiles, outfile):
    input_mts = []
    for infile in infiles:
        mt_path = infile.resolve()
        mt = hl.read_matrix_table(str(mt_path))
        input_mts.append(mt)
    outfile = outfile.resolve()
    merged_mt = hl.MatrixTable.union_rows(*input_mts)
    merged_mt.write(str(outfile), overwrite=True)


@cli.command("restrict-to-bed")
@click.argument("mt_path", type=ClickPathlibPath())
@click.argument("bed_path", type=ClickPathlibPath())
@click.argument("out_mt_path", type=ClickPathlibPath())
def restrict_to_bed(mt_path, bed_path, out_mt_path):
    interval_table = hl.import_bed(str(bed_path.resolve()), reference_genome='GRCh38')
    mt = hl.read_matrix_table(str(mt_path.resolve()))
    mt = mt.filter_rows(hl.is_defined(interval_table[mt.locus]))
    mt.write(str(out_mt_path.resolve()), overwrite=True)


@cli.command("transform-tpm-table")
@click.argument("isoform_tpm_tsv_path", type=ClickPathlibPath())
@click.argument("tx_summary_ht_path", type=ClickPathlibPath())
@click.argument("gene_maximums_ht_path", type=ClickPathlibPath())
def transform_tpm_table(isoform_tpm_tsv_path, tx_summary_ht_path, gene_maximums_ht_path):
    tx_annotation.get_gtex_summary(str(isoform_tpm_tsv_path), str(tx_summary_ht_path))
    tx_annotation.get_gene_expression(str(tx_summary_ht_path), str(gene_maximums_ht_path))


@cli.command("annotate-pext")
@click.argument("mt_path", type=ClickPathlibPath())
@click.argument("tx_summary_ht_path", type=ClickPathlibPath())
@click.argument("gene_maximums_ht_path", type=ClickPathlibPath())
def annotate_pext(mt_path, tx_summary_ht_path, gene_maximums_ht_path):
    mt_path = mt_path.resolve()
    tx_summary_ht_path = tx_summary_ht_path.resolve()
    gene_maximums_ht_path = gene_maximums_ht_path.resolve()
    mt = hl.read_matrix_table(str(mt_path))
    tx_summary_ht = hl.read_table(str(tx_summary_ht_path))
    mt = tx_annotation.tx_annotate_mt(mt, tx_summary_ht, gene_maximums_ht_path)
    outpath = replace_annotation_in_filename(mt_path, "pext_annotated")
    mt.write(str(outpath), overwrite=True)


@cli.command("filter-hi-pext")
@click.argument("mt_path", type=ClickPathlibPath())
def filter_hi_pext(mt_path):
    mt_path = mt_path.resolve()
    mt = hl.read_matrix_table(str(mt_path))
    mt = mt.filter_rows(mt.tx_annotation.any(lambda x: x.mean_proportion > 0.9))
    outpath = replace_annotation_in_filename(mt, mt_path, "pext_filtered")
    mt.write(str(outpath), overwrite=True)


@cli.command("filter-hi-pext-lof")
@click.argument("mt_path", type=ClickPathlibPath())
def filter_hi_pext_lof(mt_path):
    mt_path = mt_path.resolve()
    mt = hl.read_matrix_table(str(mt_path))
    mt = mt.filter_rows(mt.tx_annotation.any(lambda x: (x.lof == hl.literal("HC")) & (x.mean_proportion > 0.9)))
    outpath = replace_annotation_in_filename(mt_path, "pext_LOF_filtered")
    mt.write(str(outpath), overwrite=True)


if __name__ == "__main__":
    cli()
