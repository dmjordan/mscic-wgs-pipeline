# The contents of this file are stolen from the https://github.com/macarthur-lab/tx_annotation
# and https://github.com/broadinstitute/gnomad_methods repositories

import hail as hl

# from gnomad_methods/gnomad/utils/vep.py

# Note that this is the current as of v81 with some included for backwards compatibility (VEP <= 75)
CSQ_CODING_HIGH_IMPACT = [
    "transcript_ablation",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
    "stop_lost",
]

CSQ_CODING_MEDIUM_IMPACT = [
    "start_lost",  # new in v81
    "initiator_codon_variant",  # deprecated
    "transcript_amplification",
    "inframe_insertion",
    "inframe_deletion",
    "missense_variant",
    "protein_altering_variant",  # new in v79
    "splice_region_variant",
]

CSQ_CODING_LOW_IMPACT = [
    "incomplete_terminal_codon_variant",
    "start_retained_variant",  # new in v92
    "stop_retained_variant",
    "synonymous_variant",
    "coding_sequence_variant",
]


CSQ_NON_CODING = [
    "mature_miRNA_variant",
    "5_prime_UTR_variant",
    "3_prime_UTR_variant",
    "non_coding_transcript_exon_variant",
    "non_coding_exon_variant",  # deprecated
    "intron_variant",
    "NMD_transcript_variant",
    "non_coding_transcript_variant",
    "nc_transcript_variant",  # deprecated
    "upstream_gene_variant",
    "downstream_gene_variant",
    "TFBS_ablation",
    "TFBS_amplification",
    "TF_binding_site_variant",
    "regulatory_region_ablation",
    "regulatory_region_amplification",
    "feature_elongation",
    "regulatory_region_variant",
    "feature_truncation",
    "intergenic_variant",
]

CSQ_ORDER = (
    CSQ_CODING_HIGH_IMPACT
    + CSQ_CODING_MEDIUM_IMPACT
    + CSQ_CODING_LOW_IMPACT
    + CSQ_NON_CODING
)


def add_most_severe_consequence_to_consequence(
    tc: hl.expr.StructExpression,
) -> hl.expr.StructExpression:
    """
    Add most_severe_consequence annotation to transcript consequences.
    This is for a given transcript, as there are often multiple annotations for a single transcript:
    e.g. splice_region_variant&intron_variant -> splice_region_variant
    """
    csqs = hl.literal(CSQ_ORDER)

    return tc.annotate(
        most_severe_consequence=csqs.find(lambda c: tc.consequence_terms.contains(c))
    )



# from tx_annotation/tx_annotation_resources.py

# CSQ terms
lof_csqs = ["stop_gained","splice_donor_variant", "splice_acceptor_variant","frameshift_variant"]
missense_csqs = ["missense_variant", "inframe_insertion", "inframe_deletion"]
syn_csqs = ["synonymous_variant"]
all_coding_csqs = CSQ_CODING_HIGH_IMPACT + CSQ_CODING_MEDIUM_IMPACT + CSQ_CODING_LOW_IMPACT


def get_gtex_summary(gtex_rsem_path, gtex_tx_summary_out_path, get_medians=True, make_per_tissue_file=None):
    """
    Get GTEx RSEM table with ENSTs and ENSGs as rows and GTEx samples as columns (e.g. Muscle-Skeletal.12,
    Adipose.27 etc.) and write out a table with same rows, and tissues as columns (Muscle-Skeletal, Adipose etc.)
    with cells representing summary expression of transcripts across tissues (ie. mean or median).
    :param str gtex_rsem_path: Output of RSEM quantifications from GTEx
    Example: "gs://gnomad-berylc/tx-annotation/reheadered.031216.GTEx_Analysis_2016-09-07_RSEMv1.2.22_transcript_tpm.txt.bgz"
    :param str gtex_tx_summary_out_path: Path to write out.
    Example: "gs://gnomad-berylc/tx-annotation/hail2/GTEx.V7.tx_medians.030818.mt"
    :param bool get_medians: Default True. If False, returns mean transcript expression per tissue
    :return: Writes out summarized GTEx transcript expression as Table.
    :rtype: None
    """
    gtex = hl.import_matrix_table(gtex_rsem_path, row_key='transcript_id',
                                  row_fields={'transcript_id': hl.tstr, 'gene_id': hl.tstr}, entry_type=hl.tfloat64,
                                  force_bgz=True)

    gtex = gtex.annotate_cols(tissue=gtex.col_id.split("\\.")[0])

    # Will come in handy when it's time to splat
    tissue_ids = sorted(gtex.aggregate_cols(hl.agg.collect_as_set(gtex.tissue)))
    d = {tiss: i for i, tiss in enumerate(tissue_ids)}

    if get_medians:
        gtex = gtex.group_cols_by(gtex.tissue).aggregate(tx_expr_summary=hl.median(hl.agg.collect(gtex.x)))
    else:
        gtex = gtex.group_cols_by(gtex.tissue).aggregate(tx_expr_summary=hl.mean(hl.agg.collect(gtex.x)))

    # Make a new column as an array of the values across tissues (per transcript)
    # Added sorting to ensure match up, a bit more wordy then before but oh well.
    gtex_table = gtex.entries()
    gtex_table = gtex_table.key_by(gtex_table.transcript_id, gtex_table.gene_id)
    gtex_table = gtex_table.collect_by_key()
    gtex_table = gtex_table.annotate(agg_expression=hl.sorted(gtex_table.values, key=lambda x: x.tissue).
                                     map(lambda x: x.tx_expr_summary))

    # Strip version numbers
    gtex_table = gtex_table.key_by(transcript_id=gtex_table.transcript_id.split("\\.")[0],
                                   gene_id=gtex_table.gene_id.split("\\.")[0])

    # Write out ht
    gtex_table.write(gtex_tx_summary_out_path, overwrite=True)

    if make_per_tissue_file:
        gtex_table_modified = gtex_table.annotate(**{
            tissue_id.replace("-", "_").replace(" ", "_").replace("(", "_").replace(")", "_"):
                gtex_table.agg_expression[d[tissue_id]] for tissue_id in tissue_ids})

        gtex_table_modified = gtex_table_modified.drop(gtex_table_modified.values)

        gtex_table_modified.export(make_per_tissue_file)



def get_gene_expression(gtex_tx_summary_path, gene_expression_out):
    gtex = hl.read_table(gtex_tx_summary_path)
    gene_expression = gtex.group_by(ensg=gtex.gene_id).aggregate(gene_expression=hl.agg.array_sum(gtex.agg_expression))
    tissue_ids = sorted([y.tissue for y in gtex.values.take(1)[0]])
    d = {tiss: i for i, tiss in enumerate(tissue_ids)}
    gene_expression = (gene_expression.annotate(
        gene_maximum_dict={tissue_id.replace("-", "_").replace(" ", "_").replace("(", "_").replace(")", "_"):
                               gene_expression.gene_expression[d[tissue_id]] for tissue_id in tissue_ids}))

    gene_expression.show(10)
    gene_expression.write(gene_expression_out, overwrite=True)






# from tx_annotation/tx_annotation.py

def filter_table_to_gene_list(mt_kt, genes, gene_column_in_mt_kt):
    """Take a matrix table and return a table filtered down to a set of genes
    :param Table mt_kt:
    :param list of str or set of str genes: Genes of interest to which to filter table
    :param str gene_column_in_mt_kt: Column in matrix table that contains gene information within
    vep.transcript_consequences. often ["gene_id", "gene_symbol"]
    :return: Filtered table
    :rtype: Table
    """
    gene_names = hl.literal(genes)

    mt_kt = mt_kt.annotate(
        in_gene_of_interest=gene_names.find(lambda x: mt_kt.vep.transcript_consequences[gene_column_in_mt_kt] == x))

    mt_kt = mt_kt.filter(mt_kt.in_gene_of_interest != "NA")

    return mt_kt


def filter_table_to_csqs(mt_kt, csqs):
    """Take a matrix table and return a table filtered down to a set of CSQs
    :param Table mt_kt:
    :param list of str or set of str csqs: CSQs of interest to which to filter table
    :return: Filtered matrix table
    :rtype: Table
    """
    csqs = hl.literal(csqs)

    mt_kt = mt_kt.annotate(
        in_csq_of_interest=csqs.find(lambda x: mt_kt.vep.transcript_consequences.most_severe_consequence == x))
    mt_kt = mt_kt.filter(mt_kt.in_csq_of_interest != "NA")
    return mt_kt


def read_tx_annotation_tables(mt_path, gtex_tx_summary_path, mt_type="mt"):
    if mt_type == "mt":
        mt = hl.read_matrix_table(mt_path)

    elif mt_type == "ht":
        mt = hl.read_table(mt_path)
        mt = hl.MatrixTable.from_rows_table(mt)

    gtex = hl.read_table(gtex_tx_summary_path)
    return mt, gtex


def tx_annotate_mt(mt, gtex, gene_maximums_ht_path,
                   tx_annotation_type = "proportion", tissues_to_filter = None,
                   filter_to_csqs=all_coding_csqs, filter_to_genes=None, gene_column_in_mt=None, filter_to_homs=False,
                   out_tx_annotation_tsv=None, out_tx_annotation_ht=None):


    """
    Annotate variants in the input MatrixTable with transcript-based expression values accross GTEx. Returns Table.
    :param MatrixTable mt: Input variant file
    :param MatrixTable gtex: Input GTEx summary MatrixTable, must have transcript_id column to key by
    :param str tx_annotation_type: One of ["expression", "proportion"]. Select proportion if you'd like the
    tx_annotation values to be normalized by max expression of the gene
    :param None or list filter_to_csqs: Default None. If you'd like to filter the mt before annotating
    (decreases time) feed in a list or set of consequence terms.
    :param str gene_column_in_mt: Must be set if filter_to_genes != None.
    Column in matrix table that contains gene information within vep.transcript_consequences.
    often ["gene_id", "gene_symbol"]
    :param None or list filter_to_csqs: Default None. If you'd like to filter the mt before annotating
    (decreases time) feed in a list or set of consequence terms.
    Example = ["stop_gained","splice_donor_variant", "splice_acceptor_variant","frameshift_variant"]
    :param None or str out_tx_annotation_tsv: Default None.
    If you'd like to write out the results table as a tsv, provide a tsv path
    :param None or str out_tx_annotation_ht: Default None.
    If you'd like to write out the results table as a Hail 0.2 table, provide a .ht path
    :param bool filter_to_homs: Default False
    If True, filter to variants with at least one homozygote in dataset
    :return: Table with columns: variant, worst_csq, ensg, LOFTEE LOF, LOFTEE LOF Flag, transcript-aware expression
    by GTEx Tissue
    :rtype: Table with variants annotated with transcript-aware tissue expression
    """

    #check_inputs(**locals())

    gtex_table = gtex.key_by("transcript_id")

    #mt = process_consequences(mt, penalize_flags=False)
    mt_exploded = mt.distinct_by_row()
    mt_exploded = mt_exploded.annotate_rows(vep=mt_exploded.vep.annotate(
        transcript_consequences=mt_exploded.vep.transcript_consequences.map(add_most_severe_consequence_to_consequence)))

    # Explode the mt for the transcript consequences to be able to key by transcript ID
    mt_exploded = mt_exploded.explode_rows(mt_exploded.vep.transcript_consequences)

    mt_kt = mt_exploded.rows()
    # Currently testing removal of protein coding transcripts
    mt_kt = mt_kt.filter(mt_kt.vep.transcript_consequences.biotype == "protein_coding")

    if filter_to_genes:
        print("Filtering to genes of interest")
        mt_kt = filter_table_to_gene_list(mt_kt, filter_to_genes, gene_column_in_mt)

    if filter_to_csqs:
        print("Filtering to csqs in %s" % (",".join(filter_to_csqs)))
        mt_kt = filter_table_to_csqs(mt_kt, filter_to_csqs)

    if filter_to_homs:
        print("Filtering to variants with at least 1 homozygote sample in dataset")
        #mt_kt = mt_kt.filter(mt_kt.info.Hom[mt_kt.a_index - 1] > 0)
        idx = mt_kt.globals.freq_index_dict['gnomad']
        mt_kt = mt_kt.filter(mt_kt.freq[idx].homozygote_count >= 1)

    # Annotate mt with the gtex values (ie. join them)
    mt_kt = mt_kt.annotate(tx_data=gtex_table[mt_kt.vep.transcript_consequences.transcript_id])

    # Group by gene, worst_csq and variant, and do a pairwise-sum
    grouped_table = (
        mt_kt.group_by(csq=mt_kt.vep.transcript_consequences.most_severe_consequence,
                       ensg=mt_kt.vep.transcript_consequences.gene_id,
                       symbol=mt_kt.vep.transcript_consequences.gene_symbol,
                       locus=mt_kt.locus,
                       alleles=mt_kt.alleles,
                       lof=mt_kt.vep.transcript_consequences.lof,
                       lof_flag=mt_kt.vep.transcript_consequences.lof_flags).aggregate(tx_annotation=hl.agg.array_sum(mt_kt.tx_data.agg_expression)))

    # Expand the columns from the arrays and add tissues as headers
    #tissue_ids = gtex.tissue.collect()
    # Since gtex no longer has .tissue just a new way to do this, i probably want to save it as a global at some point
    tissue_ids = sorted([y.tissue for y in gtex.values.take(1)[0]])
    d = {tiss: i for i, tiss in enumerate(tissue_ids)}

    tx_annotation_table = grouped_table.annotate(
        **{tissue_id.replace("-", "_").replace(" ", "_").replace("(", "_").replace(")", "_"):
               grouped_table.tx_annotation[d[tissue_id]] for tissue_id in tissue_ids})

    tx_annotation_table = tx_annotation_table.drop(tx_annotation_table.tx_annotation)

    # First of all do you want proportions or expression?
    if tx_annotation_type == "proportion":
        print("Returning expression proportion")
        gene_maximums_ht = hl.read_table(gene_maximums_ht_path)
        tx_annotation_table = get_expression_proportion(tx_annotation_table, tissues_to_filter, gene_maximums_ht)

    #You can write the output that is exploded by variants-ensg-csq-symbol-LOFTEE-LOFTEEflag
    # and has a value for each tissue as column, either as a TSV or a KT

    if out_tx_annotation_tsv:
        print("Writing tsv file to %s" %out_tx_annotation_tsv)
        tx_annotation_table.export(out_tx_annotation_tsv)

    if out_tx_annotation_ht:
        print("Writing Table to %s" % out_tx_annotation_ht)
        tx_annotation_table.write(out_tx_annotation_ht)

    tx_annotation_table = tx_annotation_table.key_by(tx_annotation_table.locus, tx_annotation_table.alleles)
    tx_annotation_table = tx_annotation_table.collect_by_key('tx_annotation')
    mt = mt.annotate_rows(**tx_annotation_table[mt.locus, mt.alleles])

    return mt

def get_expression_proportion(tx_table, tissues_to_filter, gene_maximum_ht):

    if tissues_to_filter:
        print("Filtering tissues:", tissues_to_filter)
        tx_table = tx_table.drop(*tissues_to_filter)

    remaining_tissue_columns = list(
        set(tx_table.row) - {'locus', 'alleles','csq', 'ensg', 'symbol','lof', 'lof_flag'})

    tx_table = tx_table.annotate(
        tx_expression=
        {tissue_id: tx_table[tissue_id] for tissue_id in remaining_tissue_columns})

    tx_table = tx_table.key_by('ensg').join(gene_maximum_ht.key_by("ensg"))

    expression_proportion_table = tx_table.annotate(
        expression_proportion_dict=
        {tissue_id: tx_table.tx_expression[tissue_id] / tx_table.gene_maximum_dict[tissue_id] for tissue_id in remaining_tissue_columns})

    columns_to_drop = list(set(expression_proportion_table.row) - {'locus', 'alleles','csq', 'ensg',
                                                                   'symbol','lof', 'lof_flag','expression_proportion_dict'})

    expression_proportion_table = expression_proportion_table.drop(*columns_to_drop)

    expression_proportion_table = expression_proportion_table.annotate(
        **{tissue_id: expression_proportion_table.expression_proportion_dict[tissue_id] for tissue_id in remaining_tissue_columns})

    expression_proportion_table = expression_proportion_table.annotate(
        mean_proportion=hl.mean( hl.filter(
            lambda e : ~hl.is_nan(e),  [expression_proportion_table[tissue_id] for tissue_id in remaining_tissue_columns]), filter_missing=True))

    expression_proportion_table = expression_proportion_table.drop(
        expression_proportion_table.expression_proportion_dict).key_by(
        'locus', 'alleles', 'ensg')

    return expression_proportion_table
