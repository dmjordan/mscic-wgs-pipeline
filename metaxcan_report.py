import pandas as pd
import typing
import os
from tqdm import tqdm
import statsmodels.api as sm

if typing.TYPE_CHECKING:
    from snakemake.script import Snakemake
    snakemake: Snakemake

columns_to_extract = ["gene", "gene_name", "pvalue", "fdr_tissue", "bonferroni_tissue", "phenotype", "tissue"]
tables = []
for infile in tqdm(snakemake.input):
    phenotype, tissue, _csv = os.path.basename(infile).split(".")
    if tissue == "smultixcan":
        tissue = "Meta"
    table = pd.read_csv(infile, sep=None, engine="python")
    table["phenotype"] = phenotype
    table["tissue"] = tissue

    table["fdr_tissue"] = table["fdr"]
    table["bonferroni_tissue"] = table["bonferroni"]
    table = table[columns_to_extract]
    tables.append(table)

full_table = pd.concat(tables)
for phenotype in full_table["phenotype"].unique():
    phenotype_select = (full_table["phenotype"] == phenotype) & (full_table["tissue"] != "Meta")
    fdr = sm.stats.multipletests(full_table.loc[phenotype_select, "pvalue"], method='fdr_bh')[1]
    bonferroni = sm.stats.multipletests(full_table.loc[phenotype_select, "pvalue"], method='bonferroni')[1]
    full_table.loc[phenotype_select, "fdr_phenotype"] = fdr
    full_table.loc[phenotype_select, "bonferroni_phenotype"] = bonferroni

full_table.to_csv(snakemake.output[0], sep="\t", header=True, index=False)



