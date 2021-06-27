import pandas as pd
import typing
import os
from tqdm import tqdm

if typing.TYPE_CHECKING:
    from snakemake.script import Snakemake
    snakemake: Snakemake

columns_to_extract = ["gene", "gene_name", "pvalue", "fdr", "bonferroni", "phenotype", "tissue"]
tables = []
for infile in tqdm(snakemake.input):
    phenotype, tissue, _csv = os.path.basename(infile).split(".")
    if tissue == "smultixcan":
        tissue = "Meta"
    table = pd.read_csv(infile, sep=None, engine="python")
    table["phenotype"] = phenotype
    table["tissue"] = tissue
    table = table.loc[table.fdr < 0.1, columns_to_extract]
    tables.append(table)

pd.concat(tables).to_csv(snakemake.output[0], sep="\t", header=True, index=False)

