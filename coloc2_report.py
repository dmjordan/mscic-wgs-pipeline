import pandas as pd
import typing
import os
from tqdm import tqdm

if typing.TYPE_CHECKING:
    from snakemake.script import Snakemake
    snakemake: Snakemake

tables = []
for infile in tqdm(snakemake.input):
    phenotype, tissue, _full_table, _csv = os.path.basename(infile).split(".")
    table = pd.read_csv(infile, sep="\t")
    table["phenotype"] = phenotype
    table["tissue"] = tissue
    table = table.loc[(table["PP.H4.abf"] > 0.2) & (table["min.pval.biom"] < 1e-5)]
    tables.append(table)
pd.concat(tables).to_csv(snakemake.output[0], sep="\t", header=True, index=False)
