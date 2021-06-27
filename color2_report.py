import pandas as pd
import typing
import os

if typing.TYPE_CHECKING:
    from snakemake.script import Snakemake
    snakemake: Snakemake

tables = []
for infile in snakemake.input:
    phenotype, tissue, _full_table, _csv = os.basename(infile).split(".")
    table = pd.read_csv(infile, sep="\t")
    table["phenotype"] = phenotype
    table["tissue"] = tissue
    table = table.loc[(table["PP.H4.abf"] > 0.9) & (table["min.pval.biom"] < 0.05)]
    tables.append(table)
pd.concat(tables).to_csv(snakemake.output[0], sep="\t", header=True, index=False)
