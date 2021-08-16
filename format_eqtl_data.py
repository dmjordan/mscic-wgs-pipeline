import pandas as pd
import sys
import typing

import logging

logging.basicConfig(level="INFO")

# declaring the snakemake object so my IDE stops yelling at me
if typing.TYPE_CHECKING:
    from snakemake.script import Snakemake
    snakemake: Snakemake

try:
    bim_path, eqtl_path = snakemake.input
    out_path = snakemake.output[0]
except NameError:
    bim_path, eqtl_path, out_path = sys.argv[1:]

logging.info(f"loading alleles from {bim_path}")
bimfile = pd.read_csv(bim_path, sep="\t", header=None,
                      names=["chrom", "snpid", "cm", "pos", "A1", "A2"],
                      dtype={"chrom": str, "snpid": str,
                             "cm": float, "pos": int,
                             "A1": str, "A2": str},
                      index_col=["chrom", "pos"])
logging.info(f"loaded {len(bimfile)} alleles")
bimfile = bimfile.loc[bimfile.index.drop_duplicates(keep=False)]  # drop multiallelic alleles
logging.info(f"{len(bimfile)} alleles remain after filtering to biallelic")

logging.info(f"loading eqtls from {eqtl_path}")
eqtl_file = pd.read_csv(eqtl_path, sep=" ", header=None,
                         names=["pid", "chrom", "start", "end", "strand", "NVariants",
                                "distToTopVar", "IDTopVar", "chromTopVar", "startTopVar", "endTopVar",
                                "nominalPVal", "regressionSlope", "flag"],
                         dtype={"pid": str, "chrom": str, "start": int, "end": int, "strand": str,
                                "NVariants": int, "distToTopVar": int, "IDTopVar": str, "chromTopVar": str,
                                "startTopVar": int, "endTopVar": int, "nominalPVal": float, "regressionSlope": float,
                                "flag": int})

merged_eqtls = eqtl_file.merge(bimfile, left_on=("chromTopVar", "startTopVar"), right_index=True,
                               suffixes=(None, "_eqtl"))
logging.debug(f"matched {len(merged_eqtls)}/{len(eqtl_file)}")
merged_eqtls = merged_eqtls.rename(columns={  "pid": "ProbeID",
                                                "snpid": "SNPID",
                                                "chromTopVar": "CHR",
                                                "startTopVar": "POS",
                                                "regressionSlope": "BETA",
                                                "nominalPVal": "PVAL"})
merged_eqtls["CHR"] = merged_eqtls.CHR.str.add_prefix("chr")
merged_eqtls = merged_eqtls[["ProbeID", "SNPID", "CHR", "POS", "A1", "A2", "BETA", "PVAL"]]

logging.info(f"writing to {out_path}")
merged_eqtls.to_csv(out_path, header=True, index=False, sep=" ")
logging.info(f"done!")
