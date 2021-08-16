import pandas as pd
import sys
import typing

from tqdm import tqdm
import logging

logging.basicConfig(level="DEBUG")

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

logging.info(f"streaming eqtls from {eqtl_path}")
chunks = []
row_count = 0
matched_count = 0
for chunk in tqdm(pd.read_csv(eqtl_path, sep=" ", header=None,
                         names=["pid", "chrom", "start", "end", "strand", "NVariants",
                                "distToTopVar", "IDTopVar", "chromTopVar", "startTopVar", "endTopVar",
                                "nominalPVal", "regressionSlope", "flag"],
                         dtype={"pid": str, "chrom": str, "start": int, "end": int, "strand": str,
                                "NVariants": int, "distToTopVar": int, "IDTopVar": str, "chromTopVar": str,
                                "startTopVar": int, "endTopVar": int, "nominalPVal": float, "regressionSlope": float,
                                "flag": int},
                         chunksize=10000), desc="streaming eqtls"):
    chunk = typing.cast(pd.DataFrame, chunk)
    row_count += len(chunk)
    merged_chunk = chunk.merge(bimfile, left_on=("chromTopVar", "startTopVar"), right_index=True,
                               suffixes=(None, "_eqtl"))
    matched_count += len(merged_chunk)
    logging.debug(f"matched {len(merged_chunk)}/{len(chunk)} in chunk, cumulative {matched_count}/{row_count}")
    merged_chunk = merged_chunk.reset_index().rename(columns={  "pid": "ProbeID",
                                                                "snpid": "SNPID",
                                                                "chrom": "CHR",
                                                                "pos": "POS",
                                                                "regressionSlope": "BETA",
                                                                "nominalPVal": "PVAL"})[
        ["ProbeID", "SNPID", "CHR", "POS", "A1", "A2", "BETA", "PVAL"]
    ]
    chunks.append(merged_chunk)
logging.info(f"read {row_count} records in {len(chunks)} chunks")

result = pd.concat(chunks, ignore_index=True)
logging.info(f"kept {len(result)} records after filters")
logging.info(f"writing to {out_path}")
result.to_csv(out_path)
logging.info(f"done!")