import subprocess
import os
import typing
from tempfile import NamedTemporaryFile
import shlex

import click

if typing.TYPE_CHECKING:
    from snakemake.script import Snakemake
    snakemake: Snakemake

def main(files, outfile):
    with NamedTemporaryFile(mode='w', dir='/sc/arion/scratch/jordad05', delete=False) as metalfile:

        print("""
SCHEME STDERR
MARKER variant.id
ALLELE effect.allele other.allele
PVALUE Score.pval
EFFECT Est
STDERR Est.SE
WEIGHT n.obs
""", file=metalfile)
        for file in files:
            print(f"PROCESS {file}", file=metalfile)
        out_prefix, out_suffix = os.path.splitext(outfile)
        print(f"OUTFILE {out_prefix[:-1]} {out_suffix}", file=metalfile)
        print("ANALYZE", file=metalfile)
        print("QUIT", file=metalfile)
    metalfile_name = shlex.quote(metalfile.name)
    subprocess.run(f"ml metal && metal {metalfile_name}",
                   text=True, check=True, shell=True)
    os.unlink(metalfile.name)

click_main = click.Command("run_metal", callback=main,
                           params=[click.Argument(["files"], nargs=-1),
                                   click.Option(["--outfile", "-o"],
                                                type=click.Path(writable=True, dir_okay=False),
                                                default="METAANALYSIS.TBL")])

if __name__ == "__main__":
    try:
        main(snakemake.input, snakemake.output[0])
    except NameError:
        click_main()
