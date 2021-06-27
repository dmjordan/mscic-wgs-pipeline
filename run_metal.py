import subprocess
import click
import os
from io import StringIO

@click.command()
@click.argument("files", nargs=-1)
@click.option("--outfile", "-o", type=click.Path(writable=True, dir_okay=False), default="METAANALYSIS.TBL")
def main(files, outfile):
    metal_script_buffer = StringIO("""
SCHEME SE
MARKER variant.id
ALLELE effect.allele other.allele
PVALUE Score.pval
EFFECT Est
STDERR Est.SE
""")
    for file in files:
        print(f"PROCESS {file}", file=metal_script_buffer)
    out_prefix, out_suffix = os.path.splitext(outfile)
    print(f"OUTFILE {out_prefix} {out_suffix}", file=metal_script_buffer)
    print("ANALYZE", file=metal_script_buffer)
    print("QUIT", file=metal_script_buffer)
    subprocess.run("/hpc/packages/minerva-centos7/metal/2018-08-28/bin/metal",
                   text=True, check=True, stdin=metal_script_buffer)

if __name__ == "__main__":
    main()