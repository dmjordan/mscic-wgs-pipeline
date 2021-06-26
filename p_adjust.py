import click
import pandas as pd
import statsmodels.api as sm


@click.command()
@click.argument("infile", type=click.File('r'))
@click.option("--delimiter", "-d", default=None)
@click.option("--outfile", "-o", type=click.File('w'), default=None)
def main(infile, delimiter, outfile):
    table = pd.read_csv(infile, sep=delimiter, engine='python')
    table['fdr'] = sm.stats.multipletests(table.pvalue, method='fdr_bh')[1]
    table['bonferroni'] = sm.stats.multipletests(table.pvalue, method='bonferroni')[1]
    infile.close()

    if outfile is None:
        outfile = open(infile.name, 'w')
    if delimiter is None:
        if outfile.name.endswith(".csv"):
            delimiter = ","
        else:
            delimiter = "\t"
    table.to_csv(outfile, sep=delimiter, header=True, index=False)


if __name__ == "__main__":
    main()
