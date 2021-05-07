from utils import *
import click
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

@click.command()
@click.option("--savefig", type=click.Path(dir_okay=False, writable=True), default=None)
@click.option("--covariate-file", "-c", type=click.Path(exists=True, dir_okay=False), 
        default=CLINICAL_COVARIATES_PATH)
@click.option("--n-pcs", "-n", type=int, default=4)
@click.argument("pcafile", type=click.File('r'))
def main(savefig, covariate_file, n_pcs, pcafile):
    pca_table = load_plink_table(pcafile, double_id=True)
    covariate_table = load_covariates_file_simple_index(covariate_file, index=pca_table.index)
    pca_table["Race"] = covariate_table.Race_From_Consent.where(covariate_table.Ethnicity_From_Consent == "Hispanic or Latino", "Hispanic or Latino")
    
    g = sns.pairplot(pca_table, hue="Race", vars=[f"U{i}" for i in range(1, n_pcs+1)])
    if savefig:
        g.savefig(savefig)
    else:
        plt.show()

if __name__ == "__main__":
    main()
