import click, pathlib
from sklearn.cluster import KMeans
import pandas as pd
from utils import *

@click.command()
@click.option("--covariate-file", "-c", type=click.Path(exists=True, dir_okay=False), 
        default=CLINICAL_COVARIATES_PATH)
@click.option("--n-clusters", "-k", type=int, default=4)
@click.option("--n-pcs", "-n", type=int, default=3)
@click.option("--output-prefix", "-o", type=click.Path(writable=True), default="pca_race_cluster")
@click.argument("pcafile", type=click.File('r'))
def main(covariate_file, n_clusters, n_pcs, output_prefix, pcafile):
    pca_table = load_plink_table(pcafile, double_id=True)
    covariate_table = load_covariates_file_simple_index(covariate_file, index=pca_table.index)

    kmeans = KMeans(n_clusters).fit(pca_table[[f"U{i}" for i in range(1, n_pcs+1)]])
    
    breakpoint()
    race_clusters = pd.get_dummies(covariate_table.Race_From_Consent)\
                      .groupby(kmeans.labels_).sum().idxmax(axis=1)
    ethnicity_clusters = pd.get_dummies(covariate_table.Ethnicity_From_Consent)\
                      .groupby(kmeans.labels_).sum().idxmax(axis=1)
    cluster_table = pd.DataFrame({'race': race_clusters, 'ethnicity': ethnicity_clusters})

    cluster_table["name"] = cluster_table.race.str.split(" ").str.get(0)
    cluster_table.loc[cluster_table.ethnicity == "Hispanic or Latino", "name"] = "HispanicLatino"
    assert cluster_table.name.is_unique

    for i, index in pca_table.groupby(kmeans.labels_).groups.items():
        outfile = f"{output_prefix}.{cluster_table.loc[i,'name']}.list"
        index.to_series(name='IID').to_csv(outfile, sep="\t", header=True, index_label="FID")
        
if __name__ == "__main__":
    main()
