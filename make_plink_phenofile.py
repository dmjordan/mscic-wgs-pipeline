from utils import *
import click
import pandas as pd

@click.command()
@click.option("--fam-file", "-f", type=click.File("r"), default="-")
@click.option("--covariate-file", "-c", type=click.Path(exists=True, dir_okay=False), default=CLINICAL_COVARIATES_PATH)
@click.option("--pca-file", "-p", type=click.Path(exists=True, dir_okay=False), default=None)
@click.option("--id-column", "-I", type=str, default="Corrected_Blood_Sample")
@click.option("--covariates/--outcomes", default=False)
@click.option("--output-file", "-o", type=click.File('w'), default='-')
def main(fam_file, covariate_file, id_column, pca_file, covariates, output_file):
    fam_table = load_plink_table(fam_file, header=None, names=["FID", "IID", "Paternal_ID", "Maternal_ID", "Sex", "Phenotype"])
    covariate_table = load_covariates_file_time_index(covariate_file)

    index_column = grouped_unconflicting_values(\
            covariate_table.loc[covariate_table[id_column].isin(fam_table.index), id_column], 
            raise_on_conflict=True)

    sex_dtype = pd.api.types.CategoricalDtype(categories=["Male", "Female"])
    sex_raw = grouped_unconflicting_values(covariate_table["Sex"]).astype(sex_dtype)
    sex_formatted = sex_raw.cat.codes + 1
    age = grouped_unconflicting_values(covariate_table["Age"]).fillna(-9).astype(int)
    
    grouped = covariate_table.groupby(level=0)
    ever_positive = grouped_ever_true(covariate_table["COVID19_Ever_Positive"])
    min_antibody_titer = grouped["COVID19_Antibody_Titer"].min()
    max_cn1 =  grouped["Viral_Load_CN1_per_mL"].max()
    max_cn2 =  grouped["Viral_Load_CN2_per_mL"].max()
    max_cn3 =  grouped["Viral_Load_CN3_per_mL"].max()
    max_sofa_score = grouped["SOFA_Score"].max()
    ever_icu = grouped_ever_true(covariate_table["ICU_Status"])
    covid_ever_icu = ever_icu.mask(ever_positive != 2, 0)
    
    severity_dtype = pd.api.types.CategoricalDtype(categories=["Moderate COVID-19", "Severe COVID-19", "Severe COVID-19 with EOD"], ordered=True)
    severity = covariate_table.Severity.astype(severity_dtype)
    ever_severe = grouped_ever_true(severity >= "Severe COVID-19")
    ever_eod = grouped_ever_true(severity == "Severe COVID-19 with EOD")
    covid_severe = ever_severe.mask(ever_positive != 2, 0)
    covid_eod = ever_eod.mask(ever_positive != 2, 0)

    if covariates:
        aggregated_table =  pd.DataFrame({ "sample_id": index_column,
                                          "sex": sex_formatted,
                                          "age": age })
        aggregated_table = aggregated_table.dropna(subset=["sample_id"]).set_index("sample_id")
        if pca_file is not None:
            pca_table = load_plink_table(pca_file)
            aggregated_table = pd.concat([aggregated_table, pca_table], axis=1)

    else:        
        aggregated_table = pd.DataFrame({ "sample_id": index_column,
                                          "fatid": 0,
                                          "matid": 0,
                                          "sex": sex_formatted,
                                          "age": age,
                                          "covid_positive": ever_positive,
                                          "min_antibody_titer": min_antibody_titer,
                                          "max_cn1": max_cn1,
                                          "max_cn2": max_cn2,
                                          "max_cn3": max_cn3,
                                          "max_sofa": max_sofa_score,
                                          "ever_icu": ever_icu,
                                          "ever_icu_positive_only": covid_ever_icu,
                                          "ever_severe": ever_severe,
                                          "ever_severe_positive_only": covid_severe,
                                          "ever_eod": ever_eod,
                                          "ever_eod_positive_only": covid_eod })
                                          
        aggregated_table = aggregated_table.dropna(subset=["sample_id"]).set_index("sample_id")
    write_plink_table(aggregated_table, output_file)


if __name__ == "__main__":
    main()
