import numpy as np
import pandas as pd

import typing

# declaring the snakemake object so my IDE stops yelling at me
if typing.TYPE_CHECKING:
    from snakemake.script import Snakemake
    snakemake: Snakemake

def ever_decreased(series):
    """Calculates whether the score ever decreased during the series."""
    return series.pct_change().lt(0).any()


def ever_increased(series):
    """Calculates whether the score ever increased during the series"""
    return series.pct_change().gt(0).any()

raw_table = pd.read_csv(snakemake.input.table, parse_dates=["Encounter_Start_Date", "Encounter_End_Date"],
                        dtype = { 'Severity': pd.api.types.CategoricalDtype(['Moderate COVID-19',
                                                                                               'Severe COVID-19',
                                                                                               'Severe COVID-19 with EOD'],
                                                                                              ordered=True),
                                  'COVID19_Positive': 'boolean'})

raw_table["Encounter_Length"] = raw_table.apply(lambda x: x.Encounter_End_Date - x.Encounter_Start_Date, axis="columns")

clinical_table = raw_table.groupby("REGENID").agg(
    age=("AGE", "max"),
    sex=("GENDER", "first"),
    max_who=("WHO_Ordinal", "max"),
    max_who_minus_death=("WHO_Ordinal", lambda x: x.loc[x < 8].max()),
    who_ever_decreased=("WHO_Ordinal", ever_decreased),
    who_ever_increased=("WHO_Ordinal", ever_increased),
    severity_ever_decreased=("Severity", lambda x: ever_decreased(x.cat.codes)),
    severity_ever_increased=("Severity", lambda x: ever_increased(x.cat.codes)),
    severity_ever_severe=("Severity", lambda x: x.eq("Severe COVID-19").any()),
    severity_ever_eod=("Severity", lambda x: x.eq("Severe COVID-19 with EOD").any()),
    max_severity_moderate=("Severity", lambda x: x.max() == "Moderate COVID-19")
)

clinical_table["covid_encounter_days"] = raw_table.loc[raw_table.COVID19_Positive].groupby("REGENID").apply(lambda x: x.Encounter_End_Date.min() - x.Encounter_Start_Date.min())
clinical_table["covid_encounter_days"] = clinical_table.covid_encounter_days / np.timedelta64(1, "D")
clinical_table["covid_encounter_days"] = clinical_table.covid_encounter_days.where(clinical_table.covid_encounter_days <= 15)
clinical_table["covid_encounter_days_log"] = clinical_table.covid_encounter_days.transform("log1p")

# invert sense of "ever_decreased" so that cases remain bad
clinical_table["severity_ever_decreased"] = ~(clinical_table.severity_ever_decreased.astype("boolean"))
clinical_table["who_ever_decreased"] = ~(clinical_table.who_ever_decreased.astype("boolean"))

clinical_table["sex"] = clinical_table.sex.str.get(0)
clinical_table["age_sex"] = clinical_table["age"] * (clinical_table["sex"] == "F")
clinical_table["age_squared"] = clinical_table["age"] ** 2

excluded_samples = pd.read_csv(snakemake.input.excluded_ids)
clinical_table = clinical_table.drop(clinical_table.index.intersection(excluded_samples))

pca_table = pd.read_csv(snakemake.input.pcair, delim_whitespace=True,
                        index_col="Subject_ID")
pca_table.columns = [f"pc{i}" for i in range(1, len(pca_table.columns) + 1)]
pca_table = pca_table.iloc[:, :10]

clinical_table = clinical_table.join(pca_table, how="inner")

clinical_table.to_csv(snakemake.output[0])
