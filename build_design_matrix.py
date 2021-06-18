import matplotlib
matplotlib.use("agg")
import pandas as pd
import numpy as np
from scipy import stats
import seaborn as sns
import matplotlib.pyplot as plt

binary_phenotypes = ["ever_covid", "covid_encounter", "ever_icu", "deceased", "recovered", "discharged",
                     "recovered_not_deceased", "deceased_vs_recovered", "deceased_vs_discharged",
                     "sofa_ever_decreased", "sofa_ever_increased", "sofa_ever_decreased_icu_only",
                     "sofa_ever_increased_icu_only", "who_ever_decreased", "who_ever_increased",
                     "severity_ever_moderate", "severity_ever_severe", "severity_ever_eod", "max_severity_moderate",
                     "severity_ever_decreased", "severity_ever_increased", "severity_ever_decreased_counting_discharge",
                     "severity_ever_increased_counting_death",
                     "blood_viral_load_detected", "blood_viral_load_detected_severe_only"]
continuous_phenotypes = ["max_sofa", "max_sofa_icu_only", "max_who", "max_who_minus_death", "highest_titer",
                         "highest_titer_agm_only", "max_viral_load_cn1", "max_viral_load_cn2", "max_viral_load_cn3",
                         "days_onset_to_encounter", "covid_encounter_days", "days_onset_to_icu",
                         "days_encounter_to_icu", "icu_hours_cumulative",
                         "blood_viral_load_bitscore", "blood_viral_load_bitscore_severe_only"]
all_phenotypes = binary_phenotypes + continuous_phenotypes + \
                [f"{phenotype}_irnt" for phenotype in continuous_phenotypes] + \
                [f"{phenotype}_log" for phenotype in continuous_phenotypes] + \
                [f"{phenotype}_percentile" for phenotype in continuous_phenotypes]



def piecewise_cummax(series):
    """Calculates cumulative hours spent in ICU from the ICU_Hours_Since_Entry field."""
    cummax = series.cummax()
    value = 0
    while series.lt(cummax).any():
        split_index = np.argmax(series.lt(cummax))
        value += cummax.iloc[split_index]
        series = series[split_index:]
        cummax = series.cummax()
        value += series.max()
    return value


def ever_decreased(series):
    """Calculates whether the score ever decreased during the series"""
    return series.pct_change().lt(0).any()


def ever_increased(series):
    """Calculates whether the score ever increased during the series"""
    return series.pct_change().gt(0).any()


def build_design_matrix(covariates_path, design_matrix_path):
    all_phenotypes = []
    raw_table = pd.read_csv(covariates_path, dtype={'Severity': pd.api.types.CategoricalDtype(['Moderate COVID-19',
                                                                                               'Severe COVID-19',
                                                                                               'Severe COVID-19 with EOD'],
                                                                                              ordered=True),
                                                    'Event_Date': 'period[D]',
                                                    'Date_Of_First_Sample': 'period[D]',
                                                    'COVID19_Encounter_Start_Date': 'period[D]',
                                                    'COVID19_Encounter_End_Date': 'period[D]',
                                                    'COVID19_Symptoms_Reported_Onset_Date': 'period[D]',
                                                    'COVID19_Within_Encounter': 'boolean',
                                                    'COVID19_Deceased_Within_Encounter': 'boolean',
                                                    'COVID19_Ever_Positive': 'boolean',
                                                    'ICU_Status': 'boolean'})

    clinical_table = raw_table.groupby("Subject_ID").agg(
        age=("Age", "max"),
        sex=("Sex", "first"),
        race=("Race_From_Consent", "first"),
        ethnicity=("Ethnicity_From_Consent", "first"),
        ever_covid=("COVID19_Ever_Positive", lambda x: x.any()),
        covid_encounter=("COVID19_Within_Encounter", lambda x: x.any()),
        recruitment_date=("Date_Of_First_Sample", "first")
    )
    covid_only_table = raw_table.loc[raw_table.COVID19_Within_Encounter].groupby("Subject_ID").agg(
        ever_icu=("ICU_Status", lambda x: x.any()),
        icu_hours_cumulative=("ICU_Hours_Since_Entry", piecewise_cummax),
        max_sofa=("SOFA_Score", "max"),
        max_who=("WHO_Ordinal", "max"),
        max_who_minus_death=("WHO_Ordinal", lambda x: x.loc[x < 8].max()),
        deceased=("COVID19_Deceased_Within_Encounter", lambda x: x.any()),
        recovered=("SARSCoV2_Recovery_Status", lambda x: x.eq("Recovered").any()),
        discharged=("Encounter_Discharge_Date", lambda x: x.notnull().any()),
        highest_titer=("Serology_Highest_Titer", "max"),
        highest_titer_agm_only=("Serology_Highest_Titer_AGM_Only", "max"),
        sofa_ever_decreased=("SOFA_Score", ever_decreased),
        sofa_ever_increased=("SOFA_Score", ever_increased),
        who_ever_decreased=("WHO_Ordinal", ever_decreased),
        who_ever_increased=("WHO_Ordinal", ever_increased),
        severity_ever_decreased=("Severity", lambda x: ever_decreased(x.cat.codes)),
        severity_ever_increased=("Severity", lambda x: ever_increased(x.cat.codes)),
        severity_ever_moderate=("Severity", lambda x: x.eq("Moderate COVID-19").any()),
        severity_ever_severe=("Severity", lambda x: x.eq("Severe COVID-19").any()),
        severity_ever_eod=("Severity", lambda x: x.eq("Severe COVID-19 with EOD").any()),
        max_severity_moderate=("Severity", lambda x: x.max() == "Moderate COVID-19"),
        max_viral_load_cn1=("Viral_Load_CN1_per_mL", "max"),
        max_viral_load_cn2=("Viral_Load_CN2_per_mL", "max"),
        max_viral_load_cn3=("Viral_Load_CN3_per_mL", "max")
    )
    icu_only_table = raw_table.loc[raw_table.ICU_Status].groupby("Subject_ID").agg(
        max_sofa_icu_only=("SOFA_Score", "max"),
        sofa_ever_decreased_icu_only=("SOFA_Score", ever_decreased),
        sofa_ever_increased_icu_only=("SOFA_Score", ever_increased)
    )
    clinical_table["any_comorbidity"] = raw_table[
    raw_table.columns[raw_table.columns.str.startswith("COMORBID") |
                      raw_table.columns.str.startswith("ENCOUNTER_COMORBID")]]\
    .any(axis=1).groupby(raw_table.Subject_ID).any().astype("boolean")

    clinical_table = pd.concat([clinical_table, covid_only_table, icu_only_table], axis=1)
    clinical_table["deceased"] = clinical_table.deceased.astype("boolean")
    clinical_table["recovered"] = clinical_table.recovered.astype("boolean")
    clinical_table["discharged"] = clinical_table.discharged.astype("boolean") & ~clinical_table.deceased
    clinical_table["recovered_not_deceased"] = clinical_table.recovered & ~clinical_table.deceased
    clinical_table.loc[clinical_table.recovered, "deceased_vs_recovered"] = False
    clinical_table.loc[clinical_table.deceased, "deceased_vs_recovered"] = True
    clinical_table["deceased_vs_recovered"] = clinical_table.deceased_vs_recovered.astype("boolean")
    clinical_table.loc[clinical_table.discharged, "deceased_vs_discharged"] = False
    clinical_table.loc[clinical_table.deceased, "deceased_vs_discharged"] = True
    clinical_table["deceased_vs_discharged"] = clinical_table.deceased_vs_discharged.astype("boolean")

    clinical_table.loc[clinical_table.deceased, 'max_who'] = 8


    bvl_table = pd.read_csv("MSCIC_blood_viral_load_predictions.csv")
    bvl_table = bvl_table.drop_duplicates(subset="Subject_ID").set_index("Subject_ID")
    clinical_table.loc[clinical_table.covid_encounter, "blood_viral_load_bitscore"] = 0.0
    clinical_table["blood_viral_load_bitscore"] = bvl_table.bitscore
    clinical_table["blood_viral_load_bitscore_severe_only"] = clinical_table.blood_viral_load_bitscore.where(clinical_table.severity_ever_severe | clinical_table.severity_ever_eod)
    clinical_table["blood_viral_load_detected"] = clinical_table.blood_viral_load_bitscore.notnull().where(clinical_table.covid_encounter).astype("boolean")
    clinical_table["blood_viral_load_detected_severe_only"] = clinical_table.blood_viral_load_detected.where(clinical_table.severity_ever_severe | clinical_table.severity_ever_eod)

    clinical_table["severity_ever_increased_counting_death"] = clinical_table.severity_ever_increased | clinical_table.deceased
    clinical_table["severity_ever_decreased_counting_discharge"] = clinical_table.severity_ever_decreased | clinical_table.discharged
    for phenotype in binary_phenotypes:
        clinical_table[phenotype] = clinical_table[phenotype].astype("boolean")
    # clinical_table["sofa_ever_increased"] = clinical_table.sofa_ever_increased.astype("boolean")
    # clinical_table["sofa_ever_decreased"] = clinical_table.sofa_ever_decreased.astype("boolean")
    # clinical_table["sofa_ever_increased_icu_only"] = clinical_table.sofa_ever_increased_icu_only.astype("boolean")
    # clinical_table["sofa_ever_decreased_icu_only"] = clinical_table.sofa_ever_decreased_icu_only.astype("boolean")
    # clinical_table["who_ever_increased"] = clinical_table.who_ever_increased.astype("boolean")
    # clinical_table["who_ever_decreased"] = clinical_table.who_ever_decreased.astype("boolean")
    # clinical_table["severity_ever_increased"] = clinical_table.severity_ever_increased.astype("boolean")
    # clinical_table["severity_ever_decreased"] = clinical_table.severity_ever_decreased.astype("boolean")
    # clinical_table["severity_ever_moderate"] = clinical_table.severity_ever_moderate.astype("boolean")

    onset_date = raw_table.groupby("Subject_ID").COVID19_Symptoms_Reported_Onset_Date.first()
    encounter_start_date = raw_table.groupby("Subject_ID").COVID19_Encounter_Start_Date.first()
    encounter_end_date = raw_table.groupby("Subject_ID").COVID19_Encounter_End_Date.first()
    first_icu_date = raw_table.groupby("Subject_ID").apply(
        lambda x: x.loc[x.ICU_Status & x.COVID19_Within_Encounter, "Event_Date"]).groupby("Subject_ID").first()

    clinical_table["days_onset_to_encounter"] = encounter_start_date.sub(onset_date).div(np.timedelta64(1, "D"))
    clinical_table["covid_encounter_days"] = encounter_end_date.sub(encounter_start_date).div(np.timedelta64(1, "D"))
    clinical_table["covid_encounter_days"].mask(clinical_table.days_onset_to_encounter < 0, inplace=True)
    clinical_table["days_onset_to_encounter"].mask(
        (clinical_table.days_onset_to_encounter < 0) | (clinical_table.days_onset_to_encounter > 30), inplace=True)
    clinical_table["days_onset_to_icu"] = first_icu_date.sub(onset_date).div(np.timedelta64(1, "D"))
    clinical_table["days_onset_to_icu"].mask(clinical_table.days_onset_to_icu < 0, inplace=True)
    clinical_table["days_encounter_to_icu"] = first_icu_date.sub(encounter_start_date).div(np.timedelta64(1, "D"))
    clinical_table["days_encounter_to_icu"].mask(
        (clinical_table.days_encounter_to_icu < 0) | (clinical_table.days_encounter_to_icu > 60), inplace=True)
    clinical_table["icu_hours_cumulative"].mask(clinical_table.icu_hours_cumulative == 0, inplace=True)

    clinical_table["recruitment_date"] = (clinical_table.recruitment_date - clinical_table.recruitment_date.min()) / np.timedelta64(1, "D")
    clinical_table["sex"] = clinical_table.sex.str.get(0)
    clinical_table["age_sex"] = clinical_table["age"] * (clinical_table["sex"] == "F")
    clinical_table["age_squared"] = clinical_table["age"] ** 2
    clinical_table["is_hispanic"] = (clinical_table["ethnicity"] == "Hispanic or Latino").astype("boolean")
    clinical_table.loc[(clinical_table.ethnicity == "Unknown / Not Reported") | clinical_table.ethnicity.isnull(), "is_hispanic"] = pd.NA


    for column in binary_phenotypes:
        clinical_table[column] = clinical_table[column].astype("Int8")
        sns.countplot(x=clinical_table[column].astype('category'))
        plt.savefig(f"{column}.dist.png")
        plt.clf()
    for column in continuous_phenotypes:
            clinical_table[column + "_percentile"] = clinical_table[column].rank(pct=True)
            clinical_table[column + "_irnt"] = stats.norm.ppf(clinical_table[column + "_percentile"])
            clinical_table[column + "_log"] = clinical_table[column].transform("log1p")
            sns.histplot(x=clinical_table[column])
            plt.savefig(f"{column}.dist.png")
            plt.clf()
            sns.histplot(x=clinical_table[f"{column}_log"])
            plt.savefig(f"{column}_log.dist.png")
            plt.clf()

    clinical_table["race_factor"] = clinical_table.race.astype("category")
    # race = clinical_table.race.str.get_dummies()
    # race.columns = "race_" + race.columns
    # clinical_table = clinical_table.drop(["race", "ethnicity"], axis=1)

    flowcells = pd.read_csv("flowcells.csv", header=None, index_col=0, squeeze=True)
    flowcells = flowcells.drop("PICR7147T1a")
    flowcells.index = flowcells.index.str.split("T").str.get(0)
    clinical_table["multi_batch"] = flowcells.str.contains(r"\+")
    # clinical_table["flowcell_factor"] = flowcells.str.split("+").str.get(0).astype("category")
    flowcells = flowcells.str.get_dummies(sep="+")
    flowcells.columns = "flowcell_" + flowcells.columns

    # batches = pd.read_csv("batches.csv", header=None, index_col=0, squeeze=True)
    # batches = batches.drop("PICR7147T1a")
    # batches = batches.str.get_dummies(sep="+")
    # batches.columns = "batch_" + batches.columns
    # batches.index = batches.index.str.split("T").str.get(0)

    pca_table = pd.read_csv("625_Samples.cohort.QC_filtered.sample_matched.PCAir.txt", delim_whitespace=True, index_col="Subject_ID")
    pca_table.columns = [f"pc{i}" for i in range(1, len(pca_table.columns) + 1)]
    pca_table = pca_table.iloc[:, :10]

    clinical_table = clinical_table.join(flowcells, how="inner").join(pca_table,
                                                                      how="inner")  # .join(race, how="inner")
    clinical_table.to_csv(design_matrix_path)
    return {'phenotypes': all_phenotypes}


if __name__ == "__main__":
    build_design_matrix(snakemake.input[0], snakemake.output[0])
