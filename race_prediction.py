import pandas as pd
import numpy as np
import logging
from sklearn.calibration import CalibratedClassifierCV
from sklearn.neighbors import KNeighborsClassifier
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import GridSearchCV, cross_val_score, cross_val_predict

import seaborn as sns
import matplotlib.pyplot as plt


def predict_races():
    pc_table = pd.read_csv("625_Samples.cohort.QC_filtered.sample_matched.PCAir.txt", sep=" ", quotechar='"', index_col=False)
    pc_table.columns = ["Subject_ID"] + [f"PC{i}" for i in range(1,len(pc_table.columns))]
    pc_table = pc_table.set_index("Subject_ID")

    clinical_table = pd.read_csv("/sc/arion/projects/mscic1/data/covariates/clinical_data_deidentified_allsamples/Biobank_clinical_data_table_by_blood_sample_deidentified_UNCONSENTED.csv.gz")

    raw_races = clinical_table.groupby("Subject_ID")["Race_From_EHR"].first().reindex(pc_table.index)
    raw_ethnicities = clinical_table.groupby("Subject_ID")["Ethnicity_From_EHR"].first().reindex(pc_table.index)
    consent_races = clinical_table.groupby("Subject_ID")["Race_From_Consent"].first().reindex(pc_table.index)
    consent_ethnicities = clinical_table.groupby("Subject_ID")["Ethnicity_From_Consent"].first().reindex(pc_table.index)

    white = raw_races == "WHITE"
    black = np.isin(raw_races, ["AFRICAN-AMERICAN", "GHANAIAN", "GUINEAN", "HAITIAN", "JAMAICAN", "NIGERIAN", "OTHER: WEST AFRICAN", "OTHER: NORTH AFRICAN", "TRINIDADIAN"])
    asian = np.isin(raw_races, ["ASIAN INDIAN", "BANGLADESHI", "CHINESE", "FILIPINO", "JAPANESE", "KOREAN", "LAOTIAN", "INDONESIAN", "SINGAPOREAN"])
    hispanic = raw_ethnicities == "HISPANIC"
    race_categories = np.full(len(raw_races), "OTHER/UNKNOWN")
    race_categories[white] = "WHITE"
    race_categories[black] = "BLACK"
    race_categories[asian] = "ASIAN"
    race_categories[hispanic] = "HISPANIC"

    enc = LabelEncoder().fit(race_categories)
    model = GridSearchCV(KNeighborsClassifier(), {'n_neighbors': range(2,11)})
    score = cross_val_score(model, pc_table, race_categories)
    median_score = np.median(score)
    std_score = np.std(score)
    logging.info(f"Cross-validation accuracy of race prediction: {median_score} \xb1 {1.96*std_score}")
    
    predictions = cross_val_predict(model, pc_table, race_categories)
    pc_table.insert(0, "predicted_race_category", predictions)
    pc_table.insert(1, "race_category_from_EHR", race_categories)
    pc_table.insert(2, "raw_race_from_EHR", raw_races)
    pc_table.insert(3, "raw_ethnicity_from_EHR", raw_ethnicities)
    pc_table.insert(4, "raw_race_from_consent", consent_races)
    pc_table.insert(5, "raw_ethnicity_from_consent", consent_ethnicities)

    for race in "WHITE", "BLACK", "ASIAN", "HISPANIC":
        pc_table.loc[pc_table.predicted_race_category == race].index.to_frame().to_csv(f"{race}.indiv_list.txt", sep="\t", index=False)


    pc_table.to_csv("625_Samples.cohort.QC_filtered.sample_matched.race_and_PCA.csv")

    for i in range(1,10):
        for j in range(i+1,11):
            plt.clf()
            sns.scatterplot(x=pc_table[f"PC{i}"], y=pc_table[f"PC{j}"], hue=predictions)
            plt.savefig(f"625_Samples.cohort.QC_filtered.sample_matched.PC{i}v{j}.predicted_races.pdf")
            plt.clf()


if __name__ == "__main__":
    predict_races()
