import pandas as pd
import os
import numpy as np
import xgboost as xgb
import sys
from collections import defaultdict

model_path = sys.argv[1]
feature_folder = sys.argv[2]
past_features = sys.argv[3]
outdir = sys.argv[4]

somatic_all = pd.DataFrame()
rmd_all = pd.DataFrame()
sbs_all = pd.DataFrame()
rna_all = pd.DataFrame()

TISSUE_DICT = {0: 'BRCA',
               1: 'COADREAD',
               2: 'LUAD',
               3: 'LGG',
               4: 'HNSC',
               5: 'LUSC',
               6: 'UCEC',
               7: 'BLCA',
               8: 'STAD',
               9: 'LIHC',
               10: 'KIRC',
               11: 'PRAD',
               12: 'OV',
               13: 'CESC',
               14: 'KIRP',
               15: 'SARC',
               16: 'PAAD',
               17: 'THCA',
               18: 'GBM',
               19: 'ESCC',
               20: 'ESA',
               21: 'TGCT',
               22: 'UVM',
               23: 'MESO',
               24: 'PCPG',
               25: 'UCS',
               26: 'ACC',
               27: 'KICH',
               28: 'DLBC',
               29: 'CHOL',
               30: 'THYM'}

# According to: https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations
TISSUE_DESCRIPTION = {"LAML": "Acute Myeloid Leukemia",
                      "ACC": "Adrenocortical carcinoma",
                      "BLCA": "Bladder Urothelial Carcinoma",
                      "LGG": "Brain Lower Grade Glioma",
                      "BRCA": "Breast invasive carcinoma",
                      "CESC": "Cervical squamous cell carcinoma and endocervical adenocarcinoma",
                      "CHOL": "Cholangiocarcinoma",
                      "LCML": "Chronic Myelogenous Leukemia",
                      "COAD": "Colon adenocarcinoma",
                      "CNTL": "Controls",
                      "ESCA": "Esophageal carcinoma",
                      "FPPP": "FFPE Pilot Phase II",
                      "GBM": "Glioblastoma multiforme",
                      "HNSC": "Head and Neck squamous cell carcinoma",
                      "KICH": "Kidney Chromophobe",
                      "KIRC": "Kidney renal clear cell carcinoma",
                      "KIRP": "Kidney renal papillary cell carcinoma",
                      "LIHC": "Liver hepatocellular carcinoma",
                      "LUAD": "Lung adenocarcinoma",
                      "LUSC": "Lung squamous cell carcinoma",
                      "DLBC": "Lymphoid Neoplasm Diffuse Large B-cell Lymphoma",
                      "MESO": "Mesothelioma",
                      "MISC": "Miscellaneous",
                      "OV": "Ovarian serous cystadenocarcinoma",
                      "PAAD": "Pancreatic adenocarcinoma",
                      "PCPG": "Pheochromocytoma and Paraganglioma",
                      "PRAD": "Prostate adenocarcinoma",
                      "READ": "Rectum adenocarcinoma",
                      "SARC": "Sarcoma",
                      "SKCM": "Skin Cutaneous Melanoma",
                      "STAD": "Stomach adenocarcinoma",
                      "TGCT": "Testicular Germ Cell Tumors",
                      "THYM": "Thymoma",
                      "THCA": "Thyroid carcinoma",
                      "UCS": "Uterine Carcinosarcoma",
                      "UCEC": "Uterine Corpus Endometrial Carcinoma",
                      "UVM": "Uveal Melanoma"}


for input_checking_file in os.listdir(feature_folder):

    if input_checking_file.endswith("_somatic.csv"):

        sample_name = input_checking_file.replace("-DNA_somatic.csv", "")

        somatic_file = f'{feature_folder}/{sample_name}-DNA_somatic.csv'
        rmd_file = f'{feature_folder}/{sample_name}-DNA_rmd_all_tcga.csv'
        sbs_file = f'{feature_folder}/{sample_name}-DNA_sbs_all_tcga.csv'
        rna_feature_counts_file = f'{feature_folder}/{sample_name}-DNA_filtered_VEP_rna_feature_generation.txt'

        if (os.path.isfile(somatic_file) and os.path.isfile(rmd_file) and os.path.isfile(sbs_file) and os.path.isfile(rna_feature_counts_file)):

            if 'somatic.csv' in somatic_file:
                samp_df = pd.read_csv(somatic_file)
                samp_df.iloc[0, 0] = sample_name
                samp_df.set_index('Unnamed: 0', inplace=True)
                samp_df.rename(
                    index={samp_df.index[0]: sample_name}, inplace=True)
                samp_df = samp_df.add_suffix('_mut')
                samp_df['TMB'] = samp_df.sum(axis=1) / 2.210733
                somatic_all = pd.concat([somatic_all, samp_df])

            if 'rmd_all_tcga.csv' in rmd_file:
                samp_df = pd.read_csv(rmd_file)
                samp_df = samp_df.T
                samp_df.columns = samp_df.iloc[0]
                samp_df.drop(samp_df.index[0], inplace=True, axis=0)
                samp_df.rename(
                    index={samp_df.index[0]: sample_name}, inplace=True)
                rmd_all = pd.concat([rmd_all, samp_df])

            if 'sbs_all_tcga.csv' in sbs_file:
                samp_df = pd.read_csv(sbs_file)
                samp_df = samp_df.T
                samp_df.columns = samp_df.iloc[0]
                samp_df.drop(samp_df.index[0], inplace=True, axis=0)
                samp_df.rename(
                    index={samp_df.index[0]: sample_name}, inplace=True)
                sbs_all = pd.concat([sbs_all, samp_df])

            if '_VEP_rna_feature_generation.txt' in rna_feature_counts_file:
                samp_df = pd.read_csv(rna_feature_counts_file)
                samp_df = samp_df[['hgnc_symbol', 'FPKM']]
                samp_df['hgnc_symbol'] = samp_df['hgnc_symbol'].apply(
                    lambda x: x + '_FPKM')
                samp_df = samp_df.T
                samp_df.columns = samp_df.iloc[0]
                samp_df.drop(samp_df.index[0], inplace=True, axis=0)
                samp_df.rename(
                    index={samp_df.index[0]: sample_name}, inplace=True)
                rna_all = pd.concat([rna_all, samp_df])

print(f"somatic_all tables\n")
print(somatic_all)
print(f"\nrmd_all tables\n")
print(rmd_all)
print(f"\nsbs_all tables\n")
print(sbs_all)
print(f"\nrna_all tables\n")
print(rna_all)
print(somatic_all.index.value_counts())
print(rmd_all.index.value_counts())
print(sbs_all.index.value_counts())
print(rna_all.index.value_counts())
feat_mtx = pd.concat([somatic_all, sbs_all, rmd_all, rna_all], axis=1)
print(f"\nfeat_mtx tables\n")
print(feat_mtx)
print(feat_mtx.index.value_counts().to_list())

print(f"\nmapping_df table\n")
mapping_df = pd.read_csv(past_features)
print(mapping_df)

mapping_df.set_index('index', inplace=True)
mapping_df.drop('Weight', axis=1, inplace=True)

print(f"\nmapping_df final table\n")
print(mapping_df)

feat_new = pd.DataFrame(columns=mapping_df.columns.tolist())
print(f"\nfeat_new table\n")
print(feat_new)
feat_new = pd.concat([feat_new, feat_mtx])[feat_new.columns]
feat_new = feat_new.fillna(
    {col: 0 for col in feat_new.columns if feat_new[col].isna().all()})
feat_new.fillna(0, inplace=True)

print(f"\nfeat_new table before 2run\n")
print(feat_new)

to_run = np.array(feat_new)
print(f"\nfeat_new table after 2run\n")
print(feat_new)
print(feat_new.index.value_counts())

# imputer = KNNImputer(n_neighbors=2, weights="uniform")
# to_run = imputer.fit_transform(to_run)

print(f"\nto_run table\n")
print(to_run)

params = {'subsample': 0.8,
          'reg_lambda': 1,
          'reg_alpha': 1,
          'objective': 'multi:softprob',
          'n_estimators': 500,
          'min_child_weight': 1,
          'max_depth': 100,
          'learning_rate': 0.01,
          'gamma': 0.5,
          'colsample_bytree': 0.5,
          'booster': 'gbtree',
          'base_score': 1.2}

xgb_fit = xgb.XGBClassifier(**params)
xgb_fit.load_model(model_path)

print(f"\nresults before table\n")
results = xgb_fit.predict(to_run)
print(results)

probs = xgb_fit.predict_proba(to_run)

print(f"\nresults table\n")
print(results)

print(f"\nprobs table\n")
print(probs)

print(f"\nxgb_fit\n")
print(xgb_fit)

#Label Conversion
convert_num_label = np.vectorize(lambda x: TISSUE_DICT.get(x, "unknown"))
convert_label_full = np.vectorize(lambda x: TISSUE_DESCRIPTION.get(x, "unknown"))
results_code = convert_num_label(results)
results_full = convert_label_full(results_code)

#Probability Conversion
prob_predicted = []
for i in range(len(probs)):
    prob_predicted.append(probs[i][results[i]])
    
results_compiled = pd.DataFrame({'Sample': feat_mtx.index.tolist(), 'Tissue Abbrev': results_code, 
                                 'Model Prediction': results_full, 'Prediction Confidence': prob_predicted})

probabilities = defaultdict(dict)
for j in range(len(feat_mtx.index.tolist())):
    one_patient = {i: val for i, val in enumerate(probs[j])}
    one_patient = {TISSUE_DICT.get(k, k): v for k, v in one_patient.items()}
    probabilities[feat_mtx.index.tolist()[j]] = one_patient

regular_dict = dict(probabilities)

# Convert to MultiIndex DataFrame
probs_df = pd.DataFrame.from_dict(regular_dict, orient='index')
probs_df = probs_df.stack().to_frame(name='value').reset_index()
probs_df.columns = ['Sample Name', 'Tumor Type', 'Probability of Prediction']
probs_df = probs_df.set_index(['Sample Name', 'Tumor Type'])

results_compiled.to_csv(outdir + 'model_results.tsv', index=False, sep="\t")
probs_df.to_csv(outdir + 'model_probs_df.tsv', index = False, sep = "\t")