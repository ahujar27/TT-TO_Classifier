import pandas as pd
import os
import numpy as np
import xgboost as xgb
from sklearn.impute import KNNImputer

model_path = sys.argv[1]
feature_folder = sys.argv[2]
past_features = sys.argv[3]
outdir = sys.argv[4]

somatic_all = pd.DataFrame()
rmd_all = pd.DataFrame()
sbs_all = pd.DataFrame()
rna_all = pd.DataFrame()

for file in os.listdir(feature_folder):
    if 'somatic.csv' in file:
        samp_df = pd.read_csv(feature_folder + file)
        samp_df.iloc[0,0] = file
        samp_df.set_index('Unnamed: 0', inplace = True)
        samp_df.rename(index={samp_df.index[0]: file[0:8]}, inplace = True)
        samp_df = samp_df.add_suffix('_mut')
        samp_df['TMB'] = samp_df.sum(axis = 1) / 2.210733
        somatic_all = pd.concat([somatic_all, samp_df])
    if 'rmd_all_tcga.csv' in file:
        samp_df = pd.read_csv(feature_folder + file)
        samp_df = samp_df.T
        samp_df.columns = samp_df.iloc[0]
        samp_df.drop(samp_df.index[0], inplace = True, axis = 0)
        samp_df.rename(index={samp_df.index[0]: file[0:8]}, inplace = True)
        rmd_all = pd.concat([rmd_all, samp_df])
    if 'sbs_all_tcga.csv' in file:
        samp_df = pd.read_csv(feature_folder + file)
        samp_df = samp_df.T
        samp_df.columns = samp_df.iloc[0]
        samp_df.drop(samp_df.index[0], inplace = True, axis = 0)
        samp_df.rename(index={samp_df.index[0]: file[0:8]}, inplace = True)
        sbs_all = pd.concat([sbs_all, samp_df])
    if 'rna' in file:
        samp_df = pd.read_csv(feature_folder + file)
        samp_df = samp_df[['hgnc_symbol', 'FPKM']]
        samp_df['hgnc_symbol'] = samp_df['hgnc_symbol'].apply(lambda x: x + '_FPKM')
        samp_df = samp_df.T
        samp_df.columns = samp_df.iloc[0]
        samp_df.drop(samp_df.index[0], inplace = True, axis = 0)
        samp_df.rename(index={samp_df.index[0]: file[0:8]}, inplace = True)
        rna_all = pd.concat([rna_all, samp_df])

feat_mtx = pd.concat([somatic_all, sbs_all, rmd_all, rna_all], axis = 1)

mapping_df = pd.read_csv(past_features)
mapping_df.set_index('index', inplace = True)
mapping_df.drop('Weight', axis = 1, inplace = True)

feat_new = pd.DataFrame(columns = mapping_df.columns.tolist())
feat_new = pd.concat([feat_new, feat_mtx])[feat_new.columns]
feat_new = feat_new.fillna({col: 0 for col in feat_new.columns if feat_new[col].isna().all()})

to_run = np.array(feat_new)
imputer = KNNImputer(n_neighbors=2, weights="uniform")
to_run = imputer.fit_transform(to_run)

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
results = xgb_fit.predict(to_run)
probs = xgb_fit.predict_proba(to_run)

np.savetxt(outdir, results, fmt='%d')
