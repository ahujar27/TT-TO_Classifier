import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.colors
import os
import sys

import sklearn
from sklearn import metrics
from sklearn.svm import LinearSVC, SVC
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split, cross_val_score, KFold,StratifiedKFold,GridSearchCV,RandomizedSearchCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectFromModel,RFE
from sklearn.decomposition import PCA, NMF
from sklearn.preprocessing import MinMaxScaler, label_binarize
from sklearn.pipeline import Pipeline
from sklearn.metrics import classification_report,confusion_matrix,accuracy_score,roc_curve, auc, top_k_accuracy_score
import xgboost as xgb
from itertools import cycle
from xgboost import XGBClassifier
from sklearn.multiclass import OneVsRestClassifier
from sklearn.utils import class_weight
import random

model_path = sys.argv[0]
gene_list = sys.argv[1]
rmd_path = sys.argv[2]
cnv_path = sys.argv[3]
snv_path = sys.argv[4]
labels_path = sys.argv[5]
sbs_path = sys.argv[6]
past_features = sys.argv[7]
outdir = sys.argv[8]

xgb_model = xgb.XGBClassifier(colsample_bytree=0.3,gamma=0, learning_rate=0.1,
              max_depth=100, min_child_weight=5, n_estimators=500,subsample = 0.8,objective='multi:softprob')
xgb_model.load_model('full.model')

genelist = pd.read_csv(gene_list,names=["Genes"], encoding='latin-1')
genelist = genelist.Genes.tolist() 
rmd = pd.read_csv(rmd_path,index_col=0)
rmd_df = rmd.T

CNV_mtx = pd.read_csv(cnv_path,index_col=0)

somatic_mtx = pd.read_csv(snv_path, index_col=0)

labels = pd.read_csv(labels_path,index_col=0)
sbs = pd.read_csv(sbs_path,index_col=0)
sbs = sbs.T

CNV_mtx = CNV_mtx.div(CNV_mtx.max().max(),axis=0)
shared_samples = list(set(sbs.index.tolist()) & set(somatic_mtx.index.tolist()))
somatic_mtx = somatic_mtx.loc[shared_samples]
CNV_mtx = CNV_mtx.loc[shared_samples]
rmd_df = rmd_df.loc[shared_samples]
sbs = sbs.loc[shared_samples]

somatic_mtx = somatic_mtx.add_suffix('_mut')
CNV_mtx = CNV_mtx.add_suffix('_CNV')
feat_mtx = pd.concat([somatic_mtx,CNV_mtx,rmd_df,sbs], axis=1)

past_feat = pd.read_csv(past_features, index_col = 0)

for column in feat_mtx.columns:
    if column not in past_feat.columns:
        feat_mtx = feat_mtx.drop(column, axis=1)

x = np.array(feat_mtx)

predictions = xgb_model.predict(x)

predictions.to_csv(outdir + 'predictions.csv')
