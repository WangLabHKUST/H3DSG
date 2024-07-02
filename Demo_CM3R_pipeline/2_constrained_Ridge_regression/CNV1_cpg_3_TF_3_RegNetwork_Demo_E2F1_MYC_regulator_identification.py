#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 10:09:54 2023

@author: shanghaixia
"""

import os
import sys
os.chdir('/Users/shanghaixia/Dropbox/SP_multi-omics_data-for_figures/SP_multi_omics_integretion_Figure_code/Fig4_part/Fig4_CMR3_pipeline_code/')

#os.mkdir('Up_gene_regression_without_TF') 
#os.mkdir('paramer_optimization/lasso_L1_paramter_grid_search') 

import numpy as np
np.random.seed(42)
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
from sklearn.linear_model import Lasso
import re
#from sklearn.linear_model import ElasticNet
import pandas as pd
from sklearn.model_selection import LeaveOneOut
from sklearn.model_selection import cross_val_score
from scipy.stats.stats import pearsonr
cv=LeaveOneOut()
from sklearn.linear_model import Ridge
import re
from scipy.spatial import distance

pth_path1='./1_multi_omics_data_preparation/Prepared_data/'
store_regress_path='./2_constrained_Ridge_regression/Regression_results/'
#os.mkdir(store_regress_path)

file_dir = os.listdir(pth_path1)
file_dir.remove('.DS_Store')

def sigmoid(z):
    return(1/(1+np.exp(-z)))

alpha_list=[0.000001,0.00001,0.0001,0.01,0.011,0.012,0.013,
            0.014,0.015,0.016,0.017,0.018,0.019,
            0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]

for i in file_dir:
    
    gene=re.sub("__info_min_max_score.csv","",i)
    dat =pd.read_csv(pth_path1+i,header=0)
        
    name=dat.columns.values.tolist()
    y=sigmoid(dat.iloc[:,1].to_numpy())
    nameii=name[2:]
    cpg_loc=[]
    for oo in range(0,len(nameii)):
        if 'cg' in nameii[oo]:
            cpg_loc.append(oo)
    
    X1=dat.iloc[:,2:].to_numpy()
    tf=X1[:,(max(cpg_loc)+1):]
    print(gene)
    for kkk in range(0,tf.shape[1]):
        pcc_tf_exp,sig =pearsonr(dat.iloc[:,1],tf[:,kkk])
        print(pcc_tf_exp)
        if pcc_tf_exp <0:
            X1[:,max(cpg_loc)+1+kkk]=(-1)*X1[:,max(cpg_loc)+1+kkk]
    X1[:,1:(max(cpg_loc)+1)]=X1[:,1:(max(cpg_loc)+1)]*(-1)
    X11=sigmoid(X1)
    dts_all=[]
    num_nonzero=[]
    coeff=[]
    inte=[]
        
    for hh in range(0,len(alpha_list)):
        lin_all = Ridge(alpha=alpha_list[hh],max_iter=1000,positive=True,random_state=9999)#0.000 #, selection='random'precompute=True,
        scores_all = cross_val_score(lin_all, X11, y, scoring='neg_mean_absolute_error',
                         cv=cv, n_jobs=-1)
        y_pred_lasso_all=lin_all.fit(X11,y).predict(X11)
        lin_all.coef_
        num_nonzero.append(np.count_nonzero(np.greater(lin_all.coef_,0)))
        coeff.append(lin_all.coef_)
        lin_all.intercept_
        inte.append(lin_all.intercept_)
        dst=distance.euclidean(y,y_pred_lasso_all)
        dts_all.append(dst)
        
    store2_all=pd.concat([pd.concat([pd.DataFrame(dts_all,columns=['ecludiean distance']),
        pd.concat([pd.DataFrame(alpha_list,columns=['beta_value']),
        pd.DataFrame(coeff,columns=nameii)],axis=1)],axis=1),
                      pd.DataFrame(inte,columns=['intercept'])],axis=1)
    tar_dir=store_regress_path
    store2_all.to_csv(tar_dir+gene+'_min_max_neg_cpg_cnv_tf_sigmoidtrans_Ridge_regression.csv',
                  index=False)
    
    
    
    
    