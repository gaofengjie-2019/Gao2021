# -*- coding: utf-8 -*-
"""
Created on Fri Jun  4 10:48:37 2021

@author: A
"""


import argparse
import numpy as np
import pickle as pkl
import pandas as pd
from scipy.integrate import RK45, solve_ivp
from scipy.stats import pearsonr
from scipy.special import logsumexp
import matplotlib.pyplot as plt
import pickle as pkl
import seaborn as sns
from sklearn import preprocessing
from generalized_lotka_volterra import GeneralizedLotkaVolterra
from stein_plot_parameters import plot_heatmaps,swap_denom,choose_perturb_denom
from stein_lotka_volterra_comparison import adjust_concentrations
from compositional_lotka_volterra import CompositionalLotkaVolterra, \
                                         estimate_relative_abundances, \
                                         choose_denom, \
                                         construct_alr, \
                                         compute_rel_abun, \
                                         ridge_regression_clv
##loading data
meta = pd.read_csv('metadata.csv', index_col=0)
species = pd.read_csv('species_r.csv', index_col=0)   
meta1 = meta.loc[meta.CompleteSeries=='CS',:]
meta1 = meta1.loc[-meta1.TP.isin([3,4])]
species = species.T
species = species.loc[meta1.index, :]

##obtain the four inducers
inducer = pd.read_csv('dif_inducer.csv', index_col = 0)
species_inducer = pd.DataFrame()
for i in species.columns:
    if i in inducer.index:
        species_inducer.loc[:,i] = species.loc[:,i]
    else:
        continue

##obtain microbes that present in 20% of the samples
bac = pd.read_csv('resultIESR.csv', index_col=0)
bac1 = bac.loc[bac.P<0.05,:]
bac1.index = bac1.loc[:,'X']
species_bac = pd.DataFrame()
for i in species.columns:
    if i in bac1.index:
        species_bac.loc[:,i] = species.loc[:,i]
    else:
        continue
species_bac['id'] = meta1.loc[:,'id']
##除去在20个以上样本中都为0的菌，在所有样本中都存在的菌为0，只有一个大肠埃希菌
##除去20%以上样本都为0的菌,只有一个大肠杆菌
##除去50%以上都为0的菌, 只有一个大肠杆菌
##80%为阈值，得10个菌
aa=pd.DataFrame()
bb=pd.read_csv('microbes.csv', index_col=0)
for i in species_bac.columns:
    if i not in bb.index:
        continue   
    else:
        aa.loc[:,i] = (species_bac.loc[:,i])
aa.loc[:,'id'] = meta1.loc[:,'id']

##put species in Y1
species1 = aa.copy()
species1 = species1.replace(0,1e-8)
Y1=[]
tp = 2
bac = 10
for i in np.unique(species1.id):
    current_data=species1.loc[species1.id==i, :]
    a = np.zeros((tp,bac))
    for ia in range(tp):
        for ja in range(bac):
            a[ia,ja]=current_data.iloc[ia,ja]
    Y1.append(a)

##put perturbation in U1        
U1 = []
species1 = species_inducer.copy()
species1.loc[:,'id']=meta1.loc[:,'id']
tp = 2
bac = 4
for i in np.unique(species1.id):
    current_data=species1.loc[species1.id==i, :]
    a = np.zeros((tp,bac))
    for ia in range(tp):
        for ja in range(bac):
            a[ia,ja]=current_data.iloc[ia,ja]
    U1.append(a)


T1 = []
tt=np.array([0.,14.])
for i in np.unique(meta1.id):
    T1.append(tt)

##compute A,g,B
Y = Y1.copy()
U = U1.copy()
T = T1.copy()

Y = adjust_concentrations(Y)
# estimated previously
r_A = 1
r_g = 4
r_B = 0.5
     
clv = CompositionalLotkaVolterra(Y, T, U, pseudo_count=1e-8)
clv.r_A = r_A
clv.r_g = r_g
clv.r_B = r_B
clv.train_ridge()
A_clv, g_clv, B_clv = clv.get_params()

col_names = np.array(aa.drop('id',axis=1).columns)
ntaxa = Y[0].shape[1]

old_denom = clv.denom
taxon_row_namesa = col_names[np.array([i for i in range(ntaxa) if i != old_denom])]
antibiotics = ['B.eggerthii', 'Intestinibacter bartlettii DSM 16795',
       'Stenotrophobacter_uncultured Acidobacteria bacterium',
       'E.hallii group_uncu.'] #the four inducers
plot_heatmaps(A_clv, B_clv, taxon_row_namesa, col_names, antibiotics, antibiotics, "clv-A1", "A")

perturb_denom = choose_perturb_denom(Y)
taxon_row_namesb = col_names[np.array([i for i in range(ntaxa) if i != perturb_denom])]
A_new, g_new, B_new, antibiotic_row_names = swap_denom(A_clv, g_clv, B_clv, perturb_denom, col_names)
plot_heatmaps(A_new, B_new, taxon_row_namesb, col_names, antibiotic_row_names, antibiotics, "clv-B1", "B")

##obtain the parameters of inducer and microbes interation for visualization in cytoscape
A_clv1 = pd.DataFrame(A_clv, columns=col_names, index = taxon_row_namesa)
A_clv1.to_csv("A_clv.csv")
A = pd.DataFrame(columns=['source','target', 'value'], index=range(90))
k=0
for i in range(10):
    for j in range(9):
        A.iloc[k,0] = A_clv1.columns[i]
        A.iloc[k,1] = A_clv1.index[j]
        A.iloc[k,2] = A_clv1.iloc[j,i]
        k=k+1

B_new1 = pd.DataFrame(B_new, columns=antibiotics, index = taxon_row_namesb)
B_new1.to_csv("B_clv.csv")
B = pd.DataFrame(columns=['source','target', 'value'], index=range(36))
k=0
for i in range(4):
    for j in range(9):
        B.iloc[k,0] = B_new1.columns[i]
        B.iloc[k,1] = B_new1.index[j]
        B.iloc[k,2] = B_new1.iloc[j,i]
        k=k+1 
edge = pd.concat([A,B], axis=0)
edge.loc[:,'abs'] = abs(edge.value)
edge['level'] = 1
for i in range(edge.shape[0]):
    if edge.iloc[i,2] <0 :
        edge.iloc[i,4] = -1
    else:
        edge.iloc[i,4] = 1
#edge = edge.loc[edge.value!=0,]
#remove the edges that its source and target are the same
edge['tong'] = True
edge.index = edge.loc[:,'target']
for j in range(edge.shape[0]):
    if edge.index[j] == edge.iloc[j,0]:
        edge.iloc[j,5] = False
    else:
        continue
edge = edge.loc[edge.tong == True,: ]    
#standard the perturbation
edge.index = edge.loc[:,'source']
edge_pu = edge.loc[ ['B.eggerthii', 'Intestinibacter bartlettii DSM 16795',
       'Stenotrophobacter_uncultured Acidobacteria bacterium',
       'E.hallii group_uncu.'],:]
edge1 = np.array(edge_pu.loc[:,'value']).reshape(-1,1)
max_abs_scaler_data = preprocessing.MaxAbsScaler().fit_transform(edge1)
edge_pu['standard'] = max_abs_scaler_data
edge_pu['abs1'] = abs(edge_pu.loc[:,'standard'])
edge_pu1 = edge_pu
#standard the interaction of microbes
edge_pu = edge.drop( ['B.eggerthii', 'Intestinibacter bartlettii DSM 16795',
       'Stenotrophobacter_uncultured Acidobacteria bacterium',
       'E.hallii group_uncu.'],axis=0)
edge1 = np.array(edge_pu.loc[:,'value']).reshape(-1,1)
max_abs_scaler_data = preprocessing.MaxAbsScaler().fit_transform(edge1)
edge_pu['standard'] = max_abs_scaler_data
edge_pu['abs1'] = abs(edge_pu.loc[:,'standard'])
edge_pu2 = edge_pu

edge3 = pd.concat([edge_pu1, edge_pu2], axis = 0)
edge3 = edge3.sort_values(by = 'abs1', ascending = False)        
edge3.to_csv("edge12.csv")

