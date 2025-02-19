#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 11:38:33 2019

@author: amohamme
"""

import re
import pandas as pd 
import warnings
warnings.filterwarnings("ignore")
import markov_clustering as mc
import numpy as np
from scipy.sparse import csc_matrix
from sklearn import preprocessing

def datanormalization(df,normtype,dist,cnt):
    if (normtype == "median"):
        newdf = df.groupby('Distance')[['Count']].median().reset_index()
        newdf.columns = ['Distance','medianCount']
        newdf['medianCount'] = newdf['medianCount'].map(int) 
        df = df.merge(newdf, how='left')
        df['Normalization'] = df['Count']/df['medianCount']
    if (normtype == "minmax"):
        min_max_scaler = preprocessing.MinMaxScaler()
        ncount = df[['Count']].as_matrix()
        df['Normalization'] = min_max_scaler.fit_transform(ncount)*100
        df['Normalization'] = df['Normalization'].round().map(int)
        df = df[['id_x','id_y','Normalization','Count','Distance']]
    dfi = df.copy()
    df['id_x'] = df['id_x'].map(int)
    df['id_y'] = df['id_y'].map(int)
    df['Normalization'] = df['Normalization']*100
    df['Normalization'] = df['Normalization'].map(int)    
    df = df[df['Distance'] <= dist]
    df = df[['id_x','id_y','Normalization']]
    cnmat = np.zeros((cnt, cnt), dtype=int)
    val = df.values
    cnmat[val[:,0].astype(int),val[:,1].astype(int)] = val[:,2]
    cnmat[val[:,1].astype(int),val[:,0].astype(int)] = val[:,2]
    cscmat = csc_matrix(cnmat, dtype=int)
    return(cscmat,dfi)

def clusterminmax(cinit,wdfinit):
    nwlist = list()
    allss = list()
    for clust in cinit:
        tlist = list()
        plist = list()
        if (len(clust) > 1):
            mn = wdfinit.loc[wdfinit['id'] == min(clust),'start'].tolist()[0]
            mx = wdfinit.loc[wdfinit['id'] == max(clust),'start'].tolist()[0]
            tlist.append(mn)
            tlist.append(mx)
            nwlist.append(tuple(tlist))
        for item in clust:
            mm = wdfinit.loc[wdfinit['id'] == item,'start'].tolist()[0]
            plist.append(mm)
        allss.append(tuple(plist))
    return(nwlist,allss)

def clusterpred(cscmat, wdfi, clg, inflation):
    result = mc.run_mcl(cscmat, inflation=inflation)
    clusters = mc.get_clusters(result)
    finclust = [items for items in clusters if (len(items) >= clg)] 
    nwlist,allss = clusterminmax(finclust,wdfi)
    return(clusters,finclust,nwlist,allss)

def combinationsdf(wdf):
    chrlist = pd.unique(wdf['chrom']).tolist()
    lf = []
    for i in chrlist:
        rdf = wdf[wdf['chrom'] == i]
        for i in rdf['id'].tolist():
            for j in rdf['id'].tolist():
                lf.append([i,j])
    fd = pd.DataFrame(list(lf),columns=['id_x','id_y'])
    return (fd)


def sorted_nicely(l): 
    convert = lambda text: int(text) if text.isdigit() else text 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)
    
