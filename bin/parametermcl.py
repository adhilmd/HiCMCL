#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 13:40:57 2019

@author: amohamme
"""

import argparse

def argparser():
    parser = argparse.ArgumentParser(prog = '''\n.parametermcl.py''', description='''\n-----------Parameter MCL-------- \n
    \n[Date: 7th May 2018], \n[help: python parametermcl.py -h]\n''', usage = 'parametermcl.py *args')
    parser.add_argument('-findb','--findb', type=str, dest='findb', help="Comma seperated database file (Single file can also be provided) (Mandatory)", action = 'store', required = True)
    parser.add_argument('-chrl','--chrl', type=str, dest='chrl', help="Chromosome length file containing two columns (chromosomenumber, length) without header (Mandatory)", action = 'store', required = True)
    parser.add_argument('-icore','--icore', type=int, dest='icore', help="Core resolution for clustering, in multiples (default=10)", action = 'store', default=10)
    parser.add_argument('-nres','--nres', type=int, dest='nres', help="Base length for each random location (default=10000000)", action = 'store', default=10000000)
    parser.add_argument('-cpath','--cpath', type=str, dest='cpath', help="Main path were the codes are present for modules", action='store', required=True)
    parser.add_argument('-norm','--norm', type=str, dest='norm', help="Normalization type 'median' or 'minmax' (default=minmax)", action = 'store', default="minmax")
    parser.add_argument('-rand','--rand', type=int, dest='rand', help="Number of random location for the calculation (default=12)", action = 'store', default=12)
    parser.add_argument('-th','--th', type=int, dest='threads', help="Number of threads (default=8)", action = 'store', default=8)
    parser.add_argument('-tag','--tag', type=str, dest='tag', help="Comma seperated file tags for pickle file (single tag can also be provided), the tags should match the -bedpe (Mandatory)", action = 'store', required = True)
    parser.add_argument('-pref','--pref', type=str, dest='pref', help="Prefix for the output files (Mandatory)", action = 'store', required = True)
    parser.add_argument('-odir','--odir', type=str, dest='odir', help="Outdir (Mandatory)", action = 'store', required = True)
    args = parser.parse_args()
    return(args)

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
args = argparser()
import sys
sys.path.append(args.cpath)
from clusteringmodules import datanormalization 
import os
import sqlite3
import pandas as pd 
import random
from pybedtools import BedTool
from multiprocessing import Pool
from functools import partial
import warnings
warnings.filterwarnings("ignore")
import markov_clustering as mc
import time
import pickle

def parameterestimte(order,dfest,findbase,snamef,binres,dc,normtype):
    featuredf = pd.DataFrame(columns=['Sample','Order','bins','DistCutoff','Inflation','ModularityScore','ClusterSize'])
    chrom = dfest.loc[order]['chrom']
    irands = dfest.loc[order]['start']
    irande = dfest.loc[order]['end']
    ordr =  dfest.loc[order]['order']
    connex = sqlite3.connect(findbase)
    cur = connex.cursor()
    sql = "SELECT chr1,loc1,loc2 FROM " + snamef + "_intra WHERE chr1 == '"+ str(chrom) + "' AND loc1 >= "+ str(irands) + " AND loc2 <= "+ str(irande) +";"
    df = pd.read_sql(sql, connex)
    cur.close()
    connex.close()
    if (df.shape[0] != 0):
        for bins in binres:
            dfestbins = dfest.copy()
            dfestbins = dfestbins.loc[order:order]
            dfestbins['start'] = dfestbins['start']/bins
            dfestbins['start'] = dfestbins['start'].astype(float).round().map(int)*bins
            dfestbins['end'] = dfestbins['end']/bins
            dfestbins['end'] = dfestbins['end'].astype(float).round().map(int)*bins
            a = BedTool()
            wbed = BedTool.from_dataframe(dfestbins)
            windows = a.window_maker(b=wbed, w=bins, s=bins)
            wdfall = pd.read_table(windows.fn, names=['chrom', 'start', 'end'])
            wdfall['chrom'] = wdfall['chrom'].map(str)
            wdfall['id'] = [i+1 for i in range(0,wdfall.shape[0])]
            pdf = df.copy()
            pdf['loc1'] = pdf['loc1']/bins
            pdf['loc1'] = pdf['loc1'].astype(float).round().map(int)*bins
            pdf['loc2'] = pdf['loc2']/bins
            pdf['loc2'] = pdf['loc2'].astype(float).round().map(int)*bins
            pdf = pdf.groupby(['chr1','loc1','loc2']).size().reset_index()
            pdf.columns = ['chrom','loc1','loc2','Count']
            cnt = wdfall.shape[0]+1
            pdf = pdf.merge(wdfall[['chrom','start','id']], left_on=['chrom','loc1'], right_on=['chrom','start'], how = 'left')
            pdf = pdf.dropna()
            pdf = pdf.merge(wdfall[['chrom','start','id']], left_on=['chrom','loc2'], right_on=['chrom','start'], how = 'left')
            pdf = pdf.dropna()
            pdf = pdf[['chrom','loc1','loc2','Count','id_x','id_y']]
            pdf['Distance'] = pdf['loc2'] - pdf['loc1']
            normtype="minmax"
            cscmat,dfi = datanormalization(pdf,normtype,dc*bins,cnt)
            for inflation in [i/10 for i in range(15,30,3)]:
                result = mc.run_mcl(cscmat, inflation=inflation)
                clusters = mc.get_clusters(result)
                Q = mc.modularity(matrix=result, clusters=clusters)
                featuredf = featuredf.append({'Sample':sname[fname],'Order':ordr,'bins':bins,'DistCutoff':dc,'Inflation':inflation,'ModularityScore':Q, 'ClusterSize':len(clusters)}, ignore_index=True)
    else:
        featuredf = featuredf.append({'Sample':sname[fname],'Order':ordr,'bins':'NA','DistCutoff':'NA','Inflation':'NA','ModularityScore':'NA','ClusterSize':'NA'}, ignore_index=True)
    return(featuredf)       


if __name__ == "__main__":
    outdir = args.odir
    findb = args.findb.split(',')
    sname = args.tag.split(',')
    chrl = args.chrl
    normtype = args.norm
    nres = args.nres
    totrand = args.rand
    procs = args.threads
    dc = args.icore
    pref = args.pref
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    if (len(findb) != len(sname)):
        print ("Error: Please provide same number of -findb and -tag files")
        quit()
    #chrl = "/storage/home/amohamme/data/AnnotationFiles/GRCh38/wholechr.tsv"
    wc = pd.read_csv(chrl, sep="\t", header=None)
    wc.columns = ['chrom','end']
    wc['start'] = 0
    wc['order'] = [i+1 for i in range(0,wc.shape[0])]
    dfest = pd.DataFrame(columns=['chrom','start','end'])
    for t in range(totrand):
        i = random.randint(0,wc.shape[0]-1)
        epos = wc['end'][i]
        if (epos >= nres+10000):
            irands = random.randint(10000,epos-nres)
            irande = irands+nres
            chrr = wc['chrom'][i]
            dfest = dfest.append({'chrom':chrr, 'start':irands, 'end':irande}, ignore_index=True)
        else:
            dfest = dfest.append({'chrom':wc['chrom'][i], 'start':wc['start'][i], 'end':wc['end'][i]}, ignore_index=True)   
    dfest['order'] = [i+1 for i in range(0,dfest.shape[0])]
    if (sum(wc['end']) <= 100000000):
        binres = [(i+2)*1000 for i in range(0,20,2)]
    elif ((sum(wc['end']) > 100000000) & (sum(wc['end']) <= 1000000000)):
        binres = [(i+5)*1000 for i in range(0,50,5)]
    elif (sum(wc['end']) > 1000000000):
        binres = [(i+10)*1000 for i in range(0,100,10)]
    #distcutoff = [5,10,15,20]
    inflation = [i/10 for i in range(15,30,3)]
    #findb = ["/storage/home/amohamme/data/HiC/human/outdir/amcl/SRR5494758_fin.sqlite","/storage/home/amohamme/data/HiC/human/outdir/amcl/SRR5494759_fin.sqlite","/storage/home/amohamme/data/HiC/human/outdir/amcl/SRR5494760_fin.sqlite"]
    #sname = ["SRR5494758","SRR5494759","SRR5494760"]
    progstarts = time.time() 
    findat = pd.DataFrame()
    procs = 12
    orderlist=dfest.index.tolist()
    for fname in range(len(sname)):
        pool = Pool(processes=procs)
        mrtf=pool.map(partial(parameterestimte,dfest=dfest,findbase=findb[fname],snamef=sname[fname],binres=binres,dc=dc,normtype=normtype), orderlist)
        pool.close()
        pool.join()    
        fnb = pd.concat(mrtf) 
        findat = pd.concat([findat,fnb])
    nowtime = time.time()    
    runtime = nowtime-progstarts
    print (runtime)
    findat['ModularityScore'] = findat['ModularityScore'].map(str)     
    findat = findat[findat['ModularityScore'] !='NA']
    findat['ModularityScore']=findat['ModularityScore'].map(float)
    findat['ClusterSize']=findat['ClusterSize'].map(int)
    findat['bins']=findat['bins'].map(int)
    #bins
    finsampleall = findat.groupby(['bins','Sample'])[['ModularityScore']].median().reset_index()
    finsampleall = finsampleall.sort_values(by=['bins'])
    col = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9']
    cn=0
    fig, ax = plt.subplots()
    for items in sname:
        finind = finsampleall[finsampleall['Sample'] == items]
        cols = col[cn]+'o--'
        ax.plot(finind.bins.tolist(), finind.ModularityScore.round(2).map(float).tolist(), cols, label=items)
        cn=cn+1
    legend = ax.legend(loc='upper right', fontsize='small')
    #legend.get_frame().set_facecolor('C0')
    plt.xlabel('Binning Resolution')
    plt.ylabel('Modularity Score')
    plt.axis([min(binres),max(binres),0.5,1])
    plt.savefig(outdir+'/'+pref+'_allsamplebins_'+str(nres)+'.png')
    plt.close()
    #averagebins
    finsampleall = findat.groupby(['bins'])[['ModularityScore']].median().reset_index()
    finsampleall = finsampleall.sort_values(by=['bins'])
    fins = finsampleall.sort_values(by=['ModularityScore'],ascending=False)
    finlist = fins.iloc[0:3,0].tolist()
    plt.plot(finsampleall.bins.tolist(), finsampleall.ModularityScore.round(2).map(float).tolist(), 'o--')
    plt.xlabel('Binning Resolution')
    plt.ylabel('Modularity Score')
    plt.axis([min(binres),max(binres),0.5,1])
    plt.savefig(outdir+'/'+pref+'_averagebins_'+str(nres)+'.png')
    plt.close()
    #inflation
    inflationall = findat[findat['bins'].isin(finlist)]
    inflatesample = inflationall.groupby(['Inflation','Sample'])[['ModularityScore']].median().reset_index()
    inflatesample = inflatesample.sort_values(by=['Inflation'])
    col = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9']
    cn=0
    fig, ax = plt.subplots()
    for items in sname:
        finind = inflatesample[inflatesample['Sample'] == items]
        cols = col[cn]+'o--'
        ax.plot(finind.Inflation.tolist(), finind.ModularityScore.round(2).map(float).tolist(), cols, label=items)
        cn=cn+1
    legend = ax.legend(loc='upper right', fontsize='small')
    #legend.get_frame().set_facecolor('C0')
    plt.xlabel('Inflation')
    plt.ylabel('Modularity Score')
    plt.axis([min(inflation),max(inflation),0.5,1])
    plt.savefig(outdir+'/'+pref+'_allsampleinflation_'+str(nres)+'.png')
    plt.close()
    #averageinflation
    finsampleall = findat[findat['bins'].isin(finlist)]
    finsampleall = finsampleall.groupby(['Inflation'])[['ModularityScore']].median().reset_index()
    finsampleall = finsampleall.sort_values(by=['Inflation'])
    plt.plot(finsampleall.Inflation.tolist(), finsampleall.ModularityScore.round(2).map(float).tolist(),'o--')
    plt.xlabel('Inflation')
    plt.ylabel('Modularity Score')
    plt.axis([min(inflation),max(inflation),0.5,1])
    plt.savefig(outdir+'/'+pref+'_averageinflation_'+str(nres)+'.png')
    plt.close()
    finaldict = dict()
    finaldict['Inflation'] = finsampleall[finsampleall['ModularityScore'] == max(finsampleall['ModularityScore'])]['Inflation'].tolist()[0]
    finaldict['Bins'] = finlist
    #clusterSize
    finsampleall = findat.groupby(['bins','Sample'])[['ClusterSize']].median().reset_index()
    finsampleall = finsampleall.sort_values(by=['bins'])
    col = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9']
    cn=0
    fig, ax = plt.subplots()
    for items in sname:
        finind = finsampleall[finsampleall['Sample'] == items]
        cols = col[cn]+'o--'
        ax.plot(finind.bins.tolist(), finind.ClusterSize.round(2).map(float).tolist(), cols, label=items)
        cn=cn+1
    legend = ax.legend(loc='upper right', fontsize='small')
    #legend.get_frame().set_facecolor('C0')
    plt.xlabel('Binning Resolution')
    plt.ylabel('ClusterSize')
    plt.axis([min(binres),max(binres),min(finsampleall.ClusterSize),max(finsampleall.ClusterSize)])
    plt.savefig(outdir+'/'+pref+'_allsamplebins_clustersize_'+str(nres)+'.png')
    plt.close()
    #averageClusterSize
    finsampleall = findat.groupby(['bins'])[['ClusterSize']].median().reset_index()
    finsampleall = finsampleall.sort_values(by=['bins'])
    plt.plot(finsampleall.bins.tolist(), finsampleall.ClusterSize.round(2).map(float).tolist(), 'o--')
    plt.xlabel('Binning Resolution')
    plt.ylabel('ClusterSize')
    plt.axis([min(binres),max(binres),min(finsampleall.ClusterSize),max(finsampleall.ClusterSize)])
    plt.savefig(outdir+'/'+pref+'_averagebins_clustersize_'+str(nres)+'.png')
    plt.close()
    with open(args.odir+'/finalparameters_'+str(nres)+'_'+args.pref+'.txt', 'w') as f:
        print(finaldict, file=f)
    pickle.dump(findat, open(outdir+"/"+pref+'_'+str(nres)+"_parametermcl.pickle", "wb"))