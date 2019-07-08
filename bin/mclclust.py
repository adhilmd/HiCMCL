#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 22 15:36:30 2019

@author: amohamme
"""

import argparse

def argparser():
    parser = argparse.ArgumentParser(prog = '''\nmclclust.py''', description='''\n----------MCL clustering--------- \n
    \n[Date: 7th May 2018], \n[help: python mclclust.py -h]\n''', usage = 'mclclust.py *args')
    parser.add_argument('-findb','--findb', type=str, dest='findb', help="Comma seperated database file (Single file can also be provided) (Mandatory)", action = 'store', required = True)
    parser.add_argument('-chrl','--chrl', type=str, dest='chrl', help="Chromosome length file containing two columns (chromosomenumber, length) without header (Mandatory)", action = 'store', required = True)
    parser.add_argument('-icore','--icore', type=int, dest='icore', help="Core resolution for clustering, in multiples (default=10)", action = 'store', default=10)
    parser.add_argument('-norm','--norm', type=str, dest='norm', help="Normalization type 'median' or 'minmax' (default=minmax)", action = 'store', default="minmax")
    parser.add_argument('-cpath','--cpath', type=str, dest='cpath', help="Main path were the codes are present for modules", action='store', required=True)
    parser.add_argument('-pfile','--pfile', type=str, dest='pfile', help="parameter file containing resolution and inflation information from parametermcl script", action='store', required=True)
    parser.add_argument('-th','--th', type=int, dest='threads', help="Number of threads (default=8)", action = 'store', default=8)
    parser.add_argument('-tag','--tag', type=str, dest='tag', help="Comma seperated file tags for pickle file (single tag can also be provided), the tags should match the -bedpe (Mandatory)", action = 'store', required = True)
    parser.add_argument('-odir','--odir', type=str, dest='odir', help="Outdir (Mandatory)", action = 'store', required = True)
    args = parser.parse_args()
    return(args)

args = argparser()
import sys
sys.path.append(args.cpath)
from clusteringmodules import datanormalization 
from clusteringmodules import clusterpred
from clusteringmodules import sorted_nicely
import os
import sqlite3
import pandas as pd 
from multiprocessing import Pool
from functools import partial
import warnings
warnings.filterwarnings("ignore")
from pybedtools import BedTool
from collections import defaultdict
import pickle


def biningreq(chromnum,findb,sname,initsres,chrbed,clg,normtype,inflation):
    connex = sqlite3.connect(findb)
    cur = connex.cursor()
    sql = "SELECT chr1,((loc1 + " + str(initsres-1) + ")" + "/" + str(initsres) + ")*" +str(initsres) + " as loc1 " + ",chr2,((loc2 + " + str(initsres-1) + ")" + "/" + str(initsres) + ")*" +str(initsres) + " as loc2 " + ",COUNT(*) as Count FROM " + sname + "_intra WHERE chr1 = '" +chromnum +"' GROUP BY " + "chr1,((loc1 + " + str(initsres-1) + ")" + "/" + str(initsres) + ")*" +str(initsres) + ",chr2,((loc2 + " + str(initsres-1) + ")" + "/" + str(initsres) + ")*" +str(initsres) + ";"
    df = pd.read_sql(sql, connex)
    cur.close()
    connex.close()
    a = BedTool()
    windows = a.window_maker(b=chrbed, w=initsres, s=initsres)
    wdfall = pd.read_table(windows.fn, names=['chrom', 'start', 'end'])
    wdfall['chrom'] = wdfall['chrom'].map(str)
    wdfall['id'] = wdfall.groupby(['chrom']).cumcount()
    wdf = wdfall[wdfall['chrom'] == str(chromnum)]
    cnt = wdf.shape[0]
    df = df.merge(wdf[['chrom','start','id']], left_on=['chr1','loc1'], right_on=['chrom','start'], how = 'left')
    df = df.dropna()
    df = df.merge(wdf[['chrom','start','id']], left_on=['chr2','loc2'], right_on=['chrom','start'], how = 'left')
    df = df.dropna()
    df = df[['chr1','loc1','chr2','loc2','Count','id_x','id_y']]
    df['Distance'] = df['loc2'] - df['loc1']
    dist = initsres*10
    cscmat,dfi = datanormalization(df,normtype,dist,cnt)
    clusters,finclust,nwlist,allss = clusterpred(cscmat, wdf, clg, inflation)
    return(dfi,chromnum,wdf,finclust,allss,nwlist)


if __name__ == "__main__":
    args = argparser()
    outdir = args.odir
    findblist = args.findb.split(',')
    taglist = args.tag.split(',')
    chrl = args.chrl
    normtype = args.norm
    cpath = args.cpath
    pfile = args.pfile
    procs = args.threads
    clg = 3 
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    if (len(findblist) != len(taglist)):
        print ("Error: Please provide same number of -findb and -tag files")
        quit()
    te = eval(open(pfile, 'r').read())
    wc = pd.read_csv(chrl, sep="\t", header=None)
    wc.columns = ['chrom','end']
    wc['start'] = 0
    wc['order'] = [i+1 for i in range(0,wc.shape[0])]
    wc = wc[['chrom','start','end','order']]
    wc1 = wc[['chrom','start','end']]
    wbed = BedTool.from_dataframe(wc1)
    chromlist = sorted_nicely(list(set(wc['chrom'].tolist())))
    lsizes = te['Bins']
    inflation = te['Inflation']
    biningdict = defaultdict(dict)
    wdfdict = defaultdict(dict)
    fcdict = defaultdict(dict)
    pdict = defaultdict(dict)
    mmdict = defaultdict(dict)
    lit = 0
    for lit in range(0,len(findblist)):
        for i in lsizes:
            if i != 0:
                pool = Pool(processes=procs)
                mrt = pool.map(partial(biningreq,findb=findblist[lit],sname=taglist[lit],initsres=i,chrbed=wbed,clg=clg,normtype=normtype,inflation=inflation),chromlist)
                pool.close()
                pool.join()
                binmrt = [seq[0] for seq in mrt]
                chrmrt = [seq[1] for seq in mrt]
                wdfmrt = [seq[2] for seq in mrt]
                finclustmrt = [seq[3] for seq in mrt]
                allcmrt = [seq[4] for seq in mrt]
                mmmrt = [seq[5] for seq in mrt]
                for j in range(0,len(binmrt)):
                    chrval = chrmrt[j]
                    biningdict[i][chrval] = binmrt[j]
                    wdfdict[i][chrval] = wdfmrt[j]
                    fcdict[i][chrval] = finclustmrt[j]
                    pdict[i][chrval] = allcmrt[j]
                    mmdict[i][chrval] = mmmrt[j]
        alldataclust = pd.DataFrame()
        for lsize in lsizes:
            for chrom in chromlist:
                temp = pd.DataFrame(mmdict[lsize][chrom])
                temp.columns = ['start','end']
                temp['chrom'] = chrom
                temp['bsize'] = lsize
                alldataclust = pd.concat([alldataclust,temp])
    alldataclust = alldataclust[['chrom','start','end','bsize']]
    rsize = alldataclust[alldataclust['bsize']==lsizes[0]]
    rsize = rsize[['chrom','start','end','bsize']]
    rsize['sizebins'] = rsize['end']-rsize['start']
    rsize['id'] = list(range(rsize.shape[0]))
    x = BedTool.from_dataframe(rsize)
    clustalldict = dict()
    cnt = 0
    for lsize in lsizes[1:]:
        tclust=alldataclust[alldataclust['bsize']==lsize]
        tclust['id'] = list(range(tclust.shape[0]))
        tclust = tclust[['chrom','start','end','bsize','id']]
        y = BedTool.from_dataframe(tclust)
        mpydat = x.intersect(y,wao=True)
        mdf = pd.read_table(mpydat.fn)
        mdf = mdf.iloc[:,0:6]
        mdf.columns=['chr','start','end','bsize','length','id']
        mdf=mdf.groupby(mdf.columns.tolist()).agg({"id": "count"})
        mdf.columns=[lsize]
        mdf=mdf.reset_index()
        clustalldict[lsize] = mdf
        if (cnt == 0):
            finmerge = mdf.copy()
        else:
            finmerge = pd.merge(finmerge,mdf.copy())
        cnt = cnt+1
    #breakpoints
    bclust=alldataclust.copy()
    bclust = bclust[['chrom','start','end','bsize']]
    bclust['start'] = bclust['end']-bclust['bsize']
    bclust = bclust.sort_values(by=['chrom','start'])
    y = BedTool.from_dataframe(bclust)
    mp = y.merge(c='4,4',o='distinct,count_distinct')
    mdf = pd.read_table(mp.fn)
    mdf.columns=['chr','start','end','bnsize','bncount']
    mdf['ids'] = list(range(mdf.shape[0]))
    ydatf = mdf.sort_values(by=['bncount'],ascending=False)
    x1 = BedTool.from_dataframe(ydatf)
    bds = bclust[bclust['bsize'] == lsizes[len(lsizes)-1]]
    z = BedTool.from_dataframe(bds)
    mp = z.intersect(x1,wo=True)
    zdatf = pd.read_table(mp.fn)
    zdatf.columns=['chrx','startx','endx','bnsizex','chr','start','end','bnsize','bncount','id','overlap'] 
    zdatf = zdatf.sort_values(by=['bncount'],ascending=False)
    findata = dict()
    findata['allclust'] = alldataclust
    findata['wdfbins'] = wdfdict
    findata['zdatf'] = zdatf
    findata['finmerge'] = finmerge
    findata['allbins'] = biningdict
    pickle.dump(findata, open(outdir+"/"+taglist[lit]+"_cluster.pickle", "wb"))
