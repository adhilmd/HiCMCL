#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 22 13:39:54 2019

@author: adhil
"""

import argparse

def argparser():
    parser = argparse.ArgumentParser(prog = '''\nsql_dump.py''', description='''\n----------SQL dump of bedpe file--------- \n
    \n[Date: 7th May 2018], \n[help: python sql_dump.py -h]\n''', usage = 'sql_dump.py *args')
    parser.add_argument('-bedpe','--bedpe', type=str, dest='bedpe', help="comma seperated multiple bed file paths (single file can also be provided) containing paired interaction without header (Mandatory)", action = 'store', required = True)
    parser.add_argument('-chrl','--chrl', type=str, dest='chrl', help="chromosome length file containing two columns (chromosomenumber, length) without header (Mandatory)", action = 'store', required = True)
    parser.add_argument('-cks','--cks', type=int, dest='cks', help="Chunk size for dumping to sql (for single iteration) (default=10000000)", action = 'store', default=10000000)
    parser.add_argument('-mds','--mds', type=int, dest='mds', help="minimum intra chromosome distance (default=1000)", action = 'store', default=1000)
    parser.add_argument('-tag','--tag', type=str, dest='tag', help="comma seperated file tags (single tag can also be provided), the tags should match the -bedpe (Mandatory)", action = 'store', required = True)
    parser.add_argument('-odir','--odir', type=str, dest='odir', help="outdir (Mandatory)", action = 'store', required = True)
    args = parser.parse_args()
    return(args)

import sqlite3
import os
import pandas as pd 
import warnings
warnings.filterwarnings("ignore")

def sqldump(filepath,db,chunksize,mindis,sname,wc):
    connex = sqlite3.connect(db)
    for chunk in pd.read_csv(filepath, chunksize=chunksize, header=None, sep="\t"):
        chunk = chunk[chunk[0].isin(wc['chrom'].tolist())]
        chunk = chunk[chunk[3].isin(wc['chrom'].tolist())]
        chunk['p1'] = chunk[1] + (chunk[2]-chunk[1])/2
        chunk['p2'] = chunk[4] + (chunk[5]-chunk[4])/2
        chunk = chunk[[0,'p1',3,'p2']]
        chunk['p1'] = chunk['p1'].round().map(int)
        chunk['p2'] = chunk['p2'].round().map(int)
        chunk.columns=['chr1','loc1','chr2','loc2']
        intra = chunk[chunk['chr1']==chunk['chr2']]
        inter = chunk[chunk['chr1']!=chunk['chr2']]
        cond = intra[intra.loc1 > intra.loc2].index
        if len(cond) >= 1:
            intra.loc[cond, ['loc1', 'loc2']] = intra.loc[cond, ['loc2', 'loc1']]    
        intra['dis']=intra['loc2']-intra['loc1']
        intra = intra[intra['dis']>=mindis]
        intra.to_sql(name=sname+"_intra", con=connex, if_exists="append", index=False)
        inter.to_sql(name=sname+"_inter", con=connex, if_exists="append", index=False)
    connex.close()
    return 0

if __name__ == "__main__":
    args = argparser()
    bedpelist = args.bedpe.split(',')
    taglist = args.tag.split(',')
    if not os.path.exists(args.odir):
        os.makedirs(args.odir)
    if (len(bedpelist) != len(taglist)):
        print ("Error: Please provide same number of -bedpe and -tag files")
        quit()
    wc = pd.read_csv(args.chrl, sep="\t", header=None)
    wc.columns = ['chrom','end']
    wc['start'] = 0
    wc['order'] = [i+1 for i in range(0,wc.shape[0])]
    wc = wc[['chrom','start','end','order']]
    for i in range(0,len(bedpelist)):
        findb = args.odir+'/'+taglist[i]+'_fin.sqlite'
        rt = sqldump(bedpelist[i],findb,args.cks,args.mds,taglist[i],wc)