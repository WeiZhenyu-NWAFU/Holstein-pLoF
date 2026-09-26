"""Autosomal genotype-count homozygote-depletion screen.

Core statistics extracted unchanged from the study analysis; CLI and input
validation are packaging additions. ALT must be the pLoF allele for pLoF rows.
"""
import argparse
import math
from pathlib import Path
import numpy as np
import pandas as pd

def bh(p):
 p=np.asarray(p,float);order=np.argsort(p);out=np.ones(len(p));out[order]=np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1].clip(0,1);return out

def stats(r,h,o):
 n=r+h+o
 if n<=0:return [np.nan]*5
 af=(h+2*o)/(2*n);e=n*af*af
 def ll(k,p):return 0 if k==0 else k*math.log(p) if p>0 else -math.inf
 l0=ll(o,af*af)+ll(n-o,1-af*af);l1=ll(o,o/n)+ll(n-o,1-o/n)
 lr=max(0.,2*(l1-l0));return af,e,lr,math.erfc(math.sqrt(lr/2)),o/e if e>0 else np.nan

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',required=True,help='TSV or TSV.GZ genotype counts')
    parser.add_argument('--output',required=True,help='New output TSV; existing files are not overwritten')
    args=parser.parse_args()
    output=Path(args.output)
    if output.exists(): raise ValueError('Output already exists')
    d=pd.read_csv(args.input,sep='\t',dtype={'chrom':str})
    needed=['variant','chrom','cohort','variant_class','n_refhom','n_het','n_althom','n_called','call_rate']
    if not set(needed).issubset(d.columns): raise ValueError('Required columns: '+', '.join(needed))
    if d[['variant','chrom','cohort','variant_class']].isna().any().any(): raise ValueError('Missing keys')
    if d.duplicated(['variant','cohort','variant_class']).any(): raise ValueError('Duplicate variant/cohort/class records')
    countcols=['n_refhom','n_het','n_althom','n_called']
    for c in countcols+['call_rate']: d[c]=pd.to_numeric(d[c],errors='raise')
    counts=d[countcols].to_numpy()
    if not np.isfinite(counts).all() or (counts<0).any() or (counts%1!=0).any(): raise ValueError('Invalid genotype counts')
    if not d.call_rate.between(0,1).all(): raise ValueError('Invalid call rate')
    if not (d.n_refhom+d.n_het+d.n_althom==d.n_called).all(): raise ValueError('Inconsistent genotype counts')
    d=d[d.chrom.isin([str(i) for i in range(1,30)]) & d.n_called.gt(0)].copy()
    if d.empty: raise ValueError('No callable autosomal records')
    d[['ALT_AF','expected_HOM','LRT','P','O_E']]=np.array([stats(*x) for x in d[['n_refhom','n_het','n_althom']].itertuples(index=False,name=None)])
    d['FDR']=d.groupby(['cohort','variant_class']).P.transform(lambda s:bh(s.to_numpy()))
    d['observed_HOM']=d.n_althom
    d['depleted']=d.observed_HOM.lt(d.expected_HOM)&d.FDR.lt(.05)
    d['informative_absence']=d.observed_HOM.eq(0)&d.expected_HOM.ge(3)&d.call_rate.ge(.95)
    d['significant_informative_absence']=d.informative_absence&d.depleted
    d.to_csv(output,sep='\t',index=False,na_rep='NA')

if __name__=='__main__':main()
