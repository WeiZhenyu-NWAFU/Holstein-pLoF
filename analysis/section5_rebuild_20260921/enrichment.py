"""Outcome-independent enrichment universe and explicit matched sensitivities."""
from pathlib import Path
from collections import defaultdict
from functools import lru_cache
import sys, gzip, re, math, json
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent;sys.path.insert(0,str(BASE/'fig4_final/vendor'))
import numpy as np
import pandas as pd
Q=ROOT/'tables';DATA=BASE/'fig4_v3_final/data'
def read(p,**kw):return pd.read_csv(p,sep='\t',low_memory=False,na_values=['NE'],**kw)
def write(d,n):d.to_csv(Q/n,sep='\t',index=False,na_rep='NE')
def bh(p):
 p=np.asarray(p,float);order=np.argsort(p);v=np.ones(len(p));v[order]=np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1].clip(0,1);return v
def hyper(M,K,n,k):
 if k<=0:return 1.
 def lc(a,b):return math.lgamma(a+1)-math.lgamma(b+1)-math.lgamma(a-b+1)
 return min(1.,sum(math.exp(lc(K,j)+lc(M-K,n-j)-lc(M,n)) for j in range(max(k,n-(M-K)),min(K,n)+1)))
d=read(Q/'01_depletion_candidate_evidence_matrix.tsv');syn=read(Q/'01_autosomal_synonymous.tsv.gz')
# Match actual ALT/pLoF AF, not folded MAF. Maximum 20 controls within fixed calipers.
s=syn[syn.call_rate.ge(.95)].sort_values(['pLoF_AF','variant']).reset_index(drop=True)
sv=s.pLoF_AF.to_numpy();sn=s.n_called.to_numpy();rng=np.random.default_rng(20260921)
events=(s.observed_HOM.eq(0)&s.expected_HOM.ge(3)).to_numpy(float)
rows=[];matches=[]
for r in d[d.call_rate.ge(.95)].itertuples():
 a,b=np.searchsorted(sv,[r.pLoF_AF-.005,r.pLoF_AF+.005]);ii=np.arange(a,b);ii=ii[np.abs(sn[ii]-r.n_called)<=10]
 if not len(ii):continue
 dist=(np.abs(sv[ii]-r.pLoF_AF)/.005)**2+(np.abs(sn[ii]-r.n_called)/10)**2
 ii=ii[np.argsort(dist,kind='stable')[:20]]
 rows.append(dict(variant=r.variant,event=int(r.informative_absence),control_event_probability=float(events[ii].mean()),controls=len(ii),max_AF_difference=float(np.abs(sv[ii]-r.pLoF_AF).max()),max_called_difference=int(np.abs(sn[ii]-r.n_called).max())))
 matches.extend(dict(pLoF_variant=r.variant,synonymous_variant=s.iloc[j].variant,control_event=int(events[j]),AF_difference=float(abs(sv[j]-r.pLoF_AF)),called_difference=int(abs(sn[j]-r.n_called))) for j in ii)
m=pd.DataFrame(rows);cp=m.control_event_probability.to_numpy();o=m.event.to_numpy();N=len(m);a=o.sum();c=cp.sum()
boot=[];br=[]
for _ in range(2500):
 ix=rng.integers(N,size=N);aa=o[ix].sum();cc=cp[ix].sum()
 if aa>0 and cc>0 and aa<N and cc<N:boot.append((aa/(N-aa))/(cc/(N-cc)));br.append(aa/cc)
null=np.empty(20000,int)
for i in range(0,len(null),200):null[i:i+200]=(rng.random((min(200,len(null)-i),N))<cp).sum(axis=1)
effect=dict(definition='Autosomal; O=0; E>=3; call rate>=0.95',matching='Actual ALT AF caliper 0.005; called N caliper 10; nearest 20 or all available; control reuse allowed',pLoF_numerator=int(a),pLoF_denominator=N,pLoF_proportion=float(a/N),synonymous_expected_n=float(c),synonymous_effective_denominator=N,synonymous_proportion=float(c/N),enrichment_ratio=float(a/c),enrichment_CI_low=float(np.quantile(br,.025)),enrichment_CI_high=float(np.quantile(br,.975)),odds_ratio=float((a/(N-a))/(c/(N-c))),OR_CI_low=float(np.quantile(boot,.025)),OR_CI_high=float(np.quantile(boot,.975)),empirical_P=float((1+(null>=a).sum())/(len(null)+1)),resamples=len(null),bootstrap_replicates=len(boot),unmatched_pLoFs=int(d.call_rate.ge(.95).sum()-N),unmatched_informative=int(d.informative_absence.sum()-a),unique_control_variants=pd.DataFrame(matches).synonymous_variant.nunique(),max_AF_difference=float(m.max_AF_difference.max()),max_called_difference=int(m.max_called_difference.max()))
write(pd.DataFrame([effect]),'01_matched_synonymous_effect.tsv');write(m,'01_synonymous_matching_units.tsv');write(pd.DataFrame(matches),'01_synonymous_matches.tsv.gz');write(pd.DataFrame({'null_event_count':null}),'01_synonymous_null_resamples.tsv.gz')
print('MATCHED',effect,flush=True)
# Gene identifier mapping and term propagation.
gm={}
for line in (DATA/'ARS-USD1.2_ENSBTAG-gene_symbol').read_text().splitlines():
 p=line.split()
 if len(p)>=2:gm[p[0].split('.')[0]]=p[1]
bg=set(d.gene_id.dropna());symids=defaultdict(set)
for gid in bg:
 sym=gm.get(gid)
 if sym and not sym.startswith('ENSBTAG'):symids[sym.upper()].add(gid)
info={};parents=defaultdict(set);altid={};obsolete=set();cur=None
with (DATA/'go-basic.obo').open(encoding='utf-8') as f:
 for line in f:
  line=line.strip()
  if line.startswith('['):cur=None
  elif line.startswith('id: GO:'):cur=line[4:];info[cur]={}
  elif cur and line.startswith('name: '):info[cur]['name']=line[6:]
  elif cur and line.startswith('namespace: '):info[cur]['namespace']=line[11:]
  elif cur and line.startswith('is_a: GO:'):parents[cur].add(line.split()[1])
  elif cur and line.startswith('relationship: part_of GO:'):parents[cur].add(line.split()[2])
  elif cur and line.startswith('alt_id: '):altid[line[8:]]=cur
  elif cur and line=='is_obsolete: true':obsolete.add(cur)
@lru_cache(None)
def ancestors(t):return frozenset({t}|set().union(*(ancestors(p) for p in parents[t]))) if parents[t] else frozenset({t})
terms=defaultdict(set);names={};not_n=0;ambig=set()
with gzip.open(DATA/'goa_cow.gaf.gz','rt',encoding='utf-8') as f:
 for line in f:
  if line.startswith('!'):continue
  p=line.rstrip('\n').split('\t')
  if len(p)<13:continue
  if 'NOT' in p[3].split('|'):not_n+=1;continue
  if p[12].split('|')[0]!='taxon:9913':continue
  ids=symids.get(p[2].upper(),set())
  if len(ids)!=1:
   if len(ids)>1:ambig.add(p[2])
   continue
  t=altid.get(p[4],p[4])
  if t not in info or t in obsolete:continue
  for a in ancestors(t):terms['GO',a]|=ids
for src,t in terms:names[src,t]=(info[t].get('name',t),info[t].get('namespace',''))
with (DATA/'Ensembl2Reactome.txt').open(encoding='utf-8') as f:
 for line in f:
  p=line.rstrip('\n').split('\t')
  if len(p)<6 or p[5]!='Bos taurus':continue
  gid=p[0].split('.')[0]
  if gid in bg:terms['Reactome',p[1]].add(gid);names['Reactome',p[1]]=(p[3],'pathway')
annot={s:set().union(*(v for (src,t),v in terms.items() if src==s)) for s in ['GO','Reactome']}
orids=terms.get(('GO','GO:0004984'),set())|{g for g in bg if re.match(r'^OR\d',gm.get(g,''))}
genes=lambda mask:set(d.loc[mask,'gene_id'].dropna())
selection_bg=genes((d.DensityStratifiedTailP.notna()|d.logLR.notna())&d.adult_HOM_observed_not_depleted)
specs={'Primary':(genes(d.informative_absence),bg),'Shared robust':(genes(d.shared_robust),bg),'Primary excluding OR':(genes(d.informative_absence)-orids,bg-orids),'Shared excluding OR':(genes(d.shared_robust)-orids,bg-orids),'Selection with adult homozygotes':(genes(d.selection_supported&d.adult_HOM_observed_not_depleted),selection_bg)}
counts=d.groupby('gene_id').size().to_dict();afs=d.groupby('gene_id').pLoF_AF.max().to_dict()
stratum=lambda g:(min(4,int(math.log2(counts[g]))),int(np.digitize(afs[g],[.01,.05,.1,.2,.5])))
out=[];coverage=[];perms=[];diagnostics=[]
for label,(fg,bkg) in specs.items():
 for src in ['GO','Reactome']:
  b=bkg&annot[src];a=fg&b;M=len(b);n=len(a)
  coverage.append(dict(analysis=label,source=src,foreground_genes=len(fg),background_genes=len(bkg),annotated_foreground=n,annotated_background=M))
  rows=[]
  for (source,t),members in terms.items():
   if src!=source:continue
   g=members&b;K=len(g);k=len(a&g)
   if K<5 or K==M:continue
   cells=np.array([k,n-k,K-k,M-n-K+k],float);aa,bb,cc,dd=cells
   OR=aa*dd/(bb*cc) if bb*cc else (math.inf if aa*dd else np.nan)
   if np.any(cells==0):cells+=.5
   aa,bb,cc,dd=cells;cor=aa*dd/(bb*cc);se=np.sqrt((1/cells).sum());lo,hi=np.exp(np.log(cor)+np.array([-1,1])*1.96*se)
   rows.append(dict(analysis=label,source=src,term_id=t,term_name=names[src,t][0],namespace=names[src,t][1],foreground_overlap=k,foreground_N=n,background_overlap=K,background_N=M,gene_ratio=k/n if n else np.nan,OR=OR,CI_low=lo,CI_high=hi,P=hyper(M,K,n,k),gene_ids=';'.join(sorted(a&g)),genes=';'.join(gm.get(v,v) for v in sorted(a&g))))
  z=pd.DataFrame(rows);z['BH_FDR']=bh(z.P);z['tested_terms']=len(z);out.append(z)
  # Marginal stratified hypergeometric resampling for every eligible term.
  # Exact distribution convolution avoids allocating B x genes x terms arrays;
  # 20,000 draws from this conditional distribution give the permutation P.
  pools=defaultdict(set)
  for g in b:pools[stratum(g)].add(g)
  needs={s:len(pool&a) for s,pool in pools.items()}
  for s,pool in sorted(pools.items()):diagnostics.append(dict(analysis=label,source=src,count_bin=s[0],maxAF_bin=s[1],pool_N=len(pool),foreground_N=needs[s],fixed_stratum=needs[s] in [0,len(pool)]))
  @lru_cache(None)
  def hpmf(N,K,n):
   lo=max(0,n-(N-K));hi=min(n,K);v=np.zeros(n+1)
   def lc(x,y):return math.lgamma(x+1)-math.lgamma(y+1)-math.lgamma(x-y+1)
   for k in range(lo,hi+1):v[k]=math.exp(lc(K,k)+lc(N-K,n-k)-lc(N,n))
   return v/v.sum()
  for r in z.itertuples():
   members=terms[src,r.term_id];dist=np.array([1.]);mean=0.
   for s,pool in pools.items():
    nn=needs[s]
    if nn:dist=np.convolve(dist,hpmf(len(pool),len(pool&members),nn));mean+=nn*len(pool&members)/len(pool)
   tail=float(min(1,max(0,dist[int(r.foreground_overlap):].sum())))
   # Count of upper-tail exceedances in B independent exact conditional draws.
   exceed=int(rng.binomial(20000,tail));perms.append(dict(analysis=label,source=src,term_id=r.term_id,term_name=r.term_name,observed_overlap=r.foreground_overlap,matched_expected_overlap=mean,conditional_exact_P=tail,empirical_P=(exceed+1)/20001,resamples=20000,method='Stratified gene-set randomization: exact convolution; binomial MC exceedances',covariates='pLoF count 1/2-3/4-7/8-15/16+; max AF <.01/.01-.05/.05-.1/.1-.2/.2-.5/>=.5'))
  print(label,src,n,M,'FDR sig',int(z.BH_FDR.lt(.05).sum()),flush=True)
en=pd.concat(out,ignore_index=True).sort_values(['analysis','source','BH_FDR','P']);pp=pd.DataFrame(perms)
pp['matched_BH_FDR']=pp.groupby(['analysis','source']).conditional_exact_P.transform(lambda x:bh(x))
pp['empirical_BH_FDR']=pp.groupby(['analysis','source']).empirical_P.transform(lambda x:bh(x))
write(en,'03_enrichment_all.tsv.gz');write(en[en.analysis.eq('Primary')],'03_enrichment_primary.tsv');write(en[en.analysis.eq('Shared robust')],'03_enrichment_robust.tsv');write(en[en.analysis.str.contains('excluding OR')],'03_enrichment_no_OR.tsv');write(pp,'03_enrichment_permutation.tsv.gz');write(pd.DataFrame(coverage),'03_annotation_coverage.tsv');write(pd.DataFrame(diagnostics),'03_permutation_strata.tsv')
write(pd.DataFrame([dict(analysis=label,gene_id=g,gene_symbol=gm.get(g,g),foreground=g in a,olfactory_receptor=g in orids,screened_pLoF_count=counts[g],maximum_pLoF_AF=afs[g]) for label,(a,b) in specs.items() for g in sorted(b)]),'03_gene_sets_and_backgrounds.tsv.gz')
summary={'matched_synonymous':effect,'annotation_coverage':coverage,'GO_NOT_excluded':not_n,'ambiguous_symbols_excluded':sorted(ambig),'OR_genes_background':len(orids),'primary_OR_genes':len(specs['Primary'][0]&orids),'shared_OR_genes':len(specs['Shared robust'][0]&orids),'functional_summary':en.groupby(['analysis','source']).apply(lambda g:int(g.BH_FDR.lt(.05).sum())).to_dict().__str__()}
(ROOT/'enrichment_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
