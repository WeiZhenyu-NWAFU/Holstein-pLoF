"""Autosomal reanalysis from immutable, allele-specific genotype-count sources."""
from pathlib import Path
import sys, math, json, hashlib
ROOT=Path(__file__).resolve().parent; BASE=ROOT.parent
sys.path.insert(0,str(BASE/'fig4_final/vendor'))
import numpy as np
import pandas as pd
Q=ROOT/'tables';Q.mkdir(parents=True,exist_ok=True)
def read(p,**kw):return pd.read_csv(p,sep='\t',low_memory=False,**kw)
def write(d,n):d.to_csv(Q/n,sep='\t',index=False,na_rep='NE')
def bh(p):
 p=np.asarray(p,float);order=np.argsort(p);out=np.ones(len(p));out[order]=np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1].clip(0,1);return out
def stats(r,h,o):
 n=r+h+o
 if n<=0:return [np.nan]*5
 af=(h+2*o)/(2*n);e=n*af*af
 def ll(k,p):return 0 if k==0 else k*math.log(p) if p>0 else -math.inf
 l0=ll(o,af*af)+ll(n-o,1-af*af);l1=ll(o,o/n)+ll(n-o,1-o/n)
 lr=max(0.,2*(l1-l0));return af,e,lr,math.erfc(math.sqrt(lr/2)),o/e if e>0 else np.nan
def key(d):return d.chrom.astype(str)+':'+d.pos.astype(str)+':'+d.ref+':'+d.alt
def flag(s):return s.astype(str).str.lower().isin(['true','1'])
allv=read(BASE/'fig4_final/data/variant_level.charlier_lrt.tsv.gz',dtype={'chrom':str},na_values=['.'])
allv['variant']=key(allv);allv['site']=allv.chrom+'_'+allv.pos.astype(str)
old=allv[(allv.variant_class=='pLoF')].copy()
audit=[]
for c,g in old.groupby('cohort'):
 dep=(g.direction=='depletion')&(g.fdr_bh_within_class_cohort<.05);event=(g.observed_althom==0)&(g.expected_althom_hwe>=3)&(g.call_rate>=.95)
 audit.append(dict(cohort=c,old_tested=len(g),X_tested=int(g.chrom.eq('X').sum()),old_depleted=int(dep.sum()),X_depleted=int((dep&g.chrom.eq('X')).sum()),old_informative_absence=int(event.sum()),X_informative_absence=int((event&g.chrom.eq('X')).sum())))
auto=allv[allv.chrom.isin([str(i) for i in range(1,30)]) & allv.n_called.gt(0)].copy()
assert (auto.n_refhom+auto.n_het+auto.n_althom==auto.n_called).all()
calc=np.array([stats(r,h,o) for r,h,o in auto[['n_refhom','n_het','n_althom']].itertuples(index=False,name=None)])
max_lrt_difference=float(np.max(np.abs(calc[:,2]-auto.lrt_charlier)))
auto['pLoF_AF']=calc[:,0];auto['expected_HOM']=calc[:,1];auto['LRT']=calc[:,2];auto['P']=calc[:,3];auto['O_E']=calc[:,4]
auto['FDR']=auto.groupby(['cohort','variant_class']).P.transform(lambda s:bh(s.to_numpy()))
auto['observed_HOM']=auto.n_althom;auto['depleted']=auto.observed_HOM.lt(auto.expected_HOM)&auto.FDR.lt(.05)
auto['informative_absence']=auto.observed_HOM.eq(0)&auto.expected_HOM.ge(3)&auto.call_rate.ge(.95)
auto['significant_informative_absence']=auto.informative_absence&auto.depleted
master=read(BASE/'manuscript_supporting_information_final/phase1/variant_master.tsv',dtype={'chromosome':str})
assert master.variant_id.is_unique and (master.pLoF_allele==master.ALT).all()
meta=['variant_id','gene_id','gene_symbol','canonical_transcript','consequence','HGVS_c','HGVS_p','annotation_sources','annotation_support','EVA','1000_Bull_Genomes','RNA_support','three_tissue_RNA_support','milk_RNA_support','NMD_class','stringent_priority14']
plof=auto[auto.variant_class.eq('pLoF')].merge(master[meta],left_on='variant',right_on='variant_id',validate='many_to_one')
gm={}
with (BASE/'fig4_v3_final/data/ARS-USD1.2_ENSBTAG-gene_symbol').open() as f:
 for line in f:
  p=line.split()
  if len(p)>=2:gm[p[0].split('.')[0]]=p[1]
plof['gene_id']=plof.gene_id.str.strip().str.split('.').str[0]
plof['gene_symbol']=plof.gene_id.map(gm).fillna(plof.gene_symbol).str.strip()
full=plof[plof.cohort.eq('full574')].copy();un=plof[plof.cohort.eq('unrelated_PIHAT_lt_0.125')].copy()
uncols=['variant','n_called','pLoF_AF','observed_HOM','expected_HOM','O_E','P','FDR','call_rate','depleted','informative_absence']
full=full.merge(un[uncols].rename(columns={c:'unrelated_'+c for c in uncols if c!='variant'}),on='variant',validate='one_to_one')
full['shared_robust']=full.informative_absence&full.unrelated_informative_absence
full['pLoF_allele']=full.alt;full['RNA_supported']=flag(full.RNA_support)
full['historical14']=flag(full.stringent_priority14)
full['adult_HOM_observed_not_depleted']=full.observed_HOM.gt(0)&~full.depleted
full['fate']=np.select([full.depleted&full.observed_HOM.eq(0),full.depleted,full.observed_HOM.gt(0)],['Depleted; no homozygotes','Partial depletion','Homozygotes observed; not depleted'],default='No homozygotes; not depleted')
# Existing external summary: exact allele matching. Do not treat missing records as REF/REF.
raw=read(BASE/'fig4_v3_final/data/larger_holstein_replication_1148.raw.tsv',dtype={'chrom':str},na_values=['.'])
raw['variant']=key(raw);assert raw.variant.is_unique
for c in ['NS','AN','AC','AF','AC_Hom','AC_Het','F_MISSING']:raw[c]=pd.to_numeric(raw[c],errors='coerce')
raw['raw_external_record']=True
ext=full.merge(raw[['variant','NS','AN','AC','AF','AC_Hom','AC_Het','F_MISSING','raw_external_record']],on='variant',how='left',validate='one_to_one')
ext['raw_external_record']=ext.raw_external_record.fillna(False).astype(bool)
ext['larger_observed_HOM']=ext.AC_Hom/2;ext['larger_het']=ext.AC_Het;ext['larger_refhom']=ext.NS-ext.larger_het-ext.larger_observed_HOM
ext['larger_counts_valid']=ext.NS.gt(0)&ext.AN.eq(2*ext.NS)&ext.AC.eq(ext.AC_Hom+ext.AC_Het)&ext.larger_refhom.ge(0)&ext.AC_Hom.mod(2).eq(0)
ext['larger_call_rate']=ext.NS/1148
ext['larger_evaluable']=ext.larger_counts_valid&ext.larger_call_rate.ge(.95)
for col in ['larger_AF','larger_expected_HOM','larger_LRT','larger_P','larger_O_E']:ext[col]=np.nan
ix=ext.larger_counts_valid
ext.loc[ix,['larger_AF','larger_expected_HOM','larger_LRT','larger_P','larger_O_E']]=np.array([stats(*r) for r in ext.loc[ix,['larger_refhom','larger_het','larger_observed_HOM']].itertuples(index=False,name=None)])
ext['larger_status']=np.select([~ext.raw_external_record.fillna(False),~ext.larger_counts_valid,ext.larger_call_rate.lt(.95),ext.larger_observed_HOM.gt(0),ext.larger_expected_HOM.ge(3)],['Not evaluable','Not evaluable','Not evaluable','Homozygotes observed','Informative absence'],default='Underpowered absence')
ext['larger_NE_reason']=np.select([~ext.raw_external_record.fillna(False),~ext.larger_counts_valid,ext.larger_call_rate.lt(.95)],['No exact allele record in existing extract','Non-diploid/inconsistent or missing counts','Call rate below 0.95'],default='')
ext['larger_FDR_targeted']=np.nan
ix=ext.unrelated_informative_absence&ext.larger_evaluable
ext.loc[ix,'larger_FDR_targeted']=bh(ext.loc[ix,'larger_P'])
ext['larger_significant_partial_depletion']=ext.larger_observed_HOM.gt(0)&ext.larger_O_E.lt(1)&ext.larger_FDR_targeted.lt(.05)&ext.larger_evaluable
full=ext
# Frozen selection definitions, with exact mapping back to the pLoF master.
ihs=read(BASE/'fig4_v2/qc/regional_iHS_100kb_pLoF_annotation.tsv');clues=read(BASE/'fig4_v2/qc/clues_all_evaluable_logLR_definition.tsv')
assert full.site.is_unique
full=full.merge(ihs[['site','DensityStratifiedTailP','Frac_Extreme','Top1pct','supported_all_three_scales','site_iHS']],on='site',how='left',validate='one_to_one')
full=full.merge(clues[['site','logLR','s','clues_logLR_supported']],on='site',how='left',validate='one_to_one')
full['regional_supported']=full.Top1pct.eq(1);full['CLUES_supported']=flag(full.clues_logLR_supported);full['selection_supported']=full.regional_supported|full.CLUES_supported
full['selection_class']=np.select([full.regional_supported&full.CLUES_supported,full.regional_supported,full.CLUES_supported],['Both','Regional iHS only','CLUES only'],default='Neither/NE')
# Do not infer unmeasured technical quality. Existing genotype-level QC covers only historical14.
qc=read(BASE/'fig4_v3_final/qc/stringent14_complete_candidate_annotation.tsv')
qc_cols=['variant_key','site_QUAL','site_FILTER','mapping_quality_MQ','median_sample_DP','median_sample_GQ','heterozygote_AB_median','heterozygote_AB_q25','heterozygote_AB_q75','VCF_level_QC_status','manual_BAM_IGV_status','repetitive_low_complexity_status']
full=full.merge(qc[qc_cols].rename(columns={'variant_key':'variant'}),on='variant',how='left',validate='one_to_one')
full['technical_QC_scope']=np.where(full.historical14,'Existing VCF QC; no new BAM validation','Not audited at read level; no inferred pass')
for c in qc_cols[1:]:
 if full[c].dtype==object:full[c]=full[c].fillna('Not assessed')
write(full,'01_depletion_candidate_evidence_matrix.tsv')
write(plof,'01_autosomal_pLoF_all_cohorts.tsv.gz')
write(full[full.informative_absence],'01_primary_informative_absence.tsv')
write(full[full.unrelated_informative_absence],'01_unrelated_informative_absence.tsv')
write(full[full.shared_robust],'01_shared_robust.tsv')
write(full[full.depleted&full.observed_HOM.gt(0)].sort_values('P').head(30),'01_partial_depletion_top30.tsv')
write(full[full.unrelated_informative_absence],'02_larger_Holstein_evaluation.tsv')
larges=[]
for label,mask in [('Unrelated candidates',full.unrelated_informative_absence),('Shared robust',full.shared_robust),('Primary candidates',full.informative_absence)]:
 for status,g in full[mask].groupby('larger_status'):larges.append(dict(set=label,status=status,n=len(g),partial_FDR_lt_005=int(g.larger_significant_partial_depletion.sum())))
write(pd.DataFrame(larges),'02_larger_Holstein_summary.tsv')
write(ihs,'04_regional_iHS_frozen.tsv');write(clues,'04_CLUES_frozen.tsv')
write(full,'Supplementary_Table_S5.tsv')
syn=auto[auto.cohort.eq('full574')&auto.variant_class.eq('synonymous')].copy()
write(syn,'01_autosomal_synonymous.tsv.gz')
summary=dict(old_audit=audit,tested=len(full),depleted=int(full.depleted.sum()),depleted_O0=int((full.depleted&full.observed_HOM.eq(0)).sum()),partial_depleted=int((full.depleted&full.observed_HOM.gt(0)).sum()),informative_absence=int(full.informative_absence.sum()),informative_absence_FDR_significant=int(full.significant_informative_absence.sum()),unrelated_informative_absence=int(full.unrelated_informative_absence.sum()),shared=int(full.shared_robust.sum()),primary_only=int((full.informative_absence&~full.unrelated_informative_absence).sum()),unrelated_only=int((~full.informative_absence&full.unrelated_informative_absence).sum()),high_AF_gt_05=int(full.pLoF_AF.gt(.5).sum()),background_genes=full.gene_id.nunique(),primary_genes=full.loc[full.informative_absence,'gene_id'].nunique(),shared_genes=full.loc[full.shared_robust,'gene_id'].nunique(),larger_summary=larges,larger_raw_rows=len(raw),max_old_LRT_abs_difference=max_lrt_difference,regional_evaluable=int(ihs.DensityStratifiedTailP.notna().sum()),regional_supported=int(ihs.Top1pct.eq(1).sum()),regional_multiscale=int(flag(ihs.supported_all_three_scales).sum()),CLUES_evaluable=len(clues),CLUES_supported=int(flag(clues.clues_logLR_supported).sum()),selection_union=int(full.selection_supported.sum()),selection_both=int((full.regional_supported&full.CLUES_supported).sum()),selection_adult_observed=int((full.selection_supported&full.observed_HOM.gt(0)).sum()),selection_adult_observed_not_depleted=int((full.selection_supported&full.adult_HOM_observed_not_depleted).sum()))
summary['larger_callrate_failures']=int((full.unrelated_informative_absence&full.larger_counts_valid&~full.larger_evaluable).sum())
summary['unrelated_old_no_call_filter_autosome']=int((un.observed_HOM.eq(0)&un.expected_HOM.ge(3)&un.depleted).sum())
(ROOT/'summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps(summary,indent=2,ensure_ascii=False))
