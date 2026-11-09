from pathlib import Path
import sys,json,hashlib,re,csv
import pandas as pd
import numpy as np
sys.dont_write_bytecode=True
A=Path('/Users/adrianstier/sbc-donham-kelp-mpa-cascade');D=Path('/Users/adrianstier/sbc-donham-kelp-mpa-data');O=Path(__file__).parent
examined=[]
def read(p,**kw):
 p=Path(p);examined.append(str(p));return pd.read_excel(p,**kw) if p.suffix=='.xlsx' else pd.read_csv(p,**kw)
def save(d,n):d.to_csv(O/n,index=False)
base=D/'data/source_archive/kfm_size_feasibility_2026-10-29'
s=read(D/'data/ALL_sizefreq_2024.csv',low_memory=False)
l=read(base/'legacy/Invertebrate Size Frequency.csv')
p=read(A/'outputs/biomass_support_20261102/private/pool_members.csv',low_memory=False)
a=read(A/'outputs/biomass_resolution_20261102/review_annual_input.csv')
aff=read(A/'outputs/biomass_exclusion_reconciliation_20261102/affected_histories.csv')
cor=read(A/'outputs/biomass_exclusion_reconciliation_20261102/annual_donor_corrections.csv')
sp={11005:'Strongylocentrotus franciscanus',11006:'Strongylocentrotus purpuratus'}
h=[]
for site,ids in [('a-k-01',[1,21]),('a-k-14',[37])]:
 for code,name in sp.items():
  ss=s[s.campus.eq('KFM')&s.site.eq(site)&s.classcode.eq(name)]
  ll=l[l.Site.isin(ids)&l.Species.eq(code)]
  for y in sorted(set(ss.year)&set(ll.Year)):
   x=ss[ss.year.eq(y)].groupby('size')['count'].sum();z=ll[ll.Year.eq(y)].groupby((ll[ll.Year.eq(y)].Size*10).round(6)).size()
   err=x.subtract(z,fill_value=0).abs().max()
   h.append(dict(site=site,species_code=code,year=int(y),exact=err<1e-9,animals=int(x.sum()),max_error=err))
save(pd.DataFrame(h),'legacy_identity_checks.csv');assert all(v['exact'] for v in h)
rows=[]
for r in aff.itertuples():
 q=p[p.series_id.eq(r.series_id)&p.campus.eq('KFM')&p.site.eq(r.delivered_site)]
 yy=sorted(q.year.unique());post=sorted(a.loc[a.series_id.eq(r.series_id)&a.t.gt(0),'year'].unique());fixed=sorted(cor.loc[cor.series_id.eq(r.series_id),'year'].unique());code=11005 if 'M. franciscanus' in r.series_id else 11006
 verified=sorted(v['year'] for v in h if v['site']==r.delivered_site and v['species_code']==code and v['year'] in yy)
 unresolved=sorted(set(yy)-set(verified));required=sorted(set(post)&set(unresolved))
 assert len(q)==r.donor_record_uses and int(q['count'].sum())==r.animals
 def ys(v):return ';'.join(str(int(y)) for y in v)
 rows.append(dict(series_id=r.series_id,donor_role='reference' if r.delivered_site=='a-k-01' else 'mpa',delivered_site=r.delivered_site,donor_use_years=ys(yy),verified_full_histogram_years=ys(verified),post_years=ys(post),diagnostically_corrected_post_years=ys(fixed),unresolved_exact_measurement_post_years=ys(required),n_post=len(post),n_corrected=len(fixed),n_unresolved=len(required),identity_only_2015=int(2015 in unresolved),density_retained=True))
save(pd.DataFrame(rows),'history_years.csv')
# Independently compare historical physical count identities.
c=read(D/'data/MBON/SBCMBON_kelp_forest_integrated_quad_and_swath_20231022.csv',low_memory=False)
ld=read(base/'legacy/Benthic Density Data.csv');cc=[]
c['year']=pd.to_datetime(c.date).dt.year
for name,code in [('Mesocentrotus franciscanus',11005),('Strongylocentrotus purpuratus',11006)]:
 for year in range(2005,2012):
  z=c[c.site_id.eq('a-k-14')&c.taxon_name.eq(name)&c.year.eq(year)]
  val=(pd.to_numeric(z['count'])/z.area).mean();expected=ld[ld.Site.eq(14)&ld.Species.eq(code)&ld.Year.eq(year)].DensityMean.iloc[0]
  cc.append(dict(year=year,species_code=code,mbon_density=val,legacy14_density=expected,error=abs(val-expected)))
assert max(v['error'] for v in cc)<1e-9
save(pd.DataFrame(cc),'count_identity_checks.csv')
# Aggregated 2015 bins establish identity only; never derive individual sizes.
x=read(base/'2015/qryAppendixI_NHSF_2015.xlsx');checks=[]
for site,ids in [('a-k-01',[1,21]),('a-k-14',[37])]:
 for name in sp.values():
  z=x[x.Site_Number.isin(ids)&x.Scientific_Name.eq(name)];actual=s[s.campus.eq('KFM')&s.site.eq(site)&s.year.eq(2015)&s.classcode.eq(name)]
  ex=(z.Percent*z.Cases_N/100).groupby(z.Bin_size).sum();errs=[]
  for label,n in ex.items():
   ns=list(map(int,re.findall(r'\d+',label)));mask=actual['size'].lt(ns[0]) if label.startswith('<') else actual['size'].gt(ns[0]) if label.startswith('>') else actual['size'].between(*ns)
   errs.append(abs(actual.loc[mask,'count'].sum()-n))
  checks.append(dict(site=site,species=name,max_bin_error=max(errs),animals=int(actual['count'].sum())))
assert max(v['max_bin_error'] for v in checks)<1e-6
save(pd.DataFrame(checks),'bins2015_checks.csv')
# Source guard behavior with canonical taxa and source-specific exemption.
g=A/'code/tools/biomass_source_guard.py';examined.append(str(g));ns={};exec(compile(g.read_text(),str(g),'exec'),ns)
guard=[]
for site,taxon,source,expect in [('a-k-01','M. franciscanus','ALL_sizefreq_2024.csv',True),('a-k-14','S. purpuratus','ALL_sizefreq_2024.csv',True),('a-k-21','M. franciscanus','original',True),('a-k-37','S. purpuratus','original',True),('a-k-01','M. pyrifera','original_kelp',False),('a-k-14','M. franciscanus','original_site_specific',False)]:
 issue=ns['source_issue'](site,taxon,source);assert bool(issue)==expect;guard.append(dict(site=site,taxon=taxon,source=source,blocked=bool(issue)))
save(pd.DataFrame(guard),'guard_checks.csv')
# Verify no exclusions propagated to the six density counterparts.
d=read(A/'outputs/dataset_inventory_20261101/comparison_decisions.csv');support=read(A/'outputs/dataset_inventory_20261101/measurement_support_current.csv')
for sid in aff.series_id:
 den=d[d.series_id.eq(sid.replace('|Bio','|Den'))];assert len(den)==1 and not den['exclude'].iloc[0]
 assert len(support[support.series_id.eq(sid)])==1
# Check full actual donor guard and the one uncorrected historical zero year.
blocked={sid for sid,site,taxon in zip(p.series_id,p.site,p.Taxa) if ns['source_issue'](site,taxon)}
assert blocked==set(aff.series_id)
inv=read(A/'outputs/biomass_support_20261102/animal_pool_inventory.csv')
z=inv[inv.series_id.eq('Harris Point SMR|S. purpuratus|KFM|Bio')&inv.year.eq(2005)&inv.status.eq('reference')]
assert len(z)==1 and z.value.iloc[0]==0
zc=c[c.site_id.eq('a-k-01')&c.year.eq(2005)&c.taxon_name.eq('Strongylocentrotus purpuratus')]
assert len(zc)>0 and pd.to_numeric(zc['count']).eq(0).all()
summary_zero={'raw_count_rows':len(zc),'all_zero':True,'explanation':'2005 Harris KFM purple reference needs no size donor'}
(O/'zero_year_check.json').write_text(json.dumps(summary_zero,indent=2))
# Inventory current files, byte hashes only (no copied private measurements).
fi=read(A/'outputs/kfm_holdings_20261102/file_inventory.csv');hold=[]
for r in fi.itertuples():
 f=Path(r.path);hold.append(dict(path=str(f),exists=f.exists(),matches_prior_hash=f.exists() and hashlib.sha256(f.read_bytes()).hexdigest()==r.sha256))
save(pd.DataFrame(hold),'holdings_recheck.csv')
summary=dict(legacy_alias_histograms=len(h),legacy_exact=sum(v['exact'] for v in h),physical_count_checks=len(cc),grouped2015_checks=len(checks),guard_checks=len(guard),histories=len(rows),post_comparison_years=sum(v['n_post'] for v in rows),corrected_post_comparison_years=sum(v['n_corrected'] for v in rows),unresolved_post_comparison_years=sum(v['n_unresolved'] for v in rows),held_paths=len(hold),held_paths_unchanged=sum(v['matches_prior_hash'] for v in hold))
(O/'checks.json').write_text(json.dumps(summary,indent=2,default=int));(O/'script_examined_files.txt').write_text('\n'.join(examined)+'\n');print(json.dumps(summary,indent=2,default=int));print(pd.DataFrame(rows).to_string(index=False))
