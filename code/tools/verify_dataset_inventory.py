"""Independent reconciliation of coverage tables, keys and unchanged source artifacts."""
from pathlib import Path
import pandas as pd,hashlib,json
O=Path('outputs/dataset_inventory_20261101')
def read(n):return pd.read_csv(O/(n+'.csv'),dtype={'site_id':str},low_memory=False)
def yrs(s):return {int(float(x)) for x in str(s).replace(',',';').split(';') if x and x!='nan'}
site=read('site_response_coverage');pairs=read('candidate_pair_coverage');comp=read('comparison_eligibility');raw=read('raw_station_response_scan');year=read('raw_station_year_audit');ds=read('dataset_inventory')
for d,k in [(site,'site_record_id'),(pairs,'pair_id'),(comp,'comparison_id'),(ds,'dataset_id')]:assert not d[k].duplicated().any(),k
assert set(site.dataset_id)-{'unmatched'}<=set(ds.dataset_id)
assert set(pairs.protected_record_id)|set(pairs.reference_record_id)<=set(site.site_record_id)
for d in [site,pairs,comp]:
 for r in d.itertuples():
  ys=yrs(r.observed_years);assert len(ys)==r.n_years
  if ys:
   gaps=set(range(min(ys),max(ys)+1))-ys;assert gaps==yrs(r.missing_years);assert len(gaps)==r.n_missing;assert max(ys)-min(ys)==r.elapsed_years;assert max([b-a-1 for a,b in zip(sorted(ys),sorted(ys)[1:])],default=0)==r.longest_internal_gap
lookup=site.set_index('site_record_id')
for r in pairs.itertuples():
 a=lookup.loc[r.protected_record_id];b=lookup.loc[r.reference_record_id]
 assert yrs(r.observed_years)==yrs(a.observed_years)&yrs(b.observed_years)
 assert a.role in ['mpa','protected'] and b.role=='reference'
 assert a.MPA==b.MPA==r.MPA and a.dataset_id==b.dataset_id==r.dataset_id
 assert r.years_nearest_dates_within60<=r.years_nearest_dates_within90<=r.years_with_exact_dates_both<=r.n_years
sel=pd.read_csv('outputs/researcher_scope_20261101/selection_decisions.csv');assert set(sel.series_id)==set(comp.series_id)
assert comp.set_index('series_id').exclude.to_dict()==sel.set_index('series_id').exclude.to_dict()
canonical=pd.read_csv('data/sumstats_final.csv');keys=set(canonical[['MPA','Taxa','Source','Resp']].agg('|'.join,axis=1));assert keys==set(comp[comp.in_current_manuscript].series_id)
# Annual metadata partition agrees with station/source records, including missing full dates.
rr=raw.set_index(['dataset_id','program','site_id','taxon','method'])
for k,g in year.groupby(['dataset_id','program','site_id','taxon','method']):
 r=rr.loc[k];assert set(g.year)==yrs(r.observed_years),k;assert g.rows.sum()<=r.raw_record_rows
# Hash-pin every data input and independently preserve the frozen manuscript/results baseline.
checked={}
for name in ['raw_source_hashes.json','table_input_hashes.json','additional_input_hashes.json','site_resolution_input_hashes.json','measurement_source_hashes.json']:
 for f,h in json.loads((O/name).read_text()).items():
  assert hashlib.sha256(Path(f).read_bytes()).hexdigest()==h,f;checked[f]=h
protected=json.loads(Path('outputs/researcher_scope_20261101/protected_hashes.json').read_text())
for f,h in protected.items():assert hashlib.sha256(Path(f).read_bytes()).hexdigest()==h,f
# Verify source date/gap arithmetic independently for every dedicated LTER lobster station.
p=Path('/Users/adrianstier/Donham-Stier-CA-MPA-Data-2026/data/LTER/Lobster_Abundance_All_Years_20230831.csv');x=pd.read_csv(p);dates=pd.to_datetime(x.DATE)
for name,g in x.groupby('SITE'):
 r=raw[(raw.dataset_id=='lter_lobster')&(raw.site_id==name)].iloc[0];assert yrs(r.observed_years)==set(pd.to_datetime(g.DATE).dt.year);assert r.raw_record_rows==len(g);assert r.negative_count_rows==(g.COUNT<0).sum()
sat=read('satellite_full_delivery');assert len(sat)==72 and sat.site_number.nunique()==24 and sat.empty_mask.sum()==14
exc=read('station_link_exceptions');assert len(exc)==64 and exc.record_link_status.str.startswith('Cross-program').sum()==43
groups=read('site_group_resolution');assert len(groups)==21 and (groups.n_children>0).sum()==19
assert site.inventory_record_type.eq('Grouped-site label; not a sampled station').sum()==19
assert site.inventory_record_type.eq('Unassigned map label; no records recovered').sum()==2
assert site.loc[site.inventory_record_type.str.startswith('Grouped'),'n_years'].eq(0).all()
for r in comp.itertuples():
 ss=site[(site.MPA==r.MPA)&(site.program==r.Source)&(site.inventory_record_type=='Observed source/response record')]
 for col,role in [('mapped_protected_sites','mpa'),('mapped_reference_sites','reference')]:
  value=getattr(r,col);actual=set() if pd.isna(value) else set(value.split(';'));expected=set(ss[ss.role==role].site_id)
  assert actual==expected,(r.series_id,col)
res=dict(unique_dataset_ids=len(ds),site_rows_checked=len(site),pair_intersections_checked=len(pairs),comparison_keys_checked=len(comp),current_paper_keys_preserved=len(keys),raw_annual_metadata_rows=len(year),hashed_inputs_unchanged=len(checked),protected_files_unchanged=len(protected),models_fitted=0,checks='Keys, joins, exact year intersections, gaps/spans, duration nesting, all candidate scope decisions, current paper membership, raw annual partitions and source/protected hashes passed')
(O/'verification.json').write_text(json.dumps(res,indent=2)+'\n');print(res)

# Measurement-audit completeness and independent source checks.
measure=read('measurement_design_scan');audit=read('dataset_audit')
assert set(measure.dataset_id)==set(audit.dataset_id)==set(ds.dataset_id)
assert len(audit)==10 and audit.audit_id.nunique()==10
assert (measure.exact_repeat_rows<=measure.screened_rows).all()
assert site.measurement_audit_id.notna().sum()==site.dataset_id.isin(ds.dataset_id).sum()
D=Path('/Users/adrianstier/Donham-Stier-CA-MPA-Data-2026/data')
f=pd.read_csv(D/'LTER/Annual_fish_comb_20240307.csv');f=f[f.SP_CODE.eq('SPUL')]
assert set(f.AREA)=={80} and ((f.COUNT>0)&(f.SIZE.isna()|f.SIZE.le(0))).sum()==0
l=pd.read_csv(D/'LTER/Lobster_Abundance_All_Years_20230831.csv')
assert set(l.AREA)=={300} and ((l.COUNT>0)&(l.SIZE_MM.isna()|l.SIZE_MM.le(0))).sum()==6
assert measure.set_index('dataset_id').loc['lter_lobster','positive_count_missing_or_nonpositive_size']==6
a=pd.read_csv(O/'pisco_fish_level_coverage.csv');assert len(a)==1985 and a.unequal_level_coverage.sum()==85
assert len(pd.read_csv(O/'pisco_fish_revisited_units.csv'))==181
assert 'F7' in set(pd.read_csv(O/'freeze_readiness.csv').item_id)
res.update(measurement_datasets_checked=10,measurement_review_rows=10,additional_analysis_checks=len(read('additional_analysis_checks')),measurement_checks='Complete dataset linkage; direct LTER area/length checks; PISCO survey-support counts; measurement source hashes passed')
(O/'verification.json').write_text(json.dumps(res,indent=2)+'\n');print('Measurement/design checks passed')
