"""Build linked, review-only dataset/site/comparison inventory from raw scan and evidence."""
from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,itertools,html
O=Path('outputs/dataset_inventory_20261101');A=Path.cwd();D=Path('/Users/adrianstier/Donham-Stier-CA-MPA-Data-2026/data');M=Path('/Users/adrianstier/sbc-donham-kelp-mpa-manuscript/correspondence/coauthor_reviews/2026-10-31_candidate_sites')
inputs={}
def rd(p):
 p=Path(p);inputs[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest();return pd.read_csv(p)
def save(d,n):d.to_csv(O/(n+'.csv'),index=False)
def uid(*v):return hashlib.sha256(json.dumps(v,ensure_ascii=False).encode()).hexdigest()[:16]
def years(v):return sorted({int(float(t)) for t in str(v).replace(',',';').split(';') if t and t!='nan'})
def cov(y,start=np.nan):
 y=sorted(set(y));gap=sorted(set(range(min(y),max(y)+1))-set(y)) if y else []
 return dict(first_year=min(y) if y else np.nan,last_year=max(y) if y else np.nan,observed_years=';'.join(map(str,y)),n_years=len(y),elapsed_years=max(y)-min(y) if y else np.nan,missing_years=';'.join(map(str,gap)),n_missing=len(gap) if y else np.nan,longest_internal_gap=max([b-a-1 for a,b in zip(y,y[1:])],default=0) if y else np.nan,n_strict_before=sum(t<start for t in y) if pd.notna(start) else np.nan,n_designation_year=sum(t==start for t in y) if pd.notna(start) else np.nan,n_strict_after=sum(t>start for t in y) if pd.notna(start) else np.nan)
scan=rd(O/'raw_station_response_scan.csv');datasets=rd(O/'raw_dataset_scan.csv');assign=rd(M/'candidate_station_decisions.csv').rename(columns={'status':'historical_assignment_status','reason':'historical_assignment_reason'});meta=rd(M/'candidate_mpa_metadata.csv');oldmeta=rd(M/'candidate_station_coverage_metadata.csv')
assert not assign.duplicated(['program','MPA','site_id','role']).any()
meta=meta.rename(columns={'current_inclusion':'historical_inclusion_20261031','primary_history_count':'historical_history_count_20261031','selection_reason':'historical_selection_reason'})
take=rd(D/'mpa_take_regulations.csv').rename(columns={'CA_MPA_Name_Short':'MPA'})
meta=meta.merge(take,on='MPA',how='left',validate='one_to_one')
res=rd(O/'protection_metadata_resolutions.csv')
for _,r in res.iterrows():
 for col in ['protection_class','take_restrictions','restriction_source','restriction_check_status']:meta.loc[meta.MPA==r.MPA,col]=r[col]
save(meta,'mpa_metadata')
sel=rd('outputs/researcher_scope_20261101/selection_decisions.csv');annual=rd('outputs/upstream_sites_20261101/annual_input.csv');current=rd('data/sumstats_final.csv');registry=rd('data/harmonized/reference_primary_keys.csv');current['series_id']=current[['MPA','Taxa','Source','Resp']].astype(str).agg('|'.join,axis=1)
size=rd('tables/table_s_primary_size_inventory.csv');season=rd('tables/table_s_seasonal_paired_years.csv');foot=rd('outputs/seasonal_zone_audit_20261101/other_program_summary.csv');recon=rd('outputs/seasonal_zone_audit_20261101/reconstruction_checks.csv');fixed=rd('outputs/seasonal_zone_audit_20261101/eligibility.csv');coverage=rd('outputs/full_followup_recovery_20261101/coverage.csv')
# satellite delivery: annual polygons, not dive stations; keep known empty groups separate.
sat=rd(D/'Landsat/MPA_Runs_new.csv');ys=[c for c in sat if c.isdigit()];satrows=[]
empty=rd(D/'landsat_empty_reference_corrections.csv')
for _,r in sat.iterrows():
 y=[int(c) for c in ys if pd.notna(r[c])];satrows.append(dict(dataset_id='landsat_annual',program='Landsat',MPA=r.MPA,site_id=f'{r.MPA}|{r.status}|{r.rep}',taxon='M. pyrifera',method='Annual maximum canopy biomass',role=r.status,observed_years=';'.join(map(str,y)),months='Annual maximum of quarterly observations',raw_record_rows=len(ys),zero_count_rows=sum(r[c]==0 for c in ys),missing_count_rows=sum(pd.isna(r[c]) for c in ys),negative_count_rows=sum(r[c]<0 for c in ys),coverage_basis='Delivered annual polygon values; empty-mask registry takes precedence over numeric placeholders'))
satrows=pd.DataFrame(satrows);save(satrows,'satellite_polygon_coverage');save(empty,'satellite_empty_groups')
ss=[]
for name,g in sat.groupby('MPA'):
 good=g.merge(empty[['MPA','status','rep']].assign(empty_slot=True),on=['MPA','status','rep'],how='left');good=good[good.empty_slot!=True]
 py={int(y) for y in ys if good[good.status=='mpa'][y].notna().any()};ry={int(y) for y in ys if good[good.status=='reference'][y].notna().any()}
 ss.append(dict(MPA=name,delivered_rows=len(g),nonempty_reference_rows=int((good.status=='reference').sum()),paired_years=len(py&ry),candidate_reconstruction_present=((sel.MPA==name)&(sel.Source=='Landsat')).any(),status='Raw polygon coverage; researcher scope and model eligibility separate' if py&ry else 'No paired nonempty reference coverage'))
save(pd.DataFrame(ss),'satellite_candidate_coverage')
satentry=dict(dataset_id='landsat_annual',source_file='Landsat/MPA_Runs_new.csv',source_sha256=inputs[str(D/'Landsat/MPA_Runs_new.csv')],source_rows=len(sat),scanned_rows=len(sat),programs='Landsat',first_date='1984 (annual)',last_date='2021 (annual)',raw_sites=len(sat),scope='Annual MPA/reference polygon arrays; not independent pixels',release_label='2022 extraction / 2023 delivery; masks corrected October 2026',depth_field_present=False,size_field_present=False)
datasets=pd.concat([datasets,pd.DataFrame([satentry])],ignore_index=True)
notes={
'pisco_fish':('Non-canopy fish transects; campus/level retained in raw source','Fish counts and total lengths; density and inferred biomass','count / total length cm; area applied by processing','Explicit blank counts and missing lengths are not zeros; absent code requires survey-effort evidence','PISCO release methods; code/R/04_pisco_processing.R'),
'pisco_swath':('Benthic swaths across recorded depth strata','Kelp and invertebrate counts; kelp frond information','count; raw size interpretation is taxon-specific','Missing count/fronds remain unknown; absent taxon can be zero only within a verified sampled event','PISCO release methods; code/R/04_pisco_processing.R'),
'mbon_fish':('Integrated visual and other fish survey methods, kept separate','Focal sheephead count rows; no lengths in delivery','count and recorded area','Observed zeros retained; methods/proj_taxon_id distinguish sampling types','MBON 20231022; code/R/05_kfm_processing.R'),
'mbon_benthic':('Integrated quadrats and swaths, kept separate','Focal kelp and invertebrate count rows','count and recorded area','Negative/blank values not inferred as absence; kelp stage codes differ','MBON 20231022; code/R/05_kfm_processing.R'),
'lter_fish':('Annual LTER fish transects','Sheephead counts and size','count; size cm; area m2','-99999 is missing, not zero','code/R/06_lter_processing.R'),
'lter_kelp':('Annual LTER kelp quadrats','Fronds and holdfast diameter','fronds; diameter; area m2','-99999 is missing; same plant identity not inferred from repeated values','SBC 127.5 Table 2 Model 1 for conversion'),
'lter_lobster':('Annual LTER lobster transects','Counts, carapace length and antennae-only observations','count; SIZE_MM mm; area m2','Five missing lengths handled with flagged same-site/date mean mass in canonical pipeline','docs/LOBSTER_MISSING_LENGTHS_2026-10-30.md'),
'invertebrate_sizes':('Separate size-frequency surveys; source program/site must be crosswalked','Lobster and urchin lengths/diameters and frequencies','Raw size units require taxon-specific mapping; no global conversion','Missing day is not missing year; count frequencies are not independent dives','docs/PRIMARY_SIZE_POOL_INVENTORY_2026-10-30.md'),
'kfm_kelp_sizes':('KFM plant stipe/diameter records','Stipes per plant and diameter; used with independent plant counts','stipes and diameter cm','Do not equate plant quadrat counts with stipe counts','outputs/anacapa_comparability_20261101/README.md'),
'landsat_annual':('Fixed saved polygon masks; annual max of quarterly sums','Surface canopy wet biomass, distinct from diver standing biomass','kg wet biomass per polygon','Empty reference masks are unavailable, not biological zeros; >25% missing pixels invalidates quarter; >1 missing quarter invalidates year','docs/LANDSAT_V2_VERIFICATION_2026-10-21.md')}
for field,j in [('sampling_design',0),('measured_response',1),('units',2),('missing_zero_convention',3),('method_evidence',4)]:datasets[field]=datasets.dataset_id.map(lambda k:notes[k][j])
datasets['data_use_status']=datasets.programs.map(lambda p:'KFM permission/redistribution awaiting confirmation; no raw publication' if 'KFM' in str(p) else 'Source-specific terms must be retained; no new redistribution approval inferred')
datasets['scan_encoding']='Latin-1 lossless decoding for raw CSV text; scientific numeric fields unchanged'
save(datasets,'dataset_inventory')
# Link all mapped assignments to raw site/response/method records; no source station is silently deleted.
rows=[]
for _,a in assign.iterrows():
 q=scan[(scan.program==a.program)&(scan.site_id.astype(str)==str(a.site_id))]
 m=meta[meta.MPA==a.MPA];start=m.iloc[0].recorded_analysis_year if len(m) else np.nan
 cm=oldmeta[(oldmeta.program==a.program)&(oldmeta.MPA==a.MPA)&(oldmeta.site_id.astype(str)==str(a.site_id))]
 loc={}
 for field in ['station_latitude','station_longitude','station_coordinate_status','habitat_note','exposure_note']:
  vals=cm[field].dropna().unique();loc[field]=vals[0] if len(vals)==1 else 'Unresolved / multiple' if len(vals)>1 else ''
 for _,s in (q if len(q) else pd.DataFrame([dict(dataset_id='unmatched',taxon='Unknown',method='Unknown',observed_years='')])).iterrows():
  base={**a.to_dict(),**s.to_dict(),**loc,'designation_year':start};base['site_record_id']=uid(a.program,a.MPA,a.site_id,a.role,s.dataset_id,s.taxon,s.method)
  base.update(cov(years(s.observed_years),start));base['response_support']='Survey opportunity only' if s.taxon in ['Other / survey opportunity','Unknown'] else 'Measurement records; biomass availability not inferred from count coverage'
  ss=sel[(sel.MPA==a.MPA)&(sel.Source==a.program)]
  base['review_scope']='All reconstructed responses excluded by researcher/ecological rule' if len(ss) and ss.exclude.all() else 'Some candidate responses retained; station eligibility remains response-specific' if len(ss) else 'No reconstructed candidate history; not proof of no data'
  rows.append(base)
site=pd.DataFrame(rows);save(site,'site_response_coverage')
mapped=set(zip(assign.program,assign.site_id.astype(str)));unmapped=scan[~scan.apply(lambda r:(r.program,str(r.site_id)) in mapped,axis=1)].copy();unmapped['interpretation']='Outside supplied assignment map or uses a different source identifier; no automatic exclusion/addition';save(unmapped,'unmapped_source_records')
# Comparison-level inventory includes every reconstructed candidate, including researcher exclusions.
cr=[]
for _,r in sel.iterrows():
 h=annual[annual.series_id==r.series_id];sy=sorted(set(h.year));start=h.MPA_Start.iloc[0] if len(h) else np.nan
 rr={**r.to_dict(),**cov(sy,start),'comparison_id':uid(r.series_id),'designation_year':start,'in_current_manuscript':r.series_id in set(current.series_id),'review_scope_status':'Excluded by researcher/ecological rule' if r.exclude else 'Candidate retained; model eligibility separate'}
 mm=meta[meta.MPA==r.MPA]
 for col in ['protection_or_network_date','date_meaning','protection_class','take_restrictions','restriction_check_status','earlier_protection_note']:
  rr[col]=mm.iloc[0][col] if len(mm) else 'Unresolved'
 st=site[(site.MPA==r.MPA)&(site.program==r.Source)]
 rr['mapped_protected_sites']=';'.join(sorted(st[st.role.isin(['mpa','protected'])].site_id.astype(str).unique()));rr['mapped_reference_sites']=';'.join(sorted(st[st.role=='reference'].site_id.astype(str).unique()))
 rr['station_list_scope']='All supplied mapped candidates, including excluded stations; not fitted membership'
 post=[y for y in sy if y>start];rr['proposed_recovery_duration_pass']=len(post)>=5 and max(post,default=0)-min(post,default=0)>=5;rr['three_strict_prior_years']=sum(y<start for y in sy)>=3;rr['four_strict_prior_years']=sum(y<start for y in sy)>=4
 rr['eligibility_interpretation']='Coverage screen only; station continuity, model assumptions and matching still required; designation year separate'
 zz=season[(season.MPA==r.MPA)&(season.Source==r.Source)&(season.Taxa==r.Taxa)&(season.Resp==r.Resp)]
 rr['previous_season_audit_years']=len(zz);rr['previous_no_shared_month_years']=int((zz.n_common_months==0).sum()) if len(zz) else np.nan;rr['season_scope']='Previous canonical cohort only; not certification of expanded candidate years'
 sz=size[size.series_id==r.series_id];rr['size_inventory_side_years']=len(sz);rr['size_positive_sparse_lt10']=int((sz.positive&sz.sparse_lt10).sum()) if len(sz) else np.nan;rr['size_singletons']=int(sz.singleton.sum()) if len(sz) else np.nan;rr['size_shared_side_years']=int(sz.shared_with_other_history.sum()) if len(sz) else np.nan;rr['size_scope']='Audited canonical lobster/purple-urchin pools only' if len(sz) else 'Not covered by this size-pool inventory; no adequacy inferred'
 ff=foot[foot.series_id==r.series_id];rr['footprint_status']='; '.join(ff.scope.unique()) if len(ff) else 'See PISCO footprint audit / not applicable to satellite';rr['physical_depth_status']='Measured depths available; station/year matching requires checking' if r.Source=='PISCO' else 'No compatible survey-depth field in inspected count files' if r.Source in ['KFM','LTER'] else 'Surface canopy pixels; diver depth matching not applicable'
 cv=coverage[(coverage.series_id==r.series_id)&(coverage.scenario=='Full')];rr['previous_full_fit_status']='; '.join(cv.metric.astype(str)+':'+cv.eligible.astype(str))
 rr['reference_evidence']='Inherited mapping; no complete station-specific habitat justification recovered'
 if r.Source=='PISCO' and r.MPA=='Naples SMCA':rr['reference_evidence']='Jenn: diver expertise and community nMDS; no close match; ordination not recovered'
 if r.Source=='PISCO' and r.MPA=='Harris Point SMR':rr['reference_evidence']='Scott Hamilton: Cuyler nearby/similar habitat; Crook/Tyler differ; short reference records separately documented'
 if r.Source=='LTER':rr['reference_evidence']='Emily: consultation with Bob/Li; individual station rationales not recovered'
 if r.Source=='Landsat':rr['reference_evidence']='Tom saved masks verified; approximate 1–3 km reference band; current equal-available reference averaging; split history partly reconstructed'
 cr.append(rr)
comp=pd.DataFrame(cr);save(comp,'comparison_eligibility')
# Raw candidate pair coverage on same response/method; opportunity only, not final inclusion.
ps=[]
for key,g in site[site.dataset_id!='unmatched'].groupby(['program','MPA','dataset_id','taxon','method']):
 if key[3]=='Other / survey opportunity':continue
 for (_,p),(_,r) in itertools.product(g[g.role.isin(['mpa','protected'])].iterrows(),g[g.role=='reference'].iterrows()):
  common=set(years(p.observed_years))&set(years(r.observed_years));ps.append(dict(pair_id=uid(*key,p.site_id,r.site_id),program=key[0],MPA=key[1],dataset_id=key[2],taxon=key[3],method=key[4],protected_record_id=p.site_record_id,reference_record_id=r.site_record_id,protected_site=p.site_id,reference_site=r.site_id,**cov(common,p.designation_year),interpretation='Same-year raw recorded response/method; not same-visit/depth matched, no biomass certification or researcher exclusions applied'))
pairs=pd.DataFrame(ps);save(pairs,'candidate_pair_coverage')
# Compact view, one response history per row, preserves excluded candidates.
summary=comp[['MPA','Source','Taxa','Resp','in_current_manuscript','review_scope_status','protection_or_network_date','protection_class','first_year','last_year','n_years','n_missing','longest_internal_gap','n_strict_before','n_designation_year','n_strict_after','three_strict_prior_years','four_strict_prior_years','proposed_recovery_duration_pass','physical_depth_status','reference_evidence','size_positive_sparse_lt10','reason']].copy();save(summary,'supplement_table_preview')
checks={'dataset_rows':len(datasets),'source_rows_scanned':int(datasets.source_rows.sum()),'mapped_assignment_keys':len(assign),'site_response_method_rows':len(site),'raw_candidate_pairs':len(pairs),'candidate_histories':len(comp),'current_manuscript_histories':len(current),'current_histories_in_candidate_inventory':int(comp.in_current_manuscript.sum()),'excluded_candidate_histories':int(comp.exclude.sum()),'candidate_histories_retained':int((~comp.exclude).sum()),'mpa_metadata_rows':len(meta),'unmapped_source_groups':len(unmapped),'models_fitted':0}
for df,key in [(datasets,'dataset_id'),(site,'site_record_id'),(comp,'comparison_id'),(pairs,'pair_id')]:assert not df[key].duplicated().any()
assert set(pairs.protected_record_id)|set(pairs.reference_record_id)<=set(site.site_record_id)
assert set(site.dataset_id)-{'unmatched'}<=set(datasets.dataset_id)
assert set(current.series_id)<=set(comp.series_id)
for df in [site,comp,pairs]:
 v=df[df.n_years>0];assert ((v.n_years+v.n_missing)==v.last_year-v.first_year+1).all()
assert len(comp)==len(sel)==206
save(pd.DataFrame([checks]),'audit_counts');(O/'build_checks.json').write_text(json.dumps(checks,indent=2)+'\n');(O/'table_input_hashes.json').write_text(json.dumps(inputs,indent=2)+'\n');print(checks)
