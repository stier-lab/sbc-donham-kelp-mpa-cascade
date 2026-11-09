"""Trace seven missing visits and quantify review-only visit retention alternatives."""
from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,zipfile
A=Path(__file__).resolve().parents[2];D=A.parent/'sbc-donham-kelp-mpa-data';O=A/'outputs/kfm_missing_size_recovery_20261103';O.mkdir(exist_ok=True)
P=D/'data/source_archive/kfm_size_feasibility_2026-10-29/legacy';paths=[A/'outputs/kfm_sheephead_biomass_20261103/lost_comparison_years.csv',D/'data/source_archive/fish_lter_reconciliation_20261102/chis_kelp_forest_fish_size.csv',P/'Fish size frequency data.csv',P/'RDFC data.csv',A/'data/MBON/SBCMBON_kelp_forest_integrated_fish_20231022.csv',D/'config/kfm_physical_sites.csv',A/'data/ALL_sizefreq_2024.csv',A/'outputs/kfm_sheephead_biomass_20261103/visit_support.csv',A/'outputs/fish_lter_reconciliation_20261102/fish_count_size_matches.csv']
lost=pd.read_csv(paths[0]);s=pd.read_csv(paths[1]);mapping=pd.read_csv(paths[5]).set_index('source_name').site;s['site_id']=s.siteName.map(mapping);s['year']=s.startDate.str[:4].astype(int);assert s.site_id.notna().all()
s['sheephead']=s.scientificName.str.contains('pulcher',case=False,na=False)
old=pd.read_csv(paths[2],encoding='cp1252');old['site_id']=old.Site.map(lambda x:f'a-k-{x:02}');old['date']=pd.to_datetime(old.Date,format='%d%b%Y').dt.strftime('%Y-%m-%d');old['sheephead']=old.ScientificName.str.contains('pulcher',case=False,na=False)
rdf=pd.read_csv(paths[3],encoding='cp1252');rdf['site_id']=rdf.Site.map(lambda x:f'a-k-{x:02}');rdf['date']=pd.to_datetime(rdf.Date,format='%d%b%Y').dt.strftime('%Y-%m-%d');rdf=rdf[rdf.ScientificName.str.contains('pulcher',case=False,na=False)]
raw=pd.read_csv(paths[4],low_memory=False);raw=raw[raw.data_source.eq('kfm')&raw.taxon_name.str.contains('pulcher',case=False,na=False)].copy();raw['count']=pd.to_numeric(raw['count'],errors='coerce');raw['year']=raw.date.str[:4].astype(int)
evidence=[];source_events=[]
for r in lost.itertuples():
 st=s[s.site_id.eq(r.site_id)&s.year.eq(r.year)];target=st[st.startDate.eq(r.count_date)];ss=st[st.sheephead];oo=old[old.site_id.eq(r.site_id)&old.Year.eq(r.year)];ot=oo[oo.date.eq(r.count_date)];rd=rdf[rdf.site_id.eq(r.site_id)&rdf.date.eq(r.count_date)];rr=raw[raw.site_id.eq(r.site_id)&raw.date.eq(r.count_date)];cv=rr[rr.sample_method.eq('visualfish')]
 nearest=int((pd.to_datetime(ss.startDate)-pd.Timestamp(r.count_date)).dt.days.abs().min()) if len(ss) else None
 category='Size survey recorded other fish but no sheephead sizes' if len(target) else 'No size survey records for any fish on count date'
 evidence.append(dict(MPA=r.MPA,year=r.year,site_id=r.site_id,site_name=next((x for x,v in mapping.items() if v==r.site_id),r.site_id),count_date=r.count_date,count_rows=len(cv),observed_count=cv['count'].sum(),count_density=r.count_density,scott_all_fish_rows_on_date=len(target),scott_sheephead_rows_on_date=int(target.sheephead.sum()),scott_sheephead_dates=';'.join(sorted(ss.startDate.unique())),nearest_sheephead_gap_days=nearest,legacy_all_fish_rows_on_date=len(ot) if r.year<=2011 else np.nan,legacy_sheephead_rows_on_date=int(ot.sheephead.sum()) if r.year<=2011 else np.nan,legacy_rdfc_sheephead_positive_records=int(rd.Count.gt(0).sum()) if r.year<=2011 else np.nan,classification=category,recovered_measurements=0))
 for date,g in st.groupby('startDate'):source_events.append(dict(site_id=r.site_id,year=r.year,date=date,all_fish_rows=len(g),sheephead_rows=int(g.sheephead.sum()),sheephead_fish=g.loc[g.sheephead,'speciesCount'].sum()))
pd.DataFrame(evidence).to_csv(O/'missing_visit_evidence.csv',index=False);pd.DataFrame(source_events).drop_duplicates().to_csv(O/'available_size_events.csv',index=False)
# The combined historical size pool contains KFM benthos/kelp, not KFM fish.
combined=pd.read_csv(paths[6],low_memory=False);k=combined[combined.campus.eq('KFM')];assert not k.classcode.str.contains('pulcher|sheephead',case=False,na=False).any()
# Reconcile physical visit subsets, with the original site weights retained.
v=pd.read_csv(paths[7]);roster=pd.read_csv(paths[8]);selected=[];coverage=[];lossrows=[]
for (mpa,year),rr in roster[roster.year.ge(2007)].groupby(['MPA','year']):
 ids=set(rr.site_id);vv=v[v.site_id.isin(ids)&v.year.eq(year)].copy();vv=vv.merge(rr[['site_id','status']],on='site_id',validate='many_to_one');vv['supported']=vv['mode'].ne('unsupported');complete=bool(vv.supported.all());subset=all(vv.groupby('site_id').supported.any())
 original_n=len(vv); kept=vv[vv.supported];counts=raw[raw.sample_method.eq('visualfish')&raw.site_id.isin(ids)&raw.year.eq(year)]
 # Match source records to visits and preserve counts, not just means.
 c=counts.merge(vv[['site_id','count_date','supported']],left_on=['site_id','date'],right_on=['site_id','count_date'],validate='many_to_one')
 coverage.append(dict(MPA=mpa,year=year,all_visits_supported=complete,matched_visit_subset_possible=subset,n_original_sites=len(ids),n_original_visits=original_n,n_supported_visits=len(kept),n_original_count_rows=len(c),n_retained_count_rows=int(c.supported.sum()),original_fish_count=c['count'].sum(),retained_fish_count=c.loc[c.supported,'count'].sum()))
 if not subset:continue
 for role in ['mpa','reference']:
  role_sites=rr[rr.status.eq(role)].site_id;g=kept[kept.status.eq(role)]
  # Renormalize within site after dropping visits: explicitly a different sample.
  site_density=[];site_biomass=[];site_days=[]
  for site,vs in g.groupby('site_id'):
   w=vs.annual_weight.to_numpy();w=w/w.sum();site_density.append(np.average(vs.count_density,weights=w));site_biomass.append(np.average(vs.biomass_g_m2,weights=w));site_days.append(np.average(pd.to_datetime(vs.count_date).dt.dayofyear,weights=w))
   for j,z in enumerate(vs.itertuples()):lossrows.append(dict(MPA=mpa,year=year,status=role,site_id=site,count_date=z.count_date,new_weight=w[j],original_weight=z.annual_weight))
  dates=pd.to_datetime(g.count_date);selected.append(dict(MPA=mpa,year=year,status=role,n_sites=len(role_sites),n_visits=len(g),first_date=dates.min().date().isoformat(),last_date=dates.max().date().isoformat(),weighted_mean_day_of_year=float(np.mean(site_days)),subset_density=np.mean(site_density),subset_biomass_g_m2=np.mean(site_biomass),original_density=rr[rr.status.eq(role)].count_density.mean()))
coverage=pd.DataFrame(coverage);coverage.to_csv(O/'all_year_coverage_options.csv',index=False);pd.DataFrame(selected).to_csv(O/'matched_visit_subset_responses.csv',index=False);pd.DataFrame(lossrows).to_csv(O/'matched_visit_subset_weights.csv',index=False)
focus=coverage.merge(lost[['MPA','year']].drop_duplicates(),on=['MPA','year']);focus.to_csv(O/'five_year_options.csv',index=False)
# Hash available original Drive copies, without assuming filenames prove independence.
drive=Path('/Users/adrianstier/Library/CloudStorage/GoogleDrive-<account>/My Drive/Stier Lab/People/Emily Donham/Projects/Kelp MPA');copies=[]
for prefix in [drive,drive/'Github/MPA_Kelp']:
 for rel in ['data/ALL_sizefreq_2024.csv','data/MBON/SBCMBON_kelp_forest_integrated_fish_20231022.csv']:
  p=prefix/rel
  if p.exists():paths.append(p);copies.append(dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),matches_current=hashlib.sha256(p.read_bytes()).hexdigest()==hashlib.sha256((A/rel).read_bytes()).hexdigest()))
archive=P.parent/'KushnerEtAl2013-KFM-Data.zip';paths.append(archive)
z=zipfile.ZipFile(archive)
assert all(z.read(name)==(P/name).read_bytes() for name in ['Fish size frequency data.csv','RDFC data.csv'])
(O/'input_manifest.json').write_text(json.dumps([dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in paths],indent=2))
pd.DataFrame(copies).to_csv(O/'drive_copies.csv',index=False)
result=dict(missing_visits=len(lost),affected_comparison_years=len(focus),raw_measurements_recovered=0,existing_all_visit_years=int(coverage.all_visits_supported.sum()),matched_subset_years=int(coverage.matched_visit_subset_possible.sum()),kfm_combined_size_rows=len(k),kfm_combined_sheephead_rows=0,primary_changed=False)
(O/'summary.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2));print(focus.to_string(index=False))
# Compact reader-facing table; detailed records remain in linked CSVs.
show=focus[['MPA','year']].copy();show['Visits retained']=focus.n_supported_visits.astype(str)+' / '+focus.n_original_visits.astype(str);show['Count records retained']=focus.n_retained_count_rows.astype(str)+' / '+focus.n_original_count_rows.astype(str);show['All sites retained']=np.where(focus.matched_visit_subset_possible,'Yes, by omitting unmatched visits','No; a site lacks sizes')
(O/'index.html').write_text('<!doctype html><meta charset="utf-8"><title>KFM missing-size recovery audit</title><style>body{font:17px system-ui;max-width:1050px;margin:30px auto;padding:20px;line-height:1.5}table{border-collapse:collapse;font-size:15px}td,th{padding:10px;border-bottom:1px solid #bbb;text-align:left}</style><h1>KFM missing-size recovery audit</h1><p>No additional sheephead lengths were recovered from the checked source files. The current 52 comparison-years and paper remain unchanged.</p>'+show.to_html(index=False)+'<p>Anacapa 2007 and Scorpion 2007 and 2022 could support a separate matched-visit sensitivity: 55 comparison-years total. This drops count visits and changes seasonal coverage; it does not recover the original annual observations. Harris Point 2010 and Scorpion 2023 would require dropping a site or borrowing sizes.</p><p>At Scorpion 2022, retained protected visits are June–August and references are September. Keep the current 52-year review; test the alternative separately if useful. No fitted results were changed.</p><p>Scott’s export agrees with the older individual-size archive for the historical gaps. Emily’s project and GitHub-folder input copies are byte-identical to ours. Separate historical roving observations show that missing sheephead lengths cannot simply be read as fish absence.</p><p><a href="README.md">Full findings and recommendation</a> · <a href="missing_visit_evidence.csv">Visit-by-visit source evidence</a> · <a href="five_year_options.csv">Full coverage table</a> · <a href="matched_visit_subset_responses.csv">Candidate observations and timing</a> · <a href="verification.json">Verification</a></p>')
