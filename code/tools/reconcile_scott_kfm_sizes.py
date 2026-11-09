"""Reconcile original Scott size frequencies and refit isolated review histories."""
from pathlib import Path
import ast,hashlib,json,html,subprocess,sys
import pandas as pd
import numpy as np
from scipy.stats import t as student
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=Path('outputs/kfm_scott_reconciliation_20261102');O.mkdir(exist_ok=True)
S=Path('data/external/scott_gabara_2026-11-02/KFM Nat Hab Red Purp 2012-2022.xlsx')
assert hashlib.sha256(S.read_bytes()).hexdigest()=='c1bf620b67543bbb0f92b74c6247d3949e03fbc1e7a9d998f1e83b1836fa088e'
# Read original cells; trailing formatted blank rows are not observations.
from openpyxl import load_workbook
wb=load_workbook(S,read_only=True,data_only=True);it=wb['Sheet1'].iter_rows(values_only=True);header=next(it);rows=[]
for i,row in enumerate(it,2):
 if any(v is not None for v in row):rows.append(dict(zip(header,row),excel_row=i))
wb.close();src=pd.DataFrame(rows);assert len(src)==55142
src['year']=pd.to_datetime(src.startDate).dt.year;src['mm']=pd.to_numeric(src.invertebrateSizeInMillimeters);src['count']=pd.to_numeric(src['count']);assert src.mm.gt(0).all() and src['count'].gt(0).all()
p=pd.read_csv('outputs/biomass_support_20261102/private/pool_members.csv',low_memory=False)
h=pd.read_csv('outputs/biomass_next_steps_20261102/kfm_year_holds.csv');need=h[h.biomass_status.eq('original physical-site sizes required')];assert len(need)==54
inv=pd.read_csv('outputs/biomass_support_20261102/animal_pool_inventory.csv').set_index(['series_id','year','status'])
off=pd.read_csv('outputs/all_program_station_refits_20261102/offsets.csv').set_index('series_id').offset
a=pd.read_csv('outputs/biomass_resolution_20261102/review_annual_input.csv');old=pd.read_csv('outputs/biomass_exclusion_reconciliation_20261102/annual_donor_corrections.csv');assert len(old)==45
counts=pd.read_csv('data/MBON/SBCMBON_kelp_forest_integrated_quad_and_swath_20231022.csv',usecols=['data_source','date','site_id','taxon_name']);counts=counts[counts.data_source.eq('kfm')];counts['year']=pd.to_datetime(counts.date).dt.year
changes=[];lineage=[]
for r in need.itertuples():
 sid=r.series_id;red='M. franciscanus' in sid;power=2.917 if red else 2.870;sp='Strongylocentrotus franciscanus' if red else 'Strongylocentrotus purpuratus';harris=r.combined_label=='a-k-01';allowed='Wyckoff Ledge' if harris else 'SE Sea Lion Rookery';hyp=['Wyckoff Ledge','Miracle Mile'] if harris else ['Southeast Reef'];site='a-k-01' if harris else 'a-k-14'
 pool=p[p.series_id.eq(sid)&p.year.eq(r.year)&p.status.eq(r.donor_role)];target=pool[pool.campus.eq('KFM')&pool.site.eq(site)];s=src[src.year.eq(r.year)&src.scientificName.eq(sp)]
 if '|PISCO|' in sid:s=s[s.mm.ge(25)]
 suspect=s[s.siteName.isin(hyp)];repl=s[s.siteName.eq(allowed)];assert len(repl)>0 and not repl.siteName.isin(['Miracle Mile','Southeast Reef']).any()
 delta=target.groupby('size')['count'].sum().subtract(suspect.groupby('mm')['count'].sum(),fill_value=0);assert delta.abs().sum()==0
 keep=pool[~(pool.campus.eq('KFM')&pool.site.eq(site))];mass=lambda x:.00059*x**power
 oldmean=np.average(mass(pool['size']),weights=pool['count']);newmean=((mass(keep['size'])*keep['count']).sum()+(mass(repl.mm)*repl['count']).sum())/(keep['count'].sum()+repl['count'].sum())
 vv={st:inv.loc[(sid,r.year,st),'value'] for st in ['mpa','reference']};vv[r.donor_role]*=newmean/oldmean;new=np.log((vv['mpa']+off[sid])/(vv['reference']+off[sid]));original=a.loc[a.series_id.eq(sid)&a.year.eq(r.year),'lnDiff'];assert len(original)==1
 original_v={st:inv.loc[(sid,r.year,st),'value'] for st in ['mpa','reference']};assert abs(np.log((original_v['mpa']+off[sid])/(original_v['reference']+off[sid]))-original.iloc[0])<1e-10
 cd=pd.to_datetime(counts.loc[counts.site_id.eq(site)&counts.year.eq(r.year)&counts.taxon_name.eq('Mesocentrotus franciscanus' if red else sp),'date'].unique());sd=pd.to_datetime(repl.startDate.unique());gap=max(min(abs((d-c).days) for c in cd) for d in sd);assert gap<=1
 changes.append(dict(series_id=sid,year=r.year,status=r.donor_role,delivered_site=site,allowed_source_site=allowed,excluded_source_site='Miracle Mile' if harris else 'Southeast Reef',histogram_absolute_difference=float(delta.abs().sum()),removed_original_records=int(target['count'].sum()),allowed_replacement_records=int(repl['count'].sum()),old_mean_mass=oldmean,new_mean_mass=newmean,mass_factor=newmean/oldmean,old_lnDiff=original.iloc[0],new_lnDiff=new,source_dates=';'.join(sorted(str(x)[:10] for x in repl.startDate.unique())),max_source_count_gap_days=gap,cross_program_borrowing='|PISCO|' in sid,sparse_lt10=repl['count'].sum()<10))
 for rr in repl.itertuples():lineage.append(dict(series_id=sid,year=r.year,source_excel_row=rr.excel_row,site=allowed,size_mm=rr.mm,frequency=rr.count))
c=pd.DataFrame(changes);c.to_csv(O/'recent_corrections.csv',index=False);pd.DataFrame(lineage).to_csv(S.parent/'replacement_lineage.csv',index=False)
allc=pd.concat([old,c],ignore_index=True);assert len(allc)==99 and not allc.duplicated(['series_id','year']).any();allc.to_csv(O/'all_corrections.csv',index=False)
tree=ast.parse(Path('code/tools/audit_station_coverage_refits.py').read_text());exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in ['fit','slope']],type_ignores=[]),'fit','exec'))
effects=[]
for label,patch in [('original_review',old.iloc[:0]),('historical_only',old),('all_reconciled',allc)]:
 v=a.copy()
 for r in patch.itertuples():v.loc[v.series_id.eq(r.series_id)&v.year.eq(r.year),'lnDiff']=r.new_lnDiff
 if label=='all_reconciled':
  untouched=~pd.MultiIndex.from_frame(v[['series_id','year']]).isin(pd.MultiIndex.from_frame(allc[['series_id','year']]))
  pd.testing.assert_frame_equal(v.loc[untouched],a.loc[untouched]);assert v.drop(columns='lnDiff').equals(a.drop(columns='lnDiff'))
  v.to_csv(O/'review_annual_input.csv',index=False)
 for sid in need.series_id.unique():effects.append(dict(series_id=sid,scenario=label,**slope(v[v.series_id.eq(sid)&v.t.gt(0)],'lnDiff')))
e=pd.DataFrame(effects);e.to_csv(O/'effects.csv',index=False)
base=e[e.scenario.eq('original_review')].set_index('series_id');partial=e[e.scenario.eq('historical_only')].set_index('series_id');final=e[e.scenario.eq('all_reconciled')].set_index('series_id');summary=final[['estimate','low','high','n_years']].rename(columns={'estimate':'corrected_slope'});summary['original_slope']=base.estimate;summary['historical_only_slope']=partial.estimate;summary['recent_increment']=summary.corrected_slope-partial.estimate;summary['total_change']=summary.corrected_slope-base.estimate;summary.to_csv(O/'comparison.csv')
h['original_biomass_status']=h.biomass_status;h['alias_hold_cleared']=True;h['primary_adopted']=False;h['biomass_status']='Source identity resolved; corrected review refit complete; primary promotion pending';h['remaining_support_limits']='Fixed allometry; shared size pools; conditional temporal uncertainty';h.loc[h.series_id.str.contains('PISCO'),'remaining_support_limits']+='; KFM sizes borrowed for PISCO counts';h.loc[h.original_biomass_status.str.contains('verified zero'),'biomass_status']='Verified zero reference; no donor mass needed; unchanged';h.to_csv(O/'kfm_current_holds.csv',index=False)
fig,axes=plt.subplots(2,3,figsize=(11,6),layout='constrained')
for ax,(sid,g) in zip(axes.flat,e.groupby('series_id',sort=False)):
 for j,label in enumerate(['original_review','historical_only','all_reconciled']):
  z=g[g.scenario.eq(label)].iloc[0];ax.errorbar(z.estimate,j,xerr=[[z.estimate-z.low],[z.high-z.estimate]],fmt='o',color=['#999999','#e69f00','#0072b2'][j],capsize=3)
 ax.axvline(0,color='#bbbbbb',lw=.8);ax.set_yticks([]);ax.set_ylabel(sid.replace('|Bio','').replace('|','\n'));ax.set_xlabel('Annual change in log biomass ratio');ax.set_ylim(-.5,2.5)
fig.savefig(O/'effects.png',dpi=180);plt.close(fig)
subprocess.run([sys.executable,'code/tools/check_scott_reconciliation.py'],check=True)
checks={path:hashlib.sha256(Path(path).read_bytes()).hexdigest()==value for path,value in json.loads(Path('outputs/researcher_scope_20261101/protected_hashes.json').read_text()).items()};assert all(checks.values())
report='''# KFM original-size reconciliation\n\nAll 54 held records reproduce Scott’s exact size-frequency distributions: Wyckoff Ledge had been combined with excluded Miracle Mile; Southeast Reef had been labeled as Southeast Sea Lion. Replacements use only Wyckoff Ledge or SE Sea Lion Rookery. All 54 have replacements. Together with 45 previously verified historical corrections and one verified zero, this resolves the donor-identity issue for all 100 affected post-protection records.\n\nCounts, comparison years, offsets, other donors, model structure and temporal uncertainty method are unchanged. These are isolated review results; the manuscript and primary remain unchanged.\n\n## What changes\n\n'''+('```text\n'+summary.round(5).to_string()+'\n```')+'''\n\n![Comparison of biomass slopes](effects.png)\n\nFigure key: gray = original review; orange = historical corrections only; blue = all verified replacements. Bars are conditional 95% calendar-year HAC intervals. Points describe recovery through observed post-protection years, not before–after causal effects.\n\n## Source and assumptions\n\nExact millimeter measurements and frequency counts are retained. Duplicate-looking source rows are not removed: the file lacks transect IDs, and retaining them reproduces the delivered frequencies exactly. Date-coerced size-bin fields are ignored. Red and purple urchin mass uses 0.00059 × size(mm)^2.917 and ^2.870 respectively, applying the transformation to individual sizes before averaging. PISCO retains its existing ≥25 mm filter. Source rows shared by histories remain the same physical observations.\n\nKFM size/count dates match or differ by one day; both are retained under the approved survey-timing approach. This check concerns the KFM source counts, not temporal equivalence to PISCO recipients. KFM counts are quadrats while sizes are band-transect samples. Two PISCO histories still borrow KFM sizes. Santa Barbara red urchins in 2021 have only seven allowed measured animals. These support and uncertainty limitations remain; correcting site identity does not solve total biomass uncertainty.\n\n## Reproduce\n\nRun `python3 code/tools/reconcile_scott_kfm_sizes.py`, then `Rscript code/tools/verify_scott_kfm_refits.R`. Original workbook is preserved under `data/external/scott_gabara_2026-11-02/`; its source hash is recorded in verification.json. Row-level replacement lineage is alongside that restricted source. Public tables contain aggregate corrections. Rebuild the consolidated guide with `python3 code/tools/report_dataset_inventory.py`.\n'''
(O/'README.md').write_text(report);(O/'index.html').write_text('<meta charset="utf-8"><style>body{font:16px system-ui;max-width:1200px;margin:30px auto}table{border-collapse:collapse}td,th{padding:7px;border:1px solid #ccc}img{max-width:100%}</style><h1>KFM size reconciliation</h1><p>All 54 recent donor records resolved. Excluded physical sites removed. Paper unchanged.</p>'+summary.round(5).to_html()+'<img src="effects.png"><p>Gray: original; orange: historical corrections; blue: all corrections. Conditional 95% temporal intervals; shared-size and allometric uncertainty remain.</p><p><a href="README.md">Methods, limitations and reproducibility</a></p>')
(O/'verification.json').write_text(json.dumps(dict(source_sha256=hashlib.sha256(S.read_bytes()).hexdigest(),exact_histograms=54,historical_corrections=45,verified_zero=1,protected_files_unchanged=len(checks),output_hashes={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in O.iterdir() if f.name!='verification.json'}),indent=2))
print(summary.to_string())
