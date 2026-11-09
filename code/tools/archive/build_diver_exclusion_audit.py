from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,shutil
root=Path('/Users/adrianstier/sbc-donham-kelp-mpa-cascade');old=root/'outputs/pooled_common_estimand_20261101';out=root/'outputs/diver_exclusion_refit_20261101';out.mkdir(exist_ok=True)
audit=Path('/Users/adrianstier/sbc-donham-kelp-mpa-manuscript/correspondence/coauthor_reviews/2026-10-31_selection_rule_check')
d=pd.read_csv(old/'annual_input.csv');h=pd.read_csv(root/'data/harmonized/harmonized_response_ratios.csv').rename(columns={'CA_MPA_Name_Short':'MPA','source':'Source','resp':'Resp'})
lookup=d[['y','Taxa']].drop_duplicates();assert not lookup.y.duplicated().any();h=h.merge(lookup,on='y',validate='many_to_one');h['series_id']=h[['MPA','Taxa','Source','Resp']].agg('|'.join,axis=1)
o=pd.read_csv(audit/'response_omission_decisions.csv');o['series_id']=o[['MPA','Taxa','Source','Resp']].agg('|'.join,axis=1)
new=h[h.series_id.isin(o.series_id)].copy();assert new.series_id.nunique()==22;assert not set(new.series_id)&set(d.series_id)
m=pd.read_csv(root/'data/harmonized/harmonized_site_metadata.csv').rename(columns={'CA_MPA_Name_Short':'MPA'})[['MPA','MPA_Start']].drop_duplicates();new=new.merge(m,on='MPA',validate='many_to_one').merge(pd.read_csv(old/'reference_clusters.csv'),on='MPA',validate='many_to_one');new['t']=new.year-new.MPA_Start;new['baseline_design_verified']=new.Source!='PISCO'
r=pd.read_csv(root/'data/harmonized/harmonized_raw_responses.csv').rename(columns={'CA_MPA_Name_Short':'MPA','taxon_name':'y','source':'Source','resp':'Resp'});r=r.merge(lookup,on='y',validate='many_to_one');r['series_id']=r[['MPA','Taxa','Source','Resp']].agg('|'.join,axis=1)
assert not new.duplicated(['series_id','year']).any();assert np.isfinite(new[['lnDiff','mpa','reference','MPA_Start']]).all().all();assert (new[['mpa','reference']]>0).all().all();assert np.allclose(new.lnDiff,np.log(new.mpa/new.reference))
# Harris omitted ratios escaped the current raw-data protocol/calibration update.
hid='Harris Point SMR|M. pyrifera|KFM|Bio'
rr=r[r.series_id==hid].copy();mx=rr.value.max();floor=rr.value[rr.value>0].min()/2
rr['prop']=(rr.value+floor)/mx
wide=rr.pivot(index='year',columns='status',values='prop').dropna().reset_index()
base=new[new.series_id==hid].drop_duplicates('year').set_index('year')
rebuilt=base.reindex(wide.year).reset_index();assert rebuilt.MPA.notna().all()
rebuilt['mpa']=wide.mpa.to_numpy();rebuilt['reference']=wide.reference.to_numpy();rebuilt['Diff']=rebuilt.mpa/rebuilt.reference;rebuilt['lnDiff']=np.log(rebuilt.Diff)
new=pd.concat([new[new.series_id!=hid],rebuilt],ignore_index=True)
rawchecks=[]
for id,z in new.groupby('series_id'):
 rr=r[r.series_id==id];assert not rr.duplicated(['year','status']).any();paired=rr.pivot(index='year',columns='status',values='value').reindex(z.year);assert np.isfinite(paired[['mpa','reference']]).all().all();assert (paired[['mpa','reference']]>=0).all().all()
 floor=rr.value[rr.value>0].min()/2; derived=np.log((paired.mpa+floor)/(paired.reference+floor));assert np.allclose(derived,z.lnDiff,atol=1e-9),(id,max(abs(derived.to_numpy()-z.lnDiff.to_numpy())))
 pre=z.t<0;post=z.t.between(1,9);both=(paired.mpa==0)&(paired.reference==0)
 rawchecks.append(dict(series_id=id,before_years=int(pre.sum()),post_years=int(post.sum()),first_year=int(z.year.min()),last_year=int(z.year.max()),both_zero_before=int(both.to_numpy()[pre.to_numpy()].sum()),raw_paired_rows=len(paired),change_eligible=bool(pre.sum()>=3 and post.sum()>=3 and z.baseline_design_verified.all()),post_eligible=bool(post.sum()>=3),decision='Restore; coverage determines estimand; no slope/sign screen',source_check='Finite nonnegative raw protected/reference values for every ratio year; current harmonized transformation retained'))
ledger=o.merge(pd.DataFrame(rawchecks),on='series_id');ledger.loc[ledger.series_id==hid,'source_check']='Rebuilt from current monthly-calibrated raw biomass and modern 1996+ protocol; removed 12 stale earlier/missing ratio years; six strict Before and nine primary After years';ledger.to_csv(out/'omission_audit.csv',index=False)
full=pd.concat([d,new[d.columns]],ignore_index=True);assert full.series_id.nunique()==169;full.to_csv(out/'annual_input.csv',index=False)
shutil.copy(old/'reference_clusters.csv',out/'reference_clusters.csv');shutil.copy(old/'protected_hashes.json',out/'protected_hashes.json')
# Preserve the entire candidate-site evidence, with explicit unresolved upstream reconstruction scope.
site=pd.read_csv(audit/'candidate_mpa_rule_checks.csv');site.to_csv(out/'upstream_candidate_evidence.csv',index=False)
plan=(old/'PLAN.md').read_text();plan=plan.replace('Scope: existing129primary diver histories plus the newly audited18satellite histories. Do not restore other historical diver omissions without their own data/design audit. This scope limitation means the present refit resolves estimand mixing, not every historical inclusion decision.','Scope: 129previous diver histories +22audited omitted diver histories +18satellite histories. The added histories have finite paired raw values and the same established source/site assignments as retained responses. Eligibility uses coverage and design, never baseline significance or effect sign. All22 enter the candidate roster; each estimand applies its existing coverage rule. Upstream sites absent from harmonized data are separately documented, not declared scientifically ineligible merely because their reconstruction is incomplete.')
plan+='\nA history with zero or numerically negligible fitted sampling SE (<=1e-12 log units) cannot receive inverse-variance weight. Retain it in the audit with an explicit non-estimable-variance reason; do not invent precision or a variance floor. This applies to every history and scenario.\n'
plan+='\nAdditional prespecified sensitivity: omit Harris Point KFM lobster density (all-zero Before values), and separately omit all Gull Island KFM histories (contested reference retained by user decision). Both use the Primary time window/model. Raw zeros remain observed zeros; their corrected log ratios are pseudocount-dependent. Harris kelp rebuilt from current raw biomass, using the same adaptive correction. Twelve stale ratio years absent from the modern-protocol raw series were removed; no new allometry fitted.\n'
(out/'PLAN.md').write_text(plan)
files=[out/'PLAN.md',out/'annual_input.csv',out/'reference_clusters.csv',out/'omission_audit.csv',root/'data/harmonized/harmonized_raw_responses.csv',root/'data/harmonized/harmonized_response_ratios.csv']
(out/'freeze_manifest.json').write_text(json.dumps({str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},indent=2))
print(ledger[['MPA','Taxa','Resp','before_years','post_years','change_eligible','post_eligible','both_zero_before']].to_string(index=False))
