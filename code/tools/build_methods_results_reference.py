"""Generate site and response tables for the current analysis page."""
from pathlib import Path
import hashlib,json
import pandas as pd
R=Path(__file__).resolve().parents[2];I=R/'outputs/dataset_inventory_20261101';J=R/'outputs/integrated_review_20261103';B=R/'outputs/biomass_policy_review_20261103';O=R/'outputs/methods_results_reference_20261104';O.mkdir(exist_ok=True)
paths=[I/'comparison_decisions.csv',I/'site_response_coverage.csv',I/'mpa_metadata.csv',I/'current_rules.csv',J/'station_assignments.csv',J/'comparison_year_roster.csv',J/'model_eligibility.csv',J/'annual_values.csv',J/'pooled_results.csv',B/'comparison_support.csv']
d,raw,meta,rules,sta,ro,elig,ann,pooled,bio=[pd.read_csv(p,low_memory=False) for p in paths]
# Response history: candidate inclusion is separate from survey support and model eligibility.
assert d.series_id.is_unique and set(elig.series_id).issubset(set(d.series_id))
mm=meta.drop_duplicates('MPA').set_index('MPA');bm=bio.set_index('series_id');ev=elig.set_index('series_id');ar=ann.groupby('series_id').agg(supported_years=('response_supported','sum'),candidate_survey_years=('year','size')).reset_index().set_index('series_id')
rows=[]
for z in d.itertuples():
 sid=z.series_id;E=ev.loc[sid] if sid in ev.index else None;A=ar.loc[sid] if sid in ar.index else None;M=mm.loc[z.MPA] if z.MPA in mm.index else None;Q=bm.loc[sid] if sid in bm.index else None
 rows.append(dict(series_id=sid,MPA=z.MPA,program=z.Source,taxon=z.Taxa,response=z.Resp,decision='excluded_by_researcher_or_ecological_scope' if z.exclude else 'retained_candidate',decision_reason=z.reason,protected_sites=z.mapped_protected_sites,reference_sites=z.mapped_reference_sites,protection_date=M.protection_or_network_date if M is not None else '',protection_class=M.protection_class if M is not None else '',before_years=int(E.before_years) if E is not None else 0,after_years=int(E.after_years) if E is not None else 0,supported_years=int(A.supported_years) if A is not None else 0,model=E.eligible_model if E is not None else 'excluded_or_no_paired_coverage',survey_years=E.years if E is not None else '',unsupported_response_years=E.missing_response_years if E is not None else '',borrowed_biomass_years=int(Q.supported_borrowed_years) if Q is not None else 0,biomass_candidate_years=int(Q.candidate_years) if Q is not None else 0,biomass_supported_years=int(Q.supported_years) if Q is not None else 0,reference_evidence=z.reference_evidence,reference_sensitivity=getattr(z,'reference_sensitivity',''),issue_ids=z.issue_ids))
h=pd.DataFrame(rows);h.to_csv(O/'response_histories.csv',index=False)
# Site or site-group records: keep the raw inventory and attach exact response-specific roster decisions.
raw['program']=raw.program.fillna('');raw['site_id']=raw.site_id.fillna('');raw['MPA']=raw.MPA.fillna('[blank]')
site_names=raw.groupby(['program','MPA','site_id']).agg(site_name=('site_name',lambda s:next((str(v) for v in s if pd.notna(v) and str(v).strip()),'')),first_raw_year=('observed_years',lambda s:min((int(y) for v in s.dropna() for y in str(v).replace(',',';').split(';') if y.strip().isdigit()),default=None)),last_raw_year=('observed_years',lambda s:max((int(y) for v in s.dropna() for y in str(v).replace(',',';').split(';') if y.strip().isdigit()),default=None)),historical_assignment=('historical_assignment_status',lambda s:'; '.join(sorted(set(map(str,s.dropna()))))),historical_reason=('historical_assignment_reason',lambda s:'; '.join(sorted(set(map(str,s.dropna()))))),size_status=('size_status',lambda s:'; '.join(sorted(set(map(str,s.dropna()))))),depth_min_m=('depth_min_m','min'),depth_max_m=('depth_max_m','max')).reset_index()
sta=sta[sta.site.notna()&sta.site.ne('')].copy();sta['MPA']=sta.series_id.str.split('|').str[0];sta['program']=sta.source;sta['site_id']=sta.site
assigned=[]
for (program,mpa,site),g in sta.groupby(['program','MPA','site_id']):
 linked=g.merge(h[['series_id','taxon','response']],on='series_id',how='left',validate='many_to_one');lab=lambda v:f'{v.taxon} {v.response}'
 inc=sorted(set(lab(v) for v in linked[linked.included].itertuples()));exc=sorted(set(lab(v) for v in linked[~linked.included].itertuples()));roles=sorted(set(g.role.dropna()))
 assigned.append(dict(program=program,MPA=mpa,site_id=site,role='; '.join(roles),roster_status='included_for_some_responses' if inc else 'excluded_from_retained_responses',included_responses='; '.join(inc),excluded_responses='; '.join(exc),included_series_ids=';'.join(sorted(set(g.loc[g.included,'series_id']))),excluded_series_ids=';'.join(sorted(set(g.loc[~g.included,'series_id']))),decision_reasons='; '.join(sorted(set(g.reason.dropna())))))
a=pd.DataFrame(assigned);sites=site_names.merge(a,on=['program','MPA','site_id'],how='outer',validate='one_to_one')
sites['roster_status']=sites.roster_status.fillna('not_in_current_comparison_roster');sites['role']=sites.role.fillna('not_assigned');sites.to_csv(O/'site_decisions.csv',index=False)
# The current analysis page now renders the writing summary with the site and response tables.
checks=dict(candidate_histories=len(h),excluded_histories=int((h.decision!='retained_candidate').sum()),retained_histories=int((h.decision=='retained_candidate').sum()),MPAs=h.MPA.nunique(),site_or_group_rows=len(sites),assigned_sites=len(a),model_eligibility_rows=len(elig),input_hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
assert checks['candidate_histories']==206 and checks['excluded_histories']==31 and checks['retained_histories']==175
assert not sites.duplicated(['program','MPA','site_id']).any();assert len(h[h.model.eq('before_after_and_post_trend')])==int(elig.eligible_model.eq('before_after_and_post_trend').sum())
(O/'verification.json').write_text(json.dumps(checks,indent=2));print(checks)
