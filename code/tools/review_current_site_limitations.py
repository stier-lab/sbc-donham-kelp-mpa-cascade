"""Summarize current support and outstanding work without changing selection or estimates."""
from pathlib import Path
import json,hashlib
import pandas as pd
R=Path(__file__).resolve().parents[2]; I=R/'outputs/dataset_inventory_20261101'; B=R/'outputs/integrated_review_20261103'; O=R/'outputs/site_limitations_current';O.mkdir(exist_ok=True)
c=pd.read_csv(I/'comparison_decisions.csv');a=pd.read_csv(B/'annual_values.csv');m=pd.read_csv(B/'model_eligibility.csv');v=pd.read_csv(B/'station_annual_values.csv')
x=v.merge(a[['series_id','year','t','response_supported']],on=['series_id','year'],validate='many_to_one').merge(m[['series_id','eligible_model']],on='series_id',validate='many_to_one')
x=x[x.response_supported & x.value.notna() & x.eligible_model.ne('insufficient_post_years') & (x.eligible_model.eq('before_after_and_post_trend')|x.t.ge(0))]
check={'a-k-05':'Rodes Reef','a-k-01':'Wyckoff Ledge','a-k-02':'Hare Rock','a-k-21':'Miracle Mile','a-k-35':"Webster’s Arch",'LECHUZA':'LECHUZA'}
rows=[]
for site,name in check.items():
 q=x[x.site.eq(site)]
 if q.empty:rows.append(dict(site=site,site_name=name,series_id='',role='',years=0,first_year=None,last_year=None,note='No supported station values entering currently eligible fits'))
 for (sid,role),g in q.groupby(['series_id','status']):rows.append(dict(site=site,site_name=name,series_id=sid,role=role,years=g.year.nunique(),first_year=g.year.min(),last_year=g.year.max(),note='Supported station values in eligible model years; not independent samples'))
pd.DataFrame(rows).to_csv(O/'named_station_contributions.csv',index=False)
summary=[]
for (mpa,source),g in c.groupby(['MPA','Source']):
 kept=g[~g.exclude];ids=set(kept.series_id);mfit=m[m.series_id.isin(ids)];fitted=mfit[mfit.eligible_model.ne('insufficient_post_years')]
 local=['Review model diagnostics, ratio offsets and downstream consistency'] if len(kept) else []
 external=[];disclose=[]
 if not len(kept):disclose.append('Existing researcher/ecological exclusions remain; no new answer required to retain them')
 else:
  local.append('Current station/year membership and support are reconstructable from saved tables')
  if source=='PISCO':
   local.append('Summarize measured depths, zones, dates and station changes on the current roster; no new seasonal cutoff')
   disclose.append('Recorded depth/zone matching does not establish equal habitat or environmental exposure')
   if mpa=='Naples SMCA':disclose.append('Jenn already supplied expert/nMDS rationale and lack of a close match; ordination files optional unless quantitative validation is claimed')
   elif 'Inherited mapping' in ' '.join(g.reference_evidence.fillna('')):disclose.append('Complete station-specific habitat rationale not recovered; do not claim validated matching')
   if mpa=='Point Dume SMCA':external.append('PISCO metadata custodian: LECHUZA coordinate identity if distance or historical boundary assignment is used')
   if mpa in ['Point Vicente SMCA','South Point SMR']:disclose.append('Unrecovered stipes/fish lengths stay missing; new field sheets required to recover those values')
   if mpa=='Point Dume SMR':disclose.append('Little Dume 2011 missing stipes do not enter current annual station estimates; do not treat recovery as a present model blocker')
  if source=='KFM':
   external.append('Scott Gabara/park: data-use permit and redistribution conditions before submission/archive')
   disclose.append('Missing physical-depth metadata; no lobster lengths; restricted roving sheephead sizes do not certify transect population equivalence')
   if mpa=='Harris Point SMR':
    local.append('Review completed Wyckoff/Miracle joint refits and matched-year decomposition; reference decision still pending Scott')
    external.append('Scott Gabara: ecological recommendation for north-side Hare Rock versus south-side Wyckoff/Miracle Mile, by response')
   if mpa=='Santa Barbara Island SMR':
    local.append('Review completed Webster joint refits; fish-scope decision still pending Scott')
    external.append('Scott Gabara: whether pinnacle-habitat concern also rules out Webster’s Arch for fish')
  if source=='LTER':
   local.append('Recheck existing fixed-transect sensitivity and current raw adult-kelp conversion')
   disclose.append('Station-specific ecological rationale and actual depths incomplete; group-mean sizes cannot recover individual size variation')
   external.append('Li/program metadata only if stronger station-habitat claims or replacement by provider biomass are needed; fixed July/August discrepancy conditional on provider substitution')
  if source=='Landsat':
   local.append('Audit quarterly coverage from cached 152-quarter NetCDF and saved masks; assess reference weighting/timing without reclustering')
   disclose.append('Saved masks verified; historical clustering and ecological equivalence not independently established; pre-2012 selection is not uniformly predesignation')
 summary.append(dict(MPA=mpa,Source=source,candidate_histories=len(g),retained_histories=len(kept),fitted_histories=len(fitted),resolvable_or_testable_from_files='; '.join(local),coauthor_input_and_trigger='; '.join(external),retain_as_disclosure_or_missing='; '.join(disclose)))
pd.DataFrame(summary).to_csv(O/'site_program_actions.csv',index=False)
paths=[I/'comparison_decisions.csv',B/'annual_values.csv',B/'model_eligibility.csv',B/'station_annual_values.csv',I/'unresolved_issues.csv']
manifest={str(p.relative_to(R)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
assert len(c)==206 and c.exclude.sum()==31 and c.MPA.nunique()==25 and len(summary)==46
(O/'verification.json').write_text(json.dumps(dict(candidate_histories=len(c),excluded=int(c.exclude.sum()),protected_areas=c.MPA.nunique(),area_program_groups=len(summary),station_year_scope='Supported annual values; before years excluded from post-only fits',inputs=manifest),indent=2)+'\n')
print('Reviewed 25 protected areas / 46 area-program groups; named station contributions recorded.')
