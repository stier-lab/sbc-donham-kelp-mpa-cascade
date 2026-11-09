"""Apply user-directed researcher/ecological exclusions before review-only fitting."""
from pathlib import Path
import pandas as pd,json,hashlib,shutil
root=Path.cwd();base=root/'outputs/upstream_sites_20261101';out=root/'outputs/researcher_scope_20261101';out.mkdir(exist_ok=True)
a=pd.read_csv(base/'annual_input.csv');roster=a[['series_id','MPA','Source','Taxa','Resp']].drop_duplicates();roster['exclude']=False;roster['reason']='Retain: source and coverage rules; no baseline-trend or direction screen'
rules=[(roster.MPA.eq('Anacapa Island SMCA'),'Original ecological scope: permitted lobster take; distinct from retained Anacapa SMR'),(roster.MPA.eq('Painted Cave SMCA'),'Original ecological scope: permitted lobster take'),(roster.MPA.eq('Point Conception SMR'),'Original ecological scope: focal predators rarely observed; satellite validity does not override scope'),(roster.MPA.eq('Gull Island SMR')&roster.Source.eq('KFM'),'Researcher advice: unsuitable/uncertain KFM reference; latest user decision supersedes earlier retention')]
for mask,reason in rules:roster.loc[mask,['exclude','reason']]=[True,reason]
roster.to_csv(out/'selection_decisions.csv',index=False)
plan='''# Researcher advice and ecological scope: review-only rules

User decision November1: preserve exclusions based on researcher advice and ecological scope, but not baseline trends or response direction. Written before fitting this cohort. This is not preregistration: preceding results have been examined.

Remove Anacapa Island SMCA, Painted Cave SMCA and Point Conception SMR across any available programs/responses. Remove Gull Island KFM only, preserving PISCO. The original scope exclusions at Arrow Point, Cat Harbor, Dana Point and access-only San Miguel remain in force; they are absent from the input roster. Retain the original selected PISCO references (including Scott Hamilton's exclusions), Jenn's revised Naples reference, and distinct Anacapa SMR identity. Do not impose a new literal all-species-no-take screen on historical retained SMCAs such as Naples/Swamis: this preserves the documented site-specific ecological decisions, including existing exceptions, rather than claiming a uniform no-take network.

Retain restored historical responses excluded for baseline trends; baseline p values and response direction never enter selection. Harris lobster suspected poor sampling from absent counts is not independently verified survey invalidity; its constant primary series remains non-estimable under the uniform variance rule. Unexplained inherited filters do not become expert advice. Point Dume sheephead and Matlahuayl valid reconstructed responses remain subject to normal coverage. Other unreconstructed candidates remain pending, not newly certified.

Judith Rock: retained because the documented overlap concerns an already-excluded San Miguel extraction; no independent researcher/ecological reason recovered. Separate omission run provided. This ambiguity is disclosed, not presented as collaborator endorsement.

All calculations preserve the preceding annual measurements, zero corrections, baseline design flags and conservative shared-reference cluster assignments. No reweighting through source reconstruction. Same primary strictly prior >=3years / paired years1–9 >=3years, excluding designation year; same flat baseline primary, trend sensitivity, full post horizon,10year baseline, original timing and HAC1/3/5 checks. Before-after change and post association remain separate; same9slot multiplicity families, MPA/program heterogeneity and CR2 uncertainty safeguards. No manuscript/SI or canonical changes.

Evidence: private selection_history/HISTORY.md; reference_assignments/email_evidence.json (Scott Gabara October12 and29, Emily October13 local time); original manuscript ecological scope; upstream_sites_20261101/README.md. See parent review recommendation and Anacapa comparability audit for context.
'''
(out/'PLAN.md').write_text(plan)
for dest,omit in [(out,False),(out/'omit_judith',True)]:
 dest.mkdir(exist_ok=True);keep=roster.loc[~roster.exclude,'series_id'];d=a[a.series_id.isin(keep)].copy()
 if omit:d=d[d.MPA!='Judith Rk SMR']
 d.to_csv(dest/'annual_input.csv',index=False);shutil.copyfile(base/'protected_hashes.json',dest/'protected_hashes.json')
 if omit:(dest/'PLAN.md').write_text(plan+'\nThis companion omits Judith Rock and changes nothing else.\n')
 files=[dest/'annual_input.csv',dest/'PLAN.md',out/'selection_decisions.csv',base/'annual_input.csv'];(dest/'freeze_manifest.json').write_text(json.dumps({str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},indent=2))
print(roster[roster.exclude].groupby(['MPA','Source']).size().to_string());print('Retained histories:',roster[~roster.exclude].shape[0])
