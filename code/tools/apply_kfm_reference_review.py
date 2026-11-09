"""Expose completed reference tests without adopting their omissions."""
from pathlib import Path
import pandas as pd

def apply(views):
 report=Path('outputs/kfm_reference_sensitivity_20261103');out=Path('outputs/dataset_inventory_20261101')
 if not (report/'verification.json').exists():return views
 c=views['Eligibility decisions'].copy();c['reference_sensitivity']=''
 h=c.Source.eq('KFM')&c.MPA.eq('Harris Point SMR')&~c.exclude
 c.loc[h,'reference_sensitivity']='Wyckoff Ledge is the sole retained reference for this response; omitting it removes the paired history. Sensitivity only; no new exclusion adopted.'
 c.loc[h&c.Taxa.eq('S. pulcher'),'reference_sensitivity']='Same 15 post years: current-reference sheephead trend -16.5%/year; Miracle Mile alone -13.8%. Omitting both references removes Harris Point. Ecological choice pending Scott; no new exclusion adopted.'
 s=c.Source.eq('KFM')&c.MPA.eq('Santa Barbara Island SMR')&~c.exclude
 c.loc[s,'reference_sensitivity']='Webster’s Arch is already excluded for benthic responses; no direct annual input change in its omission test.'
 c.loc[s&c.Taxa.eq('S. pulcher'),'reference_sensitivity']='Omitting Webster’s Arch changes sheephead trend +3.4% to -0.4%/year, same 20 years; both intervals include zero. Combined omission also changes shared covariance (AR boundary). Fish-scope decision pending Scott.'
 c.to_csv(out/'comparison_decisions.csv',index=False);views['Eligibility decisions']=c.fillna('')
 q=views['Questions'].copy();i=q.issue_id.eq('U15');q.loc[i,'finding']='Four named-reference sensitivities and two matched-main joint refits completed. Wyckoff is sole reference for six Harris responses; sheephead retains 15 post years with Miracle Mile alone. Webster affects SBI sheephead only; Rodes has no current uses. Combined omission loses seven histories and 217 supported comparison-years. Kelp/urchin post trends change little, but the flat-baseline lobster contrast grows through loss of Harris. No new omission adopted.';q.loc[i,'next_action']='Show Scott named sites and completed sensitivity results; resolve Harris ecological comparability and Webster fish scope before finalizing affected results. Draft prepared, not sent.';q.loc[i,'evidence']=str(report/'README.md');q.to_csv(out/'unresolved_issues.csv',index=False);views['Questions']=q.fillna('')
 for f in ['coverage_changes','sheephead_trends','pooled_decomposition','diagnostics']:
  pd.read_csv(report/(f+'.csv')).to_csv(out/('kfm_reference_'+f+'.csv'),index=False)
 p=out/'README.md';t=p.read_text();line='- [KFM reference-site tests and draft note](../kfm_reference_sensitivity_20261103/index.html) — completed sensitivities; no new exclusions adopted.\n'
 if line not in t:t=t.replace('## Earlier scoped audits',line+'\n## Earlier scoped audits');p.write_text(t)
 return views
