from pathlib import Path
import json,hashlib,html
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=Path('outputs/biomass_resolution_20261102');e=pd.read_csv(O/'matched_donor_effects.csv');l=pd.read_csv(O/'donor_coverage_losses.csv');r=pd.read_csv(O/'red_urchin_summary.csv');d=pd.read_csv(O/'conditional_size_draws.csv');i=pd.read_csv(O/'red_urchin_role_years.csv')
# Distinguish loss of years from changing donors on those same years.
cur=e[e.scenario.eq('current')][['series_id','estimate']].rename(columns={'estimate':'all_years_original'})
comp=e[e.scenario.isin(['same_program','same_count_site'])].merge(cur,on='series_id',validate='many_to_one')
comp['change_from_year_loss']=comp.matched_original-comp.all_years_original;comp['change_from_donors']=comp.estimate-comp.matched_original
comp.to_csv(O/'separated_changes.csv',index=False)
summary=r.merge(l[l.scenario.eq('same_count_site')].drop(columns='scenario'),on='series_id',how='left',validate='one_to_one');summary['next_action']=np.where(summary.fraction_same_count_site.eq(0),'No count-site size donors: retain density, keep biomass provisional',np.where(summary.fraction_complete_dates.lt(1),'Recover original survey dates; assess count-size timing','Check sparse pools and size coverage against count protocol'))
summary.to_csv(O/'red_urchin_decisions.csv',index=False)
z=comp[comp.scenario.eq('same_count_site') & comp.series_id.str.contains('M. franciscanus',regex=False)].sort_values('series_id').reset_index(drop=True)
fig,ax=plt.subplots(figsize=(8,6));yy=np.arange(len(z))
for k,row in z.iterrows():ax.plot([row.matched_original,row.estimate],[k,k],color='#a5a5a5',lw=1.6)
ax.scatter(z.matched_original,yy,marker='o',facecolors='white',edgecolors='#333333',s=45,zorder=3);ax.scatter(z.estimate,yy,marker='D',color='#007f83',s=33,zorder=4)
labels=[s.replace('|M. franciscanus|',' · ').replace('|Bio','') for s in z.series_id]
ax.set_yticks(yy,labels,fontsize=9);ax.axvline(0,color='#bbbbbb',lw=.8);ax.set_xlabel('Annual change in log protected/reference biomass ratio');ax.spines[['top','right']].set_visible(False);fig.tight_layout();fig.savefig(O/'red_urchin_matched_donors.png',dpi=180);fig.savefig(O/'red_urchin_matched_donors.pdf');plt.close(fig)
text='''# Biomass support: corrections and red-urchin checks

The review pipeline now withholds the Point Vicente kelp biomass comparison in 2023 because two reference-site plants lack stipe measurements. Their biomass is unknown, not zero. The whole paired year is withheld to avoid silently changing reference weights. Raw counts are retained. This produces a separate corrected review input (3,224 annual rows); the manuscript and frozen results are unchanged.

The previously checked kelp recovery slope changes from 0.243 to 0.342 log-ratio units per year; both uncertainty intervals include zero. The upstream PISCO processor now also preserves the missing-stipe flags and invalidates affected station and reference-role biomass. Its exact guard expressions were tested against the real raw and reconstructed station records. A full production rebuild remains on hold with the paper.

## Sparse sizes and borrowed sizes

We tested all 35 estimable pooled-size animal biomass recovery histories. Two additional pooled-size histories have too few post-protection years for these slope checks; two LTER lobster histories use direct lengths and are outside this pooled-size calculation. No researcher exclusion was reversed.

Shared source records receive the same random weight wherever they are used. The 500 individual-weight draws and separate 500 recorded-group draws describe how size sampling alone affects slope estimates. They are not new confidence intervals: count uncertainty, allometric uncertainty and donor representativeness are not included. Recorded group identifiers can merge multiple visits when metadata are incomplete. Existing offsets and observed years are held fixed.

A single measured animal cannot reveal its population's size variation. We therefore also show a deliberately broad stress test: halve or double singleton mean masses, with opposite changes in protected and reference pools. This is an assumption check, not a probability interval or exclusion rule. Campus Point PISCO red-urchin slopes span −0.048 to 0.040 under this test, compared with −0.005 originally. Matlahuayl red-urchin slopes span −0.392 to −0.117, compared with −0.238. These histories need caution even when a conditional bootstrap looks precise.

Keeping only donors from actual count sites retains 298 of 461 post-protection comparison-years across the 35 histories. Keeping only donors from the same program retains 300. Both restrictions permit 23 slope fits. Among those fits, the largest slope change attributable to donor choice on identical years is 0.0241 log-ratio units per year. This does not establish that omitted years or borrowed donors are harmless. The separated_changes table reports the effect of losing years separately from changing donors.

## Red urchins

All 17 red-urchin biomass histories have a dedicated role-year table and decision table. The existing equation is mass in grams = 0.00059 × diameter in millimeters^2.917. PISCO pools retain diameters of at least 25 mm; KFM and LTER use unfiltered sizes. The latter therefore represent a different measured size range. We have not removed small animals merely to make the numbers agree: any standardization must match BOTH counts and sizes.

Campus Point PISCO and LTER reuse singleton size measurements in five protected years. Those measurements do not constitute independent replication. Several island PISCO histories borrow heavily from KFM; South Point PISCO has no donors from its own count sites. The donor tables identify these explicitly instead of treating an available annual size pool as proof of local support.

Missing KFM dates remain unresolved in the combined delivery. The legacy archive contains survey dates, but linking them safely requires a verified species-code/unit crosswalk and occurrence-level matching. Matching just year, station and diameter could confuse different animals. No dates were guessed or copied into the source data.

## Recommended handling

- Use the safeguarded PISCO processor when building the next production release; use the corrected review input meanwhile.
- Retain valid density histories. Keep borrowed-size biomass provisional until donor habitat, timing and size-range compatibility are justified.
- Keep sparse-pool flags visible; do not adopt an arbitrary minimum sample cutoff or exclude by result direction.
- Carry original donor identity into joint uncertainty calculations. Complete count-and-size propagation and KFM kelp plant-sampling checks before freezing biomass uncertainty.
- For red urchins, resolve KFM date/code/unit matching and confirm that each count protocol covers the same sizes as its donor pool. The existing PISCO ≥25-mm restriction should not be imposed on other programs' size samples alone.

## Figure

Red-urchin recovery slopes on identical years: open circles use the current pooled donors; teal diamonds use only donors from actual count sites. Connecting lines show the change from donor choice. Only estimable matched-year fits appear. These are point estimates, not uncertainty intervals. Missing histories and lost years are listed in the coverage table.

## Verification

Independent R calculations reproduce 35 original slopes and uncertainty estimates, plus 46 matched-donor refits. All 487 protected manuscript/result files retain their expected hashes. The shared-donor calculation preserves 77,382 original source-row identities across 104,125 uses. The raw data were not edited. The source repository PISCO processor received the missing-stipe safeguard; its other outputs were not rebuilt. This is a review result, not a final primary analysis.
'''
(O/'README.md').write_text(text)
# Plain readable report with the full searchable tables as native HTML.
body=''.join('<p>'+html.escape(p).replace('\n','<br>')+'</p>' for p in text.split('\n\n'))
body+='<img src="red_urchin_matched_donors.png" style="max-width:100%" alt="Red urchin slopes using current and count-site donors on matched years">'
for title,file in [('Red urchin decisions','red_urchin_decisions.csv'),('Separate year losses and donor changes','separated_changes.csv'),('Size sampling only: draw ranges','conditional_size_draws.csv')]:body+=f'<h2>{title}</h2><p><a href="{file}">Download table</a></p><div style="overflow:auto">'+pd.read_csv(O/file).to_html(index=False,float_format=lambda x:f'{x:.4g}')+'</div>'
(O/'index.html').write_text('<!doctype html><meta charset="utf-8"><title>Biomass support and red urchins</title><style>body{max-width:1100px;margin:30px auto;padding:0 20px;font:17px/1.5 system-ui}td,th{padding:8px;text-align:left}table{border-collapse:collapse;font-size:13px}tr:nth-child(even){background:#eef4f4}</style>'+body)
inputs=[Path('/Users/adrianstier/Donham-Stier-CA-MPA-Data-2026/code/R/04_pisco_processing.R'),Path('code/tools/verify_missing_stipe_guard.R'),Path('code/tools/resolve_biomass_support.py'),Path('code/tools/verify_biomass_resolution.R'),Path('code/tools/report_biomass_resolution.py'),Path('outputs/biomass_support_20261102/private/pool_members.csv'),Path('outputs/biomass_support_20261102/animal_pool_inventory.csv'),Path('outputs/pisco_fish_rules_20261101/approved_annual_input.csv')]
(O/'input_hashes.json').write_text(json.dumps({str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},indent=2));(O/'manifest.json').write_text(json.dumps({str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in O.iterdir() if p.is_file() and p.name!='manifest.json'},indent=2))
