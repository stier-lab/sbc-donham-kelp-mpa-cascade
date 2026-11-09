"""Human-readable audit summary generated only from the integrated output tables."""
from common import *
import html
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
an=pd.read_csv(OUT/'annual_values.csv');el=pd.read_csv(OUT/'model_eligibility.csv');vis=pd.read_csv(OUT/'visit_biomass_support.csv');j=pd.read_csv(OUT/'joint_model_diagnostics.csv');pooled=pd.read_csv(OUT/'pooled_results.csv');cv=pd.read_csv(OUT/'pooled_coverage_sensitivities.csv');scope=pd.read_csv(ROOT/'config/integrated_comparison_scope.csv')
summary=an.groupby(['Source','Resp']).agg(comparison_years=('year','size'),supported_years=('response_supported','sum'),histories=('series_id','nunique')).reset_index();summary.to_csv(OUT/'coverage_summary.csv',index=False)
# Explicitly compare against the frozen manuscript input, not a prior audit surrogate.
old=pd.concat([read_tracked(ROOT/'data/harmonized/harmonized_response_ratios.csv'),read_tracked(ROOT/'data/harmonized/harmonized_landsat_rr.csv')],ignore_index=True)
old['Taxa']=old.y.str.strip().map(TAXA);assert old.Taxa.notna().all()
def normalized(m,t,s,r):return '|'.join(str(v).strip() for v in [m,t,s,r])
old['key']=[normalized(*v) for v in old[['CA_MPA_Name_Short','Taxa','source','resp']].itertuples(index=False,name=None)];old=old[['key','year','lnDiff']];assert not old.duplicated(['key','year']).any()
new=an.copy();new['key']=[normalized(*v) for v in new[['MPA','Taxa','Source','Resp']].itertuples(index=False,name=None)];comparison=new.merge(old,on=['key','year'],how='outer',indicator=True);comparison['ln_ratio_change']=comparison.ln_ratio-comparison.lnDiff;comparison.to_csv(OUT/'paper_input_comparison.csv',index=False)
def slope(years,values):
 q=np.isfinite(values);return np.polyfit(np.asarray(years)[q],np.asarray(values)[q],1)[0] if q.sum()>=3 else np.nan
rows=[]
for key,g in new.groupby('key'):
 meta=g.iloc[0];b=old[old.key.eq(key)&old.year.ge(meta.MPA_Start)];n=g[g.year.ge(meta.MPA_Start)&g.response_supported];shared=b.merge(n[['year','ln_ratio']],on='year');rows.append(dict(series_id=meta.series_id,Source=meta.Source,Taxa=meta.Taxa,Resp=meta.Resp,old_post_years=len(b),new_post_years=len(n),shared_post_years=len(shared),old_all_slope=slope(b.year,b.lnDiff),old_shared_slope=slope(shared.year,shared.lnDiff),new_shared_slope=slope(shared.year,shared.ln_ratio),new_all_slope=slope(n.year,n.ln_ratio)))
pd.DataFrame(rows).to_csv(OUT/'input_change_decomposition.csv',index=False)
limitations=[
 ('Observed sizes','Empirical size uncertainty is estimable within replicated samples. A singleton has no empirical within-sample variance.','Same-program/taxon sparse-size variance sensitivity; intervals remain conditional on that assumption.'),
 ('LTER urchin group means','Original measurements give a mean diameter for a size class, not every individual diameter.','Within-group variation and the bias from applying a nonlinear equation to a mean cannot be recovered from these files.'),
 ('Count sampling','Shared physical-unit resampling represents observed sampling variation. All-zero or unreplicated units cannot reveal their full sampling variance.','Residual process uncertainty remains; this is not a complete detection/count-observation likelihood.'),
 ('Count–size association','Direct measurements stay with their count unit; separate size surveys have no individual links to counted animals.','Their unobserved count–size covariance cannot be estimated; do not describe the intervals as total uncertainty.'),
 ('Borrowing','Reviewed PISCO target and calibration scope is enforced; size-sampling scatter is removed from residual prediction scatter before adding donor draws.','Prediction variance is conditional on estimated calibration and a lognormal working approximation; sparse calibration remains weak evidence.'),
 ('KFM matching','Same physical site and declared date rule; no borrowing across years or unrelated sites.','Unsupported biomass remains missing, with density retained. No new KFM lobster biomass is inferred.'),
 ('Shared references','Joint process terms and shared measurement draws connect documented physical stations.','Verified shared satellite pixels inform an overlap-weighted process kernel. Their spatial covariance and cross-program overlap remain unestimated; geometry alone is not a complete sampling covariance.'),
 ('Time','Joint models estimate AR1 dependence using actual calendar gaps; fixed-dependence sensitivities remain available.','Variance and AR-boundary fits require attention. Confidence limits condition on nuisance covariance estimates.'),
 ('Conversion equations','Point coefficients and observed monthly kelp coefficients are used consistently.','Coefficient covariance, detection bias, habitat mismatch and causal counterfactual validity are not measured by these intervals.'),
 ('Publication','This is a coordinated review build, not a manuscript release.','Pathway/resilience/figure outputs in the paper must be regenerated and checked against an approved cohort before publication; no paper files changed.')]
pd.DataFrame(limitations,columns=['component','what_is_included','remaining_limit']).to_csv(OUT/'uncertainty_limits.csv',index=False)
# No titles, subtitles or embedded keys/captions in figure artwork.
names={'M. pyrifera':'Kelp','M. franciscanus':'Red urchin','S. purpuratus':'Purple urchin','S. pulcher':'Sheephead','P. interruptus':'Lobster'}
for name,model,term,xlab in [('recovery','post_trend','recovery_slope','Annual change in protected/reference ratio (%)'),('before_after','flat_before_after','average_observed_post_change','Change relative to before protection (%)')]:
 z=pooled[pooled.Source.eq('All programs')&pooled.model.eq(model)&pooled.term.eq(term)].copy();z['label']=[names[t]+' '+('biomass' if r=='Bio' else 'density') for t,r in zip(z.Taxa,z.Resp)];z=z.sort_values('label');fig,ax=plt.subplots(figsize=(8,4.6));yy=np.arange(len(z));val=100*np.expm1(z.estimate);lo=100*np.expm1(z.low);hi=100*np.expm1(z.high);ax.errorbar(val,yy,xerr=np.vstack([val-lo,hi-val]),fmt='o',color='#246a73',ecolor='#246a73',capsize=3);ax.axvline(0,color='#999999',lw=.8);ax.set_yticks(yy,z.label);ax.set_xlabel(xlab);ax.spines[['top','right']].set_visible(False);fig.tight_layout();fig.savefig(OUT/(name+'.png'),dpi=180);fig.savefig(OUT/(name+'.pdf'));plt.close(fig)
checks=json.loads((OUT/'verification.json').read_text());mainj=j[j.analysis.eq('estimated_AR1')];unmatched=scope[~scope.exclude&~scope.series_id.isin(el.series_id)]
text=f'''# Integrated pipeline review

Updated November 3, 2026. The paper and frozen outputs are unchanged.

The integrated build starts with registered contributor/provider files, applies source repairs and response-specific exclusions, rebuilds survey units and biomass, and fits before–after changes separately from post-protection trends. Main coverage uses documented stations available each year and shared program zones. Calendar years determine time; the designation year is after protection. Three matched before years permit before–after estimation; four is a sensitivity. Flat-baseline recovery is the main before–after model, with a pre-trend sensitivity.

## What the build contains

- {len(an):,} comparison-years with both roles and shared recorded zones; {int(an.response_supported.sum()):,} have a supported response under the stated measurement rules.
- {el.series_id.nunique()} histories reach the response/eligibility table; {checks['fitted_histories']} have sufficient years for at least one fitted model.
- {int(scope.exclude.sum())} candidate histories retain their recorded researcher/ecological exclusions. Another {len(unmatched)} retained candidates lack a usable comparison in this build; see the roster rather than interpreting them as ecological exclusions.
- {int(vis.biomass_supported.sum()):,} of {len(vis):,} biomass visits have measurement support. Missing biomass never defines the density roster.
- {checks['checks']} integration checks passed; {checks['protected_files']} frozen files remain unchanged.

## How to read the results

[Coverage](coverage_summary.csv), [station assignments](station_assignments.csv), [comparison-year roster](comparison_year_roster.csv) and [model eligibility](model_eligibility.csv) show exactly which data enter each stage. [Visit biomass support](visit_biomass_support.csv), [donor attempts](donor_attempts.csv), [size pools](size_pool_inventory.csv) and [borrowing variance](borrowing_variance.csv) show how sizes are used.

[Joint models](joint_model_results.csv) preserve documented shared-reference process dependence as well as shared measurements. [Pooled results](pooled_results.csv) keep before–after changes, recovery slopes and post-only trends separate, and add MPA-level heterogeneity. Across-program summaries include source fixed effects and average the represented program means equally. [Coverage and sparse-size sensitivities](pooled_coverage_sensitivities.csv) use the same fitting procedure. These are conditional estimates, not total uncertainty.

[The paper-input comparison](paper_input_comparison.csv) compares annual responses directly with the frozen paper inputs. [The decomposition](input_change_decomposition.csv) uses one simple unweighted trend on old/all, old/matched, new/matched and new/all years to isolate data and coverage changes. Those diagnostic slopes are not the joint model estimates and do not reproduce the paper's former mixture of estimands.

## Interpretation limits

Of {len(mainj)} main joint model groups, {int(mainj.boundary.sum())} estimate a variance component at zero and {int(mainj.ar_boundary.sum())} reach the AR search boundary. These flags and the dependence sensitivities must accompany interpretation; they are not reasons to exclude a site. The [uncertainty ledger](uncertainty_limits.csv) distinguishes estimable variation from missing information.

The KFM sheephead restricted-period analysis remains a separately checked, 52-comparison-year secondary branch. KFM lobster contributes density; missing lengths have not been invented. Historical site-pool biomass is reproduced in the isolated source build but is not automatically adopted as supported biomass here.

## Rebuild

Run `python3 code/integrated_review/run.py` from the analysis repository. `--analysis-only` skips the inherited R reconstruction but still rebuilds and tests the registered source layer. The final completion record records code and policy hashes. Individual measurements remain in ignored private/staging directories; the source manifest retains access conditions.

Review these results and diagnostics before regenerating pathway/resilience analyses, final paper tables/figures, or manuscript text. Those paper outputs are not certified by this review.
'''
(OUT/'README.md').write_text(text)
links=['coverage_summary','station_assignments','comparison_year_roster','model_eligibility','visit_biomass_support','size_pool_inventory','donor_attempts','borrowing_variance','joint_model_results','joint_model_diagnostics','pooled_results','pooled_coverage_sensitivities','input_change_decomposition','uncertainty_limits','validation_checks']
body='<h1>Integrated pipeline review</h1><p>November 3, 2026 · Review only. Paper and frozen results unchanged.</p><p>'+html.escape(f'{len(an):,} comparison-years; {int(an.response_supported.sum()):,} supported responses; {checks["fitted_histories"]} fitted histories.')+'</p><p><a href="README.md">Methods, limits and rebuild instructions</a></p><h2>Coverage</h2>'+summary.to_html(index=False)+'<h2>Post-protection trends</h2><img src="recovery.png" alt="Pooled post-protection trends"><p>Conditional joint estimates and 95% intervals. Positive values mean the protected/reference ratio increased through time; this is not a before–after effect of protection.</p><h2>Before–after changes</h2><img src="before_after.png" alt="Pooled average before-after changes"><p>Flat-baseline model evaluated only across observed post-protection years. The pre-trend sensitivity is separate. Neither panel establishes an unconfounded causal effect.</p><h2>Tables</h2><ul>'+''.join('<li><a href="'+n+'.csv">'+html.escape(n.replace('_',' '))+'</a></li>' for n in links)+'</ul><h2>What remains uncertain</h2>'+pd.DataFrame(limitations,columns=['Component','Included','Limit']).to_html(index=False)
(OUT/'index.html').write_text('<!doctype html><html><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Integrated pipeline review</title><style>body{font:16px/1.5 system-ui;max-width:1100px;margin:32px auto;padding:16px}table{border-collapse:collapse;width:100%;font-size:14px}td,th{padding:8px;border-bottom:1px solid #ccc;text-align:left}img{max-width:100%}a{color:#146574}</style>'+body+'</html>')
write_inputs();print('Report and figures written')
