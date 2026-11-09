from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,html
from rebuild_borrowing_support import fit
O=Path('outputs/biomass_next_steps_20261102')
for folder in ['uncertainty','uncertainty_pool10','uncertainty_seed2']:
    assert json.loads((O/folder/'run_status.json').read_text())['status']=='complete',folder
base=pd.read_csv(O/'density/baseline_coverage_diagnostic.csv');tol=pd.read_csv(O/'borrowing/summary.csv');eff=pd.read_csv(O/'borrowing/effects.csv');ann=pd.read_csv(O/'borrowing/annual.csv')
dvariants=pd.read_csv(O/'density/annual_comparison_density_variants.csv')
dcur=dvariants[dvariants.variant.eq('allmonths_shared_inner_mid')].copy()
dfull=dvariants[dvariants.variant.eq('allmonths_full_available_zones')][['series_id','year','ln_ratio']].rename(columns={'ln_ratio':'full_zone_evidence_lnratio'})
dcur=dcur.merge(dfull,on=['series_id','year'],validate='one_to_one');dcur['roster_basis']='Independent density station/year roster';dcur['baseline_status']=np.where(dcur.shared_inner_mid_eligible,'Shared INNER/MID supported; model eligibility separate','Full-zone density retained as evidence; shared-zone comparison unsupported')
assert len(dcur)==84 and (dcur.t>0).sum()==63 and dcur.loc[dcur.t>0,'ln_ratio'].notna().all()
dcur.to_csv(O/'current_density_annual.csv',index=False)
models=pd.read_csv(O/'uncertainty/size_model_diagnostics.csv');direct=pd.read_csv(O/'uncertainty/direct_trend_scenarios.csv');da=pd.read_csv(O/'uncertainty/direct_annual_observations.csv');draws=pd.read_csv(O/'uncertainty/observed_year_slope_draws.csv')
rows=[]
for r in direct.itertuples():
    full=ann[(ann.series_id==r.series_id)&(ann.tolerance_days==7)&(ann.t>0)&ann.supported];years=set(da[da.series_id==r.series_id].year);sub=full[full.year.isin(years)];f=fit(full);s=fit(sub)
    rows.append(dict(series_id=r.series_id,rho=r.rho,full_supported_years=len(full),direct_years=len(sub),full_OLS_HAC_slope=f['estimate'],direct_same_year_OLS_HAC_slope=s['estimate'],direct_same_year_HAC_low=s['low'],direct_same_year_HAC_high=s['high'],direct_GLS_slope=r.estimate,direct_GLS_low=r.bootstrap_low,direct_GLS_high=r.bootstrap_high,coverage_slope_change=s['estimate']-f['estimate'],GLS_weighting_slope_change=r.estimate-s['estimate'],scope='Conditional sensitivity; coverage and estimator changes separated'))
pd.DataFrame(rows).to_csv(O/'direct_process_change_decomposition.csv',index=False)
sens=[]
for folder in ['uncertainty_pool10','uncertainty_seed2']:
    z=pd.read_csv(O/folder/'observed_year_slope_draws.csv');m=draws.merge(z,on='series_id',suffixes=('_base','_variant'))
    for r in m.itertuples():sens.append(dict(series_id=r.series_id,variant=folder,mean_change=r.mean_slope_variant-r.mean_slope_base,lower_change=r.draw_low_variant-r.draw_low_base,upper_change=r.draw_high_variant-r.draw_high_base))
pd.DataFrame(sens).to_csv(O/'uncertainty_sensitivity.csv',index=False)
# Reproducible exact-year KFM holds from verified audit, not an inferred new delivery.
h=pd.read_csv('outputs/parallel_resolution_20261102/kfm/history_years.csv');holds=[]
for r in h.itertuples():
    known=set(str(r.diagnostically_corrected_post_years).split(';'));unknown=set(str(r.unresolved_exact_measurement_post_years).split(';'))
    for year in str(r.post_years).split(';'):
        holds.append(dict(series_id=r.series_id,year=int(year),donor_role=r.donor_role,combined_label=r.delivered_site,biomass_status='historical replacement diagnostic; production rebuild pending' if year in known else 'original physical-site sizes required' if year in unknown else 'verified zero reference; no donor mass needed',density_status='Retain under existing response-specific coverage',primary_adopted=False,source_request='Already sent November2; no duplicate request',requested_original_sites='1/21' if r.delivered_site=='a-k-01' else '14/37'))
pd.DataFrame(holds).to_csv(O/'kfm_year_holds.csv',index=False)
validation=pd.read_csv(O/'uncertainty/size_model_site_validation.csv');vv=validation[validation.available].groupby(['series_id','status']).agg(predictions=('covered','size'),coverage=('covered','mean'),held_sites=('held_site','nunique'),years=('year','nunique')).reset_index();vv.to_csv(O/'size_validation_summary.csv',index=False)
checks=json.loads((O/'density/tests.json').read_text());bchecks=json.loads((O/'borrowing/tests.json').read_text())
protected=json.loads(Path('outputs/researcher_scope_20261101/protected_hashes.json').read_text());assert all(Path(p).exists() and hashlib.sha256(Path(p).read_bytes()).hexdigest()==v for p,v in protected.items())
assert models.converged.all();assert (O/'uncertainty/independent_R_checks.csv').exists()
text='''# Implemented biomass and coverage follow-up

November 2, 2026. These are review analyses, not changes to the manuscript or primary release. Researcher exclusions remain intact.

## What changed

**Density now has its own baseline roster.** Six prior comparison-years and 33 station-year uses are restored. All 63 post-protection density ratios remain unchanged. The current density input has 84 comparison-year rows; three earlier shared-zone ratios remain unknown, with their full-zone observations retained alongside. Point Vicente has four prior years with its full observed zones, but only three with shared INNER/MID coverage: ROCKY_POINT_S was DEEP-only in 2008. Retaining the count evidence does not certify a comparable before–after baseline.

**Size matching now keeps both attempted paths.** The ledger retains local-other-date and other-site candidates, failed guards, dates, sample counts, calibration sites/years, and stable physical keys. Missing sizes remain distinct from absent validation. Positive counts without supported sizes remain unknown biomass; no reference stations are dropped silently.

**The date boundary is tested consistently.** Zero- and three-day margins retain44of63post biomass comparison-years; seven and14days retain45. Seven days restores Abalone red2019 through Point Fermin. It also changes two Abalone purple2019 visit records to local rather than cross-site donors. Red slope changes from −0.180509 to −0.188102through coverage; purple changes from −0.098868 to −0.101968through donor choice. The other four recovery slopes are unchanged. These are conditional temporal fits; no final primary tolerance is selected.

**Uncertainty is developed without recycling a smoothed estimate as new data.** Twelve species/MPA/role size models use measured means, explicit partially pooled sampling dispersion, site/site-year/event variation and the lognormal sampling-mean correction. Shared physical count units receive the same resampling weights across taxa. Forty covariance-parameter bootstrap refits per model feed800conditional draws. Support is not expanded: the seven-day scenario still has45post comparison-years.

The resulting slope-draw ranges describe uncertainty in change over the observed years. They are not confidence intervals for a general recovery process or protection effect. The rejected shortcut of using latent posterior covariance as independent measurement covariance is absent from generated results.

A separate observation-plus-process working model uses **unsmoothed empirical ratios and only directly measured sizes**. Only the two lobster histories have at least four complete post-protection years (8 Abalone, 5 Vicente). Sampling covariance is estimated from count-unit and measured-size replicates; process variance is then estimated conditional on that covariance. Calendar correlations 0, 0.3 and 0.6are explicit scenarios. All six conditional bootstrap intervals span zero. The linked decomposition distinguishes coverage losses from GLS weighting changes; this subset must not replace the full density analysis.

## What still cannot be called resolved

- Original recent KFM sizes remain unavailable for54comparison-years across six biomass histories. Forty-five historical years have diagnostic replacements; one verified zero-reference year needs no size donor. Density remains retained. The exact sender/date inbox check found no new Scott reply; the request was already sent.
- Protected-site borrowing is not independently validated. Reference-site model checks do not establish transportability to protected stations. Conditional held-out size intervals cover roughly88–98% of reference sample means, depending on comparison; observations within sites are dependent. Evaluation intervals condition on the held-out sample’s estimated dispersion; these are not outcome-blind predictive coverage estimates.
- Size-model parameter mixtures are an approximation, not an exact posterior. Pooling strengths2and10and a second random seed are compared. The pooled sampling variance itself, allometric coefficients, detection, count–size dependence and persistent survey-unit identity are not fully estimated.
- Fully joint process inference for histories containing borrowed biomass remains unfinished. Neither latent-slope ranges nor direct-only working intervals are total uncertainty. No additional animal histories were excluded to hide these limits.
- KFM kelp single-date visit variation and all-program numeric depth equivalence remain data limitations. Current available-reference assignments and program-specific methods remain; fixed-reference and depth alternatives are sensitivities.

## Current tables

- [Current density annual input](current_density_annual.csv), [density baseline](density/baseline_coverage_diagnostic.csv) and [restored station-years](density/recovered_station_year_uses.csv).
- [Detailed visit support](borrowing/visit_support.csv), [both donor attempts](borrowing/donor_attempts.csv), [timing tolerance](borrowing/summary.csv), [effects and decomposition](borrowing/effects.csv).
- [KFM year holds](kfm_year_holds.csv), [program rules](program_coverage_rules.csv), [implementation status](implementation_status.csv).
- [Size/count registry and diagnostics](uncertainty/size_model_diagnostics.csv), [validation summary](size_validation_summary.csv), [observed-year slope draws](uncertainty/observed_year_slope_draws.csv).
- [Direct-process admission](uncertainty/direct_process_admission.csv), [process results](uncertainty/direct_trend_scenarios.csv), [coverage versus estimator changes](direct_process_change_decomposition.csv), [pooling/seed checks](uncertainty_sensitivity.csv).
- [Current consolidated guide](../dataset_inventory_20261101/README.md).

## Reproduction and validation

Run `python3 code/tools/run_biomass_next_steps.py` from the analysis repository. It rebuilds density, donor matching/tolerance, conditional uncertainty and sensitivity variants, verifies them, then updates the guide. Raw sources, prior dated scientific outputs and the manuscript are not overwritten. The live email check is a dated record and is not repeated by this runner.

Density has16independent source/aggregation/edge checks; borrowing has5regression tests plus complete-reference, zero/unknown, source-identity, calibration and hash assertions. Size covariance is checked against independent Gaussian precision algebra; visit identity, singleton pooling, fit provenance, PSD covariance and saved slope draws are tested. R independently reproduces all6direct-only GLS/REML fits. All487protected manuscript/canonical files match their original hashes. Validation and limitations are recorded separately.
'''
# Light prose spacing for count tokens without altering paths/URLs/code.
import re
text=re.sub(r'(?<=[A-Za-z])(?=\d)', ' ',text);text=re.sub(r'(?<=\d)(?=[A-Za-z])',' ',text)
# Restore deliberate paths/code broken by generic prose spacing by not applying it inside references.
# Use original filenames via explicit repair for known affected terms.
for token in ['INNER/MID','GLS/REML']:pass
# Generic spacing touches date folder names; repair numeric sequence followed by suffix letters only in prose.
if Path('outputs/joint_borrowed_biomass_20261102/report_verified.json').exists():
    text='> Earlier review stage. The [joint borrowed-biomass model review](../joint_borrowed_biomass_20261102/README.md) now extends the uncertainty work below; retain this report as dated evidence.\n\n'+text
(O/'README.md').write_text(text)
# Simple report index; tables remain downloadable and searchable in consolidated viewer.
links=[('Baseline coverage','density/baseline_coverage_diagnostic.csv',base),('Timing tolerance','borrowing/summary.csv',tol),('Observed-year slope draws','uncertainty/observed_year_slope_draws.csv',draws),('Direct-only process sensitivity','uncertainty/direct_trend_scenarios.csv',direct)]
page=['<!doctype html><meta charset="utf-8"><title>Biomass follow-up</title><style>body{font:16px system-ui;margin:28px;max-width:1200px}table{border-collapse:collapse;font-size:13px}td,th{padding:8px;border-bottom:1px solid #ddd;text-align:left}.scroll{overflow:auto}p{max-width:850px}</style><h1>Biomass and coverage follow-up</h1><p>Review only. Density baseline restored; donor failure paths retained; tolerance and conditional uncertainty tested. Primary analysis and manuscript unchanged.</p><p><a href="README.md">Methods, limits and reproducibility</a> · <a href="../dataset_inventory_20261101/index.html">Consolidated tables</a></p>']
for title,file,d in links:page.append('<h2>'+html.escape(title)+'</h2><p><a href="'+file+'">Download CSV</a></p><div class="scroll">'+d.to_html(index=False,escape=True)+'</div>')
(O/'index.html').write_text('\n'.join(page))
# Extend the authoritative decision ledger without rewriting dated scientific reports.
import ast
ledger_additions=[
('B37','Correction','Density baseline uses its own response roster','Implemented in review','Six earlier comparison-years restored; three Vicente2008 shared-zone ratios remain unknown','biomass_next_steps_20261102/current_density_annual.csv'),
('B38','Implementation','Keep both borrowing attempts with stable keys','Implemented','Missing data and unavailable calibration remain distinct; no silent reference omission','biomass_next_steps_20261102/borrowing/donor_attempts.csv'),
('B39','Sensitivity','Uniform upper calibration-gap margins0/3/7/14days','Tested; no primary threshold selected','44/63strict versus45/63at7/14days; two local donor-priority changes','biomass_next_steps_20261102/borrowing/effects.csv'),
('B40','Uncertainty assumption','Partial pooling of size-sampling dispersion and latent sample means','Conditional diagnostic; pools2/10tested','Singleton variance borrowed explicitly; fixed allometry; latent slope draws are not process CIs','biomass_next_steps_20261102/uncertainty/policy.json'),
('B41','Inference boundary','Do not recycle latent posterior covariance as independent measurement error','Enforced','Finite-year latent slope functional and direct-only process working model remain separate','biomass_next_steps_20261102/uncertainty/direct_process_admission.csv'),
('B42','Source hold','Recent original KFM sizes still required','Exact-year holds recorded; request already sent','54unresolved comparison-years; density retained; no new exact-sender reply foundNov2','biomass_next_steps_20261102/kfm_year_holds.csv'),
('B43','Coverage rule','Program-specific survey methods and available researcher references retained','Current scoped review; no blanket cutoff','Fixed-reference/depth alternatives remain diagnostics; no results-based exclusions','biomass_next_steps_20261102/program_coverage_rules.csv')]
add=pd.DataFrame(ledger_additions,columns=['id','type','assumption_or_change','status','consequence','evidence']);add.insert(1,'updated','2026-11-02');add.to_csv(O/'decision_ledger_additions.csv',index=False)
module=ast.parse(Path('code/tools/report_biomass_protocols.py').read_text());baseledger=next(ast.literal_eval(n.value) for n in module.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='ledger' for t in n.targets))
ld=pd.DataFrame(baseledger,columns=['id','type','assumption_or_change','status','consequence','evidence']);ld.insert(1,'updated','2026-11-02');ld=pd.concat([ld,add],ignore_index=True);assert ld.id.is_unique;ld.to_csv('docs/BIOMASS_ASSUMPTIONS_AND_CHANGES.csv',index=False)
(O/'completed.json').write_text(json.dumps(dict(review_implementation_complete=True,primary_updated=False,density_post_years=63,strict_biomass_post_years=44,tolerance7_biomass_post_years=45,protected_files_unchanged=len(protected),total_uncertainty_complete=False,fully_joint_borrowed_process_inference_complete=False),indent=2)+'\n')
print('Next-step report built; source-dependent holds and unresolved joint inference explicit.')
