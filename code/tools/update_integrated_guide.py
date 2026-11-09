"""Final inventory overlay. Numerical review results never edit frozen paper data."""
from pathlib import Path
import pandas as pd,json
O=Path('outputs/dataset_inventory_20261101');I=Path('outputs/integrated_review_20261103')
def apply(views):
 if not (I/'completed.json').exists():return views
 x=pd.read_parquet(I/'private/supported_measurement_uses.parquet');x=x[x.imputed_density.gt(0)]
 x[['series_id','Source','site','year','date','zone','unit','source_rows','density','biomass','imputed_density','imputed_pool']].drop_duplicates().to_csv(I/'lobster_imputation_uses.csv',index=False)
 a=pd.read_csv(I/'annual_values.csv');b=pd.read_csv(I/'sensitivities/exact_KFM_dates/annual_values.csv')
 z=a.merge(b,on=['series_id','year'],suffixes=('_main','_exact'));z=z[z.response_supported_main & ~z.response_supported_exact]
 z[['series_id','year','Source_main','Taxa_main','Resp_main','ln_ratio_main']].to_csv(I/'seven_day_recovered_years.csv',index=False)
 views['Seven-day recovered years']=z[['series_id','year','Source_main','Taxa_main','Resp_main']].fillna('')
 seeds=pd.read_csv(I/'pooled_coverage_sensitivities.csv');keys=['Taxa','Resp','model','term','Source']
 mc=seeds[seeds.coverage_variant.eq('available')].merge(seeds[seeds.coverage_variant.eq('second_seed')],on=keys,suffixes=('_main','_second'))
 mc['point_change_in_main_SE']=(mc.estimate_second-mc.estimate_main).abs()/mc.se_main
 mc['interval_zero_class_changed']=((mc.low_main>0)|(mc.high_main<0))!=((mc.low_second>0)|(mc.high_second<0))
 mc.to_csv(I/'monte_carlo_check.csv',index=False);views['Integrated seed check']=mc.fillna('')
 rules=views['Current rules'].copy()
 updates={
 'R02':('Three matched strictly earlier calendar years; designation year after. At least three observed post years; all supported observed years, no nine-year cap or consecutive-year requirement. Four-before-year sensitivity. Flat-baseline level/recovery model with pre-trend sensitivity; otherwise post-protection trend only.','Frozen for integrated review; no result-direction selection'),
 'R03':('Use shared program depth zones within each comparison/year; constant-zone sensitivity. Do not equate a shared label with demonstrated physical-depth or habitat equivalence.','Integrated review and fitted coverage sensitivities'),
 'R04':('Retain valid survey dates in all months; equal dated visits within zones, equal zones within stations and available documented stations within roles.','Integrated source-unit rebuild'),
 'R06':('Use available documented stations each year; preserve response-specific exclusions. Constant-station sensitivity and explicit shared-reference models.','Integrated review and fitted coverage sensitivities'),
 'R07':('Unknown counts and positive-count missing sizes remain unknown; preserve density independently. Same-site/date mean-mass lobster imputation retains shared donors.','Integrated source-unit and biomass rebuild'),
 'R08':('Apply exclusions by verified physical site to counts, size donors and calibration; keep KFM benthic exclusions distinct from fish exclusions.','Corrected source build completed in isolated folder'),
 'R09':('Use sizes from the count visit first. KFM urchin/kelp: nearest same-site/year sizes within seven days; exact-date sensitivity. Allow calibrated PISCO reference borrowing only within the six reviewed histories and declared donor rules. Keep density where biomass lacks size support.','Integrated donor ledger and direct-size-only refit; unsupported biomass stays missing'),
 'R11':('For reviewed PISCO borrowing, require the donor–count date gap to be no larger than the largest gap observed in that comparison’s calibration data, with zero extra days of tolerance and at least two calibration site-year groups. This is not proof of strong validation.','Integrated calibration; sparse-size sensitivity'),
 'R12':('Require support for every retained surveyed station/visit in each role-year. Do not remove an unsupported station to rescue biomass; density is separate.','Integrated yearly support and eligibility tables'),
 'R14':('Review cohort and models are frozen in configuration; paper adoption remains on hold pending review and downstream pathway/resilience/figure regeneration.','Integrated review completed; no paper promotion'),
 'R13':('Shared physical count units, frequency-aware size draws, imputation donors and calibration errors feed joint calendar-AR models with shared-reference process terms. Report unestimable uncertainty separately.','Integrated conditional uncertainty; boundary and sensitivity diagnostics retained')}
 for rid,(rule,status) in updates.items():
  ix=rules.rule_id.eq(rid);rules.loc[ix,'current_rule']=rule;rules.loc[ix,'implementation_status']=status;rules.loc[ix,'evidence_path']=str(I/'README.md');rules.loc[ix,'user_decision_status']='Approved for review; paper hold retained'
 rules.loc[rules.rule_id.eq('R13'),'limitation_or_supersedes']='Intervals are conditional on recorded sizes and fixed conversion equations. LTER within-group size variation, allometric coefficient covariance, detection and missing habitat information remain unmeasured. The 20% and 30% group-size scenarios are assumptions, not confidence bounds.'
 rules.loc[rules.rule_id.eq('R07'),'limitation_or_supersedes']='Unknown records remain unknown. Source reconstruction and recovery review completed; downstream paper analyses still require regeneration.'
 rules.loc[rules.rule_id.eq('R12'),'limitation_or_supersedes']='Earlier 44/63 totals apply only to the six-series borrowing audit; use the integrated annual roster for current coverage.'
 rules.loc[rules.rule_id.eq('R14'),'limitation_or_supersedes']='Review configuration frozen; final paper release not approved or certified.'
 scopes={'R03':'Integrated diver comparisons with recorded comparable zones; physical-depth limitations retained','R04':'Integrated diver survey records; program-specific nested unit weights','R07':'Retained density and biomass responses','R08':'All documented response-specific exclusions, including count, donor and calibration records','R09':'Integrated size-based biomass; KFM urchin/kelp short-gap rule and restricted PISCO calibrated borrowing','R12':'All integrated biomass role-years; density assessed independently','R13':'All integrated response histories; shared measurement and process uncertainty'}
 for key,value in scopes.items():rules.loc[rules.rule_id.eq(key),'scope']=value
 rules.loc[rules.rule_id.eq('R10'),'implementation_status']='Incorporated into integrated biomass within the six-history target/calibration scope'
 rules.loc[rules.rule_id.eq('R10'),'evidence_path']=str(I/'README.md')
 rules.loc[rules.rule_id.eq('R03'),'limitation_or_supersedes']='Shared recorded zones do not establish physical-depth or habitat equivalence; earlier six-history zone counts are historical.'
 rules.loc[rules.rule_id.eq('R11'),'limitation_or_supersedes']='Two calibration groups is a computational minimum, not strong validation. No extra PISCO tolerance beyond the strict calibration maximum is adopted; the KFM seven-day rule is separate.'
 review=Path('outputs/biomass_policy_review_20261103')
 assert (review/'verification.json').exists()
 extra=dict(rule_id='R15',domain='Biomass reporting',current_rule='Choose direct-size and borrowing eligibility before examining effects. Report the supported biomass years and direct-size-only sensitivity alongside any result that uses borrowed sizes; preserve density when biomass is unsupported.',scope='Current integrated biomass review; six calibrated PISCO borrowing histories',implementation_status='Coverage and direct-size-only joint refit completed for review',limitation_or_supersedes='Borrowing supports 20 PISCO comparison-years. Removing it changes fitted coverage; retained annual values and measurement draws do not change. This does not establish donor representativeness.',evidence='Biomass policy review',evidence_path=str(review/'README.md'),user_decision_status='Approved for review; paper hold retained',paper_status='Review documentation; no paper promotion')
 rules=pd.concat([rules[~rules.rule_id.eq('R15')],pd.DataFrame([extra])],ignore_index=True)
 rules.to_csv(O/'current_rules.csv',index=False);views['Current rules']=rules
 el=pd.read_csv(I/'model_eligibility.csv');dec=pd.read_csv(O/'comparison_decisions.csv');extra=['before_years','after_years','eligible_model','eligible_four_before','years','missing_response_years']
 for col in extra:dec['integrated_'+col]=dec.series_id.map(el.set_index('series_id')[col])
 dec.loc[~dec.exclude,'model_decision']=dec.loc[~dec.exclude,'integrated_eligible_model'].fillna('No supported paired history in integrated review; inspect station/zone roster')
 dec['measurement_status_scope']='Earlier scoped audit field; current support is in integrated annual roster and model eligibility'
 dec['integrated_supported_years']=dec.series_id.map(a.groupby('series_id').response_supported.sum())
 dec.to_csv(O/'comparison_decisions.csv',index=False)
 for name in ['Eligibility decisions']:
  if name in views:views[name]=dec.fillna('')
 for label,name in [('Integrated annual roster','comparison_year_roster'),('Integrated model eligibility','model_eligibility'),('Integrated biomass visits','visit_biomass_support'),('Integrated donor decisions','donor_attempts'),('Integrated uncertainty limits','uncertainty_limits'),('Integrated sensitivity results','pooled_coverage_sensitivities'),('Integrated imputed lengths','lobster_imputation_uses')]:
  t=pd.read_csv(I/(name+'.csv'));t.to_csv(O/('integrated_'+name+'.csv'),index=False);views[label]=t.fillna('')
 views['Earlier coverage summary']=views['Summary'];views['Summary']=pd.read_csv(I/'coverage_summary.csv').fillna('')
 # Earlier scoped evidence remains accessible, explicitly labeled as such.
 for name in list(views):
  if name in ['Measurement support','Biomass support','Detailed visit decisions','Detailed visit support','Visit support','Donor attempts','Current density annual','Density baseline','Decisions','Station checks','Biomass refits','Count-size matching','Joint biomass draws','Biomass draw holds','All-month rerun','All-month timing','Survey date matching','Size donor validation','Borrowing uncertainty','Borrowing validation','Joint recovery comparison','Joint recovery scenarios','Joint model diagnostics','Joint model calibration','Nonlinear comparison','Nonlinear calibration','Conditional size models','Observed-year slope draws','Direct process admission','Direct process sensitivity']:views[name+' (earlier scoped audit)']=views.pop(name)
 q=views['Questions'].copy()
 updates_q={
 'U06':('Frozen available-station/shared-zone policy has been fitted, with constant-station, constant-zone and calendar-AR/shared-reference sensitivities. Shared labels still do not certify physical depth equivalence.','Review integrated coverage consequences; remaining physical-depth metadata remains U02.'),
 'U08':('Complete registered-source reconstruction and integrated visit biomass/measurement uncertainty review have run. Shared pools, calibration, imputed lengths, count sampling and temporal/reference dependence are represented conditionally.','Review integrated results, sparse-size sensitivity and uncertainty limits before paper promotion; total uncertainty is not fully estimable.'),
 'U13':('Dates are reconstructed only when original year/month/day fields are complete. The integrated donor ledger distinguishes exact, approved short-gap, calibrated and unsupported matching.','Retain incomplete-date size records as source evidence; do not claim visit support without a justified complete-date match.'),
 'U17':('User froze three matched strictly earlier calendar years, designation year after, four-year sensitivity, flat baseline with recovery and pre-trend sensitivity. This is an explicit revision of the historical convention.','Compare fitted sensitivities and disclose the revised rule when the review is adopted; no result-direction exclusions.')}
 for issue,(finding,action) in updates_q.items():
  ix=q.issue_id.eq(issue);q.loc[ix,'finding']=finding;q.loc[ix,'next_action']=action;q.loc[ix,'evidence']=str(I/'README.md')
 q.to_csv(O/'unresolved_issues.csv',index=False);views['Questions']=q
 fr=pd.read_csv(O/'freeze_readiness.csv')
 status={
 'F1':'Frozen for review: available stations/year, shared recorded zones, all valid months, constant-station/zone sensitivities. Habitat equivalence remains a separate limit.',
 'F2':'Frozen: calendar years; designation year after; three matched before years; four-year sensitivity; flat baseline main, pre-trend sensitivity; post-only trends otherwise.',
 'F3':'Source reconstruction and integrated visit biomass/support ledgers completed. Unrecoverable sizes remain missing; density is assessed separately.',
 'F4':'Conditional count/size/shared-donor/reference and calendar-process uncertainty fitted; sparse-size and seed sensitivities checked. Unestimable components remain disclosed.',
 'F5':'Integrated recovery/before-after and pooled fits completed with coverage and uncertainty sensitivities; diagnostics require scientific review before promotion.',
 'F6':'Recovery pooling completed. Downstream pathway/resilience and final paper figure/table regeneration remain on hold; paper unchanged.',
 'F7':'Raw survey-unit reconstruction uses actual dates/areas, complete PISCO fish levels and explicit nested weights; source and missingness checks passed.'}
 for key,value in status.items():fr.loc[fr.item_id.eq(key),'latest_status']=value
 fr.to_csv(O/'freeze_readiness.csv',index=False)
 for label in list(views):
  if 'Freeze' in label:views[label]=fr.fillna('')
 gates=pd.read_csv(O/'inclusion_gates.csv');gates.loc[gates.gate_id.eq('G03'),'scope_limit']='Available stations and shared recorded zones adopted for review; constant coverage tested; physical habitat/depth equivalence remains unverified.';gates.loc[gates.gate_id.eq('G04'),'scope_limit']='Separate restricted KFM sheephead review models zero-count uncertainty with size support; integrated observed-zero point biomass is retained.';gates.loc[gates.gate_id.eq('G06'),'rule']='Three matched strictly earlier calendar years; designation year after; flat-baseline recovery main, pre-trend and four-before sensitivities; at least three observed post years; use all supported observed years without a nine-year cap or consecutive-year requirement.'
 gates.to_csv(O/'inclusion_gates.csv',index=False);views['Inclusion gates']=gates
 extra=pd.read_csv(O/'additional_analysis_checks.csv')
 for key,value in {'Q01':'Integrated physical-visit, method, area, completeness and missingness reconstruction completed.','Q02':'Nested unit/visit/zone/station weights explicit; available-station/shared-zone main and constant-coverage fits completed.','Q03':'Visit-level nonlinear/frequency-aware biomass rebuilt; shared donors and imputation propagated; group-mean and allometric limits remain.'}.items():extra.loc[extra.check_id.eq(key),'latest_status']=value
 extra.to_csv(O/'additional_analysis_checks.csv',index=False)
 for label in list(views):
  if 'Additional' in label:views[label]=extra.fillna('')
 wf=views.get('Workflow status',pd.read_csv(O/'workflow_status.csv')).copy()
 # Preserve existing schema and add the unified entry point without fabricating prerequisites.
 row=dict(workflow_id='W08',workflow='Integrated source-to-model review',status='Executed and checked; conditional review',scope='Frozen candidate scope; supported density and biomass evaluated separately',working_directory_relative_to_analysis_repo='.',commands='python3 code/integrated_review/run.py',outputs=str(I),adoption_limit='No paper promotion; uncertainty limits and diagnostics remain explicit',evidence=str(I/'README.md'),input_record=str(I/'input_hashes.json'),prerequisites='Registered local sources, companion data repository, existing Python/R environments')
 wf.loc[wf.workflow_id.eq('W02'),'status']='Corrected reconstruction completed in isolated folder'
 wf.loc[wf.workflow_id.eq('W02'),'adoption_limit']='Historical broad pooling reproduced for comparison; not the adopted integrated biomass estimator'
 wf.loc[wf.workflow_id.eq('W03'),'status']='Earlier scoped audit; matching/calibration scope incorporated into integrated review'
 wf.loc[wf.workflow_id.eq('W03'),'adoption_limit']='Historical six-history results; use W08 and integrated tables for current fitted values and support'
 wf=pd.concat([wf,pd.DataFrame([row])],ignore_index=True).drop_duplicates();wf.to_csv(O/'workflow_status.csv',index=False);views['Workflow status']=wf
 (O/'README.md').write_text('''# Current decision tables

Updated November 4, 2026. Start with the [analysis guide](../../docs/ANALYSIS_GUIDE.md) and [integrated review](../integrated_review_20261103/index.html).

The coordinated source-to-model review has run. Main coverage, baseline and matching rules are recorded in configuration; station, visit, donor and model eligibility are separate linked tables. The manuscript and frozen results are unchanged. Conditional uncertainty and model diagnostics still need review before paper promotion.

- [Plain-language findings](../integrated_review_20261103/REVIEW_NOTES.md)
- [Current analysis page](index.html) with [site](../methods_results_reference_20261104/site_decisions.csv) and [response](../methods_results_reference_20261104/response_histories.csv) decisions
- [Current decisions, reasons and supporting records](index.html)
- [Decisions and reasons (CSV)](current_decisions_and_reasons.csv)
- [Current rules](current_rules.csv) and [comparison decisions](comparison_decisions.csv)
- [Integrated comparison-year roster](integrated_comparison_year_roster.csv)
- [Integrated model eligibility](integrated_model_eligibility.csv)
- [Integrated biomass support](integrated_visit_biomass_support.csv) and [donor decisions](integrated_donor_attempts.csv)
- [Current biomass coverage and borrowing decision](../biomass_policy_review_20261103/README.md)
- [Fitted sensitivities](integrated_pooled_coverage_sensitivities.csv) and [uncertainty limits](integrated_uncertainty_limits.csv)
- [Workflow commands](workflow_status.csv)

## Earlier scoped audits

The earlier measurement-support, visit-support and donor-attempt files retain their original names for reproducibility. They are omitted from the current-decision viewer. Historical counts such as 44/63 supported years and 1,086 visit uses are not totals for the integrated run. The separately labeled 52-year KFM sheephead review remains supporting work. Current support comes from the integrated links above.

Earlier dated audit tables are retained as historical evidence, not as the latest fitted roster. The unchanged 206-row candidate inventory is not the count of fitted comparisons.
''')
 from summarize_integrated_review import write_summary
 write_summary()
 from apply_kfm_reference_review import apply as apply_reference_review
 return apply_reference_review(views)
