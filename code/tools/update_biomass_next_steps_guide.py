"""Documentation overlays for verified next-step review outputs, no primary promotion."""
from pathlib import Path
import pandas as pd
O=Path('outputs/dataset_inventory_20261101');N=Path('outputs/biomass_next_steps_20261102')
def apply(views):
    if not (N/'completed.json').exists():return views
    baseline=pd.read_csv(N/'density/baseline_coverage_diagnostic.csv')
    holds=pd.read_csv(N/'kfm_year_holds.csv')
    tol=pd.read_csv(N/'borrowing/summary.csv')
    attempts=pd.read_csv(N/'borrowing/donor_attempts.csv')
    support=pd.read_csv(N/'borrowing/visit_support.csv')
    strict=support[support.tolerance_days.eq(0)]
    rules=pd.read_csv(O/'current_rules.csv')
    rules.loc[rules.rule_id.eq('R03'),'limitation_or_supersedes']='Density baseline rebuilt independently: Vicente has four full-zone prior years but three shared INNER/MID prior years. Physical depth equivalence remains unverified.'
    rules.loc[rules.rule_id.eq('R09'),'implementation_status']='Both donor attempts and stable keys implemented in lossless review ledger; strict reference and tolerance scenarios retained'
    rules.loc[rules.rule_id.eq('R11'),'implementation_status']='Uniform 0/3/7/14-day upper-gap margins tested; strict44/63, seven/fourteen45/63; no final primary cutoff'
    rules.loc[rules.rule_id.eq('R13'),'implementation_status']='Shared count-unit registry, latent size diagnostic and separate direct-only observation/process working model implemented; total uncertainty remains incomplete'
    rules.loc[rules.rule_id.eq('R13'),'limitation_or_supersedes']='Latent finite-year slope draws are not process CIs. Direct-only process model uses unsmoothed estimates; no second likelihood on smoothed posterior covariance. Allometry and pooled sampling variance fixed.'
    rules.loc[rules.rule_id.isin(['R03','R04','R07','R08','R09','R10','R11','R12','R13']),'evidence_path']='outputs/biomass_next_steps_20261102/README.md'
    rules.to_csv(O/'current_rules.csv',index=False)
    elig=pd.read_csv(O/'comparison_decisions.csv')
    lookup=baseline.set_index('series_id')
    elig['density_full_prior_years']=elig.series_id.map(lookup.full_density_prior_years)
    elig['density_shared_prior_years']=elig.series_id.map(lookup.shared_inner_mid_prior_years)
    elig.loc[elig.series_id.isin(lookup.index),'roster_basis']='Density-specific roster rebuilt; full-zone and shared INNER/MID variants explicit; all63post ratios unchanged'
    elig.to_csv(O/'comparison_decisions.csv',index=False)
    status=pd.read_csv(O/'audit_status_by_comparison.csv')
    status.loc[status.series_id.isin(lookup.index),'time_coverage']='Own density baseline restored; all63post years unchanged. Vicente four full-zone versus three shared-zone prior years; coverage is not model eligibility.'
    status.to_csv(O/'audit_status_by_comparison.csv',index=False)
    ms=pd.read_csv(O/'measurement_support_current.csv')
    for sid,g in holds.groupby('series_id'):
        mask=ms.series_id.eq(sid)
        ms.loc[mask,'alias_unresolved_post_years']=';'.join(map(str,g[g.biomass_status.eq('original physical-site sizes required')].year))
        ms.loc[mask,'alias_historical_replacement_years']=';'.join(map(str,g[g.biomass_status.str.startswith('historical')].year))
    seven=support[support.tolerance_days.eq(7)]
    cov=pd.read_csv(N/'borrowing/annual.csv')
    if 'supported' in cov:
        seven_count=cov[(cov.tolerance_days==7)&(cov.t>0)&cov.supported].groupby('series_id').size()
        ms['seven_day_sensitivity_supported_post_years']=ms.series_id.map(seven_count)
    ms['next_steps_evidence']=ms.series_id.map(lambda s:'../biomass_next_steps_20261102/README.md' if s in set(support.series_id)|set(holds.series_id) else '')
    ms.to_csv(O/'measurement_support_current.csv',index=False)
    # Preserve original strict scientific ledger; richer ledger is linked separately.
    strict.to_csv(O/'measurement_visit_support_detailed.csv',index=False)
    attempts.to_csv(O/'donor_attempts_current.csv',index=False)
    questions=pd.read_csv(O/'unresolved_issues.csv')
    questions.loc[questions.issue_id.eq('U08'),'finding']='KFM exact-year holds documented; no new Scott reply in Nov2 exact-sender search. Lossless donor ledger and0/3/7/14day sensitivities implemented; strict44/63versus45/63at7/14days. Conditional size/count diagnostics improved; fully joint borrowed-biomass process inference and allometric uncertainty remain open.'
    questions.loc[questions.issue_id.eq('U08'),'next_action']='Receive original KFM sizes; review tolerance and hierarchical assumptions. No protected borrowing validation or total-uncertainty claim; do not combine latent posterior covariance with another process likelihood.'
    questions.loc[questions.issue_id.eq('U06'),'finding']='Six-comparison density baseline rebuilt from its own roster. Three Vicente histories have four full-zone prior years but only three shared-zone years; later recovery ratios unchanged. Broader depth/reference policy remains program-specific.'
    questions.to_csv(O/'unresolved_issues.csv',index=False)
    freeze=pd.read_csv(O/'freeze_readiness.csv')
    freeze.loc[freeze.item_id.eq('F3'),'latest_status']='Density baseline and donor attempts rebuilt; timing sensitivities run; KFM sources and fully joint borrowed-biomass process uncertainty still pending'
    freeze.loc[freeze.item_id.eq('F4'),'latest_status']='Shared count registry and conditional size model implemented; separate direct-only process diagnostic checked; total uncertainty not complete'
    freeze.to_csv(O/'freeze_readiness.csv',index=False)
    views['Freeze readiness']=freeze
    views.update({'Current rules':rules,'Eligibility decisions':elig,'Measurement support':ms,'Biomass support':ms,'Audit status':status,'Questions':questions,'Density baseline':baseline,'Current density annual':pd.read_csv(N/'current_density_annual.csv'),'Donor attempts':attempts,'Detailed visit support':strict,'Timing tolerance':tol,'KFM year holds':holds})
    for label,file in [('Conditional size models','uncertainty/size_model_diagnostics.csv'),('Observed-year slope draws','uncertainty/observed_year_slope_draws.csv'),('Direct process admission','uncertainty/direct_process_admission.csv'),('Direct process sensitivity','uncertainty/direct_trend_scenarios.csv')]:views[label]=pd.read_csv(N/file)
    guide=O/'README.md';text=guide.read_text();text=text.replace('Rebuild density-specific baseline coverage before before–after integration.','The next-step rebuild restores those density records; full-zone versus shared-zone baseline availability is now explicit. Before–after model eligibility remains separate.')
    text+='\n## Implemented follow-up\n\n[Next-step implementation and results](../biomass_next_steps_20261102/README.md): restored density baseline; lossless donor attempts; uniform timing tolerances; exact KFM year holds; conditional size/count models and a separate direct-measurement process diagnostic. Strict borrowing remains44/63;7/14day sensitivity45/63. No new primary policy or manuscript changes.\n'
    guide.write_text(text)
    views['Biomass assumptions']=pd.read_csv('docs/BIOMASS_ASSUMPTIONS_AND_CHANGES.csv')
    return views
