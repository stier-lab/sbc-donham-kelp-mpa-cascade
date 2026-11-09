"""Current decision views; overlays documentation, never scientific estimates."""
from pathlib import Path
import hashlib, json
import pandas as pd
O=Path('outputs/dataset_inventory_20261101')
B=Path('outputs/pisco_biomass_borrowing_20261102')

def apply(views):
    comp=pd.read_csv(O/'audit_master.csv')
    status=pd.read_csv(O/'audit_status_by_comparison.csv')
    bio=pd.read_csv('outputs/biomass_support_20261102/biomass_history_status.csv')
    summary=pd.read_csv(B/'summary.csv')
    visits=pd.read_csv(B/'visit_biomass_support.csv')
    aliases=pd.read_csv('outputs/biomass_exclusion_reconciliation_20261102/affected_histories.csv')
    rules=[
    ('R01','Site inclusion','Preserve researcher and ecological exclusions. Do not select sites by slope, significance or effect direction.','All comparisons','Approved; retained in review','Existing exclusion reasons retained','Researcher scope audit'),
    ('R02','Model eligibility','Separate site inclusion from eligibility for before–after estimation or post-protection trends.','All comparisons','Final primary model rules not frozen','Baseline stability tolerances remain diagnostics','pBACIPS eligibility audits'),
    ('R03','Depth coverage','Use comparable shared zones; show losses from constant zones before broader adoption.','General principle; constant INNER/MID implemented for six PISCO animal comparisons','Principle approved; six-comparison review implemented','Does not establish measured-depth equivalence across programs','PISCO all-month review'),
    ('R04','Survey timing','Include valid surveys in all months; average repeated visits equally. No protected/reference date-gap cutoff adopted.','Six PISCO animal comparisons at Abalone Cove and Point Vicente','Approved and implemented in review','Replaces May–October restriction for this review only; seasonal imbalance remains','PISCO all-month and timing audits'),
    ('R05','Fish methods','Retain CRANE, average repeat surveys and omit incomplete passes from the combined fish response.','Audited PISCO fish histories','Approved and implemented in review','Do not treat missing levels as zero counts','PISCO fish rules audit'),
    ('R06','Reference membership','Compare fixed and available station sets; retain established assignments provisionally and test influence where habitat documentation is incomplete.','Retained comparisons','Diagnostics completed; final station policy pending','Explicit researcher exclusions still apply','All-program station refits'),
    ('R07','Missing measurements','Keep valid density when sizes do not support biomass. Missing measurements are unknown, not zero.','All biomass responses','Approved; review safeguards implemented','Production rebuild remains pending for missing-stipe correction','Biomass resolution and source guards'),
    ('R08','Excluded size donors','Enforce exclusions by physical source site, including aliases for KFM sites 21 and 37.','Six affected red/purple urchin biomass histories','Historical correction tested; later source reconciliation pending','Biomass provisional; density retained','Biomass exclusion reconciliation'),
    ('R09','Size matching','Use measured same-visit sizes first; other visits or sites require calibrated prediction error.','Six PISCO animal biomass histories','Implemented as review sensitivity','Supersedes unrestricted nearest-date matching in earlier local-only review','PISCO borrowing review'),
    ('R10','Borrowing boundaries','For the six reviewed PISCO animal-biomass histories, borrow reference sizes only within the same species, MPA, reference role, zone and calendar year. Exclude the target site from cross-site donors and calibration.','Six PISCO animal biomass histories','Implemented as review sensitivity','Only reference-side borrowing has been evaluated. Matching zone and year does not prove comparable animal sizes; calibration residuals pool zones and years within each response history and role.','PISCO borrowing review'),
    ('R11','Borrowing support','Require the date gap to be no larger than the calibration maximum and at least two calibration site-year groups. Leave unsupported observations missing.','Six PISCO animal biomass histories','Implemented computational guards','Two groups is not evidence of adequate validation; no universal 90-day cutoff. Minimum gap is not checked; small upper-bound tolerance remains a proposal.','PISCO borrowing review'),
    ('R12','Complete comparisons','Require biomass support for both roles and the retained reference stations; do not silently omit unsupported stations.','Six PISCO animal biomass histories','Implemented as review sensitivity','44 of 63 post-protection comparison-years supported; 38 positive species-zone-date uses unknown','PISCO borrowing review'),
    ('R13','Uncertainty','Retain shared donor identities and shared prediction errors. Distinguish conditional uncertainty from total uncertainty.','Conditional joint draws and six-history borrowing sensitivity','Implemented components; total uncertainty incomplete','Borrowing draws hold counts, coefficients and offsets fixed; no extra donor bootstrap','Joint biomass and borrowing audits'),
    ('R14','Paper release','Keep review results separate from the manuscript and canonical analysis until reviewed and approved for integration.','All current audits','Active hold','No final primary cohort or model freeze','User review hold')]
    current=pd.DataFrame(rules,columns=['rule_id','domain','current_rule','scope','implementation_status','limitation_or_supersedes','evidence'])
    evidence_paths={
    'R01':'outputs/researcher_scope_20261101',
    'R02':'outputs/model_eligibility_20261101/README.md',
    'R03':'outputs/pisco_all_months_20261102/README.md',
    'R04':'outputs/pisco_timing_donors_20261102/README.md',
    'R05':'outputs/pisco_fish_rules_20261101',
    'R06':'outputs/all_program_station_refits_20261102/README.md',
    'R07':'outputs/biomass_resolution_20261102/README.md',
    'R08':'outputs/biomass_exclusion_reconciliation_20261102/README.md',
    **{f'R{i:02d}':'outputs/pisco_biomass_borrowing_20261102/README.md' for i in range(9,14)},
    'R14':'docs/BIOMASS_ASSUMPTIONS_AND_CHANGES.csv'}
    current['evidence_path']=current.rule_id.map(evidence_paths)
    current['user_decision_status']=current.rule_id.map({
        'R01':'Approved', 'R02':'Final model choice pending', 'R03':'Shared zones approved; broader constant coverage pending',
        'R04':'Approved for review', 'R05':'Approved for review', 'R06':'Comparison and influence tests approved; final policy pending',
        'R07':'Approved', 'R08':'Correction authorized; recent sources pending',
        'R09':'Implementation authorized for review', 'R10':'Implementation authorized for review',
        'R11':'Implementation authorized for review; guard is not a biological adequacy threshold',
        'R12':'Implementation authorized for review', 'R13':'Conditional tests authorized; final uncertainty pending', 'R14':'Active user hold'})
    current['paper_status']='Review documentation; no paper promotion in this consolidation'
    for r in summary.itertuples():
        note=f'Conditional review biomass supports {r.borrowing_post_years} of {r.original_post_years} post-protection years; total uncertainty incomplete'
        comp.loc[comp.series_id.eq(r.series_id),'latest_biomass_support']=note
        status.loc[status.series_id.eq(r.series_id),'biomass_support']=note
        bio.loc[bio.series_id.eq(r.series_id),'support_status']=note
    affected=comp.series_id.isin(aliases.series_id)
    comp.loc[affected,'latest_biomass_support']='Excluded donor aliases affect biomass; original recent sizes pending; density retained'
    status.loc[status.series_id.isin(aliases.series_id),'biomass_support']='Biomass provisional: excluded donor aliases; recent original sizes pending'
    bio.loc[bio.series_id.isin(aliases.series_id),'support_status']='Biomass provisional: excluded donor aliases; recent original sizes pending'
    six=set(summary.series_id)
    twelve=six|{s.rsplit('|',1)[0]+'|Den' for s in six}
    mask=status.series_id.isin(twelve)
    status.loc[mask,'survey_timing']='All valid months included in review; equal visits; seasonal imbalance remains; no date-gap cutoff'
    status.loc[mask,'habitat_depth']='Constant shared INNER/MID zones in six-comparison review; physical depth equivalence not established'
    status.loc[mask,'methods_effort']='Dated visits averaged equally within zone; zones and retained stations weighted equally'
    status.loc[mask,'evidence']='outputs/pisco_biomass_borrowing_20261102/README.md; outputs/pisco_all_months_20261102/README.md'
    for name in ['audit_master','comparison_eligibility']: comp.to_csv(O/(name+'.csv'),index=False)
    status.to_csv(O/'audit_status_by_comparison.csv',index=False)
    decisions=pd.read_csv(O/'user_decisions.csv')
    decisions.loc[decisions.decision.eq('Timing'),'status']='All valid months approved and implemented for six affected PISCO animal comparisons; no global date-gap cutoff'
    decisions.loc[decisions.decision.eq('Biomass'),'status']='Keep density; six-history calibrated borrowing implemented as conditional sensitivity, not primary; unresolved sizes remain missing'
    decisions.to_csv(O/'user_decisions.csv',index=False)
    elig=comp[['series_id','MPA','Source','Taxa','Resp','exclude','reason','in_current_manuscript','review_scope_status','mapped_protected_sites','mapped_reference_sites','reference_evidence','eligibility_interpretation','issue_ids']].copy()
    elig['measurement_status']=comp.latest_biomass_support
    elig['model_decision']='Not frozen; before–after eligibility and recovery trends remain separate'
    elig.loc[comp.exclude,'model_decision']='Excluded from review refits; retain recorded researcher/ecological reason'
    elig['rule_ids']='R01;R02;R06;R14'
    elig.loc[elig.Source.ne('Landsat'),'rule_ids']+=';R03'
    elig.loc[elig.series_id.isin(twelve),'rule_ids']+=';R04'
    elig.loc[elig.Source.eq('PISCO') & elig.Taxa.eq('S. pulcher') & ~elig.exclude,'rule_ids']+=';R05'
    elig.loc[elig.Resp.eq('Bio'),'rule_ids']+=';R07;R13'
    elig.loc[elig.series_id.isin(six),'rule_ids']+=';R09;R10;R11;R12'
    elig.loc[elig.series_id.isin(aliases.series_id),'rule_ids']+=';R08'
    support=bio.merge(summary[['series_id','original_post_years','local_only_post_years','borrowing_post_years','unsupported_years']],on='series_id',how='left',validate='one_to_one')
    support['evidence']='outputs/biomass_support_20261102/biomass_history_status.csv'
    support['coverage_meaning']='Original support screen; blanks do not mean zero or fully resolved'
    support.loc[support.series_id.isin(six),'coverage_meaning']='Post-protection comparison-years; earlier local-only counts used different temporal rules, not a nested treatment'
    support.loc[support.series_id.isin(six),'evidence']='outputs/pisco_biomass_borrowing_20261102/summary.csv; measurement_visit_decisions.csv'
    support.loc[support.series_id.isin(aliases.series_id),'evidence']='outputs/biomass_exclusion_reconciliation_20261102/affected_histories.csv'
    history=pd.read_csv('docs/BIOMASS_ASSUMPTIONS_AND_CHANGES.csv')
    history['current_interpretation']='Dated evidence or decision; read current rules for implementation and scope'
    superseded={'B12':'Six PISCO histories now have explicit visit matching; other programs remain unresolved (R09).','B15':'Partial uncertainty work advanced in B23/B35; total uncertainty still incomplete (R13).','B24':'All 23 zone flags resolved as filtering artifacts in B26; broader representativeness remains open.','B27':'Historical season-restricted diagnostic only; all-month review adopted (R04).','B29':'Nearest-date component superseded by calibrated other-visit support (R09–R11).','B30':'Historical local-only coverage: current conditional sensitivity supports 44 years, under different rules (R12).','B33':'No-imputation statement describes that earlier audit; conditional borrowing subsequently implemented in B35–B36 (R09–R12).'}
    history['lifecycle']='Retained evidence / scoped decision'
    for bid,note in superseded.items():
        history.loc[history.id.eq(bid),'current_interpretation']=note
        history.loc[history.id.eq(bid),'lifecycle']='Earlier stage; see current interpretation'
    follow=Path('outputs/parallel_resolution_20261102')
    if (follow/'coverage/density_scope_differences.csv').exists():
        roster=pd.read_csv(follow/'coverage/density_scope_differences.csv')
        extra=roster[['series_id','omitted_density_years','omitted_prior_station_years','omitted_post_station_years']]
        elig=elig.merge(extra,on='series_id',how='left',validate='one_to_one')
        elig['roster_basis']='Program/response-specific audit; see source evidence'
        elig.loc[elig.series_id.isin(twelve),'roster_basis']='Six-comparison review uses biomass station/year roster; complete density baseline not certified'
        status.loc[status.series_id.isin(twelve),'time_coverage']='All 63 density recovery comparison-years preserved; biomass-derived roster omits six earlier density comparison-years across four histories; see linked roster audit'
        status.to_csv(O/'audit_status_by_comparison.csv',index=False)
    tables={'current_rules':current,'comparison_decisions':elig,'measurement_support_current':support,'measurement_visit_decisions':visits,'decision_history':history}
    for name,frame in tables.items(): frame.to_csv(O/(name+'.csv'),index=False)
    issue=pd.read_csv(O/'unresolved_issues.csv')
    issue.loc[issue.issue_id.eq('U08'),'finding']='Six KFM-donor alias histories remain provisional. Six PISCO animal comparisons retain 63 density years; conditional borrowing supports 44 biomass years, with 38 positive species-zone-date uses unsupported. Total uncertainty remains incomplete.'
    issue.loc[issue.issue_id.eq('U08'),'evidence']='outputs/pisco_biomass_borrowing_20261102/README.md; outputs/biomass_exclusion_reconciliation_20261102/README.md'
    issue.loc[issue.issue_id.eq('U06'),'finding']='Six affected PISCO animal comparisons use all valid months and constant shared INNER/MID zones in review. Seasonal imbalance remains; no universal date-gap cutoff or all-program depth policy has been finalized.'
    issue.to_csv(O/'unresolved_issues.csv',index=False)
    freeze=pd.read_csv(O/'freeze_readiness.csv')
    freeze.loc[freeze.item_id.eq('F1'),'latest_status']='Six PISCO animal comparisons use all months and shared INNER/MID zones; broader depth/station policy pending'
    freeze.loc[freeze.item_id.eq('F3'),'latest_status']='Conditional borrowing implemented for six histories; six donor-alias histories and total uncertainty unresolved'
    freeze.to_csv(O/'freeze_readiness.csv',index=False)
    views.update({'Audit status':status,'Decisions':decisions,'Questions':issue,'Freeze readiness':freeze,'Biomass support':support})
    front={'Current rules':current,'Eligibility decisions':elig,'Measurement support':support,'Visit support':visits,'Decision history':history}
    (O/'README.md').write_text('''# Current decision and dataset guide

Updated November 2, 2026. This is the current entry point. These tables describe review analyses; they do not freeze or replace the manuscript analysis.

- [Current rules](current_rules.csv): scope, implementation, limits and superseded choices.
- [Comparison decisions](comparison_decisions.csv): all 206 candidates, with 175 retained review histories and 31 researcher/ecological exclusions. Site inclusion and model eligibility are separate. The paper-membership column describes the existing paper, not a new selection.
- [Measurement support](measurement_support_current.csv): all 94 biomass histories, including unresolved donor identities and the latest six-comparison borrowing review.
- [Visit-level support](measurement_visit_decisions.csv): 1,086 PISCO visit/zone/species uses; measured sizes, donor identifiers, matching gaps, calibration support and unsupported observations. These are not independent samples or comparison-years.
- [Decision history](decision_history.csv): all 43 dated entries preserved, with earlier stages explicitly identified.
- [Open questions](unresolved_issues.csv) and [freeze readiness](freeze_readiness.csv): work still required before the primary analysis can be frozen.
- [Interactive tables](index.html) and [full master table](audit_master.csv).

For the six PISCO animal comparisons, the current review retains 63 post-protection density comparison-years. Conditional borrowing supports 44 biomass comparison-years. The earlier 25 local-only years used different temporal-support rules, so 44 versus 25 is not a pure borrowing effect. Of 1,086 visit uses, 624 have same-visit measurements, 131 borrow from another site, 8 from another visit, 285 have zero counts and 38 remain unsupported. Those 38 uses represent 24 site visits; 26 uses fall after protection and block 19 comparison-years. Borrowing is reference-only; protected-site prediction has not been validated. Missing biomass is never converted to zero, and unsupported reference stations are not silently dropped.

Six other biomass histories (four KFM and two PISCO histories using KFM sizes) remain provisional because intentionally excluded KFM donors entered under alternate labels. Density remains usable. Historical corrections are diagnostic; recent original sizes are still needed. Conditional uncertainty checks do not yet represent total uncertainty.

## Reproduce and verify

From the analysis repository, run `python3 code/tools/report_dataset_inventory.py` to rebuild these documentation tables from existing audit outputs, then `python3 code/tools/verify_consolidated_inventory.py` and `python3 code/tools/check_current_decisions.py`. The renderer calls `consolidate_current_decisions.py` last so older summaries cannot overwrite current decisions. `current_decision_inputs.json` records SHA-256 hashes of its immediate evidence inputs.

The six-history scientific review has its own entry point, `python3 code/tools/run_pisco_biomass_review.py`; see [its methods and limits](../pisco_biomass_borrowing_20261102/README.md). Earlier all-program station, source and uncertainty audits retain their scripts and dated outputs. This consolidation does not rerun their analyses.

## Parallel follow-up audit

[Four-area findings and recommendations](../parallel_resolution_20261102/README.md) distinguish verified corrections from proposals that still need review. The six-comparison density review uses the biomass roster: all 63 post-protection years remain, but six earlier density comparison-years are omitted. Rebuild density-specific baseline coverage before before–after integration. Temporal borrowing validation contains eight predictions at two stations in one year, for red urchins and lobster only; no purple-urchin temporal validation is available. Recommendations do not change eligibility or models.

## Archived overviews and retained evidence

Stale inventory overview documents are preserved in `docs/archive/20261102_pre_consolidation/`. Dated scientific audit directories remain at their original paths so source links and reproducibility are preserved. Their findings describe those stages; use this guide and the current rules for present decisions. No manuscript, canonical result or raw-data files changed in this consolidation.
''')
    inputs=[B/'summary.csv',B/'visit_biomass_support.csv',Path('outputs/biomass_support_20261102/biomass_history_status.csv'),Path('outputs/biomass_exclusion_reconciliation_20261102/affected_histories.csv'),Path('docs/BIOMASS_ASSUMPTIONS_AND_CHANGES.csv')]
    (O/'current_decision_inputs.json').write_text(json.dumps({str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},indent=2)+'\n')
    if (follow/'recommended_actions.csv').exists():
        front['Audit recommendations']=pd.read_csv(follow/'recommended_actions.csv')
    from update_biomass_next_steps_guide import apply as apply_next_steps
    from update_joint_biomass_guide import apply as apply_joint
    from update_nonlinear_biomass_guide import apply as apply_nonlinear
    from update_scott_reconciliation_guide import apply as apply_scott
    return apply_scott(apply_nonlinear(apply_joint(apply_next_steps({**front,**views}))))
