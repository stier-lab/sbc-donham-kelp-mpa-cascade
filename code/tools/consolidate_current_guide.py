"""Final documentation overlay; no eligibility flags or numerical results change."""
from pathlib import Path
import os

O = Path('outputs/dataset_inventory_20261101')

def apply(views):
    rules = views['Current rules'].copy()
    ix = rules.rule_id.eq('R08')
    rules.loc[ix, 'implementation_status'] = 'Historical and recent source identities reconciled; physical-site safeguards implemented; complete review rebuild pending'
    rules.loc[ix, 'user_decision_status'] = 'Correction authorized; required source gate passes'
    rules.loc[ix, 'evidence_path'] = 'outputs/pelican_pisco_recovery_20261102/README.md'
    rules.loc[ix, 'evidence'] = 'KFM source reconciliation and accepted Pelican source'
    rules.loc[rules.rule_id.eq('R07'), 'limitation_or_supersedes'] = 'Missing-stipe source and role-year safeguards implemented; full harmonization and downstream review rebuild pending'
    rules = rules.drop(columns=['latest_kfm_source_resolution'], errors='ignore')
    rules.to_csv(O / 'current_rules.csv', index=False)
    views['Current rules'] = rules
    q = views['Questions'].copy()
    ix = q.issue_id.eq('U05')
    q.loc[ix, 'issue'] = 'Remaining missing PISCO sizes/stipes'
    q.loc[ix, 'finding'] = 'Cavern Point West 2013 red-urchin count recovered (219). Underwater Arch 2023 and Little Dume West 2011 kelp stipes remain unknown; South Point 2005 lengths for three sheephead remain unrecovered.'
    q.loc[ix, 'next_action'] = 'Preserve observed counts and response-specific missingness; verify full aggregation and biomass rebuild. Do not replace unknown sizes or stipes with zero.'
    q.loc[ix, 'evidence'] = 'outputs/pelican_pisco_recovery_20261102/README.md; outputs/pisco_source_reconciliation_20261101/README.md'
    q.loc[q.issue_id.eq('U08'), 'evidence'] = 'outputs/pelican_pisco_recovery_20261102/README.md; outputs/nonlinear_biomass_calibration_20261102/README.md'
    audit_path = 'outputs/notes_email_audit_20261102/README.md'
    q.loc[q.issue_id.eq('U07'), 'finding'] = 'Scott confirmed Nov2: KFM collects no lobster carapace lengths; fish size frequencies come from one diver on a 30-minute roving survey. A new fish-size attachment is available but its matching and coverage have not been verified in this audit.'
    q.loc[q.issue_id.eq('U07'), 'next_action'] = 'Reconcile fish-size attachment with count dates, stations and populations; retain density and account for unreplicated roving size sampling. Do not infer lobster lengths from fish data.'
    q.loc[q.issue_id.eq('U07'), 'evidence'] = audit_path
    additions = [
        ('U15', 'KFM', 'New reference-site advice', 'Scott flags Webster’s Arch offshore habitat, Rodes Reef island-side mismatch and north/south San Miguel references. Some are already excluded for benthic responses; fish scope differs.', 'Response-specific reference review', 'Existing physical-site map and Scott Gabara Nov2 reply', 'Map retained contributors by response and test influence before adopting any further exclusions.'),
        ('U16', 'LTER kelp', 'Updated SBC50 annual biomass', 'Li confirmed the annual biomass package is updated. The registered annual kelp file remains the 20240305 delivery, with biomass computed in our pipeline.', 'Provider release reconciliation', 'Li Kui Nov2 reply; SBC50 and current raw-source registry', 'Compare matched observations and conversion formulas with updated SBC50; do not assume a numerical difference or replace inputs silently.'),
        ('U17', 'Before/after comparisons', 'Historical baseline versus revised rule', 'Emily clarified three observations at the same stations, separately for fish and benthic surveys, with implementation year in before. This differs from a strict pre-implementation calendar rule.', 'Model decision / explicit historical departure', 'Emily Donham Nov2 reply; station-specific coverage tables', 'Freeze the revised convention explicitly and show matched station-year coverage; do not describe it as unchanged from Emily’s original rule.'),
    ]
    import pandas as pd
    for ident, scope, issue, finding, kind, where, action in additions:
        q = q.loc[~q.issue_id.eq(ident)]
        q = pd.concat([q, pd.DataFrame([dict(issue_id=ident, scope=scope, issue=issue, finding=finding, type=kind, who_or_where=where, next_action=action, evidence=audit_path)])], ignore_index=True)
    q.to_csv(O / 'unresolved_issues.csv', index=False)
    views['Questions'] = q
    reconciliation = O.parent / 'fish_lter_reconciliation_20261102'
    if (reconciliation / 'README.md').exists():
        evidence = 'outputs/fish_lter_reconciliation_20261102/README.md'
        q.loc[q.issue_id.eq('U07'), 'finding'] = 'Scott fish export reconciled: 2007–2025, 5658 sheephead rows/10073 fish. Current count roster size coverage rises118→282 of423site-years;279have all size dates matching count dates. No pre2003sizes; roving footprint and sparse/unreplicated sampling remain limitations. KFM collects no lobster carapace lengths.'
        q.loc[q.issue_id.eq('U07'), 'next_action'] = 'Use exact count-visit sizes first in a restricted-period review; retain explicit one-day and extra-visit diagnostics. Resolve footprint support and propagate count/size uncertainty before considering primary biomass.'
        q.loc[q.issue_id.eq('U16'), 'finding'] = 'SBC50.21 retrieved and compared. Adult source observations through2023 are unchanged. Provider MAPY combines adult stipes with juvenile density and biomass; it is not the same response as our adult foliar biomass. Adult-only values imply0.10234kg/stipe while the linked protocol says August; current pipeline uses actual-month127.5coefficients.'
        q.loc[q.issue_id.eq('U16'), 'next_action'] = 'Keep raw adult stipes and actual-month conversion; clarify provider fixed-coefficient/protocol discrepancy before adopting its aggregate. Missing provider values stay unknown. No automatic addition of2024–2026.'
        q.loc[q.issue_id.isin(['U07', 'U16']), 'evidence'] = evidence
        q.to_csv(O / 'unresolved_issues.csv', index=False)
        views['Questions'] = q
        views['Fish size reconciliation'] = pd.read_csv(reconciliation / 'fish_coverage_summary.csv').fillna('')
        views['LTER biomass comparison'] = pd.read_csv(reconciliation / 'lter_ratio_change_summary.csv').fillna('')
    kfm_review = O.parent / 'kfm_sheephead_biomass_20261103'
    if (kfm_review / 'README.md').exists():
        q.loc[q.issue_id.eq('U07'), 'finding'] = 'Scott file supports sizes in282of423site-years; visit-level matching supports52comparison-years for restricted sheephead biomass. Same-site size matches<=7days; shared size draws for reused visits; no cross-site/year borrowing. Count and size uncertainty propagated with zero-count uncertainty and sparse-size sensitivity. All main t intervals include zero; Anacapa is interval-method sensitive. Roving/transect footprint and KFM lobster sizes remain unresolved.'
        q.loc[q.issue_id.eq('U07'), 'next_action'] = 'Keep full-history density primary; review restricted biomass as supporting evidence. Resolve footprint and reference suitability; do not promote favorable conditional intervals.'
        q.loc[q.issue_id.eq('U07'), 'evidence'] = 'outputs/kfm_sheephead_biomass_20261103/README.md; outputs/fish_lter_reconciliation_20261102/README.md'
        q.to_csv(O / 'unresolved_issues.csv', index=False)
        views['Questions'] = q
        trends = pd.read_csv(kfm_review / 'trends.csv')
        views['KFM sheephead biomass'] = trends.loc[(trends.scenario == 'matched') & (trends.rho == 0)].fillna('')
        views['KFM visit losses'] = pd.read_csv(kfm_review / 'lost_comparison_years.csv').fillna('')
    recovery = O.parent / 'kfm_missing_size_recovery_20261103'
    if (recovery / 'README.md').exists():
        q.loc[q.issue_id.eq('U07'), 'finding'] += ' Five-year source audit recovered no additional lengths; matched-visit filtering could support55years but changes count visits and seasonal coverage. Current52years unchanged.'
        q.loc[q.issue_id.eq('U07'), 'next_action'] = 'Keep full-history density and current52-year biomass review. Optionally test55-year matched-visit sensitivity separately; do not drop sites or widen timing to recover years. Resolve footprint/reference questions.'
        q.loc[q.issue_id.eq('U07'), 'evidence'] += '; outputs/kfm_missing_size_recovery_20261103/README.md'
        q.to_csv(O / 'unresolved_issues.csv', index=False)
        views['Questions'] = q
        views['KFM missing-size recovery'] = pd.read_csv(recovery / 'five_year_options.csv').fillna('')
    views['Notes and email audit'] = pd.read_csv(O.parent / 'notes_email_audit_20261102/findings.csv').fillna('')
    f = views['Freeze readiness'].copy()
    f.loc[f.item_id.eq('F3'), 'latest_status'] = 'Required source gate passes; original Pelican 140 accepted, four unused KFM 2023 groups deferred, Cavern count recovered. Remaining missing sizes/stipes stay unknown; full biomass rebuild and support review pending.'
    f.to_csv(O / 'freeze_readiness.csv', index=False)
    views['Freeze readiness'] = f
    (O / 'README.md').write_text('''# Current decision tables

Updated November 3, 2026. Read the [current analysis guide](../../docs/ANALYSIS_GUIDE.md) for the pipeline, site and response selection, biomass estimation, missing data, borrowing and remaining release checks.

The required-source gate passes. Complete review harmonization and coordinated model reruns remain pending. The paper and frozen results are unchanged. Review-only alternatives must not be described as implemented across all datasets.

- [Searchable tables](index.html)
- [Workflow status and rebuild commands](workflow_status.csv)
- [Inclusion gates](inclusion_gates.csv) and [biomass calculations](biomass_calculations.csv)
- [Current rules and implementation status](current_rules.csv)
- [Comparison decisions and exclusion reasons](comparison_decisions.csv)
- [Site coverage](site_response_coverage.csv) and [protection metadata](mpa_metadata.csv)
- [Biomass support](measurement_support_current.csv)
- [Detailed visit decisions](measurement_visit_support_detailed.csv) and [donor attempts](donor_attempts_current.csv)
- [Remaining issues](unresolved_issues.csv) and [freeze checklist](freeze_readiness.csv)
- [Historical decision entries](decision_history.csv): dated evidence, not current policy

The guide defines each table's row and scope. Some support metrics describe earlier donor pools and are labeled accordingly. Neither candidate rows nor source records are a frozen primary cohort.

Superseded narrative instructions are preserved in the [documentation archive](../../docs/archive/2026-11-02-consolidation/README.md). Dated output folders retain numerical evidence and reproducible audit scripts; their old status summaries do not override this guide.
''')
    from document_workflow_tables import add_workflow_tables
    views = add_workflow_tables(views, O)
    from update_integrated_guide import apply as integrated_overlay
    return integrated_overlay(views)
