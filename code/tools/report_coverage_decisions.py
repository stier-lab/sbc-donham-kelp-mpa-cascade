from pathlib import Path
import pandas as pd,json,hashlib
O=Path('outputs/coverage_decisions_20261101')
r=pd.read_csv(O/'fish_reference_coverage.csv');r=r[(r.Resp=='Den')&(r.n_years_lost_requiring_all>0)][['MPA','n_years','n_years_all_stations','n_persistent_reference_stations']]
f=pd.read_csv(O/'fish_zone_timing_summary.csv');f=f[f.Resp=='Den']
text='''# Coverage decisions: next audit

November 1, 2026. Review only; manuscript and canonical results unchanged.

## Decisions recorded

Use shared, comparable depth zones. Check constant zones and show losses before adopting that restriction. Show seasonal coverage before selecting a timing cutoff. Compare constant reference membership with available-station averaging before choosing. Retain established references with incomplete documentation provisionally and test influence; preserve researcher exclusions. Keep supported density, but require justified size support for biomass. Neither the earlier 90-day cap nor five-year duration proposal is adopted.

## What the next screen establishes

Screened 1,048 station-pair/method/taxon combinations and 10,314 paired calendar-year records from the existing raw survey inventory, after recorded researcher scope/exclusions. These are overlapping coverage records, NOT independent biological comparisons or approved model inputs. Methods remain separate; apparent breaks across protocol names are not automatic exclusions. Satellite and separate size collections need different support checks and are not assigned diver zones.

The inspected KFM and LTER annual records have no usable zone labels or measured depths for this comparison. This is unknown matching, not evidence of a mismatch. Do not discard them automatically or borrow PISCO depth labels. Existing metadata and count/size-linkage limitations remain.

## PISCO fish: complete-pass coverage

All 17 sheephead assignments have at least one zone shared by the two roles in every retained year. Therefore a constant role-level zone set can retain all currently approved fish years at this screening stage. Five assignments change their shared zone set: Abalone Cove, Point Vicente, Santa Barbara Island, Scorpion and South Point. The others still require matching the actual zone definitions and measured depths. This test pools zone availability across stations within each role; it does NOT establish that every individual station samples that zone each year. Zone restriction could therefore change contributing stations and annual means even without removing years. Biological outcomes have not been refitted for these restrictions.

## Reference composition

Six PISCO sheephead assignments change reference-station membership through time. Requiring every reference ever used in every retained year would substantially shorten histories. The table shows density; biomass has the same pattern except South Point begins with 19 valid years and retains five all-station years.

A persistent subset retains the full time span at Anacapa (two reference stations), Point Dume SMCA (one), Scorpion (six), and South Point (one). No reference station spans all currently approved years at Abalone Cove or Santa Barbara Island. A persistent subset changes the spatial population; it is a diagnostic, not a newly approved reference assignment. The source identifiers are not collapsed into alleged aliases. Available-versus-persistent reference means on identical years are saved in fish_reference_means.csv; where no persistent set exists the alternative remains missing, not zero. This is not yet a causal-effect or pooled-model comparison. Protected-station composition also needs to be held or modeled explicitly in that next comparison.

## Timing

fish_complete_pass_zone_timing.csv gives actual protected/reference month sets and date gaps for each approved fish history-year. The summary reports the largest of the nearest date gaps across shared zones. It is an optimistic role-level availability diagnostic: a nearby reference visit does not certify every contributing station or survey campaign. No dates have been discarded. In sheephead density, Abalone has four years over 90 days under this screen; the other assignments have none. A 60-day cutoff additionally flags Anacapa, Blue Cavern, Gull, Long Point, Naples and Point Vicente. Read exact years/months before choosing; within-year proximity does not remove a change in season through time.

The broader pair-year file supplies dates/months for all inspected diver programs, but includes raw response records and separately named methods, not a final effort-normalized response rebuild. Existing five-MPA swath tests remain diagnostic and must not be presented as updated results under the latest fish rules.

## Recommendation and remaining work

Do not require every historical reference station every year: it discards much of the record and does not itself establish better habitat matching. Next compare station-aware estimates on common years with constant-membership alternatives, preserving the named reference groups. Jointly verify station-specific zone coverage before adopting constant role-level zone sets. Keep missing depth as an explicit limitation. Build the seasonal decision table from actual contributing visits, not just the nearest visit, before choosing a timing rule. Extend size support to the expanded roster. Run reference-omission influence checks; incomplete documentation alone is not an exclusion.

This pass provides coverage losses and reference means, not new biological recovery fits. It does not finish all-program station-standardized modeling, biomass-support checks or reference influence tests. No paper changes, source edits, emails or pushes.
'''
(O/'README.md').write_text(text)
page='<meta charset="utf-8"><title>Coverage decisions</title><style>body{font:16px system-ui;max-width:1250px;margin:30px}td,th{padding:8px;text-align:left;border-bottom:1px solid #ddd}table{border-collapse:collapse}p{max-width:950px}.scroll{overflow:auto}</style><h1>Coverage before choosing final rules</h1><p>Coverage screen only. No new biological models or paper changes.</p><p><a href="README.md">Findings, decisions and limitations</a> · <a href="pair_year_timing_zones.csv">All diver pair-year records</a> · <a href="fish_complete_pass_zone_timing.csv">Exact fish months and years</a> · <a href="fish_reference_means.csv">Reference means on the same years</a></p><h2>Requiring all historical reference stations: sheephead density</h2>'+r.to_html(index=False)+'<h2>Shared zones and timing: sheephead density</h2><p>Role-level availability, not certification of every contributing station. Timing columns count years flagged, not adopted exclusions.</p><div class="scroll">'+f.drop(columns=['series_id','Resp']).to_html(index=False)+'</div>'
(O/'index.html').write_text(page)
# Check data identity, coverage nesting, missing alternatives, and frozen protected files.
a=pd.read_csv('outputs/pisco_fish_rules_20261101/approved_annual_input.csv');a=a[a.Source.eq('PISCO')&a.Taxa.eq('S. pulcher')];fy=pd.read_csv(O/'fish_complete_pass_zone_timing.csv');assert set(zip(a.series_id,a.year))==set(zip(fy.series_id,fy.year))
rr=pd.read_csv(O/'fish_reference_coverage.csv');assert rr.n_years_all_stations.le(rr.n_years).all()
means=pd.read_csv(O/'fish_reference_means.csv');assert means.loc[means.n_persistent.eq(0),'persistent_reference_mean'].isna().all()
protected=json.loads(Path('outputs/researcher_scope_20261101/protected_hashes.json').read_text())
for p,h in protected.items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h
checks=dict(approved_fish_rows_checked=len(fy),histories=fy.series_id.nunique(),protected_unchanged=len(protected),coverage_nesting=True,missing_alternatives_preserved=True,models_refitted=0)
(O/'verification.json').write_text(json.dumps(checks,indent=2));print(checks)
paths=list(O.glob('*'))+[Path('code/tools/audit_coverage_decisions.py'),Path('code/tools/report_coverage_decisions.py')]
(O/'manifest.json').write_text(json.dumps({str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths if p.is_file() and p.name!='manifest.json'},indent=2))
