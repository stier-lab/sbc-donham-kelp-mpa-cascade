from pathlib import Path
import pandas as pd,json,hashlib,html
O=Path('outputs/biomass_protocols_20261102')
ledger=[
('B01','Decision','Preserve researcher exclusions','Adopted','No new site inclusion/exclusion based on results','User decision; previous review cohort'),
('B02','Decision','Keep counts when biomass size support is inadequate','Adopted','Density retained separately; no missing biomass converted to zero','User decision'),
('B03','Correction','Missing stipes imply unknown biomass','Implemented; production rebuild pending','Point Vicente kelp2023 paired year withheld in review input; upstream safeguard tested','biomass_resolution_20261102'),
('B04','Assumption','One original size record retains one identity across all uses','Implemented in conditional draws','Shared donor weights preserved; frequency counts not deduplicated','resolve_biomass_support.py'),
('B05','Sensitivity assumption','Halve/double singleton mean mass','Diagnostic only','Not an uncertainty interval or exclusion cutoff','biomass_resolution_20261102'),
('B06','Assumption','Offsets remain fixed at previously verified values','Conditional review only','Size/refit checks do not fully propagate count/allometry/offset uncertainty','resolve_biomass_support.py'),
('B07','Source convention','KFM legacy11005red;11006purple;size in cm','Confirmed','Convert to mm before comparing frequencies or applying allometry','ESA Table13A/B; metadata'),
('B08','Recovery rule','Recover date only for exact full site/year/species frequency match and unambiguous bin date','Implemented as overlay','58369 frequency records recovered; source file unchanged','reconcile_kfm_urchin_dates.py'),
('B09','Source discrepancy','Combined a-k-01 merges Wyckoff and Miracle Mile','Exact match in all51 available overlap groups','Both reference sites; does not establish same physical count footprint; no arbitrary record splitting','kfm_label_hypotheses.csv'),
('B10','Source discrepancy','Combined a-k-14 sizes match Southeast Reef37, counts match Sea Lion14','Confirmed for14 species-years2005–2011','Earlier same-site support claim superseded; original Sea Lion sizes tested as separate variant','kfm_label_hypotheses.csv;sbi_count_identity_checks.csv'),
('B11','Model change','Replace SBI protected sizes with original count-site sizes2005–2011','Sensitivity only; not primary','Same annual observations; later size identity unresolved; variant partially reconciled','sbi_reconciliation_effects.csv'),
('B12','Unresolved','Annual count-size matching implies same visit','Not assumed','Recovered dates show gaps; no timing cutoff adopted','retained_count_size_date_gaps.csv'),
('B13','Protocol','PISCO conspicuous benthic organisms approximately2.5cm and larger','Broad protocol supported; exact boundary unresolved','Preserve current>=25mm rule; do not impose on KFM/LTER sizes without compatible counts','Sea Grant technical report pp22–24'),
('B14','Unresolved','Borrowed sizes represent count-site animals','Not established','Keep affected biomass provisional; same program/year is insufficient proof','red_urchin_decisions.csv; new source identity checks'),
('B15','Unresolved','Total biomass uncertainty','Incomplete','Conditional size draws exclude count/allometry uncertainty; KFM kelp plant sampling remains','biomass_resolution_20261102'),
('B16','Hold','Paper and primary release','Unchanged','No manuscript update, final eligibility freeze or primary rebuild','User review hold'),
('B17','Source confirmation','Alias problems persist in2015','Confirmed in all4grouped-distribution checks','Do not infer individual sizes from grouped report means','biomass_exclusion_reconciliation_20261102/2015_site_identity_checks.csv'),
('B18','Implementation discrepancy','Intended excluded KFM sites21and37 survive under alternate size labels','Original and current code explicitly exclude them','Six retained histories affected; restore intended donor exclusions using original physical-site records','biomass_exclusion_reconciliation_20261102/affected_histories.csv'),
('B19','Sensitivity','Reconstruct allowed original donor pools through2011','24fits independently verified; partial only','Counts/years/offsets held fixed; recent unresolved donors unchanged; no primary promotion','biomass_exclusion_reconciliation_20261102/effects.csv'),
('B20','Source access','NPS1982–2025package2318436 includes site-specific size and survey tables','Catalog visible; files markedInternal/CanDownloadfalse','Need authorized source access; public processing script obtained; no restricted file download attempted','https://irma.nps.gov/DataStore/Reference/Profile/2318436'),
('B21','Guard','Block excluded physical sources and known misleading urchin aliases','Implemented in new review resampling; production hold','Six histories withheld from joint draws; density unchanged','biomass_joint_checks_20261102/resampling_admission.csv'),
('B22','Source check','KFM exact duplicate plant rows have identical site/date/marker keys','Verified:18413duplicates;46802uniqueplants','Do not deduplicate equal sizes from distinct plants','biomass_joint_checks_20261102/verification.json'),
('B23','Conditional sampling assumption','Spatial count weights and shared size weights propagated together','1000draws;34histories','Fixed years/offsets/calibration; draw ranges are not total confidence intervals','biomass_joint_checks_20261102/README.md'),
('B24','Matching audit','Compare donor locations/date/zone/depth against retained counts','1143donor-bearing role-years screened','No new cutoff;23PISCO zone flags; borrowed-size representativeness unresolved','biomass_joint_checks_20261102/count_size_matching.csv'),
('B25','Numerical check','Deterministic KFM kelp expected mass versus cached simulated point estimate','296station-years independently verified','Station differences up to~3percent; coordinated rebuild pending','biomass_joint_checks_20261102/kelp_expected_biomass_check.csv'),
('B26','Source resolution','All23PISCO zone flags arise from seasonal filtering','43station/species/zone groups have raw count surveys on size dates','Supersedes missing-zone hypothesis; March/April/November/December counts filtered but sizes retained','pisco_size_zone_resolution_20261102/source_zone_details.csv'),
('B27','Sensitivity','Apply existing count season to size donors on matched years','66fits verified;diagnostic only','Three paired biomass years lost; no new primary season cutoff; retain density','pisco_size_zone_resolution_20261102/summary.csv'),
('B28','User decision','Allow valid surveys in all months; equal-weight repeat visits','Adopted for review rerun','Six affected PISCO histories; preserve researcher exclusions and count support','pisco_all_months_20261102/README.md'),
('B29','Matching rule','Constant shared INNER/MID zones; raw same-visit sizes then nearest within site/zone/year','Diagnostic rerun complete; primary unchanged','624exact/24nearest/153unsupported positive uses; dates exposed; exact-only sensitivity','pisco_all_months_20261102/size_matching.csv'),
('B30','Coverage finding','Local biomass completeness differs from usable density coverage','63density comparison-years;25local biomass years','Do not turn missing local sizes into zero or silently drop a reference station','pisco_all_months_20261102/summary.csv'),
('B31','Sensitivity','Pair protected/reference dates without dropping reference stations','0/14/30/60/90/366day grid; no cutoff adopted','At90days Abalone retains2postyears,Vicente6; source campaign IDs unavailable','pisco_timing_donors_20261102/timing_coverage.csv'),
('B32','Validation','Predict held-out local mean body mass from other samples','624targets;5860predictions; cross-site donors exclude target site','Typical90day same-role errors1.33–1.48fold; sparse target samples and validation coverage explicit','pisco_timing_donors_20261102/donor_validation_summary.csv'),
('B33','Limit','Donor availability does not certify biomass recovery','120of153missing visits have same-role/zone donors within90days','No imputation; Abalone protected cross-site validation unavailable; total uncertainty pending','pisco_timing_donors_20261102/README.md'),
('B34','Implementation','Carry all-month density as count-led review response','63postcomparison-years unchanged','Six PISCO histories; paper and canonical results held','pisco_biomass_borrowing_20261102/selected_density_annual.csv'),
('B35','Uncertainty assumption','Nested outer-site donor calibration and shared site-year prediction errors','1000conditional draws; reference-only interval validation','No target-site leakage; no extra donor bootstrap; fixed counts/allometry; not totalCIs','pisco_biomass_borrowing_20261102/README.md'),
('B36','Conditional coverage','Borrowed biomass available for44of63years;38positive visits remain unknown','Review sensitivity implemented; no protected borrowing','Local-only prior25used different temporal-support rule; no silent station omissions','pisco_biomass_borrowing_20261102/summary.csv')]
df=pd.DataFrame(ledger,columns=['id','type','assumption_or_change','status','consequence','evidence']);df.insert(1,'updated','2026-11-02');
extra=Path('outputs/biomass_next_steps_20261102/decision_ledger_additions.csv')
if extra.exists(): df=pd.concat([df,pd.read_csv(extra)],ignore_index=True)
assert df.id.is_unique
df.to_csv('docs/BIOMASS_ASSUMPTIONS_AND_CHANGES.csv',index=False);df.to_csv(O/'assumptions_and_changes.csv',index=False)
readme='''# Biomass assumptions and source reconciliation

The assumptions/change ledger separates confirmed source facts, adopted corrections, diagnostic assumptions and unresolved decisions. It is maintained in docs/BIOMASS_ASSUMPTIONS_AND_CHANGES.csv. Earlier results remain archived; superseded assumptions are identified rather than silently erased.

## What the source checks resolved

The original KFM metadata confirms species11005 is red urchin,11006 is purple urchin, and legacy test diameters are in centimeters. Matching complete size-frequency distributions after conversion to millimeters gives873 exact matches among904 overlapping site/year/species groups. Dates can be recovered for58,369 frequency records representing159,963 animals. These are frequency records, not58,369 independent dives or newly collected measurements. In retained donor pools,26,015 unique original records benefit from the overlay. Raw source data remain unchanged.

Recovered dates do not prove same-visit sampling. All873 matched groups have count dates in MBON;401 have size dates that match a count date exactly. The largest nearest-date gap is219days. Among511 recovered donor-site/year/species uses in retained histories,193 have a nonzero gap. These are donor-site timing comparisons, not confirmation that borrowed sizes coincide with the target count sites. No cutoff or automatic exclusion was applied.

## Two site labels require special treatment

**Wyckoff Ledge:** all51 overlapping species/year groups labeled a-k-01 match the combined legacy distributions from Wyckoff Ledge1 and Miracle Mile21. Fifteen disagree with Wyckoff alone; the other36 predate/add no Miracle Mile measurements. Both sites are references for Harris Point. This supports a grouped source identity, not a license to assign every animal to Wyckoff or divide aggregated frequencies arbitrarily.

**Santa Barbara Island:** all14 red/purple species-year distributions labeled a-k-14 in2005–2011 exactly match Southeast Reef37. The corresponding MBON count densities match Southeast Sea Lion14 in all14checks. Thus the counts and sizes represent different sites. The previous audit correctly matched the delivered labels but its physical same-site interpretation was wrong. This finding supersedes that interpretation and qualifies the earlier count-site-only donor sensitivities. Later years have not been source-reconciled.

## What changes in an isolated same-observation test

Using original Southeast Sea Lion sizes for2005–2011, with existing offsets fixed, changes the17-year red-urchin recovery slope from−0.123 [−0.193,−0.052] to−0.006 [−0.041,0.028] log-ratio units per year. Purple changes from−0.068 [−0.202,0.066] to−0.008 [−0.160,0.145]. These intervals describe temporal model uncertainty conditional on the annual inputs, not total biomass uncertainty.

This is a partially reconciled diagnostic: only2005–2011 sizes were replaced; the remaining years retain current inputs. It is not a corrected final17-year estimate. The table also reports the seven-year matched comparison separately. The result shows that source-site identity matters and must be resolved before interpreting this biomass history. It does not justify choosing one source because of the direction of its slope. Density is unchanged.

## Protocol boundaries

The PISCO technical report supports a roughly2.5cm minimum for conspicuous benthic swath organisms. Its greater-than wording does not independently settle whether exactly25mm is included; the current>=25mm implementation is retained and its boundary remains explicit. KFM size protocols include small animals and sample across the permanent transect area. Changing only the size threshold would not standardize counts. LTER-specific count/size protocol reconciliation remains open.

## Next checks

- Resolve a-k-14 and other grouped site labels in later size releases and the original harmonization code.
- Recover recent KFM survey dates and assess count-size timing without selecting a cutoff by outcome.
- Recompute physical count-site support using verified donor identities; earlier label-based percentages are not final physical-coverage estimates.
- Complete joint count/size uncertainty and KFM kelp plant-sampling checks.

No raw data, manuscript or frozen primary result was changed in this pass. Eight new model fits were independently verified in R; the14count identities and14size identities are exact within numerical precision. The source date overlay remains private; public tables contain aggregated checks.

## Sources

[Kushner et al. KFM metadata](https://esapubs.org/archive/ecol/E094/245/metadata.php), [size variable definitions](https://esapubs.org/archive/ecol/E094/245/Table13B_invertebrate_size_frequency_variables.html), and [PISCO technical report](https://caseagrant.ucsd.edu/sites/default/files/R_MPA-43_Kelp%20Forest%20Technical%20Report%20Narrative.pdf). Downloaded source copies and record-level overlays are retained privately with checksums.
'''
(O/'README.md').write_text(readme)
body=''.join('<p>'+html.escape(x).replace('\n','<br>')+'</p>' for x in readme.split('\n\n'))
for title,file in [('Assumptions and changes','assumptions_and_changes.csv'),('SBI same-observation refits','sbi_reconciliation_effects.csv'),('Retained date support','retained_history_date_support.csv'),('Source discrepancies','kfm_source_discrepancies.csv')]:body+=f'<h2>{title}</h2><a href="{file}">Download CSV</a><div style="overflow:auto">'+pd.read_csv(O/file).to_html(index=False)+'</div>'
(O/'index.html').write_text('<!doctype html><meta charset="utf-8"><title>Biomass source reconciliation</title><style>body{max-width:1100px;margin:30px auto;padding:20px;font:17px/1.5 system-ui}td,th{padding:8px;text-align:left}table{font-size:13px;border-collapse:collapse}tr:nth-child(even){background:#eef4f4}</style>'+body)
protected=json.loads(Path('outputs/researcher_scope_20261101/protected_hashes.json').read_text())
for p,h in protected.items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h
scripts=[Path('code/tools')/n for n in ['reconcile_kfm_urchin_dates.py','check_kfm_date_timing.py','check_kfm_size_site_labels.py','refit_sbi_size_reconciliation.py','verify_sbi_size_reconciliation.R','report_biomass_protocols.py']]
(O/'manifest.json').write_text(json.dumps({str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in list(O.glob('*'))+scripts if p.is_file() and p.name!='manifest.json'},indent=2))
print(f'Ledger{len(df)} entries; source report written;487protected files unchanged.')
