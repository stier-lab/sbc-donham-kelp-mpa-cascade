#!/usr/bin/env python3
"""Build an assignment ledger without changing sites or fitting models.

Run from repository root. Email evidence remains in the private manuscript repo.
Site IDs and surveyed coverage do not establish ecological comparability.
"""
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

def read(path):
    with path.open(newline='', encoding='utf-8-sig') as f:
        return list(csv.DictReader(f))

def write(path, rows):
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator='\n'); w.writeheader(); w.writerows(rows)

def main():
    registry = ROOT / 'data/harmonized/reference_sites.csv'
    rows = read(registry)
    groups = defaultdict(list)
    shared = defaultdict(set)
    for x in rows:
        groups[x['series_id']].append(x)
        if x['included'] == 'TRUE':
            shared[x['reference_id']].add(x['MPA'])
        if x['Source'] == 'LTER' and x['site_id'] in {'AHND', 'a-l-02'}:
            assert x['excluded'] == 'TRUE' and x['included'] == 'FALSE', 'LTER exclusion alias regression'
    keys = {x['series_id'] for x in read(ROOT/'data/harmonized/reference_primary_keys.csv')}
    assert keys == set(groups), 'Primary histories differ from registry'
    out = []
    for sid, entries in sorted(groups.items()):
        x = entries[0]; source = x['Source']; mpa = x['MPA']
        if source == 'PISCO':
            rule = 'Decadal-review site mapping; inherited suitability exclusions where better references existed'
            evidence = 'Historical pBACIPS_PISCO_V10.R lines 24-28,49,570ff; current PISCO mapping and exclusion constants'
            unresolved = 'Exact habitat/depth/exposure criterion and rejected alternatives for each assignment'
            if mpa == 'Harris Point SMR':
                rule = 'Cuyler retained for nearby north-side habitat; Crook Point and Tyler Bight excluded for different south-side habitats; Bay Point and Prince Island lack sustained coverage'
                unresolved = 'Qualitative habitat rationale confirmed; quantitative habitat/depth equivalence not demonstrated'
            elif mpa == 'South Point SMR':
                unresolved = 'Ford Point has only 2007 fish records in supplied files; historical ecological selection rationale still partly tentative'
            elif mpa in {'Gull Island SMR','Scorpion SMR'}:
                rule += '; Yellowbanks excluded from local PISCO pairs because of distinctive habitat and distance from Gull; Scorpion has closer alternatives'
            elif mpa == 'Naples SMCA':
                rule = 'Arroyo Quemado references selected using expert diver judgment and community-composition comparisons; no close match to Naples identified'
                unresolved = 'Qualitative basis confirmed; underlying ordination results not independently recovered and habitat equivalence not established'
                evidence += '; docs/SELECTION_CRITERIA_AUDIT_2026-10-29.md (Naples reference basis)'
            evidence += '; docs/REFERENCE_REPLY_INTEGRATION_2026-10-29.md (confirmed investigator rationale and supplied-file coverage)'
        elif source == 'LTER':
            rule = 'MBON assignment selected through consultation with local LTER investigators; Arroyo Hondo excluded because kelp was essentially absent'
            evidence = 'Historical pBACIPS_PISCO_V10.R line 2160; current MBON mapping and both LTER exclusion lists'
            unresolved = 'Consultation confirmed, but individual retained-site reasoning not recovered; excluded habitat rule does not establish retained habitat equivalence'
        elif source == 'KFM':
            rule = 'MBON assignment; intended continuity of stations across before/after periods; exclusions apply to benthic processing, not the sheephead branch'
            evidence = 'Historical pBACIPS_PISCO_V10.R lines 18-22,1124-1128,1459ff; current MBON mapping'
            unresolved = 'Per-site exclusion reason and ecological basis of retained pairing'
            if mpa == 'Gull Island SMR':
                rule += '; Yellowbanks identified as best available KFM reference, approximately 15 miles from Gull'
                unresolved = 'Distance and unverified habitat/depth comparability remain limitations; retain primary plus dedicated omission sensitivity; detailed map review pending'
            if x['Taxa'] == 'S. pulcher' and mpa in {'Harris Point SMR','Scorpion SMR'}:
                unresolved += '; after-only fish stations confirmed; same-year persistent-site sensitivity in module 53, not promoted to primary'
            evidence += '; docs/REFERENCE_REPLY_INTEGRATION_2026-10-29.md'
        else:
            rule = 'Supplied satellite reference footprints; geometric selection documented separately'
            evidence = 'docs/SOURCE_RECONCILIATION_2026-10-17.md; data/external/tom_bell_2026-10-20/'
            unresolved = 'Geometric selection does not establish habitat equivalence; registry labels do not identify shared pixels'
        included = [r for r in entries if r['included'] == 'TRUE']
        out.append(dict(series_id=sid,MPA=mpa,Source=source,Taxa=x['Taxa'],Resp=x['Resp'],
            included_reference_codes=';'.join(sorted({r['site_id'] for r in included})),
            excluded_reference_codes=';'.join(sorted({r['site_id'] for r in entries if r['excluded']=='TRUE'})),
            mapped_without_retained_coverage=';'.join(sorted({r['site_id'] for r in entries if r['excluded']=='FALSE' and r['included']=='FALSE'})),
            other_mpas_sharing_exact_reference_ids=';'.join(sorted(set().union(*(shared[r['reference_id']] for r in included))-{mpa})) if included else '',
            assignment_basis=rule,source_locator=evidence,independently_validated_by_this_audit='FALSE',
            remaining_question=unresolved,coverage_limit=';'.join(sorted({r['contribution_audit'] for r in entries}))))
    write(ROOT/'tables/table_s_reference_assignment_audit.csv',out)
    manifest={'registry_sha256':hashlib.sha256(registry.read_bytes()).hexdigest(),
        'histories':len(out),'mpa_program_comparisons':len({(r['MPA'],r['Source']) for r in out}),
        'source_counts':{s:sum(r['Source']==s for r in out) for s in sorted({r['Source'] for r in out})},
        'rule':'One row per durable primary-series key; all inherited assignments retained; code counts are not independent reef counts.',
        'private_evidence':'Email and original draft-comment evidence is held separately in the manuscript repository.'}
    (ROOT/'tables/reference_assignment_audit_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(manifest,indent=2))

if __name__ == '__main__':
    main()
