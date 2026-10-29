"""Audit source labels and co-located survey dates without changing inputs.

Run with the data repository path; outputs are written in the analysis repo.
Same-day visits within a distance threshold are candidates, not duplicates.
"""
import argparse
from collections import Counter, defaultdict
from pathlib import Path
import hashlib
import json
import math

import pandas as pd


def distance(a, b, c, d):
    a, b, c, d = map(math.radians, (a, b, c, d))
    h = math.sin((a-c)/2)**2 + math.cos(a)*math.cos(c)*math.sin((b-d)/2)**2
    return 6371 * 2 * math.asin(min(1, math.sqrt(h)))


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('data_repo', type=Path)
    parser.add_argument('--out', type=Path, default=Path('data/external/source_reconciliation_2026-10-17'))
    parser.add_argument('--report', type=Path)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    counts, visits, manifest = Counter(), defaultdict(set), []
    for kind, filename in (
        ('fish', 'SBCMBON_kelp_forest_integrated_fish_20231022.csv'),
        ('benthos', 'SBCMBON_kelp_forest_integrated_quad_and_swath_20231022.csv'),
    ):
        path = args.data_repo / 'data/MBON' / filename
        for chunk in pd.read_csv(path, usecols=['data_source', 'sample_method', 'date', 'site_id'],
                                 dtype=str, chunksize=250000, encoding='latin-1'):
            counts.update((kind, str(source), str(method)) for source, method in
                          chunk[['data_source', 'sample_method']].itertuples(index=False, name=None))
            selected = chunk[chunk.data_source.eq('kfm')]
            if kind == 'fish':
                selected = selected[selected.sample_method.eq('visualfish')]
            for site, day in selected[['site_id', 'date']].drop_duplicates().itertuples(index=False, name=None):
                if pd.notna(pd.to_datetime(day, errors='coerce')):
                    visits[(kind, 'KFM', site)].add(day)
        manifest.append({'file': 'data/MBON/' + filename, 'sha256': sha256(path), 'bytes': path.stat().st_size})
    for kind, filename in (('fish', 'MLPA_kelpforest_fish_2024.csv'),
                           ('benthos', 'MLPA_kelpforest_swath_2024.csv')):
        path = args.data_repo / 'data/PISCO' / filename
        for chunk in pd.read_csv(path, usecols=['campus', 'method', 'year', 'month', 'day', 'site'],
                                 dtype=str, chunksize=250000, encoding='latin-1'):
            counts.update((kind, 'PISCO:' + str(source), str(method)) for source, method in
                          chunk[['campus', 'method']].itertuples(index=False, name=None))
            dates = pd.to_datetime(chunk[['year', 'month', 'day']], errors='coerce')
            for site, day in zip(chunk.site, dates):
                if pd.notna(day):
                    visits[(kind, 'PISCO', site)].add(day.strftime('%Y-%m-%d'))
        manifest.append({'file': 'data/PISCO/' + filename, 'sha256': sha256(path), 'bytes': path.stat().st_size})
    pd.DataFrame([{'kind': k, 'source': s, 'method': m, 'n_rows': n}
                  for (k, s, m), n in sorted(counts.items())]).to_csv(args.out/'survey_source_counts.csv', index=False)
    ksites = pd.read_csv(args.data_repo/'data/MBON/SBCMBON_kelp_forest_site_geolocation_20210120_KFM_LTER.csv')
    ksites = ksites[['site_id', 'data_source', 'Lat', 'Lon']].drop_duplicates()
    psites = pd.read_csv(args.data_repo/'data/PISCO/master_site_table_Emilyedit.csv', encoding='utf-8-sig')
    # Deduplicate coordinates only, not observations or site-to-MPA assignments.
    psites = psites[['site', 'latitude', 'longitude']].drop_duplicates()
    rows = []
    for k in ksites[ksites.data_source.eq('kfm')].itertuples():
        for p in psites.itertuples():
            if not all(pd.notna(x) for x in (k.Lat, k.Lon, p.latitude, p.longitude)):
                continue
            km = distance(k.Lat, k.Lon, p.latitude, p.longitude)
            if km > 1:
                continue
            for kind in ('fish', 'benthos'):
                ka = visits[(kind, 'KFM', k.site_id)]
                pa = visits[(kind, 'PISCO', p.site)]
                if not ka or not pa:
                    continue
                shared = sorted(ka & pa)
                rows.append({'kind': kind, 'kfm_site': k.site_id, 'pisco_site': p.site,
                             'distance_km': km, 'kfm_visit_dates': len(ka), 'pisco_visit_dates': len(pa),
                             'same_day_candidates': len(shared), 'candidate_dates': ';'.join(shared)})
    cols = ['kind', 'kfm_site', 'pisco_site', 'distance_km', 'kfm_visit_dates',
            'pisco_visit_dates', 'same_day_candidates', 'candidate_dates']
    pd.DataFrame(rows, columns=cols).to_csv(args.out/'survey_nearby_date_candidates.csv', index=False)
    summary = {'inputs': manifest, 'distance_threshold_km': 1,
               'nearby_site_kind_pairs': len(rows),
               'same_day_candidate_pairs': sum(x['same_day_candidates'] > 0 for x in rows),
               'scope': 'All dated KFM benthic surveys, KFM visual fish surveys, and supplied PISCO fish/swath surveys; not a primary-cohort duplicate determination.',
               'limitations': 'Nearby same-day surveys are not verified duplicates; site centres do not identify transect footprints. Distinct dates do not establish ecological independence.'}
    if args.report:
        summary['technical_report'] = {
            'filename': args.report.name, 'sha256': sha256(args.report),
            'bytes': args.report.stat().st_size,
            'url': 'https://kelpforest.msi.ucsb.edu/assets/files/1658/Kelp%20Forest%20Technical%20Report%20Narrative_v2.pdf',
            'title': 'Monitoring and Evaluation of Kelp Forest Ecosystems in the MLPA Marine Protected Area Network',
            'title_page_date': '2021-12-30', 'naples_reference_table_pdf_page': 116,
            'verified_coefficient_pdf_pages': [167, 177]}
    (args.out/'survey_audit_manifest.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps({k:v for k,v in summary.items() if k != 'inputs'}, indent=2))


if __name__ == '__main__':
    main()
