"""Reproduce the two README recoveries and audit cleanup evidence, without model runs."""
from pathlib import Path
import csv
import hashlib
import json
import subprocess
import tempfile

A = Path(__file__).resolve().parents[2]
D = A.with_name('sbc-donham-kelp-mpa-data')
R = A / 'docs/archive/2026-11-02-consolidation/recovery'

def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()

def check():
    results = {}
    steps = json.loads((R/'reconstruction_steps.json').read_text())
    assert len(steps) == 18
    for label, root in [('analysis', A), ('data', D)]:
        archive = root/'docs/archive/2026-11-02-consolidation'
        with tempfile.TemporaryDirectory() as folder:
            target = Path(folder)/'README.md'
            target.write_bytes((archive/'repository_README_at_HEAD.md').read_bytes())
            selected = [step for step in steps if step['repo'] == label]
            for step in selected:
                assert sha(target) == step['before_sha256'], step['event_line']
                # Trusted, reviewed recovery recipe. Path resolves only to the temporary README;
                # no historical shell commands, network calls or other file edits are replayed.
                exec(compile(step['code'], f"readme-recovery-{step['event_line']}", 'exec'),
                     {'Path': lambda value: target, 'target': target, 'root': target.parent, 'r': target.parent})
                assert sha(target) == step['after_sha256'], step['event_line']
            recovered = archive/'repository_README_before_cleanup.md'
            assert target.read_bytes() == recovered.read_bytes()
            excerpt = (R/(label+'_pre_cleanup_excerpt.txt')).read_text()
            assert recovered.read_text().startswith(excerpt)
            results[label] = {'successful_edits_replayed': len(selected), 'recovered_lines': len(recovered.read_text().splitlines()), 'independent_excerpt_characters': len(excerpt), 'sha256': sha(recovered)}
        rows = list(csv.DictReader((archive/'manifest.csv').open()))
        assert len({r['archive'] for r in rows}) == len(rows)
        assert not set(r['archive'] for r in rows) & {str(archive/'README.md'), str(archive/'manifest.csv')}
        for row in rows:
            assert sha(row['archive']) == row['sha256'], row['archive']
            original = Path(row['original'])
            assert original.exists(), original
            if row['handling'] == 'Replaced with archive pointer':
                assert 'Preserved original' in original.read_text(), original
                # Resolve every stub's local links, including its preserved original.
                import re
                links = re.findall(r'\]\(([^)]+)\)', original.read_text())
                assert any((original.parent/link).resolve() == Path(row['archive']).resolve() for link in links)
                assert all((original.parent/link).exists() for link in links)
        results[label]['archive_entries_verified'] = len(rows)
        deleted = subprocess.check_output(['git', '-C', str(root), 'ls-files', '--deleted']).decode().splitlines()
        assert not deleted, deleted
        results[label]['deleted_tracked_files'] = deleted
    # The 47 unaffected snapshots must match their original, pre-recovery manifests.
    unchanged = 0
    for label, root in [('analysis', A), ('data', D)]:
        for row in csv.DictReader((R/(label+'_manifest_before_recovery.csv')).open()):
            if Path(row['original']) != root/'README.md':
                assert sha(row['archive']) == row['sha256']
                unchanged += 1
    assert unchanged == 47
    results['original_archive_snapshots_unchanged'] = unchanged
    protected = json.loads((A/'outputs/researcher_scope_20261101/protected_hashes.json').read_text())
    for path, digest in protected.items():
        assert sha(path) == digest, path
    results['protected_files_unchanged'] = len(protected)
    evidence = json.loads((A/'outputs/dataset_inventory_20261101/current_decision_inputs.json').read_text())
    for path, digest in evidence.items():
        assert sha(A/path) == digest, path
    results['current_decision_evidence_hashes_verified'] = len(evidence)
    sources = json.loads((D/'config/source_manifest.json').read_text())['sources']
    for source in sources:
        assert sha(D/source['path']) == source['sha256'], source['source_id']
    results['registered_raw_sources_verified'] = len(sources)
    build = json.loads((D/'data/processed/source_build/build_manifest.json').read_text())
    for path, digest in build['outputs'].items():
        assert sha(D/'data/processed/source_build'/path) == digest, path
    for path, digest in build['configuration_sha256'].items():
        assert sha(D/path) == digest, path
    results['source_build_outputs_verified'] = len(build['outputs'])
    results['source_build_configurations_verified'] = len(build['configuration_sha256'])
    signatures = {'input_manifest_sha256': 'config/source_manifest.json', 'scope_sha256': 'config/kfm_size_scope.csv', 'site_map_sha256': 'config/kfm_physical_sites.csv', 'builder_sha256': 'code/python/raw_sources.py'}
    for field, path in signatures.items():
        assert sha(D/path) == build[field], path
    results['source_builder_registry_and_scope_hashes_verified'] = len(signatures)
    return results

if __name__ == '__main__':
    print(json.dumps(check(), indent=2))
