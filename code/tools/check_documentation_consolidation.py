"""Check current guide navigation, archive integrity and current status overlay."""
from pathlib import Path
import csv
import hashlib
import re

A = Path(__file__).resolve().parents[2]
D = A.with_name('sbc-donham-kelp-mpa-data')
O = A / 'outputs/dataset_inventory_20261101'
files = [A/'README.md', A/'docs/README.md', A/'docs/ANALYSIS_GUIDE.md', O/'README.md', D/'README.md', D/'docs/RAW_TO_BIOMASS.md']
links = 0
for p in files:
    for link in re.findall(r'\]\(([^)]+)\)', p.read_text()):
        if link.startswith(('http:', 'https:', '#')):
            continue
        assert (p.parent / link.split('#')[0]).exists(), (p, link)
        links += 1
snapshots = 0
for root in [A, D]:
    archive = root / 'docs/archive/2026-11-02-consolidation'
    rows = list(csv.DictReader((archive/'manifest.csv').open()))
    destinations = [Path(row['archive']).resolve() for row in rows]
    assert len(set(destinations)) == len(destinations), 'Duplicate archive destinations'
    assert not set(destinations) & {(archive/'README.md').resolve(), (archive/'manifest.csv').resolve()}, 'Archive index collision'
    for row in rows:
        path = Path(row['archive'])
        assert hashlib.sha256(path.read_bytes()).hexdigest() == row['sha256'], path
        snapshots += 1
rules = {r['rule_id']: r for r in csv.DictReader((O/'current_rules.csv').open())}
integrated=(A/'outputs/integrated_review_20261103/completed.json').exists()
assert ('Corrected source build completed' if integrated else 'recent source identities reconciled') in rules['R08']['implementation_status']
assert 'recent sources pending' not in rules['R08']['user_decision_status']
questions = {r['issue_id']: r for r in csv.DictReader((O/'unresolved_issues.csv').open())}
assert 'recovered (219)' in questions['U05']['finding']
assert 'three sheephead remain unrecovered' in questions['U05']['finding']
assert 'same-visit sizes first' in (A/'docs/ANALYSIS_GUIDE.md').read_text()
assert ('historical production pool' if integrated else 'inherited production pool') in (A/'docs/ANALYSIS_GUIDE.md').read_text()
assert ('downstream paper analyses still require regeneration' if integrated else 'full harmonization and downstream review rebuild pending') in rules['R07']['limitation_or_supersedes']
print(f'Current guide links ({links}), archive entries ({snapshots}), and repaired source-status descriptions verified.')

# Documentation routes must resolve to existing evidence and executable scripts.
import pandas as pd
workflow = pd.read_csv(O/'workflow_status.csv')
assert workflow.workflow_id.is_unique
for row in workflow.itertuples():
    assert (A/row.working_directory_relative_to_analysis_repo).is_dir()
    assert (A/row.input_record).exists(), row.input_record
    assert (A/row.evidence).exists(), row.evidence
    root = A/row.working_directory_relative_to_analysis_repo
    for command in row.commands.split('; '):
        if command.startswith('python3 '):
            assert (root/command.split()[1]).is_file(), command
        elif command.startswith('make '):
            assert re.search(r'^'+re.escape(command.split()[1])+r':', (root/'Makefile').read_text(), re.M), command
for filename, key in [('inclusion_gates.csv','gate_id'), ('biomass_calculations.csv','method_id')]:
    table = pd.read_csv(O/filename)
    assert table[key].is_unique and not table.isna().any().any(), filename
print('Workflow commands, input records, evidence paths and documentation-table keys verified.')
