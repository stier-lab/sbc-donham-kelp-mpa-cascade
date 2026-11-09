"""Verify the primary cohort implements the satellite baseline policy."""
import csv
from pathlib import Path
root=Path(__file__).resolve().parents[2]
def read(p):return list(csv.DictReader((root/p).open()))
screen=read('tables/table_s_satellite_baseline_screen.csv')
primary=read('data/sumstats_final.csv')
rejected={r['MPA'] for r in screen if r['baseline_eligible']=='FALSE'}
assert rejected=={'Santa Barbara Island SMR','Swamis SMCA'}
assert not any(r['Source']=='Landsat' and r['MPA'] in rejected for r in primary)
assert len(primary)==143
assert sum(r['AnalysisType']=='pBACIPS' for r in primary)==40
assert sum(r['Source']=='Landsat' for r in primary)==14
for source,taxon in [('Landsat','M. pyrifera'),('LTER','S. pulcher')]:
 assert any(r['MPA']=='Naples SMCA' and r['Source']==source and r['Taxa']==taxon and r['Resp']=='Bio' for r in primary)
print('PASS: 143 histories; 14 screened satellite histories; both Naples biomass histories retained; no failed satellite baseline routed to CI.')
