"""Independent checks of fixed membership, contrasts, and unchanged primary inputs."""
from pathlib import Path
import csv, hashlib, math
ROOT = Path(__file__).resolve().parents[2]
def read(p):
    with (ROOT / p).open() as f: return list(csv.DictReader(f))
def key(r): return '|'.join(r[k] for k in ('MPA','Taxa','Source','Resp'))
primary = read('data/sumstats_final.csv')
series = read('tables/table_s_routing_clock_series.csv')
keys = {key(r) for r in primary}
scenarios = sorted({r['scenario'] for r in series})
for name in scenarios:
    subset = [r for r in series if r['scenario'] == name]
    assert len(subset) == len(keys) == 145
    assert {key(r) for r in subset} == keys
    for r in subset:
        b = next(x for x in primary if key(x) == key(r))
        assert r['N'] == b['N']
# Independently reconstruct baseline linear/mean/step contrasts with elementary OLS.
short = {'Panulirus interruptus':'P. interruptus','Semicossyphus pulcher':'S. pulcher',
         'Strongylocentrotus purpuratus':'S. purpuratus','Mesocentrotus franciscanus':'M. franciscanus','Macrocystis pyrifera':'M. pyrifera'}
annual = read('data/harmonized/harmonized_response_ratios.csv') + read('data/harmonized/harmonized_landsat_rr.csv')
for r in annual:r['key']='|'.join((r['CA_MPA_Name_Short'],short.get(r['y'],''),r['source'],r['resp']))
checked=0
for row in primary:
    if row['Model']=='Sigmoid':continue
    d=sorted((r for r in annual if r['key']==key(row)),key=lambda x:int(x['year']))
    y=[float(r['lnDiff']) for r in d];n=len(y)
    assert n==int(row['N'])
    if row['Model']=='Mean':
        est=sum(y)/n; se=math.sqrt(sum((v-est)**2 for v in y)/(n-1)/n)
    else:
        if row['Model']=='Step':x=[float(r['BA']=='After') for r in d];h=1
        else:
            supplied=row['AnalysisType']=='CI' or row['Source']=='Landsat' or key(row)=='Campus Point SMCA|P. interruptus|LTER|Bio'
            nb=sum(r['BA']=='Before' for r in d)
            x=[float(r['time']) for r in d] if supplied else [0.]*nb+list(range(n-nb));h=11
        xm=sum(x)/n;ym=sum(y)/n;sxx=sum((v-xm)**2 for v in x)
        slope=sum((a-xm)*(b-ym) for a,b in zip(x,y))/sxx
        residual=sum((b-ym-slope*(a-xm))**2 for a,b in zip(x,y))
        est=h*slope;se=h*math.sqrt(residual/(n-2)/sxx)
    assert abs(est-float(row['Mean']))<1e-7,(key(row),est)
    assert abs(se-float(row['SE']))<1e-7,(key(row),se)
    checked+=1
for r in read('tables/table_s_routing_clock_inputs.csv'):
    assert hashlib.md5((ROOT/r['path']).read_bytes()).hexdigest()==r['md5'],r['path']
for r in read('tables/table_s_routing_clock_matched_clock_difference.csv'):
    assert r['clock_direction_changed']=='FALSE' and r['clock_bh_changed']=='FALSE'
print(f'PASS: {len(scenarios)} fixed 145-history scenarios; {checked} independent OLS reconstructions; current input hashes; matched-clock conclusions unchanged.')
