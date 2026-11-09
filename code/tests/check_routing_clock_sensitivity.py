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
    assert len(subset) == len(keys) == 143
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
            x=[float(r['time']) for r in d]
            post=[int(r['year']) for r in d if r['BA']=='After']
            first=min(post);last=max(y for y in post if y<=first+11);h=last-first
        xm=sum(x)/n;ym=sum(y)/n;sxx=sum((v-xm)**2 for v in x)
        slope=sum((a-xm)*(b-ym) for a,b in zip(x,y))/sxx
        residual=sum((b-ym-slope*(a-xm))**2 for a,b in zip(x,y))
        est=h*slope;se=h*math.sqrt(residual/(n-2)/sxx)
    assert abs(est-float(row['Mean']))<1e-7,(key(row),est)
    assert abs(se-float(row['SE']))<1e-7,(key(row),se)
    checked+=1
for r in read('tables/table_s_routing_clock_inputs.csv'):
    assert hashlib.md5((ROOT/r['path']).read_bytes()).hexdigest()==r['md5'],r['path']
changes=[]
for r in read('tables/table_s_routing_clock_matched_clock_difference.csv'):
    assert r['clock_direction_changed']=='FALSE'
    assert abs(float(r['estimate_calendar'])-float(r['base_estimate_calendar']))<1e-7
    assert abs(float(r['q_calendar'])-float(r['base_q_calendar']))<1e-6
    changed=(float(r['q_calendar'])<.05)!=(float(r['q_native'])<.05)
    assert (r['clock_bh_changed']=='TRUE')==changed
    if changed:changes.append((r['Taxa'],r['Resp'],r['uncertainty']))
# Report threshold crossings; each flag is checked independently above.
print(f'PASS: {len(scenarios)} fixed 143-history scenarios; {checked} independent OLS reconstructions; current input hashes; calendar contrasts reproduce primary; threshold crossings versus sequential-index fitting: {changes}.')
