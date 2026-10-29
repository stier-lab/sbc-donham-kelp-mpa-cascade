"""Compare a completed singleton-biomass refit with its preserved checkpoint."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


def read(path):
    with path.open(newline='') as f:
        return list(csv.DictReader(f))


def same(a,b):
    try:
        return math.isclose(float(a),float(b),rel_tol=1e-10,abs_tol=1e-10)
    except ValueError:
        return a==b


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--before',type=Path,required=True)
    p.add_argument('--out',type=Path,default=Path('docs/BIOMASS_SINGLETON_REFIT_2026-10-18.json'))
    a=p.parse_args(); report={'primary_inputs':{},'primary_histories':{},'pooled':[]}
    for name in ('harmonized_raw_responses.csv','harmonized_response_ratios.csv'):
        before=read(a.before/'before_harmonized'/name)
        after=read(Path('data/harmonized')/name)
        assert len(before)==len(after)
        taxon='taxon_name' if 'raw' in name else 'y'
        columns=('value',) if 'raw' in name else ('mpa','reference','Diff','lnDiff')
        changes=[]
        for old,new in zip(before,after):
            assert all(old[k]==new[k] for k in old if k not in columns)
            if any(not same(old[k],new[k]) for k in columns):
                assert old['resp']=='Bio' and ((old['source']=='PISCO' and old[taxon]=='Panulirus interruptus') or
                    (old['source'] in ('PISCO','LTER') and old[taxon] in ('Mesocentrotus franciscanus','Strongylocentrotus purpuratus')))
                changes.append({'MPA':old['CA_MPA_Name_Short'],'source':old['source'],'year':old['year'],
                    'taxon':old[taxon],'status':old.get('status'),'before':{k:old[k] for k in columns},'after':{k:new[k] for k in columns}})
        report['primary_inputs'][name]={'rows':len(after),'changed_rows':len(changes),'changes':changes}
    keys=('MPA','Taxa','Source','Resp')
    index=lambda rs:{tuple(r[k] for k in keys):r for r in rs}
    before=index(read(a.before/'before_sumstats_final.csv'))
    after=index(read(Path('data/sumstats_final.csv')))
    assert before.keys()==after.keys() and len(after)==145
    changed=[]
    for key,r in after.items():
        old=before[key]
        if not same(old['Mean'],r['Mean']) or not same(old['SE'],r['SE']):
            assert key[3]=='Bio' and ((key[1]=='P. interruptus' and key[2]=='PISCO') or
                (key[1] in ('M. franciscanus','S. purpuratus') and key[2] in ('PISCO','LTER')))
            changed.append({'key':key,'old_model':old['Model'],'new_model':r['Model'],
                'old_mean':old['Mean'],'mean':r['Mean'],'old_se':old['SE'],'se':r['SE']})
    report['primary_histories']={'n':len(after),'changed':changed,'unchanged':len(after)-len(changed)}
    old={(r['Taxa'],r['Response']):r for r in read(a.before/'before_tables/table_02_meta_analysis.csv')}
    for r in read(Path('tables/table_02_meta_analysis.csv')):
        o=old[(r['Taxa'],r['Response'])]
        report['pooled'].append({'taxon':r['Taxa'],'response':r['Response'],
            'old_exp_contrast':float(o['RR']),'exp_contrast':float(r['RR']),
            'old_q':float(o['pval_fdr']),'q':float(r['pval_fdr']),
            'direction_changed':(float(o['Estimate'])>0)!=(float(r['Estimate'])>0),
            'detection_changed':(float(o['pval_fdr'])<.05)!=(float(r['pval_fdr'])<.05)})
    for name in ('table_s_historical_heatwave_summary.csv','table_s_historical_heatwave_pairs.csv'):
        assert read(a.before/'before_tables'/name)==read(Path('tables')/name),name
    report['historical_kelp_unchanged']=True
    hac=read(Path('tables/table_s_primary_hac_pooled.csv'))
    assert len(hac)==36 and all(r['direction_changed']=='FALSE' for r in hac)
    report['hac_detection_changes']=[{k:r[k] for k in ('Taxa','Resp','scenario','q','q_primary')}
        for r in hac if r['bh_detection_changed']=='TRUE']
    def kelp_rows(path):return [r for r in read(path) if any('Macrocystis' in str(v) or 'M. pyrifera' in str(v) for v in r.values())]
    assert kelp_rows(a.before/'before_tables/table_resistance_recovery.csv')==kelp_rows(Path('tables/table_resistance_recovery.csv'))
    report['modern_kelp_unchanged']=True
    report['files']={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in Path('data/harmonized').glob('*.csv')}
    a.out.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('primary_inputs','files')},indent=2))


if __name__=='__main__':main()
