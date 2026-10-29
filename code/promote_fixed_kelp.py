#!/usr/bin/env python3
"""Promote only kelp biomass from a staged full harmonization; preserve other rows."""
import argparse,csv,hashlib,io,json,math
from pathlib import Path
FILES=('harmonized_raw_responses.csv','harmonized_response_ratios.csv')
SOURCE_MD5='8ddaf77458d21adfb6d0398cc1c51fda'
def md5(p):return hashlib.md5(p.read_bytes()).hexdigest()
def promote(baseline,candidate):
    before=list(csv.DictReader(io.StringIO(baseline.decode())))
    after=list(csv.DictReader(io.StringIO(candidate.decode())))
    tax='taxon_name' if 'taxon_name' in before[0] else 'y'
    keys=['CA_MPA_Name_Short','year',tax,'source','resp']+(['status'] if tax=='taxon_name' else [])
    fields=['value'] if tax=='taxon_name' else ['mpa','reference','Diff','lnDiff']
    key=lambda r:tuple(r[k] for k in keys)
    lookup={key(r):r for r in after}
    assert len(lookup)==len(after)==len(before) and set(lookup)=={key(r) for r in before},'Record identities changed'
    lines=baseline.decode().splitlines(keepends=True);assert len(lines)==len(before)+1
    changed=[];other_differences=[]
    for i,r in enumerate(before,1):
        n=lookup[key(r)];assert set(n)==set(r)
        target=r[tax]=='Macrocystis pyrifera' and r['resp']=='Bio' and r['source'] in ('KFM','PISCO','LTER')
        assert all(r[k]==n[k] for k in r if k not in fields),'Metadata changed'
        different=any(not math.isclose(float(r[f]),float(n[f]),rel_tol=1e-10,abs_tol=1e-10) for f in fields)
        if not target:
            if different:other_differences.append(key(r))
            continue
        if different:
            new=dict(r)
            for f in fields:
                assert math.isfinite(float(n[f]));new[f]=float(n[f])
            s=io.StringIO();w=csv.DictWriter(s,fieldnames=list(r),lineterminator='\n',quoting=csv.QUOTE_NONNUMERIC);w.writerow(new)
            lines[i]=s.getvalue();changed.append({'key':key(r),'before':{f:r[f] for f in fields},'after':{f:n[f] for f in fields}})
    assert not other_differences, f'Unexpected non-kelp-biomass changes: {other_differences[:8]}'
    return ''.join(lines).encode(),changed

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--baseline',type=Path,required=True);ap.add_argument('--candidate',type=Path,required=True);ap.add_argument('--target',type=Path,action='append',required=True);ap.add_argument('--calibration',type=Path,required=True);ap.add_argument('--apply',action='store_true');a=ap.parse_args()
    assert md5(a.calibration)==SOURCE_MD5
    report={'release':'kelp_fixed_sbc1275_2026_10_29','package':'knb-lter-sbc.127.5','model':1,'calibration_md5':SOURCE_MD5,'mass_basis':'dry foliar','conversion':'recorded month before annual aggregation','applied':a.apply,'files':{}}
    staged={}
    for name in FILES:
        old=(a.baseline/name).read_bytes();new,changes=promote(old,(a.candidate/name).read_bytes());staged[name]=new
        report['files'][name]={'before_md5':hashlib.md5(old).hexdigest(),'after_md5':hashlib.md5(new).hexdigest(),'rows':len(list(csv.DictReader(io.StringIO(new.decode())))),'changed_rows':len(changes),'changed_records':changes}
        for t in a.target:assert (t/name).read_bytes()==old,f'Destination changed: {t/name}'
    if a.apply:
        for t in a.target:
            for name,b in staged.items():(t/name).write_bytes(b)
            (t/'kelp_fixed_calibration.json').write_text(json.dumps(report,indent=2)+'\n')
            (t/'checksums.md5').write_text(''.join(md5(p)+' '+p.name+'\n' for p in sorted(t.glob('*.csv'))))
    print(json.dumps({k:v for k,v in report.items() if k!='files'}));print(json.dumps({n:{k:v for k,v in z.items() if k!='changed_records'} for n,z in report['files'].items()},indent=2))
if __name__=='__main__':main()
