from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,shutil
parent=Path('outputs/diver_exclusion_refit_20261101');out=parent/'consistent_kfm';out.mkdir(exist_ok=True)
d=pd.read_csv(parent/'annual_input.csv');r=pd.read_csv('data/harmonized/harmonized_raw_responses.csv');mask=(d.Source=='KFM')&(d.Taxa=='M. pyrifera')&(d.Resp=='Bio');updates=[];audit=[]
for m,z in d[mask].groupby('MPA'):
 raw=r[(r.source=='KFM')&(r.CA_MPA_Name_Short==m)&(r.taxon_name.str.strip()=='Macrocystis pyrifera')&(r.resp=='Bio')].copy();assert not raw.duplicated(['year','status']).any()
 mx=raw.value.max();floor=raw.value[raw.value>0].min()/2;assert mx>0 and floor>0
 raw['prop']=(raw.value+floor)/mx
 wide=raw.pivot(index='year',columns='status',values='prop').dropna().reset_index();assert set(wide.year)<=set(z.year)
 x=z.set_index('year').reindex(wide.year).reset_index();x['mpa']=wide.mpa.to_numpy();x['reference']=wide.reference.to_numpy();x['Diff']=x.mpa/x.reference;x['lnDiff']=np.log(x.Diff);updates.append(x)
 audit.append(dict(MPA=m,old_rows=len(z),new_rows=len(x),old_before=int((z.t<0).sum()),new_before=int((x.t<0).sum()),post=int(x.t.between(1,9).sum()),removed_years=';'.join(map(str,sorted(set(z.year)-set(x.year)))),max_common_year_difference=float(max(abs(x.set_index('year').lnDiff-z.set_index('year').loc[x.year].lnDiff))),raw_additive_floor=floor))
full=pd.concat([d[~mask],*updates],ignore_index=True);assert not full.duplicated(['series_id','year']).any();full.to_csv(out/'annual_input.csv',index=False);pd.DataFrame(audit).to_csv(out/'kfm_protocol_audit.csv',index=False)
for f in ['reference_clusters.csv','protected_hashes.json']:shutil.copy(parent/f,out/f)
(out/'PLAN.md').write_text((parent/'PLAN.md').read_text()+'\nUniform KFM kelp source reconciliation: all six KFM kelp biomass histories, both retained and restored, use current raw paired annual biomass from the modern quadrat protocol. Recompute the same half-minimum-positive correction from each full source history. This isolates a source/protocol repair from exclusion membership. Do not use the sign or significance of any new estimate. All other annual histories are unchanged.\n')
files=[out/'annual_input.csv',out/'reference_clusters.csv',out/'PLAN.md',out/'kfm_protocol_audit.csv',Path('data/harmonized/harmonized_raw_responses.csv')]
(out/'freeze_manifest.json').write_text(json.dumps({str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in files},indent=2))
print(pd.DataFrame(audit).to_string(index=False))
