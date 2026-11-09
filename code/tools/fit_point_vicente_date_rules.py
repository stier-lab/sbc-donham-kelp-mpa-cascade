from pathlib import Path
import pandas as pd,numpy as np
O=Path('outputs/point_vicente_seasons_20261101');rules=pd.read_csv(O/'date_rule_grid.csv');z=pd.read_csv(O/'grid_footprints.csv');a=pd.read_csv(O/'annual_input.csv');a=a[~a.scenario.str.contains('gap90')];local=pd.read_csv(O/'local_annual.csv');local=local[~local.scenario.str.contains('gap90')];base=a[a.scenario=='Previous_review'];offsets=pd.read_csv(O/'reconstruction.csv').set_index(['Taxa','Resp']).offset;parts=[];foot=[];role=[]
for r in rules[(rules.max_gap_rule_days==90)&(rules.n_references>0)].itertuples():
 name=f'{r.season}_{r.zones}_gap90_full11';q=z[z.scenario==f'{r.season}_{r.zones}_full11'];q=q[(q.status=='mpa')|q.site.isin(r.references.split(';'))].copy();q['scenario']=name;foot.append(q);rr=q.groupby(['site','status','year','Taxa','Resp'],as_index=False).value.mean().groupby(['status','year','Taxa','Resp'],as_index=False).value.mean();rr['scenario']=name;role.append(rr)
 for (tx,resp),g in rr.groupby(['Taxa','Resp']):
  h=base[(base.MPA=='Point Vicente SMCA')&(base.Source=='PISCO')&(base.Taxa==tx)&(base.Resp==resp)];v=g.pivot(index='year',columns='status',values='value');p=pd.DataFrame([h.iloc[0].to_dict() for _ in v.index]);p.year=v.index.to_numpy();p.t=p.year-p.MPA_Start;p.time=p.t;p.BA='After';p.lnDiff=np.log((v.mpa.to_numpy()+offsets[(tx,resp)])/(v.reference.to_numpy()+offsets[(tx,resp)]));p.scenario=name;p.baseline_design_verified=False;parts.append(p)
patch=pd.concat(parts);pd.concat([local,patch]).to_csv(O/'local_annual.csv',index=False)
for sc,p in patch.groupby('scenario'):
 q=base[~(base.series_id.isin(p.series_id)&(base.t>0))];q=pd.concat([q,p],ignore_index=True);q.scenario=sc;a=pd.concat([a,q],ignore_index=True)
a.to_csv(O/'annual_input.csv',index=False);pd.concat(foot).to_csv(O/'date_matched_footprints.csv',index=False);pd.concat(role).to_csv(O/'date_matched_role_values.csv',index=False)
