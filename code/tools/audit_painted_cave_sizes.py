from pathlib import Path
import pandas as pd, json,hashlib
O=Path('outputs/pisco_painted_cave_20261102');O.mkdir(exist_ok=True)
B=Path('/Users/adrianstier/Library/CloudStorage/GoogleDrive-<account>/My Drive/Stier Lab/People/Emily Donham/Projects/Kelp MPA/Github/MPA_Kelp')
old=pd.read_csv('data/ALL_sizefreq_2024.csv',low_memory=False);drive=pd.read_csv(B/'data/ALL_sizefreq_2024.csv',low_memory=False)
sites=pd.read_csv(B/'data/PISCO/master_site_table_Emilyedit.csv');keys=['site','site_new','campus'];meta=sites.drop_duplicates(['site','CA_MPA_Name_Short'])
x=old[old.campus.eq('UCSB')&old.method.eq('SBTL_SIZEFREQ_PISCO')&old.classcode.isin(['MESFRA','STRPUR','PANINT'])].copy();x['audit_row']=x.index+1
cross=sites[['site','campus','site_new']].drop_duplicates();assert not cross.duplicated(['site','campus']).any()
y=x.rename(columns={'site_new':'old_site_new'}).merge(cross,on=['site','campus'],how='left',validate='many_to_one');assert y.site_new.notna().all()
changed=y[y.old_site_new!=y.site_new].copy();assert len(changed)==99;assert changed.site.eq('SCI_PAINTED_CAVE_E').all();assert changed.site_new.eq('SCI_PAINTED_CAVE').all()
assert drive.loc[changed.audit_row-1,'site_new'].eq('SRI_SOUTH_POINT').all()
oldj=x.merge(meta,on=keys,how='left');newj=y.merge(meta,on=keys,how='left');cid=set(changed.audit_row)
a=oldj[oldj.audit_row.isin(cid)];b=newj[newj.audit_row.isin(cid)]
assert a.CA_MPA_Name_Short.isna().all();assert b.CA_MPA_Name_Short.eq('Painted Cave SMCA').all()
assert not b.CA_MPA_Name_Short.eq('South Point SMR').any()
# Exact retained membership table: correction can only restore excluded Painted Cave membership.
cols=['audit_row','CA_MPA_Name_Short','site_status']
def retained(z):return z[z.CA_MPA_Name_Short.notna()&~z.CA_MPA_Name_Short.eq('Painted Cave SMCA')][cols].sort_values(cols).reset_index(drop=True)
pd.testing.assert_frame_equal(retained(oldj),retained(newj))
changed['decision']='excluded_Painted_Cave_ecological_scope';changed.to_csv(O/'corrected_records.csv',index=False)
summary=changed.groupby(['site','old_site_new','site_new','year','classcode']).agg(rows=('count','size'),animals=('count','sum')).reset_index();summary.to_csv(O/'correction_summary.csv',index=False)
pd.concat([a.assign(version='legacy_join'),b.assign(version='corrected_join')])[['version','audit_row','site','site_new','CA_MPA_Name_Short','site_status']].to_csv(O/'membership_changes.csv',index=False)
report={'corrected_rows':len(changed),'animals':int(changed['count'].sum()),'legacy_matched_rows':int(a.CA_MPA_Name_Short.notna().sum()),'corrected_Painted_Cave_rows':len(b),'corrected_South_Point_rows':0,'retained_donor_memberships_changed':0,'paper_changed':False,'primary_models_rerun':False,'source_hashes':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [B/'code/pBACIPS_PISCO_V10.R',B/'data/ALL_sizefreq_2024.csv',B/'data/PISCO/master_site_table_Emilyedit.csv']}}
(O/'verification.json').write_text(json.dumps(report,indent=2)+'\n');print(summary.to_string(index=False));print(json.dumps(report,indent=2))
