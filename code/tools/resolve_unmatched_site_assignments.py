"""Resolve inventory-only unmatched assignments through explicit source-map hierarchy."""
from pathlib import Path
import pandas as pd,numpy as np,re,json,hashlib,html
O=Path('outputs/site_assignment_resolution_20261101');D=Path('/Users/adrianstier/Donham-Stier-CA-MPA-Data-2026/data');I=Path('outputs/dataset_inventory_20261101');used={}
def rd(p,**kw):
 p=Path(p);used[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest();return pd.read_csv(p,**kw)
a=rd(I/'alias_evidence.csv').fillna({'MPA':''});m=rd(D/'PISCO/master_site_table_Emilyedit.csv').fillna({'CA_MPA_Name_Short':'','site_status':''});m['site_original']=m.site;m['site']=m.site.str.strip();m['site_new']=m.site_new.str.strip()
r=rd(I/'raw_station_response_scan.csv');inventory=rd(I/'site_response_coverage.csv');canonref=rd('data/harmonized/reference_sites.csv');season=rd('tables/table_s_seasonal_site_months.csv')
u=D.parent/'code/R/01_utils.R';used[str(u)]=hashlib.sha256(u.read_bytes()).hexdigest();ut=u.read_text()
def constants(name):return set(re.findall(r'"([^"]+)"',re.search(name+r' <- c\((.*?)\)',ut,re.S).group(1)))
excluded_refs=constants('EXCLUDED_REFERENCE_SITES');excluded_mpas=constants('EXCLUDED_MPAS');excluded_kfm=constants('EXCLUDED_KFM_SITES')
# Verify underlying source rows directly, including old-name columns and complete dataset identity.
raw={}
for did,file in [('pisco_fish','PISCO/MLPA_kelpforest_fish_2024.csv'),('pisco_swath','PISCO/MLPA_kelpforest_swath_2024.csv'),('invertebrate_sizes','ALL_sizefreq_2024.csv')]:
 z=rd(D/file,usecols=['site','site_name_old','year','month','day']+(['site_new'] if did=='invertebrate_sizes' else []),encoding='latin1',low_memory=False);raw[did]=z
rows=[];links=[]
for ar in a.itertuples():
 parent=ar.assigned_site;mpa=ar.MPA
 pm=m[(m.site==parent)&(m.CA_MPA_Name_Short==mpa)&m.campus.isin(['UCSB','VRG'])];assert len(pm)==1,(parent,mpa)
 role=pm.site_status.iloc[0]
 children=m[(m.site_new==parent)&m.site.ne(parent)&m.campus.isin(['UCSB','VRG'])&(m.CA_MPA_Name_Short==mpa)]
 # Different MPA assignments of the same raw station remain distinct.
 childnames=sorted(children.site.unique());assert all(children.site_status==role)
 observed=[];ys=set();count=0
 for child in childnames:
  actual=0
  for did,z in raw.items():
   q=z[z.site==child];y=sorted(q.year.dropna().astype(int).unique());actual+=len(q)
   if did!='invertebrate_sizes':ys.update(y);count+=len(q)
   links.append(dict(MPA=mpa,group_name=parent,raw_station=child,campus='UCSB/VRG',dataset_id=did,raw_rows=len(q),first_year=min(y) if y else np.nan,last_year=max(y) if y else np.nan,years=';'.join(map(str,y)),explicit_reference_exclusion=child in excluded_refs,whole_mpa_exclusion=mpa in excluded_mpas,already_in_station_inventory=bool(((inventory.program=='PISCO')&(inventory.site_id==child)&(inventory.MPA.fillna('')==mpa)).any()),canonical_count_month_evidence=bool(((season.Source=='PISCO')&(season.site==child)&(season.MPA==mpa)).any()),mapping_evidence='Exact master.site_new = parent; master.site = raw child; same MPA/role, UCSB/VRG only'))
  if actual:observed.append(child)
 k=m[(m.site_new==parent)&(m.campus=='KFM')]
 kcodes=sorted(k.site.unique());kexc=[s for s in kcodes if s in excluded_kfm]
 direct=sum((z.site==parent).sum() for z in raw.values())
 if childnames:
  assert observed and direct==0
  identity='Confirmed grouped-site identifier; observed constituent stations'
  if mpa in ['Anacapa Island SMCA','Anacapa Island SMR 1978']:
   use='Real whole-MPA exclusion';why='Anacapa SMCA excluded by retained ecological/researcher scope; 1978 reserve outside the new-protection cohort.' if mpa=='Anacapa Island SMCA' else 'Earlier 1978 protection; excluded from the new-protection cohort. Same reference records can separately support the 2003 comparison.'
  elif parent=='SCI_YELLOWBANKS':
   use='Real reference exclusion / short remaining branch';why='CEN/W explicitly excluded; Hamilton describes poor local Gull match. E is not in the explicit list but has only two raw survey years; do not treat name resolution as permission to restore it.'
  elif parent=='SMI_TYLER_BIGHT':
   use='Real reference exclusion';why='Both E/W explicitly excluded. Hamilton confirms southern-island habitat differs from Harris/Cuyler; data are not absent or poor quality.'
  else:
   use='No whole-group exclusion found for this PISCO assignment';why='Group label resolves to existing station records. Preserve individual stations, years and response-specific eligibility; no newly recovered observations or automatic additions.'
 else:
  identity='Unassigned map-only label; no source records recovered';use='Not an eligible study assignment; reason for map entry unresolved';why='No MPA, role or coordinates in source row; no observed child station or exact raw record. Do not invent an MPA or claim a documented ecological exclusion.'
 rows.append(dict(MPA=mpa or 'Unassigned',assigned_group=parent,role=role or 'Unassigned',identity_class=identity,analysis_status=use,raw_children=';'.join(childnames),observed_children=';'.join(observed),n_children=len(childnames),raw_count_survey_rows=count,first_raw_year=min(ys) if ys else np.nan,last_raw_year=max(ys) if ys else np.nan,n_raw_years=len(ys),direct_parent_records=int(direct),linked_kfm_codes=';'.join(kcodes),linked_kfm_codes_in_exclusion_list=';'.join(kexc),reason=why,action='Link inventory group to constituent stations; do not merge or duplicate raw rows' if childnames else 'Leave unassigned; ask original mapping custodian only if these sites become relevant'))
out=pd.DataFrame(rows);lk=pd.DataFrame(links);out.to_csv(O/'assignment_resolution.csv',index=False);lk.to_csv(O/'constituent_station_evidence.csv',index=False)
# Exact/raw and normalized legacy-name absence for the two unassigned labels.
remaining=[]
for par in ['SCI_TRAILER','SCI_WILLOWS']:
 for did,z in raw.items():
  norm=lambda x:re.sub('[^A-Z0-9]','',str(x).upper())
  cols=['site','site_name_old']+(['site_new'] if 'site_new' in z else []);hits=z[z[cols].astype(str).apply(lambda c:c.str.contains('TRAILER' if par=='SCI_TRAILER' else 'WILLOW',case=False,regex=False)).any(axis=1)]
  remaining.append(dict(label=par,dataset_id=did,exact_parent_rows=int((z.site==par).sum()),name_fragment_rows=len(hits),matched_names=';'.join(sorted(hits.site.astype(str).unique()))))
pd.DataFrame(remaining).to_csv(O/'remaining_name_searches.csv',index=False)
# Do not conflate a PISCO/KFM geographic crosswalk with shared survey observations.
assert len(out)==21 and out.assigned_group.nunique()==19
assert (out.n_children>0).sum()==19 and (out.n_children==0).sum()==2
assert lk.already_in_station_inventory.all()
assert out.direct_parent_records.eq(0).all()
for file,digest in used.items():assert hashlib.sha256(Path(file).read_bytes()).hexdigest()==digest,file
protected=json.loads(Path('outputs/researcher_scope_20261101/protected_hashes.json').read_text())
for file,digest in protected.items():assert hashlib.sha256(Path(file).read_bytes()).hexdigest()==digest,file
for file,digest in json.loads((I/'manifest.json').read_text()).items():assert hashlib.sha256(Path(file).read_bytes()).hexdigest()==digest,file
checks=dict(assignments=21,unique_group_labels=19,resolved_group_assignments=19,unique_resolved_groups=17,remaining_unassigned_labels=2,real_excluded_group_assignments=int(out.analysis_status.str.startswith('Real').sum()),not_whole_group_excluded_assignments=int(out.analysis_status.str.startswith('No whole').sum()),unique_observed_pisco_child_stations=lk[lk.raw_rows>0].raw_station.nunique(),child_dataset_assignment_checks=len(lk),all_children_already_in_inventory=True,protected_files_unchanged=len(protected),previous_inventory_unchanged=True,models_fitted=0)
(O/'verification.json').write_text(json.dumps(checks,indent=2)+'\n');(O/'source_hashes.json').write_text(json.dumps(used,indent=2)+'\n')
def tab(d):return '| '+' | '.join(d.columns)+' |\n| '+' | '.join(['---']*len(d.columns))+' |\n'+''.join('| '+' | '.join(str(x).replace('|',' / ') for x in row)+' |\n' for row in d.itertuples(index=False,name=None))
text='''# Resolution of the 21 unmatched site assignments

November 1, 2026. Review only; paper, source data, models, reference choices and previous inventory remain unchanged.

## Main finding

Nineteen assignment rows (17 distinct group labels) are explicit grouped-site identifiers in the source master table. They link through `site_new` to the east/west/center/north/south station codes actually recorded in the fish and swath files. These are one-to-many group links, not interchangeable names for one physical station. They should be navigation/grouping entries, not additional sampled stations or missing datasets.

The previous audit checked direct station IDs and legacy-name fields but did not follow the master table's explicit `site_new` hierarchy. That made these 19 assignments appear unresolved. No observations were recovered or added by resolving the hierarchy. All linked PISCO constituent stations were already present in the inventory.

The other two labels, SCI_TRAILER and SCI_WILLOWS, are unassigned map entries: no MPA, status or coordinates, no constituent station mapping, and no exact or name-fragment records in the supplied fish, swath or size sources. They are unresolved historical/map labels, not confirmed ecological exclusions and not evidence that surveys never occurred.

## Identity and exclusion are separate

Of the19 mapped assignments,14 have no whole-group PISCO exclusion identified;5 are associated with real exclusions: Anacapa West Isle in the SMCA; East Fish Camp and Lighthouse Reef in the1978 reserve assignment; Yellowbanks; and Tyler Bight. East Fish Camp and Lighthouse Reef each also appear in the2003 assignment, where the same reference data have a different comparison role. Counts therefore refer to21 assignments, not21 independent sites.

Yellowbanks CEN/W and Tyler Bight E/W are explicitly excluded in source code. Yellowbanks E has observed records in only two years and is not in that explicit exclusion vector; the ecological advice still does not justify restoring it as a local Gull reference. Anacapa SMCA and the1978 reserve are whole-comparison exclusions, not missing records. Finding source data does not reverse any exclusion.

KFM codes also map to many of these geographic group names. These remain a separate network. KFM exclusions at Cavern Point, Coche Point, Little Scorpion, Potato Pasture and San Pedro Point must not be transferred automatically to PISCO constituent stations. Jenn's October19 reply explicitly distinguishes separate KFM and PISCO observations even at nearby sites/dates. Geographic links do not authorize pooling or deduplicating across programs. Shared size samples remain a separate, previously verified issue.

## Every assignment

'''+tab(out[['MPA','assigned_group','raw_children','analysis_status']])+'''
## Evidence and remaining questions

- Original master site table: exact `site_new` group links, retaining program, MPA and role. Leading/trailing whitespace is normalized only for map matching; original source fields remain unchanged.
- Raw fish, swath and size files: constituent record counts and observed years independently checked; all source hashes retained. See constituent_station_evidence.csv. These are raw years before seasonal/response filtering, not final analysis coverage.
- Source exclusion vectors: PISCO references, whole MPAs and KFM stations are recorded separately. No sign/significance rule is introduced.
- Scott Hamilton's October29 email supports retaining the Tyler Bight and Yellowbanks habitat-based exclusions despite the existence of useful data.
- Emily's October29 email describes intended before/after station continuity, with fish and benthic coverage allowed to differ. It does not explain every KFM station-specific omission.
- Jenn's October19 email says KFM and PISCO survey observations are separate. Her reply does not establish that their physical footprints are identical.
- Targeted coauthor email searches found no explanation of Trailer or Willows and no contrary instruction to collapse these group labels. Search absence is not evidence that no private discussion occurred. Original email bodies remain private; durable IDs are recorded below without copying correspondence.

Only Trailer and Willows still need an identity/history answer: what did those labels refer to, and were observations ever collected under another name? If worth pursuing, ask the PISCO mapping/data custodian (Jenn or the current data manager); do not ask coauthors to reconfirm the19 source-resolved group links. The missing count/length/frond measurements identified in other audits remain separate; this name audit did not recover them.

## Files and reproduction

- assignment_resolution.csv: all21 assignments, identity classification, separate exclusion status and recommended treatment.
- constituent_station_evidence.csv: every exact group-to-station link by source, raw years, exclusion flags and existing inventory membership.
- remaining_name_searches.csv: bounded raw-name searches for Trailer and Willows.
- verification.json / source_hashes.json: assertions, source fingerprints and preservation checks.

Run `python3 code/tools/resolve_unmatched_site_assignments.py`. The next inventory revision should represent parents as group records linked to station records, preserving response-specific coverage. No analytical joins should substitute a parent for a constituent station or count both as samples.

Private-source message IDs checked: Scott Hamilton 1a188f9a2394d9e3; Emily Donham 1a189059ab077678; Jenn Caselle 1a15475000fd8c4a. No emails sent.
'''
text=text.replace('the19','the 19').replace('assignments,14','assignments, 14').replace(';5','; 5').replace('the1978','the 1978').replace('the2003','the 2003').replace('to21','to 21').replace('not21','not 21').replace('all21','all 21').replace('October19','October 19').replace('October29','October 29')
(O/'README.md').write_text(text)
page='<!doctype html><meta charset="utf-8"><title>Site assignment resolution</title><style>body{font:16px system-ui;margin:24px;color:#172e34}table{border-collapse:collapse;font-size:14px}td,th{padding:10px;border-bottom:1px solid #ddd;text-align:left;vertical-align:top}th{background:#eef3f4}input{padding:10px;font:inherit;width:360px}.scroll{overflow:auto}p{max-width:1000px}</style><h1>21 unmatched assignments: resolved</h1><p>19 assignments are explicit group names for existing stations; 2 remain unassigned map labels. Five of the mapped assignments have separate exclusion reasons. Paper and analyses unchanged.</p><p>Group names are not additional survey stations. Keep KFM and PISCO observations separate.</p><input id="q" aria-label="Filter assignments" placeholder="Filter by site, MPA or status"><div class="scroll">'+out[['MPA','assigned_group','identity_class','raw_children','analysis_status','reason']].to_html(index=False,escape=True,border=0)+'</div><p><a href="README.md">Full evidence and limitations</a> · <a href="assignment_resolution.csv">Download all rows</a></p><script>document.getElementById("q").oninput=function(){let q=this.value.toLowerCase();for(let r of document.querySelectorAll("tbody tr"))r.hidden=!r.textContent.toLowerCase().includes(q)}</script>'
(O/'index.html').write_text(page);print(checks);print(out[['assigned_group','analysis_status']].to_string(index=False))
