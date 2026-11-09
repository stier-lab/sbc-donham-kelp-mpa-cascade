#!/usr/bin/env python3
"""Read-only input, diagnostic-only borrowing rebuild. Requires numpy/pandas/scipy.
Usage: python rebuild_borrowing_support.py --repo REPO --output SCRATCH
BORROWING_OUTPUT may replace --output. No existing analysis script is imported.
Tolerance relaxes ONLY the upper calibration gap, uniformly for both donor kinds.
"""
from pathlib import Path
import argparse, hashlib, json, os
from types import SimpleNamespace
import numpy as np
import pandas as pd
from scipy.stats import t as student

KEY = ['series_id', 'site', 'status', 'year', 'zone', 'date']
PHYSICAL = ['site', 'year', 'zone', 'date', 'Taxa']
KINDS = ['other_visit', 'other_site']
TOLS = [0, 3, 7, 14]

def digest(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1024*1024), b''): h.update(b)
    return h.hexdigest()

def key(row, columns=KEY):
    return json.dumps([str(getattr(row, k)) for k in columns], separators=(',', ':'))

def qweighted(z):
    w = 1/z.cluster.map(z.groupby('cluster').size()).to_numpy()
    ix = np.argsort(z.residual.to_numpy())
    w, r = w[ix], z.residual.to_numpy()[ix]
    return np.interp([.025, .5, .975], (np.cumsum(w)-.5*w)/w.sum(), r)

def fit(z, col='biomass_lnratio'):
    z = z[np.isfinite(z[col])].sort_values('year')
    n = len(z)
    if n < 3:
        return dict(estimate=np.nan, se=np.nan, low=np.nan, high=np.nan,
                    n_years=n, df=max(n-2, 0), fit_status='fewer_than_3_years')
    years, y = z.year.to_numpy(int), z[col].to_numpy(float)
    X = np.column_stack([np.ones(n), z.t])
    inv = np.linalg.inv(X.T@X)
    b = inv@X.T@y
    score = X*(y-X@b)[:, None]
    scores = np.zeros((years.max()-years.min()+1, 2))
    for yr in set(years): scores[yr-years.min()] = score[years==yr].sum(axis=0)
    meat = scores.T@scores
    lag = min(3, len(scores)-1)
    for h in range(1, lag+1):
        cross = scores[h:].T@scores[:-h]
        meat += (1-h/(lag+1))*(cross+cross.T)
    cov = inv@meat@inv*n/(n-2)
    se = np.sqrt(max(0, cov[1,1]))
    return dict(estimate=b[1], se=se, low=b[1]-student.ppf(.975,n-2)*se,
                high=b[1]+student.ppf(.975,n-2)*se, n_years=n, df=n-2,
                fit_status='estimated', hac_calendar_lag=lag)

def aggregate(v, base, offsets):
    """Equal visits within zones, equal zones, equal sites; NEVER skip unknowns."""
    rows=[]
    for b in base.itertuples():
        z=v[(v.series_id==b.series_id)&(v.year==b.year)]
        roles={}
        for role in ['mpa','reference']:
            sites=[]
            for _, ss in z[z.status==role].groupby('site'):
                zones=[zz.biomass_point.to_numpy().mean() for _,zz in ss.groupby('zone')]
                sites.append(np.mean(zones))
            roles[role]=np.mean(sites) if sites else np.nan
        pp, rr = roles['mpa'], roles['reference']
        ln=np.log((pp+offsets[b.series_id])/(rr+offsets[b.series_id]))
        rows.append(dict(series_id=b.series_id,year=b.year,t=b.t,
            protected_biomass=pp,reference_biomass=rr,biomass_lnratio=ln,
            local_only_lnratio=b.revised,original_lnratio=b.original,
            supported=bool(np.isfinite(ln)),n_visits=len(z),
            n_reference_sites=z[z.status=='reference'].site.nunique(),
            n_protected_sites=z[z.status=='mpa'].site.nunique(),
            unknown_visits=int(z.biomass_point.isna().sum()),
            unknown_reference_visits=int(z.loc[z.status=='reference','biomass_point'].isna().sum())))
    return pd.DataFrame(rows)

class Engine:
    def __init__(self,s,v):
        self.s=s
        self.groups={k:g for k,g in s.groupby(['series_id','status','year','zone'])}
        self.cal={}
        rows=[]
        for sid,outer,status in v[['series_id','site','status']].drop_duplicates().itertuples(index=False,name=None):
            for kind in KINDS:
                records=[]
                for t in s[(s.series_id==sid)&(s.status==status)&(s.site!=outer)].itertuples():
                    d=self.donors(t,kind,outer)
                    if d.empty: continue
                    pred=self.predict(d)
                    rec=dict(series_id=sid,outer_site=outer,status=status,kind=kind,
                        cal_site=t.site,cal_year=t.year,cal_zone=t.zone,cal_date=t.date,
                        observed=t.mean_mass,predicted=pred,residual=np.log(t.mean_mass/pred),
                        gap_days=int(d.gap.max()),donor_sites=';'.join(sorted(d.site.unique())),
                        cluster=t.site+'|'+str(t.year),calibration_target_key=key(t),
                        donor_keys=json.dumps(d.sample_key.tolist()))
                    assert outer!=t.site and not (d.site==outer).any()
                    records.append(rec); rows.append(rec)
                self.cal[sid,outer,status,kind]=pd.DataFrame(records)
        self.calibration=pd.DataFrame(rows)

    def donors(self,t,kind,outer=None):
        d=self.groups.get((t.series_id,t.status,t.year,t.zone),self.s.iloc[:0])
        d=d[d.site!=outer] if outer is not None else d
        d=d[d.site!=t.site] if kind=='other_site' else d[(d.site==t.site)&(d.date!=t.date)]
        if d.empty: return d.copy()
        d=d.copy();d['gap']=(pd.to_datetime(d.date)-pd.Timestamp(t.date)).abs().dt.days
        d=d[d.gap==d.groupby('site').gap.transform('min')]
        assert (d.series_id==t.series_id).all() and (d.status==t.status).all()
        assert (d.year==t.year).all() and (d.zone==t.zone).all()
        return d

    @staticmethod
    def predict(d):
        g=d.groupby('site').agg(total_mass=('total_mass','sum'),n=('n_sized','sum'))
        return (g.total_mass/g.n).mean()

    def attempt(self,t,kind,tol):
        d=self.donors(t,kind);c=self.cal[t.series_id,t.site,t.status,kind]
        ng=0 if c.empty else c.cluster.nunique()
        gap=np.nan if d.empty else int(d.gap.max())
        cmax=np.nan if c.empty else int(c.gap_days.max())
        reason=('no_donors' if d.empty else 'fewer_than_two_calibration_groups' if ng<2
                else 'upper_calibration_gap_exceeded' if gap>cmax+tol else 'eligible')
        q=qweighted(c) if ng>=2 else [np.nan]*3
        pred=self.predict(d) if len(d) else np.nan
        out=dict(visit_key=key(t),kind=kind,priority=1 if kind=='other_visit' else 2,
            eligible=reason=='eligible',reason=reason,predicted=pred,
            mass=pred*np.exp(q[1]) if reason=='eligible' else np.nan,
            low=pred*np.exp(q[0]) if reason=='eligible' else np.nan,
            high=pred*np.exp(q[2]) if reason=='eligible' else np.nan,
            max_gap=gap,cal_min_gap=np.nan if c.empty else int(c.gap_days.min()),
            cal_max_gap=cmax,gap_excess=np.nan if c.empty or d.empty else max(0,gap-cmax),
            gap_below_cal_min=False if c.empty or d.empty else bool(gap<c.gap_days.min()),
            cal_samples=len(c),cal_clusters=ng,
            cal_unique_sites=0 if c.empty else c.cal_site.nunique(),
            cal_unique_years=0 if c.empty else c.cal_year.nunique(),
            cal_sites='' if c.empty else ';'.join(sorted(c.cal_site.unique())),
            cal_years='' if c.empty else ';'.join(map(str,sorted(c.cal_year.unique()))),
            cal_groups='' if c.empty else ';'.join(sorted(c.cluster.unique())),
            donor_sample_n=0 if d.empty else d.n_sized.sum(),donor_sample_groups=len(d),
            donor_singleton_groups=0 if d.empty else int(d.singleton.sum()),
            donor_zero_empirical_variance_groups=0 if d.empty else int(d.zero_empirical_variance.sum()),
            donors=';'.join(sorted(d.site.unique())),donor_keys=json.dumps(d.sample_key.tolist()),
            donor_evidence=json.dumps([dict(sample_key=x.sample_key,site=x.site,year=int(x.year),
                 zone=x.zone,date=x.date,n=int(x.n_sized),gap_days=int(x.gap),mean_mass=x.mean_mass,
                 singleton=bool(x.singleton),zero_empirical_variance=bool(x.zero_empirical_variance))
                 for x in d.itertuples()]),
            cal_residual_variance=np.nan if len(c)<2 else c.residual.var(ddof=1),
            cal_zero_empirical_variance=False if c.empty else bool(c.residual.nunique()==1))
        return out

def run(repo, data, output):
    repo=repo.resolve();data=data.resolve();output=output.resolve()
    if ((output==repo or repo in output.parents) and not output.is_relative_to(repo/'outputs/biomass_next_steps_20261102')) or output==data or data in output.parents:
        raise ValueError('Output must be outside inputs or inside outputs/biomass_next_steps_20261102')
    output.mkdir(parents=True,exist_ok=True)
    inputs={}
    def read(rel, root=repo, **kw):
        p=root/rel;inputs[str(p)]=digest(p);return pd.read_csv(p,**kw)
    B='outputs/pisco_biomass_borrowing_20261102/'
    v=read(B+'visit_biomass_support.csv')
    s=read('outputs/pisco_timing_donors_20261102/size_sample_inventory.csv')
    savedcal=read(B+'outer_site_calibration.csv')
    oldannual=read(B+'biomass_annual.csv');oldeffects=read(B+'biomass_effects.csv')
    rawannual=read('outputs/pisco_all_months_20261102/annual_values.csv')
    base=rawannual[(rawannual.scenario=='allmonths_sharedzones_visit_sizes')&rawannual.series_id.str.endswith('|Bio')].copy()
    density=read(B+'selected_density_annual.csv')
    off=read('outputs/all_program_station_refits_20261102/offsets.csv').set_index('series_id').offset
    assert not v.duplicated(KEY).any() and not s.duplicated(KEY).any()
    assert np.isfinite(v.density).all() and (v.density>=0).all()
    v['visit_key']=[key(t) for t in v.itertuples()]
    s['sample_key']=[key(t) for t in s.itertuples()]
    # Frequency-aware descriptive size variation; no invented singleton variance.
    raw=read('data/ALL_sizefreq_2024.csv',root=data,low_memory=False)
    codes={'MESFRA':'M. franciscanus','STRPUR':'S. purpuratus','PANINT':'P. interruptus'}
    r=raw[raw.campus.isin(['VRG','UCSB'])&raw.classcode.isin(codes)].copy()
    r['Taxa']=r.classcode.map(codes)
    r=r[(r['count']>0)&(r['size']>0)&((r.Taxa=='P. interruptus')|(r['size']>=25))].copy()
    r['date']=pd.to_datetime(r[['year','month','day']],errors='coerce').dt.strftime('%Y-%m-%d')
    r=r.merge(s[PHYSICAL].drop_duplicates(),on=PHYSICAL,validate='many_to_one')
    r['mass']=np.select([r.Taxa=='M. franciscanus',r.Taxa=='S. purpuratus'],
                       [.00059*r['size']**2.917,.00059*r['size']**2.870],
                       default=.001352821*r['size']**2.913963113)
    stat=[]
    for kk,g in r.groupby(PHYSICAL):
        n=g['count'].sum();mu=np.average(g.mass,weights=g['count'])
        ss=np.sum(g['count']*(g.mass-mu)**2)
        stat.append(dict(zip(PHYSICAL,kk),raw_n=n,raw_mean=mu,
            empirical_population_variance=ss/n,
            empirical_sample_variance=ss/(n-1) if n>1 else np.nan,
            singleton=bool(n==1),zero_empirical_variance=bool(g.mass.nunique()==1)))
    s=s.merge(pd.DataFrame(stat),on=PHYSICAL,validate='many_to_one')
    assert np.allclose(s.n_sized,s.raw_n) and np.allclose(s.mean_mass,s.raw_mean)
    assert set(s.sample_key)==set(key(t) for t in s.itertuples())
    eng=Engine(s,v)
    ck=['series_id','outer_site','status','kind','cal_site','cal_year','cal_zone','cal_date']
    pd.testing.assert_frame_equal(eng.calibration[savedcal.columns].sort_values(ck).reset_index(drop=True),
        savedcal.sort_values(ck).reset_index(drop=True),check_dtype=False,rtol=1e-10,atol=1e-10)
    local=s.set_index('sample_key')
    # Re-evaluate actual held-out targets under every tolerance; protected evidence stays separate.
    valid=[]
    for tol in TOLS:
        for t in v[v.density>0].itertuples():
            if t.visit_key not in local.index: continue
            for kind in KINDS:
                a=eng.attempt(t,kind,tol)
                valid.append(dict(series_id=t.series_id,status=t.status,site=t.site,year=t.year,
                    zone=t.zone,date=t.date,tolerance_days=tol,kind=kind,visit_key=t.visit_key,
                    eligible=a['eligible'],reason=a['reason'],observed=local.loc[t.visit_key,'mean_mass'],
                    predicted=a['predicted'],low=a['low'],high=a['high'],
                    covered=bool(a['low']<=local.loc[t.visit_key,'mean_mass']<=a['high'])))
    validation=pd.DataFrame(valid)
    def vstatus(sid,role,kind,tol):
        z=validation[(validation.series_id==sid)&(validation.status==role)&
            (validation.kind==kind)&(validation.tolerance_days==tol)&validation.eligible]
        p=validation[(validation.series_id==sid)&(validation.status=='mpa')&
            (validation.kind==kind)&(validation.tolerance_days==tol)&validation.eligible]
        return dict(validation_samples=len(z),validation_unique_sites=z.site.nunique(),
            validation_unique_years=z.year.nunique(),
            validation_groups=len(z[['site','year']].drop_duplicates()),
            protected_validation_samples=len(p),
            protected_validation_status='available_conditional' if len(p) else 'unavailable_no_eligible_protected_holdouts')
    visits=[];attempts=[];annual=[]
    for tol in TOLS:
        scenario='strict_reference' if tol==0 else f'upper_gap_plus_{tol}d_diagnostic'
        vv=[]
        for t in v.itertuples():
            row={k:getattr(t,k) for k in KEY+['density','row_id','visit_key']}
            row.update(scenario=scenario,tolerance_days=tol,reason='',sample_n=np.nan,
                singleton=False,zero_empirical_variance=False,empirical_sample_variance=np.nan,
                empirical_population_variance=np.nan,donor_sample_n=np.nan)
            if t.density==0: row.update(method='zero_count',mass=0.,biomass_point=0.)
            elif t.visit_key in local.index:
                m=local.loc[t.visit_key]
                row.update(method='measured_same_visit',mass=m.mean_mass,biomass_point=t.density*m.mean_mass,
                    sample_n=m.n_sized,**{x:m[x] for x in ['singleton','zero_empirical_variance',
                    'empirical_sample_variance','empirical_population_variance']})
            else:
                chosen=None;aa=[]
                # Evaluate BOTH, even when one succeeds: failures cannot be overwritten.
                for kind in KINDS:
                    a=eng.attempt(t,kind,tol)
                    a.update({k:getattr(t,k) for k in KEY+['row_id']})
                    a.update(scenario=scenario,tolerance_days=tol,**vstatus(t.series_id,t.status,kind,tol))
                    a['selected']=bool(chosen is None and a['eligible'])
                    a['attempt_key']=json.dumps([scenario,t.visit_key,kind],separators=(',',':'))
                    if a['selected']: chosen=a
                    aa.append(a);attempts.append(a)
                if chosen:
                    row.update({k:x for k,x in chosen.items() if k not in ['reason','kind','selected','priority','attempt_key']})
                    row.update(method='borrowed_'+chosen['kind'],biomass_point=t.density*chosen['mass'])
                else:
                    row.update(method='unsupported',mass=np.nan,biomass_point=np.nan,
                        reason=';'.join(a['kind']+':'+a['reason'] for a in aa))
                for a in aa:
                    for field in ['reason','cal_samples','cal_clusters','cal_unique_sites','cal_unique_years',
                                  'cal_min_gap','cal_max_gap','max_gap','protected_validation_status']:
                        row[a['kind']+'_'+field]=a[field]
            vv.append(row)
        vv=pd.DataFrame(vv)
        # Explicitly retain all original columns, with inapplicable metadata missing.
        for col in v.columns:
            if col not in vv: vv[col]=np.nan
        visits.append(vv)
        a=aggregate(vv,base,off);a['scenario']=scenario;a['tolerance_days']=tol;annual.append(a)
    visits=pd.concat(visits,ignore_index=True);attempts=pd.DataFrame(attempts);annual=pd.concat(annual,ignore_index=True)
    strict=visits[visits.tolerance_days==0].set_index('visit_key').loc[v.visit_key]
    assert (strict.method.to_numpy()==v.method.to_numpy()).all()
    for col in ['mass','biomass_point']:
        assert np.allclose(strict[col],v[col],equal_nan=True,rtol=1e-10,atol=1e-10),col
    stricta=annual[annual.tolerance_days==0].sort_values(['series_id','year']).reset_index(drop=True)
    olda=oldannual.sort_values(['series_id','year']).reset_index(drop=True)
    for col in olda.columns:
        if col in ['series_id','year','t','supported']: assert stricta[col].equals(olda[col]),col
        else: assert np.allclose(stricta[col],olda[col],equal_nan=True,rtol=1e-10,atol=1e-10),col
    # Full reference safeguard, and no zero substitution for unknown positives.
    assert visits.loc[visits.method=='unsupported','biomass_point'].isna().all()
    assert (visits.loc[visits.method=='zero_count','density']==0).all()
    assert (visits.loc[visits.density>0,'biomass_point'].dropna()>0).all()
    assert (annual.supported==(annual.unknown_visits==0)).all()
    assert annual.loc[annual.unknown_reference_visits>0,'reference_biomass'].isna().all()
    assert attempts.groupby(['scenario','visit_key']).size().eq(2).all()
    assert not attempts.duplicated('attempt_key').any()
    assert attempts.loc[attempts.selected,'cal_clusters'].ge(2).all()
    assert not attempts.loc[attempts.selected,'status'].eq('mpa').any()
    for tol in TOLS:
        z=visits[visits.tolerance_days==tol]
        assert set(z.visit_key)==set(v.visit_key)
        assert np.array_equal(z.density.to_numpy(),v.density.to_numpy())
        assert set(zip(annual[annual.tolerance_days==tol].series_id,annual[annual.tolerance_days==tol].year))==set(zip(base.series_id,base.year))
    # Contrasts separate altered values on common years from added-year coverage.
    effects=[]
    for (tol,sid),g in annual.groupby(['tolerance_days','series_id']):
        g=g[(g.t>0)&g.supported];b=stricta[(stricta.series_id==sid)&(stricta.t>0)&stricta.supported]
        common=set(g.year)&set(b.year);gm=g[g.year.isin(common)];bm=b[b.year.isin(common)]
        ff,fs,bs,bf=fit(g),fit(gm),fit(bm),fit(b)
        for scope,zz,ans in [('full_supported',g,ff),('same_year',gm,fs)]:
            effects.append(dict(series_id=sid,tolerance_days=tol,scenario='strict_reference' if tol==0 else f'upper_gap_plus_{tol}d_diagnostic',
                observation_scope=scope,years=';'.join(map(str,sorted(zz.year))),**ans,
                strict_full_estimate=bf['estimate'],strict_same_year_estimate=bs['estimate'],
                value_change=fs['estimate']-bs['estimate'],
                coverage_change=(ff['estimate']-fs['estimate'])-(bf['estimate']-bs['estimate']),
                total_change=ff['estimate']-bf['estimate'],interval_scope='conditional temporal HAC; not total uncertainty'))
    effects=pd.DataFrame(effects)
    joined=effects[(effects.tolerance_days==0)&(effects.observation_scope=='full_supported')].merge(oldeffects,on='series_id',suffixes=('_new','_old'),validate='one_to_one')
    assert len(joined)==6
    for col in ['estimate','se','low','high']: assert np.allclose(joined[col+'_new'],joined[col+'_old'],rtol=1e-9,atol=1e-10)
    summary=[]
    for tol in TOLS:
        a=annual[annual.tolerance_days==tol];vv=visits[visits.tolerance_days==tol]
        common=vv.merge(visits[visits.tolerance_days==0][['visit_key','biomass_point','method']],on='visit_key',suffixes=('','_strict'))
        both=common.biomass_point.notna()&common.biomass_point_strict.notna()
        summary.append(dict(tolerance_days=tol,post_total=int((a.t>0).sum()),
            post_supported=int(((a.t>0)&a.supported).sum()),unknown_visits=int((vv.method=='unsupported').sum()),
            recovered_visits=int((common.biomass_point.notna()&common.biomass_point_strict.isna()).sum()),
            changed_common_values=int((both&~np.isclose(common.biomass_point,common.biomass_point_strict,equal_nan=True)).sum()),
            changed_methods=int((common.method!=common.method_strict).sum())))
    summary=pd.DataFrame(summary)
    # Expected recovery is tested after a uniform all-case rebuild; never a row-ID exception.
    a7=annual[(annual.tolerance_days==7)&(annual.t>0)].merge(stricta[['series_id','year','supported']],on=['series_id','year'],suffixes=('','_strict'))
    recovered=a7[a7.supported&~a7.supported_strict]
    assert len(recovered)==1 and recovered.year.iloc[0]==2019
    assert recovered.series_id.iloc[0]=='Abalone Cove SMCA|M. franciscanus|PISCO|Bio'
    assert summary.set_index('tolerance_days').post_supported.to_dict()=={0:44,3:44,7:45,14:45}
    pf=attempts[(attempts.series_id==recovered.series_id.iloc[0])&(attempts.site=='POINT_FERMIN')&(attempts.year==2019)&(attempts.kind=='other_visit')&attempts.selected&(attempts.tolerance_days==7)]
    assert len(pf)==1 and pf.max_gap.iloc[0]==294 and pf.cal_max_gap.iloc[0]==290
    assert (density.t>0).sum()==63
    assert np.array_equal(summary.post_total.to_numpy(),np.repeat(63,4))
    # Poison one supported reference cell: complete-reference annual must become unknown.
    poisoned=visits[visits.tolerance_days==0].copy()
    eligible=set(zip(stricta.loc[stricta.supported,'series_id'],stricta.loc[stricta.supported,'year']))
    ix=next(i for i,t in poisoned.iterrows() if t.status=='reference' and (t.series_id,t.year) in eligible)
    target=poisoned.loc[ix];poisoned.loc[ix,'biomass_point']=np.nan
    pa=aggregate(poisoned,base,off)
    assert not pa[(pa.series_id==target.series_id)&(pa.year==target.year)].supported.iloc[0]
    # Hash every declared protected file, before any output writes.
    protected_path=repo/'outputs/researcher_scope_20261101/protected_hashes.json'
    inputs[str(protected_path)]=digest(protected_path)
    protected=json.loads(protected_path.read_text());pc=[]
    for p,h in protected.items():
        path=Path(p);path=path if path.is_absolute() else repo/path
        actual=digest(path);pc.append(dict(path=str(path),expected_sha256=h,actual_sha256=actual,unchanged=actual==h))
    assert all(x['unchanged'] for x in pc)
    assert all(digest(Path(p))==h for p,h in inputs.items())
    # Preserve density byte-for-byte, including all months and the documented baseline scope.
    (output/'historical_density_reference.csv').write_bytes((repo/B/'selected_density_annual.csv').read_bytes())
    assert digest(output/'historical_density_reference.csv')==inputs[str(repo/B/'selected_density_annual.csv')]
    change=visits[visits.tolerance_days==7].merge(
        visits[visits.tolerance_days==0][['visit_key','method','biomass_point']],
        on='visit_key',suffixes=('','_strict'),validate='one_to_one')
    change=change[change.method!=change.method_strict]
    frames={'changed_visits_7days.csv':change,'visit_support.csv':visits,'donor_attempts.csv':attempts,'annual.csv':annual,
            'effects.csv':effects,'summary.csv':summary,'outer_site_calibration.csv':eng.calibration,
            'heldout_validation.csv':validation,'sample_diagnostics.csv':s,
            'protected_checks.csv':pd.DataFrame(pc),'recovered_post_years_7days.csv':recovered}
    for name,frame in frames.items(): frame.to_csv(output/name,index=False)
    tests=dict(status='PASS',strict_visit_identity=len(v),strict_annual_identity=len(stricta),
        strict_model_estimate_se_interval_identity=len(joined),calibration_predictions_reproduced=len(eng.calibration),
        source_samples_reproduced=len(s),all_scenario_visits=len(visits),lossless_donor_attempts=len(attempts),
        protected_files_unchanged=len(pc),input_files_unchanged=len(inputs),
        complete_reference_poison_test=True,zero_unknown_test=True,all_month_density_byte_identity=True,
        no_cross_species_mpa_role_zone_year=True,outer_site_exclusion=True,minimum_two_groups=True,
        uniform_tolerance_actual_allcases=True,singleton_measured_uses=int(((strict.method=='measured_same_visit')&strict.singleton).sum()))
    (output/'tests.json').write_text(json.dumps(tests,indent=2)+'\n')
    (output/'input_hash_manifest.json').write_text(json.dumps(inputs,indent=2)+'\n')
    report='''# Borrowing support rebuild — review sensitivities only

Implemented a standalone, read-only-input rebuild with uniform upper calibration-gap tolerances of 0, 3, 7 and 14 days. Zero reproduces every current visit mass, annual biomass and six HAC slopes/SEs/intervals. The whole current all-month biomass roster is preserved. The earlier density CSV is preserved byte-for-byte as historical_density_reference.csv for validation only. The parent current_density_annual.csv uses the restored density-specific roster.

Both borrowing attempts are retained for every positive visit lacking same-visit sizes, including eligible alternatives not selected by the priority. Priority remains same visit, local other date, then other site. Keys include comparison/species/MPA, physical site, role, year, zone and date; row IDs are retained only for compatibility. Donor evidence includes every selected nearest-date sample (ties retained), sample n and stable keys. Calibration is recomputed with the target site removed from targets and donors, compared with all saved calibration predictions, and requires at least two site-year groups. Calibration residuals continue to pool eligible years/zones within species/MPA/role and outer-site exclusion; individual donor transfers never cross year or zone.

The tolerance changes only the upper observed calibration-gap boundary. Minimum observed gap is exposed but is not a new lower restriction. This reproduces existing behavior and makes extrapolation below the minimum visible. Calibration and held-out summaries expose sites, years, groups, min/max and protected validation status. Protected borrowing remains unsupported. Empirical frequency-weighted mass variance, n=1 and zero empirical variance are explicit; singleton sample variance remains missing. Zero empirical variation does not establish zero population uncertainty.

'''+summary.to_string(index=False)+'''\n
Seven days recovers exactly one post-protection paired year: Abalone Cove red urchins, 2019, through Point Fermin's 294-day local donor versus a 290-day calibration maximum. This is obtained by testing every case uniformly, with no row-ID exception. Three days restores none; 14 days adds no further post year. The uniform 7/14-day rule also changes two already-supported Point Fermin purple-urchin visits (INNER and MID, January 4, 2019) from cross-site to higher-priority local donors. changed_visits_7days.csv records all three changed visits. Red-urchin slope changes from -0.180509 to -0.188102 through coverage; purple changes from -0.098868 to -0.101968 through values on unchanged years. The other four slopes are unchanged.

Annual means propagate any unknown visit through zone, station and role means; no reference station is silently dropped. All valid count-zero visits remain zero; unsized positive counts remain unknown. A poison-cell regression test verifies full-reference failure propagation. effects.csv reports full supported post-protection observations and common years with the strict scenario; total change is decomposed into changed values and coverage. Calendar gaps are preserved in Bartlett HAC (lag 3, finite-sample n/(n-2), Student t intervals). These are conditional temporal intervals, not total uncertainty. No transfer-error simulation or primary model adoption is implied.

Run with `python3 rebuild_borrowing_support.py --repo /path/to/analysis --data-root /path/to/data --output /path/to/scratch`. Alternatively set BORROWING_OUTPUT. Dependencies: numpy, pandas, scipy. No top-level analysis script is imported, no source output overwritten, no shared continuity edit, no sending/push. Input hashes are in input_hash_manifest.json; protected checks and actual executed assertions are in protected_checks.csv and tests.json. Repository outputs are restricted to outputs/biomass_next_steps_20261102; raw input directories are never output locations.
'''
    (output/'report.md').write_text(report)
    (output/'manifest.json').write_text(json.dumps(dict(scope='Diagnostic only; no final primary adoption',
        repo=str(repo),data_root=str(data),tolerances_days=TOLS,
        source_script_sha256=digest(Path(__file__)),inputs=inputs,
        outputs={p.name:digest(p) for p in output.iterdir() if p.is_file() and p.name!='manifest.json' and p.suffix not in ['.log']},
        dependencies=dict(pandas=pd.__version__,numpy=np.__version__)),indent=2)+'\n')
    print(summary.to_string(index=False));print(json.dumps(tests,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repo',type=Path,default=Path('/Users/adrianstier/sbc-donham-kelp-mpa-cascade'))
    p.add_argument('--data-root',type=Path,default=Path('/Users/adrianstier/sbc-donham-kelp-mpa-data'))
    p.add_argument('--output',type=Path,default=os.environ.get('BORROWING_OUTPUT'))
    args=p.parse_args()
    if args.output is None: p.error('--output or BORROWING_OUTPUT is required')
    run(args.repo,args.data_root,args.output)
