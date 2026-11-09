"""Shared input, exclusion and identity rules for the integrated review."""
from pathlib import Path
import hashlib,json,re,os
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT.parent/'sbc-donham-kelp-mpa-data'
BASE_OUT=ROOT/'outputs/integrated_review_20261103'
VARIANT=os.environ.get('INTEGRATED_REVIEW_VARIANT','')
OUT=BASE_OUT/'sensitivities'/VARIANT if VARIANT else BASE_OUT
PRIVATE=OUT/'private';PRIVATE.mkdir(parents=True,exist_ok=True)
POLICY=json.loads((ROOT/'config/integrated_review_policy.json').read_text())
if VARIANT=='exact_KFM_dates':POLICY['kfm_same_site_size_gap_days']=0
if VARIANT=='four_before':POLICY['minimum_before_years']=4
MANIFEST=json.loads((DATA/'config/source_manifest.json').read_text())
SOURCES={s['source_id']:s for s in MANIFEST['sources']}
INPUTS={}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read_source(name,**kw):
 s=SOURCES[name];p=DATA/s['path'];h=sha(p)
 if h!=s['sha256']:raise ValueError('Changed source '+name)
 INPUTS[str(p)]=h
 return pd.read_csv(p,low_memory=False,**kw)
def read_tracked(path,**kw):
 p=Path(path);INPUTS[str(p)]=sha(p);return pd.read_csv(p,low_memory=False,**kw)
def sid(source,site,*parts):return '|'.join(map(str,(source,site,*parts)))
def strict_sum(s):return s.sum() if s.notna().all() else np.nan
def strict_mean(s):return s.mean() if len(s) and s.notna().all() else np.nan
TAXA={'MESFRAAD':'M. franciscanus','MESFRA':'M. franciscanus','STRPURAD':'S. purpuratus','STRPUR':'S. purpuratus','PANINT':'P. interruptus','MACPYRAD':'M. pyrifera','SPUL':'S. pulcher','SFL':'M. franciscanus','SFS':'M. franciscanus','SPL':'S. purpuratus','SPS':'S. purpuratus','Mesocentrotus franciscanus':'M. franciscanus','Strongylocentrotus franciscanus':'M. franciscanus','Strongylocentrotus purpuratus':'S. purpuratus','Macrocystis pyrifera':'M. pyrifera','Panulirus interruptus':'P. interruptus','Semicossyphus pulcher':'S. pulcher'}
def mass(tx,size):
 if tx=='M. franciscanus':return .00059*np.asarray(size)**2.917
 if tx=='S. purpuratus':return .00059*np.asarray(size)**2.870
 if tx=='P. interruptus':return .001352821*np.asarray(size)**2.913963113
 if tx=='S. pulcher':return .0144*np.asarray(size)**3.04
 raise ValueError(tx)
utils=(DATA/'code/R/01_utils.R').read_text();INPUTS[str(DATA/'code/R/01_utils.R')]=sha(DATA/'code/R/01_utils.R')
def excluded(name):return set(re.findall(r'"([^"]+)"',re.search(name+r' <- c\((.*?)\)',utils,re.S).group(1)))
PISCO_EXCLUDED=excluded('EXCLUDED_REFERENCE_SITES');KFM_EXCLUDED=excluded('EXCLUDED_KFM_SITES');LTER_EXCLUDED={'AHND','a-l-02'}
def blocked(source,site,taxon=None):
 if source=='KFM' and taxon=='S. pulcher':return False # benthic exclusions are not fish exclusions
 return site in (PISCO_EXCLUDED if source=='PISCO' else KFM_EXCLUDED if source=='KFM' else LTER_EXCLUDED if source=='LTER' else set())
def write_inputs():
 p=OUT/'input_hashes.json';old=json.loads(p.read_text()) if p.exists() else {};old.update(INPUTS);p.write_text(json.dumps(old,indent=2))


def measurement_covariance(draws, years, source, response, taxon):
    """Respect known independent survey-year measurement draws, not MC noise.

    KFM count and size records are year-specific; approved near-date donors stay
    within year. PISCO density, direct fish and kelp measurements likewise have
    no cross-year donor error. LTER permanent-transect weights and calibrated
    PISCO borrowed biomass can span years, so retain their full covariance.
    Biological serial dependence is estimated separately in the process model.
    """
    covariance = np.cov(draws)
    independent_years = source == 'KFM' or (source == 'PISCO' and (response == 'Den' or taxon in ['S. pulcher', 'M. pyrifera']))
    if independent_years:
        years = np.asarray(years)
        covariance *= years[:, None] == years[None, :]
    return (covariance + covariance.T) / 2
