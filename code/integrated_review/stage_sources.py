"""Isolate the corrected source pipeline; never write through source symlinks."""
from pathlib import Path
import shutil,json,hashlib
A=Path(__file__).resolve().parents[2];D=A.parent/'sbc-donham-kelp-mpa-data';O=A/'outputs/integrated_review_20261103';S=O/'stage';S.mkdir(parents=True,exist_ok=True)
(O/'.gitignore').write_text('stage/\nprivate/\n')
for name in ['code','config']:shutil.copytree(D/name,S/name,dirs_exist_ok=True)
(S/'data').mkdir(exist_ok=True)
for p in (D/'data').iterdir():
 q=S/'data'/p.name
 if q.exists() or q.is_symlink():continue
 if p.name in ['cache','processed']:shutil.copytree(p,q) if p.is_dir() else shutil.copy2(p,q)
 else:q.symlink_to(p.resolve(),target_is_directory=p.is_dir())
(S/'.here').touch()
for name in ['tables','output','logs']:(S/name).mkdir(exist_ok=True)
# Preserve full objects for integration adapters, rather than reading frozen summaries.
r='''suppressPackageStartupMessages({library(plyr);library(tidyverse);library(here);library(Rmisc);library(sf);library(jsonlite);library(digest)})
here::i_am("integrated_source_run.R")
Sys.setenv(DONHAM_KELP_CALIBRATION_CSV="/Users/adrianstier/sbc-donham-kelp-mpa-cascade/data/external/sbc1275/Macrocystis_Biomass_Relationships_20261016.csv")
Sys.setenv(DONHAM_HARMONIZED_OUTPUT=file.path(here::here(),"output","harmonized"))
source("code/R/00c_analysis_constants.R")
source("code/R/01_utils.R")
source("code/R/03_data_import.R")
source("code/R/04_pisco_processing.R")
save.image("after_pisco.RData")
source("code/R/05_kfm_processing.R")
source("code/R/06_lter_processing.R")
source("code/R/06b_landsat_processing.R")
save.image("reconstruction.RData")
source("code/R/07_combine_and_export.R")
save.image("completed.RData")
'''
(S/'integrated_source_run.R').write_text(r)
files=list((D/'code').rglob('*.R'))+list((D/'code').rglob('*.py'))+list((D/'config').glob('*'))
(O/'source_code_manifest.json').write_text(json.dumps({str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in files if p.is_file()},indent=2))
print(S)
