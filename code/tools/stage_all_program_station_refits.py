from pathlib import Path
import shutil,json,hashlib
D=Path('/Users/adrianstier/Donham-Stier-CA-MPA-Data-2026');O=Path('outputs/all_program_station_refits_20261102');S=O/'stage';(S/'code/R').mkdir(parents=True,exist_ok=True);(S/'data').mkdir(exist_ok=True);(S/'.here').touch();(O/'.gitignore').write_text('stage/\n')
for p in (D/'code/R').glob('*.R'):shutil.copy2(p,S/'code/R'/p.name)
for p in (D/'data').iterdir():
 q=S/'data'/p.name
 if q.exists():continue
 if p.name=='cache':shutil.copytree(p,q)
 else:q.symlink_to(p.resolve(),target_is_directory=p.is_dir())
runner='''suppressPackageStartupMessages({library(plyr);library(tidyverse);library(here)})
here::i_am("run.R")
Sys.setenv(DONHAM_KELP_CALIBRATION_CSV=normalizePath("../../../data/external/sbc1275/Macrocystis_Biomass_Relationships_20261016.csv"))
source("code/R/00c_analysis_constants.R");source("code/R/01_utils.R");source("code/R/03_data_import.R")
source("code/R/04_pisco_processing.R")
write.csv(All_PISCO.short.sub,"../pisco_sites.csv",row.names=FALSE)
source("code/R/05_kfm_processing.R")
write.csv(Swath.KFM,"../kfm_sites.csv",row.names=FALSE);write.csv(mbon.spul.ave,"../kfm_fish_sites.csv",row.names=FALSE)
source("code/R/06_lter_processing.R")
save.image("reconstruction.RData")
for(n in c("Swath.LTER","lter.lob.site.ave","lter.macro.site.ave","lter.fish.ave"))write.csv(get(n),paste0("../",n,".csv"),row.names=FALSE)
'''
(S/'run.R').write_text(runner);(O/'source_code_hashes.json').write_text(json.dumps({str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (D/'code/R').glob('*.R')},indent=2))
