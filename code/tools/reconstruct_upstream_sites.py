from pathlib import Path
import os,shutil,json,hashlib
root=Path.cwd();dr=Path(os.environ.get('MPA_DATA_REPO',str(root.parent/'Donham-Stier-CA-MPA-Data-2026')))
out=root/'outputs/upstream_sites_20261101';stage=out/'stage';stage.mkdir(parents=True,exist_ok=True)
(out/'.gitignore').write_text('stage/\n')
for name in ['AGENTS.md','CLAUDE.md']:
 if (dr/name).exists():raise RuntimeError('Read newly present data-repository instructions before proceeding: '+str(dr/name))
(stage/'.here').touch();(stage/'code/R').mkdir(parents=True,exist_ok=True);(stage/'data').mkdir(exist_ok=True)
for f in (dr/'code/R').glob('*.R'):shutil.copy2(f,stage/'code/R'/f.name)
for f in (dr/'data').iterdir():
 dest=stage/'data'/f.name
 if f.name=='cache':
  if not dest.exists():shutil.copytree(f,dest)
 elif not dest.exists():dest.symlink_to(f.resolve(),target_is_directory=f.is_dir())
# Only remove named whole-MPA omissions. Other expert reference exclusions remain.
p=stage/'code/R/04_pisco_processing.R';s=p.read_text();s=s.replace('"Point Dume SMR"','"__no_legacy_Point_Dume_omission__"');
needle='dplyr::filter(!duplicated(cbind(site, CA_MPA_Name_Short)))'
assert s.count(needle)==1
s=s.replace(needle,needle+'\nreview_mpas <- c("Anacapa Island SMCA","Painted Cave SMCA","Point Dume SMR","Point Dume SMCA","Abalone Cove SMCA","Matlahuayl SMR")\nreview_sites <- unique(sites.short$site[sites.short$CA_MPA_Name_Short %in% review_mpas])\nSwath.edit <- subset(Swath.edit, site %in% review_sites)')
s=s.replace('Fish <- subset(Fish, campus == "UCSB" | campus == "VRG")','Fish <- subset(Fish, (campus == "UCSB" | campus == "VRG") & site %in% review_sites)')
p.write_text(s)
p=stage/'code/R/05_kfm_processing.R' ;s=p.read_text();s=s.replace('"Anacapa Island SMCA"','"__no_legacy_Anacapa_omission__"');needle='kfm.site <- kfm.site[!(kfm.site$site_id %in% EXCLUDED_KFM_SITES), ]';assert s.count(needle)==1;s=s.replace(needle,needle+'\nkfm.site <- subset(kfm.site, CA_MPA_Name_Short == "Anacapa Island SMCA")')
p.write_text('if (exists("Urchin.site")) rm(Urchin.site)\nif (exists("Macro.site")) rm(Macro.site)\n'+s)
runner='''Sys.setenv(DONHAM_KELP_CALIBRATION_CSV=normalizePath("../../../data/external/sbc1275/Macrocystis_Biomass_Relationships_20261016.csv"))
suppressPackageStartupMessages({library(plyr);library(tidyverse);library(here)})
here::i_am("run_review.R")
source("code/R/00c_analysis_constants.R")
source("code/R/01_utils.R")
EXCLUDED_MPAS<-setdiff(EXCLUDED_MPAS,c("Anacapa Island SMCA","Painted Cave SMCA","Point Dume SMR","Point Dume SMCA"))
# Keyhole was a whole-site omission without a recovered method-specific failure.
# Reconstruct from source observations, evaluate actual coverage, disclose restoration.
EXCLUDED_KFM_SITES<-setdiff(EXCLUDED_KFM_SITES,"a-k-31")
source("code/R/03_data_import.R")
source("code/R/04_pisco_processing.R")
saveRDS(list(raw=PISCO.resp,ratios=Swath.join.sub,annual=All_PISCO.mpa,sites=All_PISCO.short.sub,fish=Fish.SPUL.All),"../pisco_reconstruction.rds")
write.csv(PISCO.resp,"../pisco_raw_all.csv",row.names=FALSE)
write.csv(All_PISCO.short.sub,"../pisco_site_values.csv",row.names=FALSE)
write.csv(Swath.join.sub,"../pisco_ratios_all.csv",row.names=FALSE)
saveRDS(SizeFreq.Urch.OG,"../kfm_urchin_size_pool.rds")
source("../../../code/tools/repair_upstream_size_assignments.R")
FORCE_BOOTSTRAP <- TRUE
source("code/R/05_kfm_processing.R")
saveRDS(list(raw=KFM.resp,fish=KFM.fish.den.long,sites=Swath.KFM,fish_sites=mbon.spul.ave,macro=Macro.site,urchin=Urchin.site),"../kfm_reconstruction.rds")
write.csv(KFM.resp,"../kfm_raw.csv",row.names=FALSE)
write.csv(KFM.fish.den.long,"../kfm_fish_raw.csv",row.names=FALSE)
write.csv(Swath.KFM,"../kfm_site_values.csv",row.names=FALSE)
write.csv(mbon.spul.ave,"../kfm_fish_site_values.csv",row.names=FALSE)
'''
(stage/'run_review.R').write_text(runner)
files=list((dr/'code/R').glob('*.R'))+list((dr/'data/cache').glob('*.rds'))
(out/'source_code_cache_hashes.json').write_text(json.dumps({str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in files},indent=2))
print(stage)
