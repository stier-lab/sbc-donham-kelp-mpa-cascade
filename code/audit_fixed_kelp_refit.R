# Run after the SBC127.5 refit; baseline path is an explicit frozen audit input.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1)
baseline <- args[1]
source('code/R/00_libraries.R')
old <- read.csv(file.path(baseline,'data/sumstats_final.csv'))
new <- read.csv('data/sumstats_final.csv')
keys <- c('MPA','Taxa','Source','Resp')
id <- function(d) do.call(paste,c(d[keys],sep='|'))
stopifnot(setequal(id(old),id(new)),!anyDuplicated(id(new)))
old <- old[match(id(new),id(old)),]
changed <- abs(new$Mean-old$Mean)>1e-8 | abs(new$SE-old$SE)>1e-8
stopifnot(all(new$Taxa[changed]=='M. pyrifera'),all(new$Resp[changed]=='Bio'))
out <- cbind(new[keys],old_model=old$Model,new_model=new$Model,old_estimate=old$Mean,new_estimate=new$Mean,old_se=old$SE,new_se=new$SE)
write.csv(out,'tables/table_s_fixed_kelp_refit_comparison.csv',row.names=FALSE)
source('code/R/51_kelp_selection_sensitivity.R')
cat('Fixed-coefficient audit passed:',nrow(new),'histories;',sum(changed),'changed kelp histories;',length(switched),'changed model selections.\n')
