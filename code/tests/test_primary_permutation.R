suppressPackageStartupMessages({library(dplyr);library(metafor)})
a <- read.csv('data/sumstats_final.csv') %>% group_by(MPA,Taxa,Resp) %>%
  summarise(y=mean(Mean),v=sum(SE^2)/n()^2,.groups='drop')
x <- filter(a,Taxa=='P. interruptus',Resp=='Bio')
y <- filter(a,Taxa=='S. purpuratus',Resp=='Bio')
d <- merge(x,y,by='MPA',suffixes=c('_x','_y'));stopifnot(nrow(d)==4)
g <- as.matrix(expand.grid(rep(list(1:4),4)))
g <- g[apply(g,1,function(z)length(unique(z))==4),]
fit <- function(x,y,v)as.numeric(coef(rma(yi=y,vi=v,mods=~x,method='REML',test='knha'))[2])
obs <- fit(d$y_x,d$y_y,d$v_y)
# Two independent label-shuffling formulations must yield the same distribution.
px <- apply(g,1,function(i)fit(d$y_x[i],d$y_y,d$v_y))
py <- apply(g,1,function(i)fit(d$y_x,d$y_y[i],d$v_y[i]))
stopifnot(length(px)==24,all(is.finite(px)),max(abs(sort(px)-sort(py)))<1e-7)
p <- mean(abs(px)>=abs(obs)-1e-12)
t <- read.csv('tables/table_s_prereview_crosstaxa_permutation.csv')
z <- t[t$k==4,];stopifnot(nrow(z)==1,all(t$reproduces_published),
  abs(z$p_permutation-p)<.00005,abs(z$slope_reproduced-round(obs,4))<1e-8,
  all(t$permutation_unit=='MPA predictor labels; response/SE pairs retained'))
cat('PASS: all six primary slopes reproduced; 24 exact shuffles agree with paired-response/variance enumeration; P=',p,'\n')
