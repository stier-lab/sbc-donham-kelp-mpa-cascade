#!/usr/bin/env Rscript
# Run from repository root; isolated outputs only. See output PLAN.md.
suppressPackageStartupMessages(library(ggplot2))
source('code/R/primary_uncertainty_utils.R')
args <- commandArgs(trailingOnly=TRUE)
out <- if(length(args)) args[1] else 'outputs/satellite_baseline_comparison_20261101'
cohort <- read.csv(file.path(out,'cohort.csv'))
if(file.exists(file.path(out,'freeze_manifest.json'))) {
 f <- jsonlite::fromJSON(file.path(out,'freeze_manifest.json'))$files
 stopifnot(all(vapply(names(f),function(p) digest::digest(file=p,algo='sha256')==f[[p]],logical(1))))
}
d <- read.csv(file.path(out,'annual_input.csv'))
stopifnot(setequal(unique(d$MPA),cohort$MPA), !anyDuplicated(d[c('MPA','year')]),all(is.finite(d$lnDiff)))
forms <- list(Stable=lnDiff~post+post_time, Changing=lnDiff~t+post+post_time)
results <- list(); preds <- list(); inputs <- list(); holdout <- list(); coefficients <- list()
for (site in unique(d$MPA)) {
 z <- d[d$MPA==site,];z <- z[order(z$year),]
 z$t <- z$year-z$MPA_Start;z$post <- as.integer(z$BA=='After');z$post_time <- z$post*z$t
 stopifnot(all(z$post==(z$year>z$MPA_Start)))
 for (window in c('Full','10','5')) {
  x <- if(window=='Full') z else z[z$post==1 | z$year>unique(z$MPA_Start)-as.integer(window),]
  stopifnot(sum(x$post==0)>=3,sum(x$post==1)>=3)
  inputs[[length(inputs)+1]] <- data.frame(x,window)
  for(model in names(forms)) {
   fit <- lm(forms[[model]],data=x);X <- model.matrix(fit)
   stopifnot(fit$rank==ncol(X))
   cf <- x;cf$post <- 0;cf$post_time <- 0
   X0 <- model.matrix(delete.response(terms(fit)),cf)
   D <- X-X0
   gradients <- list(Mean_after=colMeans(D[x$post==1,,drop=FALSE]),Last_after=D[nrow(x),],Level=setNames(as.numeric(names(coef(fit))=='post'),names(coef(fit))),Slope=setNames(as.numeric(names(coef(fit))=='post_time'),names(coef(fit))))
   for(lag in c(1,3,5)) {
    V <- uncertainty_hac(fit,x$year,lag)
    for(quantity in names(gradients)) {
     g <- gradients[[quantity]];est <- sum(g*coef(fit));se <- sqrt(drop(g%*%V%*%g));crit <- qt(.975,df.residual(fit))
     results[[length(results)+1]] <- data.frame(MPA=site,window,model,lag,quantity,n=nrow(x),n_before=sum(x$post==0),n_after=sum(x$post==1),first_year=min(x$year),last_year=max(x$year),estimate=est,se,low=est-crit*se,high=est+crit*se,p=2*pt(-abs(est/se),df.residual(fit)),AICc=AIC(fit)+2*(length(coef(fit))+1)*(length(coef(fit))+2)/(nrow(x)-length(coef(fit))-2),residual_lag1=if(sum(diff(x$year)==1)>2) cor(head(resid(fit),-1)[diff(x$year)==1],tail(resid(fit),-1)[diff(x$year)==1]) else NA)
    }
    if(lag==3 && window=='Full') {
     fitted_se <- sqrt(rowSums((X%*%V)*X));counter_se <- sqrt(rowSums((X0%*%V)*X0));crit <- qt(.975,df.residual(fit))
     preds[[length(preds)+1]] <- data.frame(MPA=site,year=x$year,t=x$t,post=x$post,model,observed=x$lnDiff,fitted=as.numeric(X%*%coef(fit)),low=as.numeric(X%*%coef(fit))-crit*fitted_se,high=as.numeric(X%*%coef(fit))+crit*fitted_se,counter=as.numeric(X0%*%coef(fit)),counter_low=as.numeric(X0%*%coef(fit))-crit*counter_se,counter_high=as.numeric(X0%*%coef(fit))+crit*counter_se)
    }
   }
   coefficients[[length(coefficients)+1]] <- data.frame(MPA=site,window,model,term=names(coef(fit)),coefficient=unname(coef(fit)))
  }
 }
 b <- z[z$post==0,];train <- head(b,-3);test <- tail(b,3)
 for(model in names(forms)) {
  f <- lm(if(model=='Stable') lnDiff~1 else lnDiff~t,data=train)
  pp <- predict(f,test)
  holdout[[length(holdout)+1]] <- data.frame(MPA=site,model,train_n=nrow(train),test_years=paste(test$year,collapse=','),rmse=sqrt(mean((test$lnDiff-pp)^2)),bias=mean(pp-test$lnDiff))
 }
}
r <- do.call(rbind,results)
r$q <- ave(r$p,interaction(r$window,r$model,r$lag,r$quantity),FUN=function(x)p.adjust(x,'BH'))
p <- do.call(rbind,preds);p$model <- factor(p$model,levels=c('Stable','Changing'));i <- do.call(rbind,inputs);h <- do.call(rbind,holdout)
write.csv(r,file.path(out,'estimates.csv'),row.names=FALSE);write.csv(p,file.path(out,'predictions.csv'),row.names=FALSE)
write.csv(i,file.path(out,'model_input_years.csv'),row.names=FALSE);write.csv(h,file.path(out,'baseline_holdout.csv'),row.names=FALSE)
write.csv(do.call(rbind,coefficients),file.path(out,'coefficients.csv'),row.names=FALSE)
cols <- c(Stable='#0072B2',Changing='#D55E00')
th <- theme_minimal(base_size=12)+theme(legend.position='none',panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),strip.text=element_text(size=11),plot.margin=margin(25,25,25,30))
# Separate facets for models prevent obscuring identical post-designation fits.
trajectory <- function(sites,file,width=11,height=10) {
 pp <- p[p$MPA%in%sites,];pp$MPA <- factor(pp$MPA,levels=sites)
 g <- ggplot(pp,aes(year,observed))+geom_hline(yintercept=0,color='grey80')+geom_vline(data=unique(d[d$MPA%in%sites,c('MPA','MPA_Start')]),aes(xintercept=MPA_Start),linetype='dotted',color='grey45')+
 geom_ribbon(aes(ymin=low,ymax=high,fill=model,group=interaction(MPA,model,post)),alpha=.16)+
 geom_ribbon(data=pp[pp$post==1,],aes(ymin=counter_low,ymax=counter_high,fill=model),alpha=.08)+
 geom_line(aes(y=fitted,color=model,group=interaction(MPA,model,post)),linewidth=.65)+
 geom_line(data=pp[pp$post==1,],aes(y=counter,color=model),linetype='dashed',linewidth=.7)+
 geom_point(aes(color=model),size=1.3,alpha=.8)+facet_grid(MPA~model,scales='free_y')+scale_color_manual(values=cols)+scale_fill_manual(values=cols)+labs(x='Year',y='Kelp biomass: log(MPA / reference)')+th
 ggsave(file.path(out,paste0(file,'.png')),g,width=width,height=height,dpi=180)
 ggsave(file.path(out,paste0(file,'.pdf')),g,width=width,height=height)
}
trajectory(c('Santa Barbara Island SMR','Swamis SMCA','Naples SMCA'),'key_trajectories',11,9)
sites <- sort(unique(d$MPA))
for(j in seq_len(ceiling(length(sites)/4))) trajectory(sites[((j-1)*4+1):min(j*4,length(sites))],paste0('trajectories_',j),11,11)
if(all(c('Judith Rk SMR','Point Conception SMR') %in% sites)) trajectory(c('Judith Rk SMR','Point Conception SMR'),'new_pair_trajectories',11,7)
a <- subset(r,window=='Full' & lag==3 & quantity=='Mean_after')
a$MPA <- factor(a$MPA,levels=rev(sites))
g <- ggplot(a,aes(estimate,MPA,color=model))+geom_vline(xintercept=0,color='grey60',linetype='dashed')+geom_errorbar(aes(xmin=low,xmax=high),position=position_dodge(width=.5),width=.2,orientation='y')+geom_point(position=position_dodge(width=.5),size=2.4)+scale_color_manual(values=cols)+labs(x='Average departure from the earlier trajectory (log ratio)',y=NULL)+th
ggsave(file.path(out,'all_pair_comparison.png'),g,width=10,height=max(8,length(sites)*.48),dpi=180);ggsave(file.path(out,'all_pair_comparison.pdf'),g,width=10,height=max(8,length(sites)*.48))
a <- subset(r,lag==3 & quantity=='Mean_after');a$window <- factor(a$window,levels=c('Full','10','5'));a$MPA <- factor(a$MPA,levels=sites)
g <- ggplot(a,aes(window,estimate,color=model,group=model))+geom_hline(yintercept=0,color='grey70')+geom_line(position=position_dodge(.25),linewidth=.4)+geom_errorbar(aes(ymin=low,ymax=high),width=.15,position=position_dodge(.25))+geom_point(position=position_dodge(.25),size=1.5)+facet_wrap(~MPA,ncol=4,scales='free_y')+scale_color_manual(values=cols)+labs(x='Baseline window (years)',y='Average departure (log ratio)')+th+theme(strip.text=element_text(size=9))
ggsave(file.path(out,'baseline_windows.png'),g,width=13,height=10,dpi=180)
cat('Completed',nrow(r),'contrasts;',nrow(p),'trajectory rows\n')
