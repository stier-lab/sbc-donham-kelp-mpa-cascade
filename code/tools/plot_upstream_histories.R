suppressPackageStartupMessages(library(ggplot2))
out<-'outputs/upstream_sites_20261101'
e<-read.csv(file.path(out,'effects.csv'));old<-read.csv('outputs/diver_exclusion_refit_20261101/consistent_kfm/annual_input.csv')
x<-subset(e,scenario=='Primary'&metric=='Post_difference'&!series_id%in%old$series_id)
labels<-c('M. pyrifera'='Kelp','M. franciscanus'='Red urchin','S. purpuratus'='Purple urchin','P. interruptus'='Lobster','S. pulcher'='Sheephead')
x$label<-paste(labels[x$Taxa],ifelse(x$Resp=='Bio','biomass','density'))
x$label<-factor(x$label,rev(c('Kelp biomass','Sheephead biomass','Red urchin biomass','Purple urchin biomass','Sheephead density','Lobster density','Red urchin density','Purple urchin density')))
x$low<-x$estimate-qt(.975,x$df)*x$se;x$high<-x$estimate+qt(.975,x$df)*x$se
pd<-position_dodge(width=.45)
g<-ggplot(x,aes(estimate,label,color=Source,shape=Source))+geom_vline(xintercept=0,color='grey65',linetype='dashed')+geom_errorbar(aes(xmin=low,xmax=high),orientation='y',width=.15,position=pd,linewidth=.6)+geom_point(position=pd,size=2.5)+facet_wrap(~MPA,ncol=2,scales='free_y')+scale_color_manual(values=c(KFM='#0072B2',PISCO='#D55E00'))+scale_shape_manual(values=c(KFM=16,PISCO=17))+labs(x='Post-protection log protected/reference ratio',y=NULL)+theme_minimal(base_size=12)+theme(legend.position='none',panel.grid.minor=element_blank(),panel.grid.major.y=element_blank())
ggsave(file.path(out,'restored_histories.png'),g,width=13,height=8,dpi=180);ggsave(file.path(out,'restored_histories.pdf'),g,width=13,height=8)
write.csv(x,file.path(out,'restored_primary_effects.csv'),row.names=FALSE)
