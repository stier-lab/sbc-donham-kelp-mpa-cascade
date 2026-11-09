suppressPackageStartupMessages(library(ggplot2))
out<-'outputs/pooled_common_estimand_20261101'
e<-read.csv(file.path(out,'effects.csv'));p<-read.csv(file.path(out,'pooled.csv'))
labels<-c('M. pyrifera'='Kelp','M. franciscanus'='Red urchin','S. purpuratus'='Purple urchin','P. interruptus'='Lobster','S. pulcher'='Sheephead')
label<-function(x)paste(labels[x$Taxa],ifelse(x$Resp=='Bio','biomass','density'))
x<-subset(e,scenario=='Primary');x$label<-label(x)
y<-subset(p,scenario=='Primary'&view=='All');y$label<-label(y)
levels<-rev(c('Kelp biomass','Sheephead biomass','Lobster biomass','Red urchin biomass','Purple urchin biomass','Sheephead density','Lobster density','Red urchin density','Purple urchin density'))
x$label<-factor(x$label,levels);y$label<-factor(y$label,levels)
y$precision <- ifelse(is.finite(y$low),'Interval available','Descriptive only')
x$metric<-factor(x$metric,c('Change','Post_difference'),labels=c('Before–After change','After-protection difference'));y$metric<-factor(y$metric,c('Change','Post_difference'),labels=levels(x$metric))
th<-theme_minimal(base_size=12)+theme(legend.position='none',panel.grid.minor=element_blank(),plot.margin=margin(20,20,20,20),panel.grid.major.y=element_blank())
g<-ggplot(x,aes(estimate,label))+geom_vline(xintercept=0,color='grey65',linetype='dashed')+geom_point(aes(color=Source),position=position_jitter(width=0,height=.14,seed=68),alpha=.55,size=1.8)+geom_errorbar(data=y,aes(xmin=low,xmax=high),orientation='y',width=.22,linewidth=.65,na.rm=TRUE)+geom_point(data=y[y$n_mpa>1,],aes(shape=precision),size=3)+scale_shape_manual(values=c('Interval available'=18,'Descriptive only'=5))+facet_wrap(~metric,nrow=1)+scale_color_manual(values=c(Landsat='#D55E00',KFM='#0072B2',LTER='#009E73',PISCO='#CC79A7'))+labs(x='Log protected/reference ratio or change in that ratio',y=NULL)+th
ggsave(file.path(out,'common_estimands.png'),g,width=13,height=7.5,dpi=180);ggsave(file.path(out,'common_estimands.pdf'),g,width=13,height=7.5)
