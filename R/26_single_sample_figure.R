# Supplementary Figure S9 from the fixed single-sample application outputs.
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork)})
APP <- file.path(RESULT_DIR,'single_sample')
FIG <- file.path(APP,'figures');dir.create(FIG,recursive=TRUE,showWarnings=FALSE)
results <- fread(file.path(APP,'tables','fixed_pipeline_results.csv'))
mix <- fread(file.path(APP,'tables','synthetic_signal_mixing.csv'))
cohorts <- c('GSE33630','GSE35570','GSE60542','GSE29265')
d <- results[method=='contrast'&operating_point=='youden'&dataset%in%cohorts]
plotdat <- rbind(d[,.(dataset,metric='Sensitivity',estimate=100*sensitivity,lower=100*sensitivity_lo,upper=100*sensitivity_hi)],
                d[,.(dataset,metric='Specificity',estimate=100*specificity,lower=100*specificity_lo,upper=100*specificity_hi)])
plotdat[,dataset:=factor(dataset,levels=rev(cohorts))]
plotdat[,metric:=factor(metric,levels=c('Sensitivity','Specificity'))]
theme_pub <- theme_classic(base_size=9,base_family='Arial')+
  theme(axis.title=element_text(size=9),axis.text=element_text(size=8,color='black'),
        legend.position='bottom',legend.title=element_blank(),legend.text=element_text(size=8),
        plot.title=element_text(size=10,face='bold',margin=margin(b=9)),
        plot.tag=element_text(face='bold',size=11),plot.margin=margin(8,8,8,8))
p1 <- ggplot(plotdat,aes(x=estimate,y=dataset,color=metric,shape=metric))+
  geom_errorbar(aes(xmin=lower,xmax=upper),orientation='y',width=.13,
                position=position_dodge(width=.45),linewidth=.45)+
  geom_point(position=position_dodge(width=.45),size=2.1)+
  scale_color_manual(values=c(Sensitivity='#0072B2',Specificity='#D55E00'))+
  scale_shape_manual(values=c(Sensitivity=16,Specificity=17))+
  scale_x_continuous(limits=c(40,102),breaks=c(40,60,80,100),expand=c(0,0))+
  labs(x='Performance (%)',y=NULL,title='Fixed-threshold classification')+theme_pub
m <- mix[method=='contrast'];m[,label:=paste0(positive,'/',n_profiles)];m[,hjust:=ifelse(tumor_signal_weight%in%c(0,.75,1),1,.5)]
p2 <- ggplot(m,aes(x=100*tumor_signal_weight,y=100*positive_fraction))+
  geom_line(color='#0072B2',linewidth=.7)+geom_point(color='#0072B2',size=2.2)+
  geom_text(aes(label=label,hjust=hjust),vjust=-.9,size=2.7,family='Arial')+
  scale_x_continuous(limits=c(-10,103),breaks=c(0,25,50,75,100),expand=c(0,0))+
  scale_y_continuous(limits=c(0,108),breaks=c(0,25,50,75,100),expand=c(0,0))+
  labs(x='PTC expression-signal weight (%)',y='Positive profiles (%)',title='Matched-profile signal mixing')+
  theme_pub+theme(legend.position='none')
p <- p1+p2+plot_layout(widths=c(1.08,1))+plot_annotation(tag_levels='a')
ggsave(file.path(FIG,'Supplementary_Fig_S9.pdf'),p,width=7.2,height=3.65,device=cairo_pdf)
ggsave(file.path(FIG,'Supplementary_Fig_S9.png'),p,width=7.2,height=3.65,dpi=600,bg='white')
ggsave(file.path(FIG,'Supplementary_Fig_S9.tif'),p,width=7.2,height=3.65,dpi=600,compression='lzw',bg='white')
cat('Supplementary Figure S9 exported.\n')
