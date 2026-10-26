suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork)})
ROOT <- file.path(RESULT_DIR,'weighted_single_sample')
FIG <- file.path(ROOT,'figures')
dir.create(FIG,recursive=TRUE,showWarnings=FALSE)
style <- theme_classic(base_family='Arial',base_size=8)+
  theme(axis.text=element_text(colour='black',size=7),axis.title=element_text(size=8),
        legend.position='bottom',legend.text=element_text(size=7.5),legend.key.width=unit(13,'pt'),
        strip.background=element_blank(),strip.text=element_text(face='bold',size=8),
        plot.margin=margin(5,5,5,5),plot.tag=element_text(face='bold',size=11))
selected <- c('equal_contrast','difference__linear_svm','rank__lasso','raw__enet','raw__radial_svm','raw__random_forest')
labels <- c(equal_contrast='Equal contrast',difference__linear_svm='Relative-expression SVM',
            rank__lasso='Rank LASSO',raw__enet='Expression elastic net',
            raw__radial_svm='Expression radial SVM',raw__random_forest='Expression random forest')
d <- fread(file.path(ROOT,'tables/fixed_transfer_results.csv'))[
  operating_point=='youden'&dataset%in%c('GSE35570','GSE60542','GSE29265')&candidate%in%selected]
long <- rbindlist(list(d[,.(dataset,candidate,measure='Sensitivity',value=sensitivity,lo=sensitivity_lo,hi=sensitivity_hi)],
                       d[,.(dataset,candidate,measure='Specificity',value=specificity,lo=specificity_lo,hi=specificity_hi)]))
long[,method:=factor(unname(labels[candidate]),levels=rev(unname(labels[selected])))]
long[,dataset:=factor(dataset,levels=c('GSE35570','GSE60542','GSE29265'))]
p1 <- ggplot(long,aes(100*value,method,colour=measure,shape=measure))+
  geom_errorbar(aes(xmin=100*lo,xmax=100*hi),position=position_dodge(width=.38),width=.1,orientation='y',linewidth=.35)+
  geom_point(position=position_dodge(width=.38),size=1.7)+facet_wrap(~dataset,nrow=1)+
  scale_colour_manual(values=c(Sensitivity='#0072B2',Specificity='#D55E00'))+
  scale_shape_manual(values=c(Sensitivity=16,Specificity=17))+
  scale_x_continuous(limits=c(-2,102),breaks=c(0,50,100))+
  labs(x='Performance (%, Wilson 95% CI)',y=NULL,colour=NULL,shape=NULL,tag='a')+
  style+theme(panel.spacing=unit(11,'pt'))
m <- fread(file.path(ROOT,'tables/weighted_signal_mixing.csv'))[candidate%in%c('equal_contrast','difference__linear_svm')]
m[,method:=factor(ifelse(candidate=='equal_contrast','Equal contrast','Weighted relative-expression SVM'),
                    levels=c('Equal contrast','Weighted relative-expression SVM'))]
m[,`:=`(x=100*ptc_signal_weight,y=100*positive_fraction,label=paste0(positive,'/24'))]
mid <- m[ptc_signal_weight>0&ptc_signal_weight<1]
edge <- m[candidate=='equal_contrast'&ptc_signal_weight%in%c(0,1)]
edge[,texty:=y+5]
p2 <- ggplot(m,aes(x,y,colour=method,shape=method,linetype=method))+
  geom_line(linewidth=.55)+geom_point(size=2)+
  geom_text(data=mid[candidate=='equal_contrast'],aes(label=label),vjust=1.6,show.legend=FALSE,size=2.5,family='Arial')+
  geom_text(data=mid[candidate=='difference__linear_svm'],aes(label=label),vjust=-.75,show.legend=FALSE,size=2.5,family='Arial')+
  geom_text(data=edge,aes(x=x,y=texty,label=label),inherit.aes=FALSE,colour='black',size=2.5,family='Arial')+
  scale_colour_manual(values=c('#0072B2','#D55E00'))+scale_shape_manual(values=c(16,17))+
  scale_linetype_manual(values=c('solid','dashed'))+
  scale_x_continuous(limits=c(-4,104),breaks=c(0,25,50,75,100))+
  scale_y_continuous(limits=c(0,108),breaks=c(0,25,50,75,100))+
  labs(x='Simulated PTC expression-signal weight (%)',y='Positive profiles (%)',colour=NULL,shape=NULL,linetype=NULL,tag='b')+style
plot <- (p1/p2)+plot_layout(heights=c(1,1.03))
ggsave(file.path(FIG,'Supplementary_Fig_S10.pdf'),plot,width=7.2,height=6.5,device=cairo_pdf)
ggsave(file.path(FIG,'Supplementary_Fig_S10.png'),plot,width=7.2,height=6.5,dpi=600,bg='white')
ggsave(file.path(FIG,'Supplementary_Fig_S10.tif'),plot,width=7.2,height=6.5,dpi=600,compression='lzw',bg='white')
cat('Supplementary Figure S10 generated from frozen comparison tables.\n')
