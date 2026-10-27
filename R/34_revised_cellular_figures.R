options(warn=1)
Sys.setenv(OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2',MKL_NUM_THREADS='2')
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork);library(png);library(grid)})
F <- file.path(OUTPUT_DIR,'figures');dir.create(F,recursive=TRUE,showWarnings=FALSE)
TAB <- file.path(SOURCE_ROOT,'data/figure_tables');REV<-file.path(TAB,'revision')
SIG<-c('ABHD17C','ZDHHC6','ABHD17B','ZDHHC20')
gene_display<-function(x){
  x<-as.character(x)
  ifelse(x%in%SIG,paste0('<i>',x,'</i>'),x)
}
ST<-theme_bw(base_family='Arial',base_size=9)+theme(panel.grid.minor=element_blank(),
  panel.grid.major=element_line(colour='grey94',linewidth=.25),axis.text=element_text(colour='grey20',size=7.5),
  axis.title=element_text(size=8),plot.title=element_text(size=9,face='bold'),plot.subtitle=element_text(size=7.3),
  plot.tag=element_text(size=11,face='bold'),legend.position='bottom',legend.text=element_text(size=7),
  legend.title=element_text(size=7),plot.margin=margin(5,5,4,4))
savefig<-function(p,name,h){
  ggsave(file.path(F,paste0(name,'.png')),p,width=6.8,height=h,dpi=300,bg='white',limitsize=FALSE)
  ggsave(file.path(F,paste0(name,'.pdf')),p,width=6.8,height=h,device=cairo_pdf,bg='white',limitsize=FALSE)
  ggsave(file.path(F,paste0(name,'.tif')),p,width=6.8,height=h,dpi=600,compression='lzw',bg='white',limitsize=FALSE)
}

# Molecular associations and mutation-stratified tissue differences are contiguous.
cl<-fread(file.path(REV,'clinical_immune_fixed.csv'))
param<-fread(file.path(TAB,'frozen_tcga_parameters.csv'))
intercept<-29.269077
if(all(c('term','coefficient') %in% names(param)))intercept<-param[term=='(Intercept)',coefficient]
cl[,risk:=risk+intercept]
bx<-rbind(cl[!is.na(BRAF_V600E),.(group=ifelse(BRAF_V600E==1,'BRAF V600E','BRAF WT'),value=risk,panel='BRAF V600E')],
           cl[!is.na(RAS_MUT),.(group=ifelse(RAS_MUT==1,'RAS mutant','RAS WT'),value=risk,panel='RAS mutant')])
bx[,group:=factor(group,levels=c('BRAF WT','BRAF V600E','RAS WT','RAS mutant'))]
bx[,panel:=factor(panel,levels=c('BRAF V600E','RAS mutant'))]
nns<-bx[,.(n=.N),by=group];labsx<-setNames(paste0(nns$group,'\n(n = ',nns$n,')'),nns$group)
labsx<-setNames(gsub('\n','<br>',sub('^BRAF','<i>BRAF</i>',unname(labsx)),fixed=TRUE),names(labsx))
panel_labels<-c('BRAF V600E'='<i>BRAF</i> V600E','RAS mutant'='RAS mutant')
pa<-ggplot(bx,aes(group,value,fill=group))+geom_boxplot(width=.55,outlier.shape=NA,linewidth=.4,alpha=.85)+
  geom_point(position=position_jitter(width=.13,seed=11),size=.45,alpha=.3)+facet_wrap(~panel,scales='free_x',labeller=as_labeller(panel_labels))+
  scale_fill_manual(values=c('grey75','#D55E00','grey75','#0072B2'),guide='none')+scale_x_discrete(labels=labsx)+
  labs(x=NULL,y='Signature score')+ST+theme(strip.background=element_rect(fill='grey97'),
  strip.text=ggtext::element_markdown(size=8),axis.text.x=ggtext::element_markdown())
axisdat<-fread(file.path(REV,'molecular_axis_plot_data.csv'))
scatter<-function(v,label){
 d<-axisdat[is.finite(get(v)) & is.finite(risk)]
 ggplot(d,aes(.data[[v]],risk))+geom_point(size=.5,alpha=.35,colour='#0072B2')+
   geom_smooth(method='lm',formula=y~x,se=FALSE,colour='#D55E00',linewidth=.5)+
   labs(x=label,y='Signature score',subtitle=sprintf('Pearson r = %+.3f; n = %d',cor(d[[v]],d$risk),nrow(d)))+ST
}
mu<-fread(file.path(REV,'mutation_stratified_smd.csv'))
mu[,unit:=factor(unit,levels=rev(c(SIG,'SCORE')))];mu[,subgroup:=factor(subgroup,levels=c('BRAFmut','RASmut','driverneg'))]
# Fixed vertical subrows keep labels separate even for nearly identical effects.
mu[,plot_y:=as.numeric(unit)+c(BRAFmut=.30,RASmut=0,driverneg=-.30)[as.character(subgroup)]]
pm<-ggplot(mu,aes(d_cohen,plot_y,colour=subgroup))+geom_vline(xintercept=0,linetype=2,colour='grey60',linewidth=.35)+
  geom_point(size=1.7)+
  geom_text(aes(label=sprintf('%+.2f',d_cohen)),nudge_x=.10,hjust=0,size=2.6,show.legend=FALSE)+
  scale_colour_manual(values=c(BRAFmut='#D55E00',RASmut='#0072B2',driverneg='#009E73'),
    labels=c('<i>BRAF</i> V600E (n = 283)','RAS mutant (n = 58)','Neither recorded (n = 146)'),name=NULL)+
  scale_x_continuous(limits=c(-2.7,4.85),breaks=seq(-2,4,2))+
  scale_y_continuous(breaks=seq_along(levels(mu$unit)),labels=gene_display(levels(mu$unit)),
                     limits=c(.45,5.55),expand=expansion(mult=0))+
  labs(x='Standardised tumour-normal difference',y=NULL,title='Expression across mutation groups')+ST+
  theme(axis.text.y=ggtext::element_markdown(),legend.text=ggtext::element_markdown(size=6.8),
        legend.key.width=unit(.3,'cm'))
fig5<-pa/(scatter('BRAF_RAF_score','BRAF-RAS score')|scatter('differentiation_score','Thyroid differentiation score'))/
  (scatter('ERK_score','ERK-output score')|scatter('purity','Tumour purity'))/pm+
  plot_layout(heights=c(.85,1,1,1.45))+plot_annotation(tag_levels='a')
savefig(fig5,'Fig5',8.2)

# Tissue adjustment leads into broad localisation and donor-level lineage inference.
adj<-fread(file.path(REV,'composition_coefficients_long.csv'))
adj[,unit:=factor(unit,levels=rev(c(SIG,'SCORE')))];adj[,model:=factor(model,levels=c('Unadjusted','Adjusted'))]
adjplot<-function(dd,title){
 ggplot(dd,aes(beta,unit,colour=model))+geom_vline(xintercept=0,linetype=2,colour='grey60',linewidth=.35)+
 geom_errorbar(aes(xmin=lo,xmax=hi),orientation='y',width=.18,position=position_dodge(.45),linewidth=.4)+
 geom_point(size=1.7,position=position_dodge(.45))+
 scale_colour_manual(values=c(Unadjusted='#0072B2',Adjusted='#D55E00'),name=NULL)+
 labs(x='Tumour-normal coefficient',y=NULL,title=title)+ST+scale_y_discrete(labels=gene_display)+
 theme(axis.text.y=ggtext::element_markdown(),legend.text=element_text(size=6.8))
}
U<-fread(file.path(REV,'umap_current_annotations.csv'))
ctcol<-c(Thyrocyte='#D55E00',Endothelial='#009E73',Fibroblast='#CC79A7',`T/NK`='#0072B2',`B/Plasma`='#56B4E9',Myeloid='#7A5195')
ctlev<-names(ctcol);U[,celltype:=factor(celltype,levels=ctlev)]
pu<-ggplot(U,aes(UMAP1,UMAP2,colour=celltype))+geom_point(size=.08,alpha=.65,stroke=0)+
  scale_colour_manual(values=ctcol,name=NULL)+coord_equal()+theme_void(base_family='Arial',base_size=8)+
  labs(title='Cellular landscape')+guides(colour=guide_legend(ncol=3,override.aes=list(size=1.4,alpha=1)))+
  theme(plot.title=element_text(size=9,face='bold'),plot.tag=element_text(size=11,face='bold'),
        legend.position='bottom',legend.text=element_text(size=6.3),legend.key.height=unit(.16,'cm'),
        legend.key.width=unit(.2,'cm'),plot.margin=margin(4,4,4,4))
cs<-fread(file.path(REV,'celltype_detection.csv'));cs[,celltype:=factor(celltype,levels=rev(ctlev))];cs[,gene:=factor(gene,levels=SIG)]
pd<-ggplot(cs,aes(gene,celltype,size=pct_detect,colour=mean_expr))+geom_point()+
 scale_size_area(max_size=4.5,name='Detected (%)',breaks=c(5,15,25))+
 scale_colour_viridis_c(option='C',name='Mean expression',breaks=c(0,.1,.2))+
 labs(x=NULL,y=NULL,title='Transcript distribution')+ST+
 theme(panel.grid=element_blank(),panel.border=element_blank(),axis.ticks=element_blank(),
       axis.text=element_text(size=6.7),axis.text.x=element_text(face='italic'),legend.box='vertical',legend.spacing.y=unit(0,'pt'),legend.margin=margin(0,0,0,0))+
 guides(colour=guide_colourbar(barwidth=unit(2.5,'cm'),barheight=unit(.15,'cm'),title.position='top'),
        size=guide_legend(nrow=1,title.position='top'))
q1<-fread(file.path(REV,'Q1_results.csv'))[analysis=='Q1_primary_independent']
q2<-fread(file.path(REV,'Q2_results.csv'))[analysis=='Q2_primary_independent']
forest<-function(d,title){
 d<-copy(d);d[,gene:=factor(gene,levels=rev(SIG))]
 ggplot(d,aes(logFC,gene))+geom_vline(xintercept=0,linetype=2,colour='grey60',linewidth=.35)+
 geom_errorbar(aes(xmin=CI.L,xmax=CI.R),orientation='y',width=.17,colour='#0072B2',linewidth=.5)+
 geom_point(size=1.9,colour='#0072B2')+
 geom_text(aes(x=2.6,label=sprintf('%.3f',BH4)),hjust=0,size=2.3)+
 annotate('text',x=2.6,y=4.6,label='BH P',size=2.3,hjust=0)+
 scale_x_continuous(limits=c(-1.2,3.35),breaks=c(-1,0,1,2))+scale_y_discrete(expand=expansion(add=c(.55,.85)))+
 labs(x='log2 fold change (95% CI)',y=NULL,title=title)+ST+
 theme(panel.grid.major.y=element_blank(),axis.text.y=element_text(face='italic'))
}
lv<-fread(file.path(REV,'per_unit_expression.csv'))[analysis=='Q1_primary_independent']
donorplot<-function(g,title=g,label=FALSE){
 title_label<-if(identical(title,g))bquote(bolditalic(.(g)))else if(identical(title,'ABHD17C by donor'))bquote(bolditalic('ABHD17C')*' by donor')else title
 dd<-copy(lv);dd[,group:=factor(tissue,levels=c('NORM','PTC'),labels=c('Normal','PTC'))]
 p<-ggplot(dd,aes(group,.data[[g]],colour=group))+
   geom_point(position=position_jitter(width=.1,seed=11),size=2)+
   stat_summary(fun=median,fun.min=median,fun.max=median,geom='crossbar',width=.4,colour='black',linewidth=.4)+
   scale_colour_manual(values=c(Normal='#0072B2',PTC='#D55E00'),guide='none')+
   scale_x_discrete(labels=c('Normal\n(n = 3)','PTC\n(n = 4)'))+
   labs(x=NULL,y='log2 CPM',title=title_label)+ST+theme(panel.grid.major.x=element_blank())
 if(label)p<-p+geom_text(aes(label=sample),position=position_jitter(width=.1,seed=11),nudge_y=.12,size=2.1,show.legend=FALSE)
 p
}
fig6<-((adjplot(adj[unit!='SCORE'],'Composition-adjusted genes')|adjplot(adj[unit=='SCORE'],'Signature score'))+plot_layout(widths=c(1.65,1)))/
  (wrap_elements(full=pu)|wrap_elements(full=pd))/(forest(q1,'Epithelial-lineage differences')|donorplot('ABHD17C','ABHD17C by donor'))+
  plot_layout(heights=c(.8,1.08,.95))+plot_annotation(tag_levels='a')
savefig(fig6,'Fig6',7.0)

# Display transcript detection using retained spot coordinates on a neutral background.
SP<-fread(file.path(REV,'spatial_spot_detection.csv'));QC<-fread(file.path(REV,'spatial_coordinate_qc.csv'))
sections<-fread(file.path(REV,'spatial_section_detection.csv'))
stopifnot(nrow(SP)==52906L,nrow(QC)==16L,nrow(sections)==64L,!anyNA(SP[,.(x,y)]))
sectionpanel<-function(gsm,g,small=FALSE){
 d<-copy(SP[sample==gsm]);rr<-QC[sample==gsm]
 stopifnot(nrow(rr)==1L,nrow(d)>0L,!anyNA(d[[g]]))
 xr<-range(d$x);yr<-range(d$y);pad<-max(diff(xr),diff(yr))*.045
 section_label<-if(small)rr$label else paste0(if(rr$tissue=='NORM')'Normal' else rr$tissue,' (',rr$label,')')
 title<-if(small)bquote(bold(.(section_label))~bolditalic(.(g))) else bquote(bolditalic(.(g)))
 ggplot()+
  geom_point(data=d[!get(g)],aes(x,-y),colour=if(small)'#C8C8C8' else '#DDDDDD',size=if(small).30 else .42,stroke=0)+
  geom_point(data=d[get(g)],aes(x,-y),colour='#C44E00',size=if(small).42 else .68,stroke=0)+
  coord_fixed(xlim=xr+c(-pad,pad),ylim=c(-yr[2]-pad,-yr[1]+pad),expand=FALSE)+
  theme_void(base_family='Arial',base_size=8)+
  labs(title=title,subtitle=if(small)sprintf('Detected %.1f%%\n(%s/%s spots)',100*mean(d[[g]]),
                                  format(sum(d[[g]]),big.mark=','),format(nrow(d),big.mark=',')) else
       sprintf('%s\nDetected %.1f%% (%s/%s)',section_label,100*mean(d[[g]]),
               format(sum(d[[g]]),big.mark=','),format(nrow(d),big.mark=',')))+
  theme(plot.title=element_text(size=if(small)7 else 9,hjust=.5),
        plot.subtitle=element_text(size=if(small)6.3 else 7.5,hjust=.5),
        plot.margin=margin(5,4,3,4))
}
norm<-QC[tissue=='NORM'][order(label)]$sample[1];ptc<-QC[tissue=='PTC'][order(label)]$sample[1]
spgrid<-wrap_plots(sectionpanel(norm,'ABHD17C'),sectionpanel(norm,'ZDHHC20'),
                  sectionpanel(ptc,'ABHD17C'),sectionpanel(ptc,'ZDHHC20'),ncol=2)+
 plot_annotation(caption='Orange: detected (raw count > 0)    Grey: not detected',
                 theme=theme(plot.caption=element_text(family='Arial',size=7,hjust=.5)))
ss<-sections[gene %in% c('ABHD17C','ZDHHC20')];ss[,gene:=factor(gene,levels=c('ABHD17C','ZDHHC20'))]
ss[,group:=factor(tgrp,levels=c('NORM','PTC','LPTC','ATC'),labels=c('Normal','PTC','LPTC','ATC'))]
ps<-ggplot(ss,aes(group,det_pct,colour=group))+geom_point(position=position_jitter(width=.09,height=0,seed=11),size=1.7)+
 stat_summary(fun=median,fun.min=median,fun.max=median,geom='crossbar',width=.4,colour='black',linewidth=.4)+
 facet_wrap(~gene,ncol=1)+scale_colour_manual(values=c(Normal='#0072B2',PTC='#D55E00',LPTC='#009E73',ATC='#CC79A7'),guide='none')+
 scale_y_continuous(limits=c(0,NA),expand=expansion(mult=c(.02,.12)))+
 labs(x=NULL,y='Detected spots (%)',title='Detection by section')+ST+
 theme(panel.grid.major.x=element_blank(),strip.background=element_blank(),strip.text=element_text(face='italic',size=8))
savefig((wrap_elements(full=spgrid)|ps)+plot_layout(widths=c(1.5,1))+plot_annotation(tag_levels='a'),'Fig7',4.9)

# Complete maps use the same retained coordinates and detection states in all 16 sections.
sp_groups<-c('NORM','PTC','LPTC','ATC');all_maps<-list()
for(i in seq_along(sp_groups)){
 ids<-QC[tissue==sp_groups[i]][order(label)]$sample
 panels<-c(lapply(ids,function(s)sectionpanel(s,'ABHD17C',TRUE)),
           lapply(ids,function(s)sectionpanel(s,'ZDHHC20',TRUE)))
 pp<-wrap_plots(panels,ncol=4)+
  plot_annotation(title=c('Normal thyroid','PTC','Locally advanced PTC','ATC')[i],
                  caption='Orange: detected (raw count > 0)    Grey: not detected',
                  theme=theme(plot.title=element_text(family='Arial',size=11,face='bold'),
                              plot.caption=element_text(family='Arial',size=7,hjust=.5)))
 savefig(pp,paste0('Supplementary_Fig_S7',letters[i]),4.4);all_maps[[i]]<-pp
}
cairo_pdf(file.path(F,'Supplementary_Fig_S7.pdf'),width=6.8,height=4.4)
for(pp in all_maps)print(pp)
invisible(dev.off())

loo<-fread(file.path(REV,'Q1_leave_one_donor_out.csv'))[gene=='ABHD17C']
loo[,omitted:=factor(sub('Q1_leave_out_','',analysis),levels=rev(c('7','18','19','1','4','5','6')))]
ploo<-ggplot(loo,aes(logFC,omitted))+geom_vline(xintercept=q1[gene=='ABHD17C',logFC],linetype=2,colour='grey50',linewidth=.4)+
 geom_errorbar(aes(xmin=CI.L,xmax=CI.R),orientation='y',width=.16,linewidth=.4,colour='#0072B2')+geom_point(size=1.8,colour='#0072B2')+
 labs(x=bquote(italic('ABHD17C')*' log2 fold change (95% CI)'),y='Omitted donor',title='Leave-one-donor-out sensitivity')+ST
supp<-(donorplot('ABHD17C')|donorplot('ZDHHC6'))/(donorplot('ABHD17B')|donorplot('ZDHHC20'))/
  (forest(q2,'Malignant versus normal follicular')|ploo)+plot_annotation(tag_levels='a')
savefig(supp,'Supplementary_Fig_S6',7.2)
cat('Saved Figures 5, 6, 7 and Supplementary Figure S6\n')
