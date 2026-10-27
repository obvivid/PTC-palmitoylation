# Coordinated manuscript figures and complete cellular/spatial visualisation.
# Run as step 30 after the prerequisites in config/steps.csv.
# Curated plotting tables reproduce the published coordinates and annotations;
# raw preprocessing remains available in steps 18-21.
suppressPackageStartupMessages({library(data.table);library(Matrix);library(ggplot2);library(patchwork);library(jsonlite);library(png)})
DEST <- OUTPUT_DIR
MAIN <- file.path(DEST,"figures")
DAT <- file.path(DEST,"figure_source_data")
dir.create(MAIN,recursive=TRUE,showWarnings=FALSE)
dir.create(DAT,recursive=TRUE,showWarnings=FALSE)
SIG4 <- c("ABHD17C","ZDHHC6","ABHD17B","ZDHHC20")
SUPSRC <- file.path(SOURCE_ROOT,"data/figure_tables/revision")
U <- fread(file.path(SUPSRC,"umap_current_annotations.csv"))
cell_summary <- fread(file.path(SUPSRC,"celltype_detection.csv"))
SP <- fread(file.path(SUPSRC,"spatial_spot_detection.csv"))
QC <- fread(file.path(SUPSRC,"spatial_coordinate_qc.csv"))
sections <- fread(file.path(SUPSRC,"spatial_section_detection.csv"))
stopifnot(nrow(U)==67672L,!anyDuplicated(U$cell),nrow(SP)==52906L,nrow(QC)==16L,nrow(sections)==64L)
# Spatial rendering uses curated retained coordinates and no tissue-image input.
CHECK <- list(umap_cells=nrow(U),spatial_spots=nrow(SP),spatial_sections=nrow(QC),curated_coordinates=TRUE)
# Reconstruct plot objects using R/23_publication_figures.R;
# feature selection and model fitting are not repeated.
cat("Loading publication plot objects\n")
original_output <- OUTPUT_DIR
OUTPUT_DIR <- file.path(DEST,"revision_baseline")
source(file.path(SOURCE_ROOT,"R/23_publication_figures.R"),local=environment(),encoding="UTF-8")
expr[,gene:=factor(gene,levels=SIG4)]
px <- box(expr,"log2(CPM + 1)")+facet_wrap(~gene,nrow=1)
fig3 <- (f1a|f1b)/((px|ps)+plot_layout(widths=c(3,1)))/(f2a|f2b)+
 plot_layout(heights=c(1,.8,1.1))+plot_annotation(tag_levels="a")
export(fig3,"Fig3",8.0)
OUTPUT_DIR <- original_output
for(map in list(c(1,1),c(3,2),c(4,3),c(5,5)))for(ext in c("png","pdf","tif")){
  from <- file.path(DEST,"revision_baseline/figures",paste0("Fig",map[1],".",ext))
  to <- file.path(MAIN,paste0("Fig",map[2],".",ext))
  stopifnot(file.copy(from,to,overwrite=TRUE))
}
for(n in c("Supplementary_Fig_S1","Supplementary_Fig_S4","Supplementary_Fig_S3"))for(ext in c("png","pdf","tif"))
 file.copy(file.path(DEST,"revision_baseline/figures",paste0(n,".",ext)),file.path(MAIN,paste0(n,".",ext)),overwrite=TRUE)
STYLE <- theme_bw(base_family="Arial",base_size=9)+theme(panel.grid.minor=element_blank(),panel.grid.major=element_line(colour="grey93",linewidth=.25),axis.text=element_text(colour="grey20",size=8),plot.title=element_text(size=10,face="bold"),plot.subtitle=element_text(size=8),plot.tag=element_text(size=12,face="bold"),legend.position="bottom",legend.text=element_text(size=7),legend.title=element_text(size=7),plot.margin=margin(5,6,5,5))
save_new <- function(p,name,h,w=6.8){
  ggsave(file.path(MAIN,paste0(name,".png")),p,width=w,height=h,dpi=300,bg="white",limitsize=FALSE)
  ggsave(file.path(MAIN,paste0(name,".pdf")),p,width=w,height=h,device=cairo_pdf,bg="white",limitsize=FALSE)
  ggsave(file.path(MAIN,paste0(name,".tif")),p,width=w,height=h,dpi=600,compression="lzw",bg="white",limitsize=FALSE)
}

cat("Building functional background figure\n")
go2 <- go+STYLE+labs(title="GO enrichment")+scale_x_continuous(expand=expansion(mult=c(.04,.12)))+scale_y_discrete(labels=function(z)vapply(z,function(s)paste(strwrap(s,32),collapse="\n"),""))+theme(axis.text.y=element_text(size=7),legend.box="vertical",legend.title=element_text(size=7),legend.text=element_text(size=6.5))
net2 <- ggplot()+geom_curve(data=edge,aes(x,y,xend=xend,yend=yend),curvature=.045,linewidth=.22,alpha=.14)+
 geom_point(data=node,aes(x,y,size=degree,colour=signature))+
 geom_text(data=node,aes(label_x,label_y,label=gene,hjust=label_hjust,vjust=label_vjust,fontface=ifelse(signature,"bold","plain")),size=2.25)+
 scale_colour_manual(values=c(`FALSE`="#0072B2",`TRUE`="#D55E00"),guide="none")+scale_size_continuous(range=c(1.6,2.8),guide="none")+
 coord_equal(xlim=c(-1.7,1.7),ylim=c(-1.3,1.3),clip="off")+theme_void(base_family="Arial",base_size=8)+
 labs(title="STRING associations",subtitle="25 genes; 200 associations",caption="Orange nodes: four selected genes\nNode size: degree")+
 theme(plot.title=element_text(size=10,face="bold"),plot.subtitle=element_text(size=8),plot.caption=element_text(size=7,hjust=.5),plot.margin=margin(5,6,5,3),plot.tag=element_text(size=12,face="bold"))
save_new((wrap_elements(full=go2)|wrap_elements(full=net2))+plot_layout(widths=c(.95,1.05))+plot_annotation(tag_levels="a"),"Fig4",4.5)
save_new(kg+STYLE,"Supplementary_Fig_S2",3.7)

cat("Building cellular and spatial figures\n")
ctcol <- c(Thyrocyte="#D55E00",Endothelial="#009E73",Fibroblast="#CC79A7",`T/NK`="#0072B2",`B/Plasma`="#56B4E9",Myeloid="#7A5195")
ctlev <- c("Thyrocyte","Endothelial","Fibroblast","T/NK","B/Plasma","Myeloid")
stopifnot(setequal(unique(U$celltype),names(ctcol)))
U[,celltype:=factor(celltype,levels=ctlev)]
p_um <- ggplot(U,aes(UMAP1,UMAP2,colour=celltype))+geom_point(size=.11,alpha=.7,stroke=0)+
 scale_colour_manual(values=ctcol,name=NULL)+coord_equal()+theme_void(base_family="Arial",base_size=8)+
 labs(title="Cellular landscape")+guides(colour=guide_legend(ncol=3,override.aes=list(size=1.7,alpha=1)))+
 theme(plot.title=element_text(size=10,face="bold"),legend.position="bottom",legend.text=element_text(size=6.5),legend.key.height=unit(.16,"cm"),legend.key.width=unit(.22,"cm"),plot.margin=margin(4,4,4,4),plot.tag=element_text(size=12,face="bold"))
cell_summary[,celltype:=factor(celltype,levels=rev(ctlev))];cell_summary[,gene:=factor(gene,levels=SIG4)]
p_dot <- ggplot(cell_summary,aes(gene,celltype,size=pct_detect,colour=mean_expr))+geom_point()+
 scale_size_area(max_size=4.5,name="Detected (%)",breaks=c(5,15,25))+scale_colour_viridis_c(option="C",name="Mean expression",breaks=c(0,.1,.2))+
 labs(x=NULL,y=NULL,title="Transcript distribution")+STYLE+
 theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_blank(),axis.ticks=element_blank(),axis.text=element_text(size=7),legend.box="vertical",legend.margin=margin(0,0,0,0),legend.spacing.y=unit(0,"pt"))+
 guides(colour=guide_colourbar(barwidth=unit(2.6,"cm"),barheight=unit(.17,"cm"),title.position="top"),size=guide_legend(nrow=1,title.position="top"))
map_plot<-function(gsm,g,small=FALSE){
 d<-copy(SP[sample==gsm]);rr<-QC[sample==gsm]
 stopifnot(nrow(rr)==1L,nrow(d)>0L,!anyNA(d[[g]]))
 xr<-range(d$x);yr<-range(d$y);pad<-max(diff(xr),diff(yr))*.045
 section_label<-if(small)rr$label else paste0(if(rr$tissue=='NORM')'Normal' else rr$tissue,' (',rr$label,')')
 title<-bquote(bold(.(section_label))~bolditalic(.(g)))
 ggplot()+
  geom_point(data=d[!get(g)],aes(x,-y),colour=if(small)'#C8C8C8' else '#DDDDDD',size=if(small).30 else .42,stroke=0)+
  geom_point(data=d[get(g)],aes(x,-y),colour='#C44E00',size=if(small).42 else .68,stroke=0)+
  coord_fixed(xlim=xr+c(-pad,pad),ylim=c(-yr[2]-pad,-yr[1]+pad),expand=FALSE)+
  theme_void(base_family='Arial',base_size=8)+
  labs(title=title,subtitle=sprintf('Detected %.1f%%\n(%s/%s spots)',100*mean(d[[g]]),
                                  format(sum(d[[g]]),big.mark=','),format(nrow(d),big.mark=',')))+
  theme(plot.title=element_text(size=if(small)7 else 9,hjust=.5),
        plot.subtitle=element_text(size=if(small)6.3 else 7.5,hjust=.5),
        plot.margin=margin(5,4,3,4))
}
normal <- QC[tissue=="NORM"][order(label)]$sample[1];ptc <- QC[tissue=="PTC"][order(label)]$sample[1]
CHECK$main_spatial_sections <- c(normal,ptc)
# Intermediate cellular figure; step 34 writes Figures 6 and 7 used in the paper.
spatial_grid <- wrap_plots(map_plot(normal,"ABHD17C"),map_plot(normal,"ZDHHC20"),
                          map_plot(ptc,"ABHD17C"),map_plot(ptc,"ZDHHC20"),ncol=2)+
 plot_annotation(caption="Orange: detected (raw count > 0)    Grey: not detected",
                 theme=theme(plot.caption=element_text(family="Arial",size=7,hjust=.5)))
section_data <- copy(sections[gene %in% c("ABHD17C","ZDHHC20")])
section_data[,gene:=factor(gene,levels=c("ABHD17C","ZDHHC20"))]
section_data[,group:=factor(tgrp,levels=c("NORM","PTC","LPTC","ATC"),labels=c("Normal","PTC","LPTC","ATC"))]
p_sections <- ggplot(section_data,aes(group,det_pct,colour=group))+
 geom_point(position=position_jitter(width=.09,height=0,seed=11),size=1.7)+
 stat_summary(fun=median,fun.min=median,fun.max=median,geom="crossbar",width=.4,colour="black",linewidth=.35)+
 facet_wrap(~gene,nrow=1)+scale_colour_manual(values=c("Normal"="#0072B2","PTC"="#D55E00","LPTC"="#009E73","ATC"="#CC79A7"),guide="none")+
 scale_y_continuous(limits=c(0,NA),expand=expansion(mult=c(.02,.12)))+
 labs(x=NULL,y="Detected spots (%)",title="Detection across tissue sections")+STYLE+
 theme(panel.grid.major.x=element_blank(),strip.background=element_blank(),strip.text=element_text(face="italic",size=8),plot.title=element_text(size=9,face="bold"),axis.text=element_text(size=8))
fig6 <- (wrap_elements(full=p_um)|wrap_elements(full=p_dot))/wrap_elements(full=spatial_grid)/p_sections+
 plot_layout(heights=c(1.05,1.8,.77))+plot_annotation(tag_levels="a")
save_new(fig6,"Fig6",7.8)

groups <- c("NORM","PTC","LPTC","ATC");all_maps <- list()
for(i in seq_along(groups)){
 ids <- QC[tissue==groups[i]][order(label)]$sample
 pp <- c(lapply(ids,function(s)map_plot(s,"ABHD17C",TRUE)),lapply(ids,function(s)map_plot(s,"ZDHHC20",TRUE)))
 pp <- wrap_plots(pp,ncol=4)+plot_annotation(title=c("Normal thyroid","PTC","Locally advanced PTC","ATC")[i],
     caption="Orange: detected (raw count > 0)    Grey: not detected",
     theme=theme(plot.title=element_text(family="Arial",size=11,face="bold"),
                 plot.caption=element_text(family="Arial",size=7,hjust=.5)))
 save_new(pp,paste0("Supplementary_Fig_S7",letters[i]),4.4);all_maps[[i]]<-pp
}
cairo_pdf(file.path(MAIN,"Supplementary_Fig_S7.pdf"),width=6.8,height=4.4);for(pp in all_maps)print(pp);dev.off()
ts <- fread(file.path(SUPSRC,"Supplementary_Table_S13_thyrocyte_samples.csv"))
tl <- melt(ts,id.vars=intersect(c("sample","tissue","cells","median_lib"),names(ts)),measure.vars=SIG4,variable.name="gene",value.name="detection")
tl[,group:=factor(tissue,levels=c("NORM","PTC"),labels=c("Normal","PTC"))]
q1 <- ggplot(tl[is.finite(detection)],aes(group,detection,colour=group))+geom_point(aes(size=cells),position=position_jitter(width=.1,seed=11))+
 scale_size_continuous(name="Thyrocytes",range=c(.8,2.2),breaks=c(10,100,1000,3000))+
 stat_summary(fun=median,fun.min=median,fun.max=median,geom="crossbar",width=.35,colour="black",linewidth=.35)+facet_wrap(~gene,nrow=1)+scale_colour_manual(values=c(Normal="#0072B2",PTC="#D55E00"),guide="none")+
 labs(x=NULL,y="Thyrocytes detected (%)",title="Detection by sample")+STYLE+theme(strip.background=element_blank(),strip.text=element_text(face="italic",size=8))
sd <- copy(sections[gene %in% c("ABHD17C","ZDHHC20")]);sd[,group:=factor(tgrp,levels=groups,labels=c("Normal","PTC","LPTC","ATC"))]
q2 <- ggplot(sd,aes(group,det_pct,colour=group))+geom_point(position=position_jitter(width=.08,seed=11),size=1.5)+
 stat_summary(fun=median,fun.min=median,fun.max=median,geom="crossbar",width=.35,colour="black",linewidth=.35)+facet_wrap(~gene,nrow=1)+scale_colour_manual(values=c("#0072B2","#D55E00","#009E73","#CC79A7"),guide="none")+
 labs(x=NULL,y="Spots detected (%)",title="Detection by tissue section")+STYLE+theme(strip.background=element_blank(),strip.text=element_text(face="italic",size=8))
save_new((q1/q2)+plot_layout(heights=c(1,1))+plot_annotation(tag_levels="a"),"Supplementary_Fig_S5",4.9)

cat("Building adjustment and mutation-stratified effect figure\n")
genes_order <- rev(c(SIG4,"SCORE"))
md2 <- copy(md);md2[,unit:=factor(unit,levels=genes_order)]
md2[,subgroup:=factor(subgroup,levels=c("BRAFmut","RASmut","driverneg"))]
mp2 <- ggplot(md2,aes(d_cohen,unit,colour=subgroup))+geom_vline(xintercept=0,linetype=2,colour="grey55",linewidth=.4)+
 geom_point(position=position_dodge(width=.55),size=2)+
 geom_text(aes(label=sprintf("%+.2f",d_cohen)),position=position_dodge(width=.55),hjust=-.3,size=2.5,show.legend=FALSE)+
 scale_colour_manual(values=c(BRAFmut="#D55E00",RASmut="#0072B2",driverneg="#009E73"),labels=c("BRAF V600E (n = 283)","RAS mutant (n = 58)","No recorded BRAF/RAS (n = 146)"),name=NULL)+
 scale_x_continuous(limits=c(-2.7,5.0),breaks=seq(-2,4,2))+labs(x="Standardised mean difference",y=NULL,title="Mutation-stratified tumour-normal differences")+STYLE+
 theme(legend.text=element_text(size=6.7))
adj1 <- adjplot(al[unit!="SCORE"],"Signature genes")+STYLE
adj2 <- adjplot(al[unit=="SCORE"],"Signature score")+STYLE
save_new((adj1|adj2)/mp2+plot_layout(heights=c(1,1.25))+plot_annotation(tag_levels="a"),"Fig7",5.9)

cat("Building compact histological comparison figure\n")
atc2 <- atc+STYLE+theme(legend.position="inside",legend.position.inside=c(.98,.03),legend.justification=c(1,0),legend.text=element_text(size=6),legend.key.width=unit(.25,"cm"),plot.title=element_text(size=9))+
 guides(colour=guide_legend(ncol=1))
dist82 <- f4b+STYLE+theme(plot.subtitle=element_blank())+labs(title="Follicular tumours",x=NULL)+scale_x_discrete(labels=c(FTA="Adenoma\n(n = 25)",FTC="Carcinoma\n(n = 27)"))
roc82 <- f4c+STYLE+labs(title="Carcinoma versus adenoma",subtitle="AUC 0.708")+theme(plot.title=element_text(size=8.5))
save_new((atc2|dist82|roc82)+plot_layout(widths=c(1,1,1))+plot_annotation(tag_levels="a"),"Supplementary_Fig_S8",2.85)
fwrite(md2,file.path(DAT,"mutation_stratified_smd.csv"));fwrite(al,file.path(DAT,"composition_coefficients_long.csv"))
fwrite(data.table(dataset=names(aucs),auc=as.numeric(aucs)),file.path(DAT,"external_auc_check.csv"))
for(n in c("external_gene_distribution_plot_data.csv","external_auc_plot_data.csv","molecular_axis_plot_data.csv"))file.copy(file.path(DEST,"revision_baseline/figures",n),DAT,overwrite=TRUE)
CHECK$external_auc <- as.list(aucs)
CHECK$atc_auc <- c(GSE33630=as.numeric(auc(r30)),GSE65144=as.numeric(auc(r65)))
CHECK$ftc_fta_auc <- as.numeric(auc(ro82))
write_json(CHECK,file.path(DEST,"figure_data_checks.json"),pretty=TRUE,auto_unbox=TRUE)
cat("COMPLETED\n");print(CHECK)
