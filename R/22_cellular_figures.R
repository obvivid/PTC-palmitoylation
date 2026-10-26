# Figure 6 uses cells for distribution and samples/sections for group displays.
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork)})
SC <- file.path(RESULT_DIR,"extended","singlecell")
SP <- file.path(RESULT_DIR,"extended","spatial")
OUT <- file.path(OUTPUT_DIR,"figures")
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
SIG4 <- c("ABHD17C","ZDHHC6","ABHD17B","ZDHHC20")
d <- fread(file.path(SC,"celltype_detection.csv"))
d[,gene:=factor(gene,levels=SIG4)]
d[,celltype:=factor(celltype,levels=rev(c("Thyrocyte","Endothelial","Fibroblast","Squamous","Prolif","T/NK","B/Plasma","Myeloid","Mast")))]
theme_set(theme_classic(base_size=9,base_family="Arial")+theme(plot.title=element_text(size=10,face="bold"),axis.title=element_text(size=9),strip.background=element_blank(),strip.text=element_text(size=9,face="italic"),legend.text=element_text(size=8),legend.title=element_text(size=8)))
pa <- ggplot(d,aes(gene,celltype,size=pct_detect,colour=mean_expr))+geom_point()+scale_size_area(max_size=6,breaks=c(5,15,25),name="Detected (%)")+scale_colour_viridis_c(option="C",name="Mean expression")+labs(x=NULL,y=NULL,title="a  Transcript distribution across cell types")+theme(legend.position="right")
s <- fread(file.path(SC,"thyrocyte_by_sample.csv"))
long <- melt(s,id.vars=c("sample","tissue","cells","median_lib"),measure.vars=SIG4,variable.name="gene",value.name="pct")
long[,tissue:=factor(tissue,levels=c("NORM","PTC"),labels=c("Normal","PTC"))]
long[,gene:=factor(gene,levels=SIG4)]
pb <- ggplot(long,aes(tissue,pct,colour=tissue))+geom_point(position=position_jitter(width=0.10,height=0,seed=9),size=1.6)+stat_summary(fun=median,geom="point",shape=95,size=7,colour="black")+facet_wrap(~gene,nrow=1)+scale_colour_manual(values=c("#0072B2","#D55E00"))+scale_y_continuous(limits=c(0,105),breaks=c(0,25,50,75,100))+labs(x=NULL,y="Thyrocytes detected (%)",title="b  Detection by sample")+theme(legend.position="none",strip.text=element_text(size=7.6))
v <- fread(file.path(SP,"spatial_by_section.csv"))[gene %in% c("ABHD17C","ZDHHC20")]
v[,tgrp:=factor(tgrp,levels=c("NORM","PTC","LPTC","ATC"),labels=c("Normal","PTC","LPTC","ATC"))]
pc <- ggplot(v,aes(tgrp,det_pct,colour=tgrp))+geom_point(position=position_jitter(width=0.09,height=0,seed=9),size=1.9)+stat_summary(fun=median,geom="point",shape=95,size=7,colour="black")+facet_wrap(~gene,nrow=1)+scale_colour_manual(values=c("#0072B2","#D55E00","#009E73","#CC79A7"))+scale_y_continuous(limits=c(0,50),breaks=c(0,10,20,30,40,50))+labs(x=NULL,y="Spots detected (%)",title="c  Detection by tissue section")+theme(legend.position="none")
figure <- (pa / pb / pc)+plot_layout(heights=c(2.4,1.8,1.8))
ggsave(file.path(OUT,"Fig6.pdf"),figure,width=6.8,height=6.9,device=cairo_pdf)
ggsave(file.path(OUT,"Fig6.png"),figure,width=6.8,height=6.9,dpi=300,bg="white")
ggsave(file.path(OUT,"Fig6.tif"),figure,width=6.8,height=6.9,dpi=300,compression="lzw",bg="white")
