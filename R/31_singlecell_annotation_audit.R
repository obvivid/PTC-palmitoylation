# Reconstruct the recorded broad marker assignment and descriptive sample checks.
suppressPackageStartupMessages({library(data.table);library(Matrix);library(ggplot2)})
SC <- file.path(RESULT_DIR,"extended/singlecell")
L <- readRDS(file.path(SC,"counts_annotated_193581.rds"))
if(is.null(L$cluster))stop("Run step 18 to generate processed counts with cluster identities.")
C <- L$counts;cl <- as.character(L$cluster);ct <- as.character(L$celltype)
stopifnot(ncol(C)==length(cl),ncol(C)==length(ct),!anyNA(cl),!anyNA(ct))
MARK<-list(`T/NK`=c("CD3D","CD3E","CD2","TRAC","IL7R","NKG7","KLRD1"),`B/Plasma`=c("MS4A1","CD79A","CD79B","IGKC","MZB1","JCHAIN","IGHG1"),Myeloid=c("LYZ","CD68","CD14","AIF1","ITGAX","C1QA","FCN1"),Mast=c("TPSAB1","TPSB2","MS4A2","CPA3","KIT"),Endothelial=c("PECAM1","VWF","EGFL7","PLVAP","CLDN5"),Fibroblast=c("COL1A1","COL1A2","DCN","FAP","PDGFRA","LUM"),Thyrocyte=c("TG","TPO","SLC5A5","TSHR","PAX8","FOXE1","SLC26A4","DIO1","DIO2"),Prolif=c("MKI67","TOP2A","PCNA","CCNB1","BIRC5"),Squamous=c("KRT5","KRT14","KRT6A","KRT17","S100A2"))
g<-intersect(unique(unlist(MARK)),rownames(C));lib<-Matrix::colSums(C)
K<-as.matrix(C[g,,drop=FALSE]);N<-sweep(K,2,pmax(lib,1),"/")*1e4
clusters<-as.character(sort(as.integer(unique(cl))))
av<-sapply(clusters,function(k)rowMeans(N[,cl==k,drop=FALSE]));rownames(av)<-g
scores<-sapply(MARK,function(mk)colMeans(av[intersect(mk,g),,drop=FALSE]));rownames(scores)<-clusters
labels<-colnames(scores)[max.col(scores,ties.method="first")]
review<-rbindlist(lapply(seq_along(clusters),function(i){k<-clusters[i];ex<-unique(ct[cl==k]);v<-sort(scores[i,],decreasing=TRUE);data.table(cluster=k,cells=sum(cl==k),assigned_label=paste(ex,collapse=";"),reconstructed_label=labels[i],label_matches=length(ex)==1&&ex==labels[i],best_marker_score=v[1],second_label=names(v)[2],relative_gap=(v[1]-v[2])/pmax(abs(v[1]),1e-12))}))
stopifnot(all(review$label_matches))
fwrite(review,file.path(SC,"cluster_annotation_review.csv"))
fwrite(data.table(cluster=rownames(scores),as.data.table(scores)),file.path(SC,"cluster_marker_scores.csv"))
show<-intersect(c("CD3D","NKG7","MS4A1","JCHAIN","IGKC","IGHG1","LYZ","C1QA","TPSAB1","CPA3","PECAM1","VWF","COL1A1","DCN","TG","TPO","MKI67","TOP2A","KRT5","KRT14"),g)
pd<-rbindlist(lapply(clusters,function(k){j<-cl==k;data.table(cluster=k,gene=show,pct=100*rowMeans(K[show,j,drop=FALSE]>0),mean_log=rowMeans(log1p(N[show,j,drop=FALSE])))}))
pd[,scaled_mean:=as.numeric(scale(mean_log)),by=gene];pd[!is.finite(scaled_mean),scaled_mean:=0]
pd[,cluster_label:=paste0(cluster,"  ",labels[match(cluster,clusters)])]
fwrite(pd,file.path(SC,"cluster_marker_plot_data.csv"))
ts<-fread(file.path(SC,"thyrocyte_by_sample.csv"));SIG4<-c("ABHD17C","ZDHHC6","ABHD17B","ZDHHC20")
sens<-rbindlist(lapply(c(0L,50L,100L),function(nmin)rbindlist(lapply(SIG4,function(g){d<-ts[cells>=max(1L,nmin)];data.table(minimum_thyrocytes=nmin,gene=g,normal_samples=sum(d$tissue=="NORM"),PTC_samples=sum(d$tissue=="PTC"),normal_median=median(d[tissue=="NORM"][[g]]),PTC_median=median(d[tissue=="PTC"][[g]]))}))))
depth<-rbindlist(lapply(SIG4,function(g){z<-ts[cells>0];data.table(gene=g,samples=nrow(z),spearman_depth_detection=cor(z$median_lib,z[[g]],method="spearman"))}))
fwrite(sens,file.path(SC,"low_count_sensitivity.csv"));fwrite(depth,file.path(SC,"depth_detection_summary.csv"))
pd[,gene:=factor(gene,levels=show)];pd[,cluster_label:=factor(cluster_label,levels=rev(paste0(clusters,"  ",labels)))]
p<-ggplot(pd,aes(gene,cluster_label,size=pct,colour=scaled_mean))+geom_point()+
 scale_colour_gradient2(low="#0072B2",mid="#EEEEEE",high="#D55E00",midpoint=0,name="Scaled mean\nlog expression")+
 scale_size_continuous(range=c(.15,3.8),limits=c(0,100),breaks=c(0,25,50,75,100),name="Detected (%)")+
 labs(x=NULL,y="Cluster and assigned broad cell type")+theme_bw(base_family="Arial",base_size=8)+
 theme(panel.grid.major=element_line(colour="grey94",linewidth=.3),panel.grid.minor=element_blank(),axis.text.x=element_text(angle=55,hjust=1,size=7),axis.text.y=element_text(size=7),legend.position="bottom",legend.title=element_text(size=8))
F<-file.path(OUTPUT_DIR,"figures");dir.create(F,recursive=TRUE,showWarnings=FALSE)
ggsave(file.path(F,"Supplementary_Fig_S11.png"),p,width=6.8,height=7.2,dpi=300,bg="white")
ggsave(file.path(F,"Supplementary_Fig_S11.tif"),p,width=6.8,height=7.2,dpi=600,compression="lzw",bg="white")
ggsave(file.path(F,"Supplementary_Fig_S11.pdf"),p,width=6.8,height=7.2,device=cairo_pdf,bg="white")
print(review)
