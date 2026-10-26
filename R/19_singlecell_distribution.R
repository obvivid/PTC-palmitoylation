# Transcript distribution and descriptive detection summaries.
suppressPackageStartupMessages({library(data.table);library(Matrix)})
SC <- file.path(RESULT_DIR,"extended","singlecell")
L <- readRDS(file.path(SC,"counts_annotated_193581.rds"))
C <- L$counts; CT <- L$celltype; TT <- as.character(L$tissue)
SP <- sub("_.*$","",as.character(L$sample))
stopifnot(ncol(C)==length(CT),length(CT)==length(TT),length(SP)==ncol(C),!anyNA(CT))
SIG4 <- c("ABHD17C","ZDHHC6","ABHD17B","ZDHHC20")
stopifnot(all(SIG4 %in% rownames(C)))
lib <- Matrix::colSums(C)
N <- log1p(C[SIG4,,drop=FALSE] %*% Matrix::Diagonal(x=1e4/pmax(lib,1)))
out <- rbindlist(lapply(SIG4,function(g) data.table(gene=g,celltype=CT,det=as.numeric(C[g,])>0,v=as.numeric(N[g,]))))
out <- out[,.(pct_detect=round(100*mean(det),2),mean_expr=round(mean(v),4),n=.N),by=.(gene,celltype)]
fwrite(out,file.path(SC,"celltype_detection.csv"))
k <- CT=="Thyrocyte" & TT %in% c("NORM","PTC")
ST <- rbindlist(lapply(sort(unique(SP[TT %in% c("NORM","PTC")])),function(s){
  j <- k & SP==s
  r <- data.table(sample=s,tissue=TT[which(SP==s)[1]],cells=sum(j),median_lib=if(any(j))median(lib[j])else NA_real_)
  for(g in SIG4)r[[g]] <- if(any(j))100*mean(as.numeric(C[g,j])>0)else NA_real_
  r
}))
setorder(ST,tissue,-cells)
fwrite(ST,file.path(SC,"thyrocyte_by_sample.csv"))
# Two normal/PTC tissue pairs share donors; all tissue rows are descriptive.
summary <- rbindlist(lapply(SIG4,function(g){
  x <- ST[tissue=="NORM" & cells>0][[g]]; y <- ST[tissue=="PTC" & cells>0][[g]]
  data.table(gene=g,NORM_n=length(x),PTC_n=length(y),NORM_median=median(x),PTC_median=median(y),NORM_range=sprintf("%.1f-%.1f",min(x),max(x)),PTC_range=sprintf("%.1f-%.1f",min(y),max(y)))
}))
fwrite(summary,file.path(SC,"thyrocyte_detection_summary.csv"))
print(summary)
