DATA_DIR <- file.path(PROJECT_ROOT, "data")
# Section-level descriptive detection, with no spot-level group tests.
suppressPackageStartupMessages({library(data.table);library(Matrix)})
SC <- file.path(DATA_DIR,"spatial","GSE250521")
OUT <- file.path(RESULT_DIR,"extended","spatial")
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
files <- list.files(SC,pattern="_visium_matrix\\.mtx\\.gz$")
if(length(files)!=16)stop("Expected 16 public Visium matrices.")
SIG4 <- c("ABHD17C","ZDHHC6","ABHD17B","ZDHHC20")
rows <- list(); qc <- list()
for(f in files){
  nm <- sub("_visium_matrix\\.mtx\\.gz$","",f)
  gsm <- sub("_.*$","",nm); label <- sub("^GSM[0-9]+_","",nm)
  features <- fread(file.path(SC,sub("_matrix\\.mtx\\.gz$","_features.tsv.gz",f)),header=FALSE,sep="\t")
  symbols <- as.character(if(ncol(features)>=2 && any(grepl("^[A-Za-z]",features[[2]])))features[[2]]else features[[1]])
  m <- readMM(gzfile(file.path(SC,f))); stopifnot(nrow(m)==length(symbols))
  rownames(m) <- symbols
  keep <- Matrix::colSums(m)>=500 & Matrix::colSums(m>0)>=200
  tgrp <- if(grepl("^N-",label))"NORM" else if(grepl("^LPTC",label))"LPTC" else if(grepl("^PTC",label))"PTC" else if(grepl("^ATC",label))"ATC" else stop("Unrecognised tissue label.")
  qc[[f]] <- data.table(sample=gsm,label=label,tissue=tgrp,spots=ncol(m),kept=sum(keep))
  rows[[f]] <- rbindlist(lapply(SIG4,function(g){
    if(!g %in% symbols)stop("Signature gene absent from spatial features.")
    v <- as.numeric(m[match(g,symbols),keep])
    data.table(gene=g,tgrp=tgrp,sample=gsm,spots=sum(keep),det_pct=100*mean(v>0))
  }))
}
fwrite(rbindlist(qc),file.path(OUT,"spatial_qc.csv"))
fwrite(rbindlist(rows),file.path(OUT,"spatial_by_section.csv"))
