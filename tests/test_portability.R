# Parse every module and test privacy-safe synthetic clinical and barcode fixtures.
args <- commandArgs(trailingOnly=FALSE)
f <- sub("^--file=","",grep("^--file=",args,value=TRUE)[1])
root <- dirname(dirname(normalizePath(f,mustWork=TRUE)))
for (p in list.files(root,pattern="\\.R$",recursive=TRUE,full.names=TRUE)) parse(p)
source(file.path(root,"R/common.R"))
samples <- c("TCGA-AA-0001-01A","TCGA-AA-0001-11A","TCGA-AA-0002-06A")
m <- make_tcga_metadata(samples)
stopifnot(identical(m$group,c("Tumor","Normal",NA_character_)))
input <- tempfile(fileext=".txt")
clinical <- data.frame(PATIENT_ID=c("TCGA-AA-0001","TCGA-AA-0002"),
  ICD_O_3_HISTOLOGY=c("8260/3",NA),AGE=c(42,NA),SEX=c("FEMALE",NA),
  AJCC_PATHOLOGIC_TUMOR_STAGE=c("STAGE IVA",NA),PATH_T_STAGE=c("T3B",NA),PATH_N_STAGE=c("NX",NA))
write.table(clinical,input,sep="\t",row.names=FALSE,quote=FALSE,na="[Not Available]")
d <- read_clinical_annotation(input,samples)
stopifnot(nrow(d)==3,d$stage4[1]=="IV",d$t4[1]=="T3",is.na(d$n_pos[1]),
          d$histgrp[1]=="PTC_classic",is.na(d$histgrp[3]),is.na(d$age[3]))
unlink(input)
library(data.table)
for (filename in c("04_nested_cv.R","12_histology_sensitivity.R")) {
  expressions <- parse(file.path(root,"R",filename))
  fold_function <- Filter(function(e) is.call(e) && identical(e[[1]],as.name("<-")) &&
    identical(e[[2]],as.name("make_folds")),as.list(expressions))
  stopifnot(length(fold_function)==1)
  e <- new.env()
  e$N_FOLD <- 10
  e$meta <- data.table(pid=c(rep(sprintf("P%02d",1:20),2),sprintf("P%02d",21:40)),
                      group=c(rep("Tumor",20),rep("Normal",20),rep("Tumor",20)))
  eval(fold_function[[1]],envir=e)
  for (seed in 1001:1005) {
    f <- e$make_folds(seed)
    stopifnot(length(f)==nrow(e$meta), all(f %in% 1:10),
              all(vapply(split(f,e$meta$pid),function(x)length(unique(x))==1,logical(1))),
              identical(f,e$make_folds(seed)))
  }
}
cat("R syntax and synthetic portability checks passed.\n")
