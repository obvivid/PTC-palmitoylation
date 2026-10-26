make_tcga_metadata <- function(samples) {
  type <- substr(samples, 14, 15)
  data.table::data.table(sample=samples,
    group=ifelse(type == "01", "Tumor", ifelse(type == "11", "Normal", NA_character_)))
}

read_clinical_annotation <- function(path, samples) {
  cl <- data.table::fread(path, skip="PATIENT_ID",
    na.strings=c("", "NA", "[Not Available]", "[Not Applicable]", "[Unknown]"))
  fields <- c("PATIENT_ID", "ICD_O_3_HISTOLOGY", "AGE", "SEX",
              "AJCC_PATHOLOGIC_TUMOR_STAGE", "PATH_T_STAGE", "PATH_N_STAGE")
  stopifnot(all(fields %in% names(cl)), !anyDuplicated(cl$PATIENT_ID))
  cl <- cl[match(substr(samples,1,12), cl$PATIENT_ID)]
  h <- cl$ICD_O_3_HISTOLOGY
  hist <- ifelse(is.na(h), NA_character_,
    ifelse(h == "8260/3", "PTC_classic", ifelse(h == "8340/3", "FVPTC",
    ifelse(h == "8344/3", "Tall_columnar", ifelse(h == "8290/3", "Hurthle", "Other")))))
  st <- cl$AJCC_PATHOLOGIC_TUMOR_STAGE
  stage <- ifelse(st == "STAGE I", "I", ifelse(st == "STAGE II", "II",
    ifelse(st %in% c("STAGE III","STAGE IIIA","STAGE IIIB","STAGE IIIC"), "III",
    ifelse(grepl("^STAGE IV", st), "IV", NA_character_))))
  t <- cl$PATH_T_STAGE
  t4 <- ifelse(grepl("^T1",t),"T1",ifelse(grepl("^T2",t),"T2",
        ifelse(grepl("^T3",t),"T3",ifelse(grepl("^T4",t),"T4",NA_character_))))
  n <- cl$PATH_N_STAGE
  data.table::data.table(sample=samples, age=as.numeric(cl$AGE), sex=cl$SEX,
    stage4=stage, t4=t4, n_pos=ifelse(n %in% c("N1","N1A","N1B"),1L,
       ifelse(n == "N0",0L,NA_integer_)), n_stage=n, histgrp=hist)
}

extract_published_scores <- function(path) {
  sheets <- readxl::excel_sheets(path)
  sheet <- sheets[grepl("^THCA-TP \\(496\\)", sheets)]
  if (length(sheet) != 1) stop("Expected primary tumor worksheet missing.")
  d <- data.table::as.data.table(readxl::read_excel(path, sheet=sheet))
  fields <- c("sample", "BRAF_RAF_score", "differentiation_score", "ERK_score", "purity")
  stopifnot(all(fields %in% names(d)))
  d[, ..fields]
}


# The unweighted, complete-data duplicateCorrelation calculation with the
# common design decomposition cached across genes. Gene-wise REML fitting,
# trimming and correlation bounds follow limma/statmod unchanged.
patient_correlation <- function(E,design,block,trim=.15,progress=FALSE) {
 stopifnot(is.matrix(E),all(is.finite(E)),ncol(E)==nrow(design),length(block)==ncol(E))
 A<-factor(block);bs<-table(A);mx<-nrow(design)
 if(max(bs)<2L)return(list(consensus.correlation=0,cor=0,atanh.correlations=rep(0,nrow(E))))
 Z<-model.matrix(~0+A)
 fz<-lm.fit(design,Z);r<-fz$rank;mq<-mx-r
 if(length(levels(A))>=mx-1L || mx<=ncol(design)+2L)stop("Insufficient repeated observations")
 sv<-La.svd(fz$effects[(r+1L):mx,,drop=FALSE],nu=mq,nv=0)
 fy<-lm.fit(design,t(E))
 uy<-crossprod(sv$u,fy$effects[(r+1L):mx,,drop=FALSE])
 d<-rep(0,mq);d[seq_along(sv$d)]<-sv$d^2
 dx<-cbind(Residual=1,Block=d);dy<-uy^2
 initial<-lm.fit(dx,dy)
 do_glm<-mq>2L && sum(abs(d)>1e-15)>1L && var(d)>1e-15
 rho<-rep(NA_real_,nrow(E))
 for(i in seq_len(nrow(E))){
  cc<-initial$coefficients[,i]
  if(do_glm){
   start<-if(all(initial$fitted.values[,i]>=0))cc else c(Residual=mean(dy[,i]),Block=0)
   cc<-tryCatch(suppressWarnings(statmod::glmgam.fit(dx,dy[,i],coef.start=start,tol=1e-6,maxit=20)$coefficients),error=function(e)c(NA_real_,NA_real_))
  }
  if(!is.na(cc[1]))rho[i]<-cc[2]/sum(cc)
  if(progress && i%%5000L==0L)cat("Correlation genes",i,"of",nrow(E),"\n")
 }
 rho<-pmin(.99,pmax(1/(1-max(bs))+.01,rho));arho<-atanh(rho)
 consensus<-tanh(mean(arho,trim=trim,na.rm=TRUE))
 list(consensus.correlation=consensus,cor=consensus,atanh.correlations=arho)
}


# Evaluate both one-sided tails directly. R 4.6.0's asymptotic two-sided
# implementation may cancel the upper tail in 1 - p for large positive W.
stable_wilcox_test <- function(..., alternative="two.sided") {
  if (alternative != "two.sided") return(stats::wilcox.test(...,alternative=alternative))
  out <- stats::wilcox.test(...,alternative="two.sided")
  lo <- stats::wilcox.test(...,alternative="less")$p.value
  hi <- stats::wilcox.test(...,alternative="greater")$p.value
  out$p.value <- min(1,2*min(lo,hi))
  out
}

# HC1 clustered covariance with finite-cluster correction and t reference.
patient_cluster_coefficient <- function(fit,patient,term="grpTumor") {
  V <- sandwich::vcovCL(fit,cluster=patient,type="HC1",cadjust=TRUE)
  b <- unname(coef(fit)[term]);se <- sqrt(V[term,term]);df <- length(unique(patient))-1L
  c(beta=b,se=se,lo=b-qt(.975,df)*se,hi=b+qt(.975,df)*se,
    p=2*pt(abs(b/se),df=df,lower.tail=FALSE),patients=df+1L)
}
