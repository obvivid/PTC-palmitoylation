# Correlations with published TCGA molecular scores and purity estimates.
suppressMessages({library(data.table); library(readxl)})
clinical <- fread(file.path(RESULT_DIR,"clinical_immune_fixed.csv"), na.strings=c("","NA"))
published <- extract_published_scores(file.path(PROJECT_ROOT,"data/reference/NIHMS633017-supplement-15.xlsx"))
published[, sid := substr(sample,1,15)]
stopifnot(!anyDuplicated(published$sid))
i <- match(clinical$sid, published$sid)
out <- rbindlist(lapply(c("BRAF_RAF_score","differentiation_score","ERK_score","purity"),function(v) {
  x <- published[[v]][i]; y <- clinical$risk
  k <- is.finite(x) & is.finite(y)
  test <- cor.test(x[k],y[k],method="pearson")
  data.table(variable=v,n=sum(k),pearson_r=unname(test$estimate),
             ci_lo=test$conf.int[1],ci_hi=test$conf.int[2],p=test$p.value)
}))
fwrite(out,file.path(RESULT_DIR,"published_molecular_correlations.csv"))
print(out)
