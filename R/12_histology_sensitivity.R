# Histology exclusion and nested CV sensitivity
# Run through run_analysis.R; see docs/METHODS.md.

suppressMessages({ library(data.table); library(limma); library(glmnet)
                   library(randomForest); library(e1071); library(pROC) })

ROOT <- PROJECT_ROOT
FM   <- MATRIX_DIR
DIR  <- RESULT_DIR
sec  <- function(t){ cat("\n"); cat(strrep("=",78),"\n",sep=""); cat(t,"\n"); cat(strrep("=",78),"\n",sep="") }

GENES <- c("ABHD17C","ZDHHC6","ABHD17B","ZDHHC20")
GSET  <- fread(file.path(DIR, "palmitoylation_geneset_full.csv"))$gene
N_REP <- 5; N_FOLD <- 10

sci <- fread(file.path(DIR, "clinical_immune_fixed.csv"), na.strings=c("","NA"))
exc <- sci[histgrp %in% c("Hurthle","Other")]
clin <- fread(file.path(ROOT, "data", "tcga", "data_clinical_patient.txt"),
              skip="PATIENT_ID")
icd_col <- "ICD_O_3_HISTOLOGY"
stopifnot(icd_col %in% names(clin))
exc[, pid12 := substr(sid, 1, 12)]
exc[, icd := clin[[icd_col]][match(pid12, clin$PATIENT_ID)]]
sec("The 7 excluded cases (histgrp = Hurthle/Other)")
print(exc[, .(sid, sample, histgrp, icd)])


E    <- readRDS(file.path(FM, "expr_full_logcpm.rds"))
meta <- as.data.table(readRDS(file.path(FM, "meta_full.rds")))
mi   <- match(colnames(E), meta$sample)
E    <- E[, !is.na(mi), drop=FALSE]; meta <- meta[mi[!is.na(mi)]]
meta[, pid := substr(sample, 1, 12)]
meta <- meta[!is.na(group) & group %in% c("Tumor","Normal")]
E    <- E[, meta$sample, drop=FALSE]
cat(sprintf("\nOriginal analysis set %d cases; after excluding %d cases, %d cases (tumor %d / normal %d)\n",
            nrow(meta), nrow(exc), nrow(meta) - nrow(exc),
            sum(meta$group=="Tumor") - nrow(exc), sum(meta$group=="Normal")))
keep <- !meta$sample %in% exc$sample
E <- E[, keep, drop=FALSE]; meta <- meta[keep]

y <- factor(meta$group, levels=c("Normal","Tumor"))
X <- t(E[GENES, , drop=FALSE])
df <- data.frame(y=y, X); names(df) <- c("y", make.names(GENES))
m_new <- glm(y ~ ., data=df, family=binomial)
cf_new <- coef(m_new)
cf_old <- c(`(Intercept)`=29.269077, ABHD17C=1.843595, ZDHHC6=-6.379393,
            ABHD17B=-4.568215, ZDHHC20=2.774859)
pr_new <- predict(m_new, type="response")
auc_new <- as.numeric(auc(roc(y, pr_new, levels=c("Normal","Tumor"), direction="<", quiet=TRUE)))
sec("Refit of the final model (after exclusion, n=557)")
cmp <- data.table(term=names(cf_old), original=cf_old,
                  excluded7=as.numeric(cf_new[names(cf_old)]),
                  rel_change_pct=100*(as.numeric(cf_new[names(cf_old)])-cf_old)/abs(cf_old))
print(cmp)
cat(sprintf("  Apparent AUC: original 0.979 (full data 564) -> after exclusion %.3f (n=%d)\n", auc_new, nrow(meta)))

make_folds <- function(seed) {
  set.seed(seed)
  pu <- unique(meta[, .(pid)]); pu[, has_normal := pid %in% meta[group=="Normal"]$pid]
  f  <- integer(nrow(pu))
  for (h in c(TRUE, FALSE)) {
    idx <- which(pu$has_normal == h)
    f[idx] <- rep(1:N_FOLD, length.out=length(idx))[order(runif(length(idx)))]
  }
  pu[, fold := f]
  return(pu$fold[match(meta$pid, pu$pid)])
}
svm_rfe <- function(X, y, n_stop) {
  feats <- colnames(X)
  while (length(feats) > n_stop) {
    m <- svm(X[, feats, drop=FALSE], y, kernel="linear", scale=TRUE, cost=1)
    w <- abs(as.numeric(t(m$coefs) %*% m$SV)); names(w) <- feats
    feats <- setdiff(feats, names(sort(w))[1])
  }
  feats
}
fit_in_fold <- function(tr_idx, te_idx) {
  Etr <- E[, tr_idx, drop=FALSE]; gtr <- factor(meta$group[tr_idx], levels=c("Normal","Tumor"))
  if (min(table(gtr)) < 5) return(NULL)
  design_tr <- model.matrix(~ gtr)
  cor_tr <- patient_correlation(Etr,design_tr,meta$pid[tr_idx])
  fit <- eBayes(lmFit(Etr, design_tr, block=meta$pid[tr_idx], correlation=cor_tr$consensus.correlation), trend=TRUE, robust=TRUE)
  d   <- topTable(fit, coef=2, number=Inf, sort.by="none")
  cand <- intersect(rownames(d)[d$adj.P.Val < 0.05], GSET)
  if (length(cand) < 6) return(NULL)
  X <- t(Etr[cand, , drop=FALSE]); y <- gtr
  s_lasso <- tryCatch({
    fid <- { p <- meta$pid[tr_idx]; u <- unique(p); set.seed(11)
             asg <- setNames(rep(1:5, length.out=length(u))[order(runif(length(u)))], u)
             as.integer(asg[p]) }
    cv <- cv.glmnet(X, y, family="binomial", alpha=1, foldid=fid)
    names(which(coef(cv, s="lambda.min")[-1,1] != 0))
  }, error=function(e) character(0))
  s_rf <- tryCatch({
    rf <- randomForest(X, y, ntree=500, importance=TRUE)
    names(sort(importance(rf, type=1)[,1], decreasing=TRUE))[1:min(10, ncol(X))]
  }, error=function(e) character(0))
  s_svm <- tryCatch(svm_rfe(X, y, n_stop=max(5, round(length(cand)*0.2))), error=function(e) character(0))
  sig <- Reduce(intersect, list(s_lasso, s_rf, s_svm))
  if (length(sig) < 2) sig <- intersect(s_lasso, s_rf)
  if (length(sig) < 2) return(NULL)
  df_tr <- data.frame(X[, sig, drop=FALSE]); names(df_tr) <- make.names(sig)
  m <- glm(y ~ ., data=data.frame(y=y, df_tr), family="binomial")
  df_te <- data.frame(t(E[sig, te_idx, drop=FALSE])); names(df_te) <- make.names(sig)
  pr <- tryCatch(predict(m, newdata=df_te, type="response"), error=function(e) rep(NA_real_, length(te_idx)))
  list(pred=pr, genes=sig)
}

sec("Nested CV rerun (after exclusion, 5 x 10 folds, seeds 1001-1005)")
oof <- data.table(sample = meta$sample, pid = meta$pid, group = meta$group)
gene_sel <- list()
for (r in 1:N_REP) {
  meta$fold <- make_folds(1000 + r)
  oof[[paste0("pred_r", r)]] <- NA_real_
  for (k in 1:N_FOLD) {
    tr <- which(meta$fold != k); te <- which(meta$fold == k)
    res <- fit_in_fold(tr, te)
    if (is.null(res)) { cat(sprintf("  repeat %d fold %d skipped\n", r, k)); next }
    data.table::set(oof, i=te, j=paste0("pred_r", r), value=res$pred)
    gene_sel[[length(gene_sel)+1]] <- res$genes
  }
}
stopifnot(all(vapply(oof[, paste0("pred_r", 1:N_REP), with=FALSE], function(x) all(is.finite(x)), logical(1))))
aucs <- sapply(1:N_REP, function(r) {
  p <- oof[[paste0("pred_r", r)]]; k <- !is.na(p)
  if (sum(k) < 20 || length(unique(oof$group[k])) < 2) return(NA_real_)
  as.numeric(auc(roc(oof$group[k], p[k], levels=c("Normal","Tumor"), direction="<", quiet=TRUE)))
})
cat(sprintf("  AUC per repeat: %s\n", paste(sprintf("%.3f", aucs), collapse="  ")))
cat(sprintf("  Mean AUC = %.3f (SD %.3f)  vs  original 0.965 (SD 0.007)\n",
            mean(aucs, na.rm=TRUE), sd(aucs, na.rm=TRUE)))

sec("In-fold selection of the original 4 genes")
cat(sprintf("  %s\n",
    paste(sprintf("%s %d/%d", GENES,
                  sapply(GENES, function(g) sum(sapply(gene_sel, function(s) g %in% s))),
                  length(gene_sel)), collapse="  ")))

out <- list(histology_excluded=exc[, .(sid, sample, histgrp, icd)],
            coef_compare=cmp,
            auc_excluded_model=auc_new,
            nested_cv_auc=aucs)
saveRDS(out, file.path(DIR, "histology_exclusion_sensitivity_2026-10-03.rds"))
fwrite(cmp, file.path(DIR, "histology_exclusion_coefs_2026-10-03.csv"))
fwrite(exc[, .(sid, sample, histgrp, icd)],
       file.path(DIR, "histology_excluded7_2026-10-03.csv"))
cat(sprintf("\n  output: %s\n", DIR))
