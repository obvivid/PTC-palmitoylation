# Patient-grouped repeated nested cross-validation
# Run through run_analysis.R; see docs/METHODS.md.

suppressMessages({ library(data.table); library(limma); library(glmnet)
                   library(randomForest); library(e1071); library(pROC) })

ROOT <- PROJECT_ROOT
FM   <- MATRIX_DIR
DIR  <- RESULT_DIR
sec  <- function(t){ cat("\n"); cat(strrep("=",78),"\n",sep=""); cat(t,"\n"); cat(strrep("=",78),"\n",sep="") }

N_REP <- 5; N_FOLD <- 10
sid12 <- function(x) substr(as.character(x), 1, 12)

GSET <- fread(file.path(DIR, "palmitoylation_geneset_full.csv"))$gene
E    <- readRDS(file.path(FM, "expr_full_logcpm.rds"))
meta <- as.data.table(readRDS(file.path(FM, "meta_full.rds")))
mi   <- match(colnames(E), meta$sample)
E    <- E[, !is.na(mi), drop=FALSE]; meta <- meta[mi[!is.na(mi)]]
meta[, pid := sid12(sample)]
meta <- meta[!is.na(group) & group %in% c("Tumor","Normal")]
E    <- E[, meta$sample, drop=FALSE]

cat(sprintf("Samples %d (tumor %d / normal %d); patients %d; palmitoylation gene set %d\n",
            ncol(E), sum(meta$group=="Tumor"), sum(meta$group=="Normal"),
            uniqueN(meta$pid), length(GSET)))

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
meta[, fold := make_folds(1)]
stopifnot(all(meta[, .(nfold = uniqueN(fold)), by=pid]$nfold == 1))
cat(sprintf("Fold self-check passed: all samples of a patient share one fold; normals per fold = %s\n",
            paste(table(meta$group, meta$fold)["Normal", ], collapse="/")))

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
             as.integer(asg[p]) }          # Patient-grouped inner folds.
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
  df_te <- data.frame(t(E[sig, te_idx, drop=FALSE]))
  names(df_te) <- make.names(sig)
  pr <- tryCatch(predict(m, newdata=df_te, type="response"), error=function(e) rep(NA_real_, length(te_idx)))
  list(pred=pr, genes=sig, cand_n=length(cand))
}

sec("Repeated nested cross-validation (patient-grouped, 5 x 10 folds)")
oof <- data.table(sample = meta$sample, pid = meta$pid, group = meta$group)
gene_sel <- list(); n_cand <- c(); fold_registry <- list(); fold_selections <- list()
for (r in 1:N_REP) {
  meta$fold <- make_folds(1000 + r)
  oof[[paste0("pred_r", r)]] <- NA_real_
  for (k in 1:N_FOLD) {
    tr <- which(meta$fold != k); te <- which(meta$fold == k)
    res <- fit_in_fold(tr, te)
    if (is.null(res)) { cat(sprintf("  repeat %d fold %d skipped (insufficient samples / selection failed)\n", r, k)); next }
    data.table::set(oof, i=te, j=paste0("pred_r", r), value=res$pred)
    fold_registry[[length(fold_registry)+1L]] <- data.table(repeat_id=r,fold=k,sample=meta$sample[te],patient=meta$pid[te],group=meta$group[te])
    fold_selections[[length(fold_selections)+1L]] <- data.table(repeat_id=r,fold=k,training_tissues=length(tr),heldout_tissues=length(te),training_patients=uniqueN(meta$pid[tr]),heldout_patients=uniqueN(meta$pid[te]),candidate_genes=res$cand_n,selected_genes=paste(res$genes,collapse=";"))
    gene_sel[[length(gene_sel)+1]] <- res$genes
    n_cand <- c(n_cand, res$cand_n)
  }
}

stopifnot(all(vapply(oof[, paste0("pred_r", 1:N_REP), with=FALSE], function(x) all(is.finite(x)), logical(1))))
sec("Results")
aucs <- sapply(1:N_REP, function(r) {
  p <- oof[[paste0("pred_r", r)]]; k <- !is.na(p)
  if (sum(k) < 20 || length(unique(oof$group[k])) < 2) return(NA_real_)
  as.numeric(auc(roc(oof$group[k], p[k], levels=c("Normal","Tumor"), direction="<", quiet=TRUE)))
})
cat(sprintf("  Outer hold-out AUC per repeat: %s\n", paste(sprintf("%.3f", aucs), collapse="  ")))
cat(sprintf("  Mean AUC = %.3f (SD %.3f, range %.3f-%.3f)\n",
            mean(aucs, na.rm=TRUE), sd(aucs, na.rm=TRUE), min(aucs,na.rm=TRUE), max(aucs,na.rm=TRUE)))

set.seed(7)
p1 <- oof$pred_r1; ok <- !is.na(p1)
d1 <- oof[ok]
bs <- replicate(2000, {
  pid_s <- sample(unique(d1$pid), replace=TRUE)
  idx <- unlist(lapply(pid_s, function(p) which(d1$pid == p)))
  dd <- d1[idx]
  if (length(unique(dd$group)) < 2) return(NA_real_)
  as.numeric(auc(roc(dd$group, dd$pred_r1, levels=c("Normal","Tumor"), direction="<", quiet=TRUE)))
})
bs <- bs[!is.na(bs)]
cat(sprintf("  Patient-level bootstrap 95%% CI (repeat 1): %.3f-%.3f (n=%d)\n",
            quantile(bs,.025), quantile(bs,.975), length(bs)))

sec("Feature selection frequency (across the 50 folds)")
allg <- table(unlist(gene_sel))
allg <- sort(allg, decreasing=TRUE)
cat(sprintf("  Candidate genes per fold: median %d (range %d-%d)\n", median(n_cand), min(n_cand), max(n_cand)))
cat(sprintf("  %d genes were selected into a signature; those present in >= 25%% of folds:\n", length(allg)))
for (i in seq_along(allg)) if (allg[i] >= 0.25*length(gene_sel))
  cat(sprintf("    %-10s %2d/%d folds (%.0f%%)\n", names(allg)[i], allg[i], length(gene_sel), 100*allg[i]/length(gene_sel)))
cat(sprintf("  Published four-gene selection frequency: %s\n",
            paste(sprintf("%s %d/%d", c("ABHD17C","ZDHHC6","ABHD17B","ZDHHC20"),
                          sapply(c("ABHD17C","ZDHHC6","ABHD17B","ZDHHC20"), function(g) sum(sapply(gene_sel, function(s) g %in% s))),
                          length(gene_sel)), collapse="  ")))

fwrite(oof, file.path(DIR, "nested_cv_oof_predictions.csv"), na="NA")
fwrite(rbindlist(fold_registry),file.path(DIR,"nested_cv_fold_registry.csv"))
fwrite(rbindlist(fold_selections),file.path(DIR,"nested_cv_fold_selections.csv"))
fwrite(data.table(repeat_id=1:N_REP,auc=aucs),file.path(DIR,"nested_cv_auc_repeats.csv"))
fwrite(data.table(gene=names(allg), n_folds=as.integer(allg)),
       file.path(DIR, "nested_cv_gene_selection_freq.csv"))
cat(sprintf("\n  output: %s\n", DIR))
