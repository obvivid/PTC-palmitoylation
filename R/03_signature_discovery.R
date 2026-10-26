# Full-data signature discovery and descriptive fit
# Run through run_analysis.R; see docs/METHODS.md.

suppressMessages({ library(data.table); library(glmnet); library(randomForest); library(e1071); library(pROC) })

ROOT <- PROJECT_ROOT
FM   <- MATRIX_DIR
DIR  <- RESULT_DIR
set.seed(20261001)

sec <- function(t){ cat("\n"); cat(strrep("=",78),"\n",sep=""); cat(t,"\n"); cat(strrep("=",78),"\n",sep="") }

E    <- readRDS(file.path(FM, "expr_full_logcpm.rds"))
meta <- as.data.table(readRDS(file.path(FM, "meta_full.rds")))
ok   <- !is.na(meta$group)
E    <- E[, ok, drop=FALSE]; meta <- meta[ok]

degp <- fread(file.path(DIR, "deg_full_palmitoylation_genes.csv"))
CAND <- degp[adj.P.Val < 0.05]$gene
cat(sprintf("Candidate palmitoylation DEGs: %d\n", length(CAND)))
cat(sprintf("Samples: tumor %d / normal %d\n", sum(meta$group=="Tumor"), sum(meta$group=="Normal")))

X <- t(E[CAND, , drop=FALSE])                 # samples x genes
y <- factor(meta$group, levels=c("Normal","Tumor"))
cat(sprintf("Feature matrix: %d x %d\n", nrow(X), ncol(X)))

sec("[1] LASSO (glmnet, 10-fold CV, lambda.min)")
cvfit <- cv.glmnet(X, y, family="binomial", alpha=1, nfolds=10, type.measure="auc")
co <- coef(cvfit, s="lambda.min")
sel_lasso <- names(which(co[-1,1] != 0))
cat(sprintf("  selected %d: %s\n", length(sel_lasso), paste(sel_lasso, collapse=", ")))
cat(sprintf("  CV AUC = %.3f\n", max(cvfit$cvm)))

sec("[2] RandomForest (importance, top 10)")
rf <- randomForest(X, y, ntree=1000, importance=TRUE)
imp <- importance(rf, type=1)[,1]
sel_rf <- names(sort(imp, decreasing=TRUE))[1:min(10, length(imp))]
cat(sprintf("  top 10: %s\n", paste(sel_rf, collapse=", ")))
cat(sprintf("  OOB error rate = %.3f\n", rf$err.rate[1000,"OOB"]))

sec("[3] SVM-RFE (linear kernel, drop lowest weight each round)")
svm_rfe <- function(X, y, n_stop=10) {
  feats <- colnames(X)
  while (length(feats) > n_stop) {
    m <- svm(X[, feats, drop=FALSE], y, kernel="linear", scale=TRUE, cost=1)
    w <- t(m$coefs) %*% m$SV          # weight vector
    w <- abs(as.numeric(w)); names(w) <- feats
    drop <- names(sort(w))[1]          # drop the smallest weight
    feats <- setdiff(feats, drop)
  }
  feats
}
sel_svm <- svm_rfe(X, y, n_stop = max(5, round(length(CAND)*0.2)))
sel_svm <- sel_svm[order(match(sel_svm, colnames(X)))]
cat(sprintf("  selected %d: %s\n", length(sel_svm), paste(sel_svm, collapse=", ")))

sec("[4] Intersection of the three algorithms -> signature")
SIG <- Reduce(intersect, list(sel_lasso, sel_rf, sel_svm))
if (length(SIG) < 2) {
  cat("  WARNING: intersection < 2, falling back to LASSO ∩ RF\n")
  SIG <- intersect(sel_lasso, sel_rf)
}
cat(sprintf("  Signature genes (%d): %s\n", length(SIG), paste(SIG, collapse=", ")))

Xs <- X[, SIG, drop=FALSE]
dfall <- as.data.frame(Xs)
mfull <- glm(y ~ ., data=data.frame(y=y, dfall), family="binomial")
pfull <- predict(mfull, type="response")
rfull <- roc(y, pfull, levels=c("Normal","Tumor"), direction="<", quiet=TRUE)
ci <- ci.auc(rfull)
cat(sprintf("  Full-data AUC = %.3f (95%% CI %.3f-%.3f)\n", as.numeric(rfull$auc), ci[1], ci[3]))
cat("  Per-gene coefficients:\n")
s <- summary(mfull)$coef
for (i in 2:nrow(s)) cat(sprintf("    %-10s beta %+7.3f  P %.3g\n", rownames(s)[i], s[i,1], s[i,4]))

fwrite(data.table(gene=sel_lasso), file.path(DIR,"sel_lasso.csv"))
fwrite(data.table(gene=sel_rf),    file.path(DIR,"sel_rf.csv"))
fwrite(data.table(gene=sel_svm),   file.path(DIR,"sel_svm.csv"))
fwrite(data.table(gene=SIG),       file.path(DIR,"signature_genes.csv"))
fwrite(data.table(sample=rownames(Xs), prob=pfull, group=as.character(y), risk=as.numeric(predict(mfull,type="link"))),
       file.path(DIR,"signature_scores_tcga.csv"))
environment(mfull$formula) <- baseenv()
attr(mfull$terms, ".Environment") <- baseenv()
saveRDS(list(model=mfull, genes=SIG), file.path(DIR,"signature_model.rds"))
patient <- substr(rownames(dfall),1,12)
V <- sandwich::vcovCL(mfull,cluster=patient,type="HC1",cadjust=TRUE)
se <- sqrt(diag(V));df_cluster <- length(unique(patient))-1L
fwrite(data.table(term=names(coef(mfull)),coefficient=as.numeric(coef(mfull)),
 se=as.numeric(se),ci_lo=as.numeric(coef(mfull)-qt(.975,df_cluster)*se),
 ci_hi=as.numeric(coef(mfull)+qt(.975,df_cluster)*se),
 p=as.numeric(2*pt(abs(coef(mfull)/se),df=df_cluster,lower.tail=FALSE)),n_patients=df_cluster+1L),
 file.path(DIR,"signature_coefficients.csv"))
cat(sprintf("  output: %s\n", DIR))
