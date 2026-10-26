# Paired tumor versus normal sensitivity
# Run through run_analysis.R; see docs/METHODS.md.

suppressMessages({ library(data.table); library(pROC) })

ROOT <- PROJECT_ROOT
FM   <- MATRIX_DIR
DIR  <- RESULT_DIR
sec  <- function(t){ cat("\n"); cat(strrep("=",78),"\n",sep=""); cat(t,"\n"); cat(strrep("=",78),"\n",sep="") }

GENES <- c("ABHD17C","ZDHHC6","ABHD17B","ZDHHC20")
COEF  <- c(ABHD17C=1.843595, ZDHHC6=-6.379393, ABHD17B=-4.568215, ZDHHC20=2.774859)

E    <- readRDS(file.path(FM, "expr_full_logcpm.rds"))
meta <- as.data.table(readRDS(file.path(FM, "meta_full.rds")))
mi   <- match(colnames(E), meta$sample)
E    <- E[, !is.na(mi), drop=FALSE]; meta <- meta[mi[!is.na(mi)]]
meta[, pid := substr(sample, 1, 12)]
meta <- meta[!is.na(group) & group %in% c("Tumor","Normal")]
E    <- E[, meta$sample, drop=FALSE]
stopifnot(all(GENES %in% rownames(E)))

setorder(meta, sample)
tu <- meta[group=="Tumor",  .SD[1], by=pid]
no <- meta[group=="Normal", .SD[1], by=pid]
pairs <- merge(tu[, .(pid, tumor=sample)], no[, .(pid, normal=sample)], by="pid")
sec("Paired structure")
cat(sprintf("Tumor %d / normal %d / patients %d; patients with both sides (pairs) = %d\n",
            nrow(tu), nrow(no), uniqueN(meta$pid), nrow(pairs)))

Xtu <- E[GENES, pairs$tumor, drop=FALSE]
Xno <- E[GENES, pairs$normal, drop=FALSE]
score_tu <- as.numeric(COEF %*% Xtu)   # intercept is irrelevant to the paired difference, omitted
score_no <- as.numeric(COEF %*% Xno)

run_one <- function(a, b, label) {
  d  <- a - b
  wt <- suppressWarnings(stable_wilcox_test(a, b, paired=TRUE, exact=FALSE))
  data.table(item=label, n_pairs=length(d),
             median_tumor=median(a), median_normal=median(b),
             median_paired_diff=median(d),
             n_direction_tumor_higher=sum(d > 0),
             p_wilcoxon=wt$p.value)
}
res <- rbindlist(lapply(seq_along(GENES), function(i)
  run_one(Xtu[i,], Xno[i,], GENES[i])))
res <- rbind(res, run_one(score_tu, score_no, "four_gene_score"))
res[, p_BH := p.adjust(p_wilcoxon, method="BH")]

sec("Paired Wilcoxon result (log2(CPM+1) scale)")
print(res)

roc_p <- roc(c(rep("Tumor",nrow(pairs)), rep("Normal",nrow(pairs))),
             c(score_tu, score_no), levels=c("Normal","Tumor"), direction="<", quiet=TRUE)
cat(sprintf("\n  Reference: paired subset (%d tumor vs %d normal) score AUC = %.3f (unpaired definition)\n",
            nrow(pairs), nrow(pairs), as.numeric(auc(roc_p))))

fwrite(res, file.path(DIR, "paired_sensitivity_59pairs_2026-10-03.csv"))
cat(sprintf("\n  output: %s\n", file.path(DIR, "paired_sensitivity_59pairs_2026-10-03.csv")))
