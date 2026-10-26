# Differential expression and gene-family enrichment
# Run through run_analysis.R; see docs/METHODS.md.

suppressMessages({ library(data.table) })

ROOT <- PROJECT_ROOT
FM   <- MATRIX_DIR
DIR  <- RESULT_DIR

sec <- function(t){ cat("\n"); cat(strrep("=",78),"\n",sep=""); cat(t,"\n"); cat(strrep("=",78),"\n",sep="") }

E    <- readRDS(file.path(FM, "expr_full_logcpm.rds"))
meta <- as.data.table(readRDS(file.path(FM, "meta_full.rds")))

sec("[0] Data")
cat(sprintf("  matrix %d genes x %d samples\n", nrow(E), ncol(E)))
cat(sprintf("  group distribution: %s (NA %d)\n",
            paste(names(table(meta$group, useNA="ifany")), table(meta$group, useNA="ifany"),
                  sep="=", collapse="  "), sum(is.na(meta$group))))

keep <- !is.na(meta$group)
E    <- E[, keep, drop = FALSE]
meta <- meta[keep]
cat(sprintf("  samples with group retained: %d (tumor %d / normal %d)\n",
            ncol(E), sum(meta$group=="Tumor"), sum(meta$group=="Normal")))

sec("[1] Palmitoylation gene set (extracted from the new matrix by naming rule)")
gn <- rownames(E)
gs_parts <- list(
  ZDHHC    = grep("^ZDHHC[0-9]+$", gn, value=TRUE),
  LYPLA    = grep("^LYPLA[0-9]+$", gn, value=TRUE),
  ABHD17   = grep("^ABHD17[ABC]$",  gn, value=TRUE),
  PPT      = grep("^PPT[12]$",      gn, value=TRUE),
  ABHD10   = grep("^ABHD10$",       gn, value=TRUE)
)
for (k in names(gs_parts)) cat(sprintf("  %-8s (%2d): %s\n", k, length(gs_parts[[k]]),
                                       paste(gs_parts[[k]], collapse=", ")))
GS <- unique(unlist(gs_parts, use.names=FALSE))
cat(sprintf("  --> %d genes in total\n", length(GS)))

sec("[2] limma tumor vs normal (genome-wide)")
suppressMessages(library(limma))
grp <- factor(meta$group, levels=c("Normal","Tumor"))
design <- model.matrix(~ grp)
pid <- substr(meta$sample,1,12)
corfit <- patient_correlation(E,design,pid)
fit <- eBayes(lmFit(E, design, block=pid, correlation=corfit$consensus.correlation), trend=TRUE, robust=TRUE)
fwrite(data.table(consensus_correlation=corfit$consensus.correlation),file.path(DIR,"differential_patient_correlation.csv"))
tt  <- as.data.table(topTable(fit, coef=2, number=Inf, sort.by="none"), keep.rownames="gene")
tt[,se:=fit$stdev.unscaled[,2]*sqrt(fit$s2.post)]
tt[,`:=`(ci_lo=logFC-qt(.975,fit$df.total)*se,ci_hi=logFC+qt(.975,fit$df.total)*se)]
setorder(tt, adj.P.Val)

n_all <- nrow(tt); n_de <- tt[adj.P.Val<0.05, .N]
cat(sprintf("  total genes %d; FDR<0.05 %d (%.1f%%)\n", n_all, n_de, 100*n_de/n_all))
BASE <- mean(tt$adj.P.Val < 0.05)

sec("[3] Per-gene results for the palmitoylation genes")
sub <- tt[gene %in% GS]; setorder(sub, adj.P.Val)
cat(sprintf("  %-10s %8s %8s %10s %10s  %s\n","gene","logFC","AveExpr","P","FDR","decision"))
for (i in seq_len(nrow(sub))) {
  r <- sub[i]
  cat(sprintf("  %-10s %+8.3f %8.2f %10.3g %10.3g  %s\n",
              r$gene, r$logFC, r$AveExpr, r$P.Value, r$adj.P.Val,
              ifelse(r$adj.P.Val<0.05,"DE","ns")))
}

sec("[4] Enrichment test and direction")
nh <- nrow(sub); nhd <- sub[adj.P.Val<0.05, .N]
cat(sprintf("  DE count: %d / %d = %.1f%% (genome-wide %.1f%%)\n", nhd, nh, 100*nhd/nh, 100*BASE))
cat(sprintf("  hypergeometric P = %.4g (observed %d vs expected %.1f)\n",
            phyper(nhd-1, n_de, n_all-n_de, nh, lower.tail=FALSE), nhd, nh*BASE))
up <- sub[adj.P.Val<0.05 & logFC>0,.N]; dn <- sub[adj.P.Val<0.05 & logFC<0,.N]
cat(sprintf("  direction: up %d / down %d (concordance %.0f%%)\n", up, dn, 100*max(up,dn)/max(1,up+dn)))
big <- sub[abs(logFC)>0.5, .N]
mid <- sub[abs(logFC)>0.3, .N]

sec("[5] ZDHHC family expression direction")
z <- sub[grepl("^ZDHHC", gene)]
cat(sprintf("  ZDHHC DE: %d; up %d / down %d\n",
            z[adj.P.Val<0.05,.N], z[adj.P.Val<0.05 & logFC>0,.N], z[adj.P.Val<0.05 & logFC<0,.N]))
recovered <- intersect(c("ZDHHC5","ZDHHC2","ZDHHC20","ZDHHC8","ZDHHC9","ZDHHC16","ZDHHC1","ZDHHC12","ZDHHC21","ZDHHC15","ZDHHC23"),
                       sub$gene)
for (g in recovered) { r <- sub[gene==g]; cat(sprintf("    %-9s logFC %+7.3f  FDR %.3g  %s\n",
                                                      g, r$logFC, r$adj.P.Val, ifelse(r$adj.P.Val<0.05,"DE","ns"))) }

fwrite(tt,  file.path(DIR, "deg_full_tumor_vs_normal.csv"))
fwrite(sub, file.path(DIR, "deg_full_palmitoylation_genes.csv"))
fwrite(data.table(gene=GS), file.path(DIR, "palmitoylation_geneset_full.csv"))
cat(sprintf("  output directory: %s\n", DIR))
