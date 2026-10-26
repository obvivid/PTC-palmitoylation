# TCGA expression preprocessing
# Run through run_analysis.R; see docs/METHODS.md.

suppressMessages({ library(data.table) })

ROOT <- PROJECT_ROOT
RAW  <- file.path(ROOT, "data", "tcga", "TCGA-THCA.star_counts.tsv.gz")
HGNC <- file.path(ROOT, "data", "reference", "hgnc_complete_set.txt")
DIR  <- MATRIX_DIR; dir.create(DIR, showWarnings = FALSE)

sec <- function(t){ cat("\n"); cat(strrep("=",78),"\n",sep=""); cat(t,"\n"); cat(strrep("=",78),"\n",sep="") }

sec("[1] Read raw counts")
raw <- fread(RAW, header = TRUE)
setnames(raw, 1, "ensg_raw")
samples <- setdiff(names(raw), "ensg_raw")
cat(sprintf("  rows %d; sample columns %d\n", nrow(raw), length(samples)))

v <- as.numeric(raw[[2]])
cat(sprintf("  column 2 range %.3f ~ %.3f; all integers after 2^x-1: %s\n",
            min(v), max(v, na.rm=TRUE),
            ifelse(max(abs((2^v - 1) - round(2^v - 1)), na.rm=TRUE) < 1e-6, "yes", "no")))

ensg <- sub("\\..*$", "", raw$ensg_raw)
mat  <- as.matrix(raw[, ..samples])
rm(raw); invisible(gc())
cat(sprintf("  matrix %d x %d\n", nrow(mat), ncol(mat)))

sec("[2] Recover counts -> CPM -> log2")
cnt <- 2^mat - 1
cnt[cnt < 0] <- 0
rm(mat); invisible(gc())
lib <- colSums(cnt, na.rm = TRUE)
cat(sprintf("  median library size %.3g (expected ~tens of millions)\n", median(lib)))
cpm    <- sweep(cnt, 2, lib, "/") * 1e6
logcpm <- log2(cpm + 1)
rm(cnt, cpm); invisible(gc())

sec("[3] Filter low-expression genes")
keep <- rowSums(logcpm > 1, na.rm = TRUE) >= ceiling(0.10 * ncol(logcpm))
cat(sprintf("  %d -> %d genes (%.1f%% retained)\n", nrow(logcpm), sum(keep), 100*mean(keep)))
logcpm <- logcpm[keep, , drop = FALSE]
ensg   <- ensg[keep]

sec("[4] ENSG -> symbol (via HGNC, not from memory)")
hgnc <- fread(HGNC, select = c("symbol", "ensembl_gene_id"), quote = "", na.strings = c("", "NA"))
hgnc <- unique(hgnc[!is.na(ensembl_gene_id)])
hgnc <- hgnc[!duplicated(ensembl_gene_id)]
m    <- hgnc[match(ensg, hgnc$ensembl_gene_id)]
sym  <- m$symbol
cat(sprintf("  annotated %d / %d (%.1f%%)\n", sum(!is.na(sym)), length(sym), 100*mean(!is.na(sym))))
rownames(logcpm) <- ifelse(is.na(sym), ensg, sym)

mu <- rowMeans(logcpm)
logcpm <- logcpm[order(mu, decreasing = TRUE), , drop = FALSE]
dup <- duplicated(rownames(logcpm))
cat(sprintf("  %d duplicate symbols removed\n", sum(dup)))
logcpm <- logcpm[!dup, , drop = FALSE]

meta <- make_tcga_metadata(colnames(logcpm))
zd_new <- grep("^ZDHHC[0-9]+$", rownames(logcpm), value=TRUE)
saveRDS(logcpm, file.path(DIR, "expr_full_logcpm.rds"))
saveRDS(meta,   file.path(DIR, "meta_full.rds"))
cat(sprintf("  written: %s\n", file.path(DIR, "expr_full_logcpm.rds")))
cat(sprintf("  ZDHHC expression in the new matrix (median log2CPM):\n"))
z <- sort(rowMeans(logcpm[zd_new, , drop=FALSE]), decreasing = TRUE)
for (i in seq_along(z)) cat(sprintf("    %-9s %.2f\n", names(z)[i], z[i]))
