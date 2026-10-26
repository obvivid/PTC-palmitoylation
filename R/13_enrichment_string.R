# MSigDB enrichment and STRING network
# Run through run_analysis.R; see docs/METHODS.md.

suppressMessages({ library(data.table); library(httr); library(jsonlite) })

ROOT <- PROJECT_ROOT
FM   <- MATRIX_DIR
DIR  <- RESULT_DIR
sec  <- function(t){ cat("\n"); cat(strrep("=",78),"\n",sep=""); cat(t,"\n"); cat(strrep("=",78),"\n",sep="") }

E  <- readRDS(file.path(FM, "expr_full_logcpm.rds"))
GS <- fread(file.path(DIR, "deg_full_palmitoylation_genes.csv"))[adj.P.Val < 0.05]$gene
UNIVERSE <- rownames(E)
cat(sprintf("Gene set %d; background genes %d\n", length(GS), length(UNIVERSE)))

sec("[1] MSigDB hypergeometric enrichment")
gs_sets <- list()
try_msig <- function(expr) tryCatch(expr, error = function(e) NULL)

m1 <- try_msig(msigdbr::msigdbr(species="Homo sapiens", collection="H"))
if (!is.null(m1)) {
  gs_sets[["Hallmark"]] <- m1
} else {
  gs_sets[["Hallmark"]] <- try_msig(msigdbr::msigdbr(species="Homo sapiens", category="H"))
}
m2 <- try_msig(msigdbr::msigdbr(species="Homo sapiens", collection="C2", subcollection="CP:KEGG_MEDICUS"))
if (is.null(m2)) m2 <- try_msig(msigdbr::msigdbr(species="Homo sapiens", category="C2", subcategory="CP:KEGG"))
if (!is.null(m2)) gs_sets[["KEGG"]] <- m2
m3 <- try_msig(msigdbr::msigdbr(species="Homo sapiens", collection="C5", subcollection="GO:BP"))
if (is.null(m3)) m3 <- try_msig(msigdbr::msigdbr(species="Homo sapiens", category="C5", subcategory="GO:BP"))
if (!is.null(m3)) gs_sets[["GO_BP"]] <- m3
if (!all(c("Hallmark","KEGG","GO_BP") %in% names(gs_sets)) || any(vapply(gs_sets,is.null,logical(1))))
  stop("Required MSigDB collection unavailable; restore access before rerunning this step.")

for (nm in names(gs_sets)) cat(sprintf("  %-9s gene sets %d\n", nm, length(unique(gs_sets[[nm]]$gs_name))))

hyper <- function(df) {
  if (is.null(df)) return(NULL)
  gcol <- if ("gene_symbol" %in% names(df)) "gene_symbol" else "gene_symbol"
  d <- unique(as.data.table(df)[, .(set=gs_name, gene=get(gcol))])
  d <- d[gene %in% UNIVERSE]
  N <- length(UNIVERSE); n <- length(GS)
  res <- d[, {
    k <- length(intersect(gene, GS)); m <- .N
    .(hits = k, set_size = m,
      p = if (k == 0) NA_real_ else phyper(k-1, m, N-m, n, lower.tail=FALSE))
  }, by=set][!is.na(p)]
  res[, FDR := p.adjust(p, "BH")]
  res[, genes := sapply(set, function(s) paste(intersect(d[set==s]$gene, GS), collapse=","))]
  setorder(res, p)
  res
}
enr <- lapply(gs_sets, hyper)
for (nm in names(enr)) {
  r <- enr[[nm]]; if (is.null(r) || !nrow(r)) next
  cat(sprintf("\n  --- %s top 8 ---\n", nm))
  for (i in seq_len(min(8, nrow(r))))
    cat(sprintf("    %-52s %2d/%3d  P=%.3g  FDR=%.3g\n",
                substr(r$set[i],1,52), r$hits[i], r$set_size[i], r$p[i], r$FDR[i]))
}
for (nm in names(enr)) if (!is.null(enr[[nm]])) fwrite(enr[[nm]], file.path(DIR, paste0("enrich_", nm, ".csv")))

sec("[2] STRING PPI + degree hub")
# Frozen functional associations for the analysed 25-gene candidate set.
snapshot<-file.path(SOURCE_ROOT,"data/annotation_snapshots/string_network_25.tsv")
res<-fread(snapshot)
if(!setequal(unique(c(res$preferredName_A,res$preferredName_B)),GS))stop("STRING snapshot and analysed gene set differ; update the versioned query explicitly.")

if (!is.null(res) && nrow(res) > 0) {
  e <- unique(res[, .(a = preferredName_A, b = preferredName_B, s = as.numeric(score))])
  e <- e[a %in% GS & b %in% GS]
  deg <- rbind(e[, .(gene=a, w=s)], e[, .(gene=b, w=s)])[
         , .(degree=.N, mean_score=mean(w)), by=gene][order(-degree, -mean_score)]
  cat(sprintf("  %d interaction edges covering %d genes\n", nrow(e), nrow(deg)))
  cat("  top 10 by degree:\n")
  for (i in seq_len(min(10, nrow(deg))))
    cat(sprintf("    %-10s degree=%2d  mean confidence=%.3f\n", deg$gene[i], deg$degree[i], deg$mean_score[i]))
  fwrite(e,   file.path(DIR, "ppi_edges_string.csv"))
  fwrite(deg, file.path(DIR, "ppi_degree.csv"))
} else {
  cat("  WARNING: no STRING network retrieved (network may be restricted)\n")
}

cat(sprintf("  output directory: %s\n", DIR))
