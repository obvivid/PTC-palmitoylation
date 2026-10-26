# GO and KEGG enrichment
# Run through run_analysis.R; see docs/METHODS.md.

suppressMessages({ library(data.table); library(clusterProfiler); library(org.Hs.eg.db) })

ROOT <- PROJECT_ROOT
FM   <- MATRIX_DIR
DIR  <- RESULT_DIR
sec  <- function(t){ cat("\n"); cat(strrep("=",78),"\n",sep=""); cat(t,"\n"); cat(strrep("=",78),"\n",sep="") }

E  <- readRDS(file.path(FM, "expr_full_logcpm.rds"))
GS <- fread(file.path(DIR, "deg_full_palmitoylation_genes.csv"))[adj.P.Val < 0.05]$gene
UNIVERSE <- rownames(E)
cat(sprintf("Gene set %d; background %d\n", length(GS), length(UNIVERSE)))

sec("[1] enrichGO (BP)")
go <- tryCatch(
  enrichGO(gene = GS, universe = UNIVERSE, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
           ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
           minGSSize = 5, maxGSSize = 1000, readable = TRUE),
  error = function(e) { cat("  failed:", conditionMessage(e), "\n"); stop("Required external resource unavailable; retry this step after restoring access.") })
if (!is.null(go)) {
  r <- as.data.table(go@result)
  r[, geneID_len := lengths(strsplit(geneID, "/"))]
  setorder(r, pvalue)
  cat(sprintf("  Significant terms (FDR<0.05): %d\n", r[p.adjust < 0.05, .N]))
  cat("  top 10:\n")
  for (i in seq_len(min(10, nrow(r))))
    cat(sprintf("    %-44s %2d/%3d  P=%.3g  FDR=%.3g\n",
                substr(r$Description[i],1,44), r$geneID_len[i], r$Count[i]*0 + r$Count[i] + (r$geneID_len[i]-r$Count[i]),
                r$pvalue[i], r$p.adjust[i]))
  cat(sprintf("  Terms containing \'palmitoyl\': %d, best FDR=%.3g\n",
              r[grepl("palmitoyl", Description, ignore.case=TRUE), .N],
              r[grepl("palmitoyl", Description, ignore.case=TRUE)][1]$p.adjust))
  fwrite(r, file.path(DIR, "enrich_GO_BP_clusterProfiler.csv"))
}

sec("[2] enrichKEGG")
# Use the supplied KEGG annotation snapshot to retain the analysed mapping.
kgmap<-fread(file.path(SOURCE_ROOT,"data/annotation_snapshots/kegg_pathway_genes.tsv"),header=FALSE)
knames<-fread(file.path(SOURCE_ROOT,"data/annotation_snapshots/kegg_pathway_names.tsv"),header=FALSE)
t2g<-data.frame(term=sub("^path:","",kgmap[[2]]),gene=sub("^hsa:","",kgmap[[1]]))
t2n<-data.frame(term=knames[[1]],name=sub(" - Homo sapiens \\(human\\)$","",knames[[2]]))
kg<-enricher(gene=bitr(GS,fromType="SYMBOL",toType="ENTREZID",OrgDb=org.Hs.eg.db)$ENTREZID,
 universe=bitr(UNIVERSE,fromType="SYMBOL",toType="ENTREZID",OrgDb=org.Hs.eg.db)$ENTREZID,
 TERM2GENE=t2g,TERM2NAME=t2n,pAdjustMethod="BH",pvalueCutoff=1,qvalueCutoff=1,minGSSize=5,maxGSSize=1000)

if (!is.null(kg)) {
  rk <- as.data.table(kg@result); setorder(rk, pvalue)
  cat(sprintf("  Significant terms (FDR<0.05): %d\n", rk[p.adjust < 0.05, .N]))
  for (i in seq_len(min(8, nrow(rk))))
    cat(sprintf("    %-44s %2d  P=%.3g  FDR=%.3g\n", substr(rk$Description[i],1,44),
                rk$Count[i], rk$pvalue[i], rk$p.adjust[i]))
  fwrite(rk, file.path(DIR, "enrich_KEGG_clusterProfiler.csv"))
}
