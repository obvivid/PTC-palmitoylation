# CTD associations and published null control
# Run through run_analysis.R; see docs/METHODS.md.

suppressMessages({ library(data.table) })

DIR <- RESULT_DIR
set.seed(20261001)
sec <- function(t){ cat("\n"); cat(strrep("=",78),"\n",sep=""); cat(t,"\n"); cat(strrep("=",78),"\n",sep="") }

sec("[0] Read CTD human chemical-gene pairs")
cg <- fread(file.path(PROJECT_ROOT, "data/ctd/ctd_human_chem_gene.tsv"), header=FALSE,
            col.names=c("chem","gene"), quote="", na.strings="")
cg <- unique(cg[!is.na(chem) & !is.na(gene) & chem!="" & gene!=""])
cat(sprintf("  chemical-gene pairs %d; chemicals %d, genes %d\n",
            nrow(cg), uniqueN(cg$chem), uniqueN(cg$gene)))

GS <- fread(file.path(DIR, "deg_full_palmitoylation_genes.csv"))[adj.P.Val<0.05]$gene
SIG <- fread(file.path(DIR, "signature_genes.csv"))$gene
cat(sprintf("  palmitoylation DEGs %d (of which %d are recorded in CTD)\n",
            length(GS), sum(GS %in% cg$gene)))
cat(sprintf("  records for the 4-gene signature in CTD: %s\n",
            paste(sprintf("%s(%d)", SIG, sapply(SIG, function(g) uniqueN(cg[gene==g]$chem))), collapse="  ")))

g2c <- split(cg$chem, cg$gene)
cnt <- function(gs) { x <- unlist(g2c[gs]); if (is.null(x)) return(integer(0)); sort(table(x), decreasing=TRUE) }

sec("[1] Main analysis: chemicals interacting with the most palmitoylation DEGs")
tt <- cnt(GS)
cat(sprintf("  chemicals hit %d; those interacting with >=3 genes %d\n",
            length(tt), sum(tt >= 3)))
cat("  top 20:\n")
for (i in seq_len(min(20, length(tt)))) {
  ch <- names(tt)[i]
  gs <- cg[chem==ch & gene %in% GS]$gene
  cat(sprintf("    %-32s %2d genes: %s\n", ch, tt[i], paste(gs, collapse=", ")))
}

sec("[2] Matched null control: 1000 random gene sets of the same size")
allg <- readLines(file.path(SOURCE_ROOT,"config/ctd_null_gene_universe.txt"),encoding="UTF-8")
if (!setequal(allg, unique(cg$gene))) stop("CTD gene universe differs from the manuscript snapshot; update the analysis definition explicitly.")
N <- 1000; n <- length(GS)
rand_max <- integer(N); rand_ge3 <- integer(N); rand_any <- integer(N)
for (i in seq_len(N)) {
  rs <- sample(allg, n)
  x  <- cnt(rs)
  rand_any[i]  <- length(x)
  rand_max[i]  <- if (length(x)) as.integer(x[1]) else 0L
  rand_ge3[i]  <- sum(x >= 3)
}
obs_max  <- if (length(tt)) as.integer(tt[1]) else 0L
obs_any  <- length(tt)
obs_ge3  <- sum(tt >= 3)
p_max <- (sum(rand_max >= obs_max) + 1) / (N + 1)
p_any <- (sum(rand_any >= obs_any) + 1) / (N + 1)
p_ge3 <- (sum(rand_ge3 >= obs_ge3) + 1) / (N + 1)

cat(sprintf("  chemicals hit        %6d | %.1f (%s~%s) | %.3f\n",
            obs_any, median(rand_any), quantile(rand_any,.025), quantile(rand_any,.975), p_any))
cat(sprintf("  max single-chemical  %6d | %.1f (%s~%s) | %.3f\n",
            obs_max, median(rand_max), quantile(rand_max,.025), quantile(rand_max,.975), p_max))
cat(sprintf("  chemicals >=3 genes  %6d | %.1f (%s~%s) | %.3f\n",
            obs_ge3, median(rand_ge3), quantile(rand_ge3,.025), quantile(rand_ge3,.975), p_ge3))
cat("\n  Interpretation: if the empirical P > 0.05, the palmitoylation DEGs yield no more candidates\n")
cat("        than a random gene set, i.e. the corresponding association count is not enriched against this CTD background.\n")

out <- data.table(chem=names(tt), n_genes=as.integer(tt))
out[, genes := sapply(chem, function(ch) paste(cg[chem==ch & gene %in% GS]$gene, collapse=";"))]
fwrite(out, file.path(DIR, "predicted_drugs_CTD.csv"))
fwrite(data.table(metric=c("n_chem_hits","max_gene_per_chem","n_chem_ge3",
                           "null_n_chem_median","null_max_median","null_ge3_median",
                           "p_n_chem","p_max","p_ge3"),
                  value=c(obs_any, obs_max, obs_ge3,
                          median(rand_any), median(rand_max), median(rand_ge3),
                          p_any, p_max, p_ge3)),
       file.path(DIR, "drug_prediction_control.csv"))
cat(sprintf("  output: %s\n", DIR))

background_summary<-data.table(statistic=c("n_chem_hits","n_chem_ge3","max_gene_per_chem"),
 observed=c(obs_any,obs_ge3,obs_max),null_median=c(median(rand_any),median(rand_ge3),median(rand_max)),
 empirical_P=c(p_any,p_ge3,p_max))
background_summary[,BH_adjusted_P:=p.adjust(empirical_P,"BH")]
fwrite(background_summary,file.path(DIR,"chemical_background_statistics.csv"))
