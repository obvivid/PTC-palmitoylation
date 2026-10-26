# FTC versus FTA transfer assessment
# Run through run_analysis.R; see docs/METHODS.md.

suppressMessages({ library(data.table); library(pROC) })
ROOT <- PROJECT_ROOT
GEO  <- file.path(ROOT, "data/geo")
DIR  <- RESULT_DIR
sec  <- function(t){ cat("\n"); cat(strrep("=",78),"\n",sep=""); cat(t,"\n"); cat(strrep("=",78),"\n",sep="") }

SIG   <- fread(file.path(DIR, "signature_genes.csv"))$gene
model <- readRDS(file.path(DIR, "signature_model.rds"))
COEF  <- coef(model$model)[-1]; names(COEF) <- SIG

ann <- fread(file.path(ROOT,"data/geo/GPL570.annot.gz"), skip=27, header=TRUE, sep="\t",
             quote="", na.strings=c("","NA"), fill=TRUE); setnames(ann,1,"probe")
A <- ann[!is.na(`Gene symbol`) & `Gene symbol` %in% SIG,
         .(probe=as.character(probe), symbol=as.character(`Gene symbol`))]
cat(sprintf("Probe counts for the four genes on GPL570:\n")); print(A[, .N, by=symbol])

fp <- file.path(GEO, "GSE82208_series_matrix.txt.gz")
read_gpl570 <- function(fp) {
  hl <- grep("!series_matrix_table_begin", readLines(gzfile(fp), n=300))[1]
  M  <- fread(fp, skip=hl, header=TRUE, na.strings=c("","NA")); setnames(M,1,"probe")
  M[, probe := gsub('"', '', as.character(probe))]; M
}
M  <- read_gpl570(fp); sc <- setdiff(names(M), "probe")

tl   <- readLines(gzfile(fp), n=300)
gsm_l  <- grep("^!Sample_geo_accession", tl, value=TRUE)
tit_l  <- grep("^!Sample_title", tl, value=TRUE)
gsm  <- strsplit(gsub('"','',sub("^!Sample_geo_accession\\t?","",gsm_l)), "\t")[[1]]
tit  <- strsplit(gsub('"','',sub("^!Sample_title\\t?","",tit_l)), "\t")[[1]]
lab  <- ifelse(grepl("adenoma", tit, ignore.case=TRUE), "FTA",
        ifelse(grepl("cancer",  tit, ignore.case=TRUE), "FTC", NA))
names(lab) <- gsm; lab <- lab[sc]
sec("[0] Sample composition")
print(table(lab, useNA="ifany"))

sub <- merge(A, M, by="probe")
G   <- sub[, lapply(.SD, function(x) mean(as.numeric(x), na.rm=TRUE)), by=symbol, .SDcols=sc]
Gm  <- as.matrix(G[, ..sc]); rownames(Gm) <- G$symbol
cat(sprintf("Mapped to genes: %s (%d signature genes)\n", paste(rownames(Gm),collapse=","), nrow(Gm)))

keep <- !is.na(lab) & lab %in% c("FTA","FTC")
Gt   <- Gm[, keep, drop=FALSE]; labt <- lab[keep]
sec("[1] Score after z-standardisation within the target samples")
Z  <- t(scale(t(Gt)))
lp <- as.numeric(t(Z) %*% COEF[rownames(Gt)])
ro <- roc(labt, lp, levels=c("FTA","FTC"), direction="<", quiet=TRUE)
ci <- ci.auc(ro, method="delong")
cat(sprintf("Four-gene weighted score  FTA(n=%d) vs FTC(n=%d):  AUC=%.4f (DeLong %.3f-%.3f)\n",
            sum(labt=="FTA"), sum(labt=="FTC"), as.numeric(ro$auc), ci[1], ci[3]))

sec("[3] Single-gene difference between FTC and FTA (Wilcoxon, descriptive only)")
deg <- rbindlist(lapply(rownames(Gt), function(g){
  w <- stable_wilcox_test(Gt[g, labt=="FTC"], Gt[g, labt=="FTA"])
  data.table(gene=g, median_FTC=median(Gt[g,labt=="FTC"]), median_FTA=median(Gt[g,labt=="FTA"]), p=w$p.value)
}))
deg[, p_BH := p.adjust(p, method="BH")]; print(deg)

sec("[4] Score distribution (for describing where FTA falls in the text)")
cat(sprintf("  median FTA score %.2f (IQR %.2f-%.2f)\n", median(lp[labt=="FTA"]), quantile(lp[labt=="FTA"],.25), quantile(lp[labt=="FTA"],.75)))
cat(sprintf("  median FTC score %.2f (IQR %.2f-%.2f)\n", median(lp[labt=="FTC"]), quantile(lp[labt=="FTC"],.25), quantile(lp[labt=="FTC"],.75)))

out <- data.table(dataset="GSE82208", contrast="FTC_vs_FTA",
                  n=length(labt), n_mal=sum(labt=="FTC"), n_benign=sum(labt=="FTA"),
                  auc=as.numeric(ro$auc), delong_lo=ci[1], delong_hi=ci[3])
fwrite(out, file.path(DIR, "benign_gse82208_summary.csv"))
fwrite(deg, file.path(DIR, "benign_gse82208_per_gene.csv"))
cat(sprintf("\nOutput written to %s\n", DIR))
