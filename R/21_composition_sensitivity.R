# Descriptive composition proxies, OLS coefficients and mutation-stratified contrasts.
DATA_DIR <- file.path(PROJECT_ROOT, "data")
suppressPackageStartupMessages({
  library(data.table); library(estimate); library(Matrix)
})

WS  <- OUTPUT_DIR
OUT <- file.path(RESULT_DIR,"extended","composition")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

SIG_GENES <- c("ABHD17C", "ZDHHC6", "ABHD17B", "ZDHHC20")

E <- as.matrix(readRDS(file.path(MATRIX_DIR,"expr_full_logcpm.rds")))
m <- as.data.table(readRDS(file.path(MATRIX_DIR,"meta_full.rds")))
ci <- fread(file.path(RESULT_DIR,"clinical_immune_fixed.csv"))   # 505 tumours; includes the study ssGSEA annotations.
cat(sprintf("Expression matrix %d genes x %d samples | Tumour %d Normal %d Unlabelled %d\n",
            nrow(E), ncol(E), sum(m$group=="Tumor", na.rm=TRUE),
            sum(m$group=="Normal", na.rm=TRUE), sum(is.na(m$group))))

ok <- !is.na(m$group)
E  <- E[, ok]; grp <- m$group[ok]
E  <- E[rowSums(E > 0) > 10, ]
tsamp <- colnames(E)[grp == "Tumor"]
cat(sprintf("Filtered %d x %d (Tumour %d)\n", nrow(E), ncol(E), length(tsamp)))

res_list <- list()

ti <- tempfile(fileext = ".txt")
fwrite(data.table(GENE = rownames(E), as.data.table(E)), ti, sep = "\t")
filterCommonGenes(input.f = ti, output.f = paste0(ti, ".f"), id = "GeneSymbol")
estimateScore(paste0(ti, ".f"), paste0(ti, ".est"), platform = "illumina")
d <- fread(paste0(ti, ".est"), skip = 2, header = TRUE, sep = "\t")
d <- d[, -2]                       # Remove the Description column.
setnames(d, 1, "metric")
mm <- t(as.matrix(d[, -1])); colnames(mm) <- d$metric
ES <- data.table(sample = rownames(mm), as.data.table(mm))
ES[, sample := gsub("\\.", "-", sample)]          # Restore hyphens in sample identifiers.
for (cc in setdiff(names(ES), "sample")) ES[[cc]] <- as.numeric(ES[[cc]])
fwrite(ES, file.path(OUT, "estimate_scores.csv"))
cat("\n=== [A] ESTIMATE ===\n"); print(head(ES, 3))
es_cols <- intersect(c("StromalScore","ImmuneScore","ESTIMATEScore","TumorPurity"), names(ES))
ES[, grp := ifelse(sample %in% tsamp, "Tumor", "Normal")]
cat("\n  Tumour vs Adjacent normal (median):\n")
print(ES[, lapply(.SD, function(x) round(median(x, na.rm=TRUE),2)), by = grp, .SDcols = es_cols])

MC <- file.path(DATA_DIR,"reference","MCPcounter","genes.txt")
sg <- fread(MC, sep = "\t"); setnames(sg, 1:2, c("gene","pop"))
sg <- sg[!is.na(gene) & nzchar(gene)]
cat(sprintf("\n=== [B] MCP-counter:%d cell populations / %d marker genes ===\n",
            uniqueN(sg$pop), nrow(sg)))
MCPL <- rbindlist(lapply(unique(sg$pop), function(p) {
  gs <- intersect(sg[pop == p]$gene, rownames(E))
  if (length(gs) < 3) { cat(sprintf("  Skip %-24s (only %d markers present in the matrix)\n", p, length(gs))); return(NULL) }
  v <- colMeans(pmax(E[gs, , drop = FALSE], 0))     # Mean non-negative log-scale marker expression.
  data.table(sample = colnames(E), pop = p, n_marker = length(gs), score = as.numeric(v))
}))
MCPA <- dcast(MCPL, sample ~ pop, value.var = "score")
setnames(MCPA, setdiff(names(MCPA), "sample"),
         paste0("MCP_", setdiff(names(MCPA), "sample")))
fwrite(MCPA, file.path(OUT, "mcpcounter_scores.csv"))
S <- MCPL[sample %in% tsamp, .(mean_score = round(mean(score),3), n_marker = n_marker[1]),
          by = pop][order(-mean_score)]
print(S)


suppressPackageStartupMessages({ library(data.table); library(Matrix) })
WS  <- OUTPUT_DIR
IM  <- file.path(RESULT_DIR,"extended","composition")
OUT <- file.path(RESULT_DIR,"extended","composition")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
SIG4 <- c("ABHD17C","ZDHHC6","ABHD17B","ZDHHC20")

# Column contract for the two published CSVs (Supplementary Data 1, Tables S10–S11
# sources). The export below is the only step that writes them.
M1_COLS <- c("unit","n","n_patients","beta_raw","lo_raw","hi_raw","beta_adj",
             "lo_adj","hi_adj","pct_beta_retained","p_adj")
M3_COLS <- c("unit","subgroup","n_tumor","n_normal","d_cohen","p","beta","ci_lo",
             "ci_hi","n_patients","wilcoxon_p_recomputed","statistic","p_BH")

E <- as.matrix(readRDS(file.path(MATRIX_DIR,"expr_full_logcpm.rds")))
m <- as.data.table(readRDS(file.path(MATRIX_DIR,"meta_full.rds")))
ok <- !is.na(m$group); E <- E[, ok]; mp <- m[ok]
grp <- factor(mp$group, levels=c("Normal","Tumor"))
cat(sprintf("bulk %d genes x %d samples (Tumor %d / Normal %d)\n",
            nrow(E), ncol(E), sum(grp=="Tumor"), sum(grp=="Normal")))

es  <- fread(file.path(IM, "estimate_scores.csv"))
mcp <- fread(file.path(IM, "mcpcounter_scores.csv"))
ss  <- fread(file.path(RESULT_DIR,"signature_scores_tcga.csv"))[, .(sample, risk)]
THY_MARK <- c("TG","TPO","TSHR","SLC5A5","PAX8","FOXE1","NKX2-1","SLC26A4","DIO1","DIO2")
gm  <- intersect(THY_MARK, rownames(E))
gth <- as.numeric(rowMeans(scale(t(E[gm, , drop=FALSE]))))

COMP <- data.table(sample = colnames(E), THY = gth)
COMP[, STR := es$StromalScore[match(sample, es$sample)]]
COMP[, IMM := es$ImmuneScore [match(sample, es$sample)]]
COMP[, FIB := mcp$`MCP_Fibroblasts`      [match(sample, mcp$sample)]]
COMP[, END := mcp$`MCP_Endothelial cells`[match(sample, mcp$sample)]]
COMP[, TCL := mcp$`MCP_T cells`          [match(sample, mcp$sample)]]
COMP[, risk := ss$risk[match(sample, ss$sample)]]
cat("Thyroid marker genes used:", paste(gm, collapse=", "), "\n")
cat("Covariates: THY (mean of ten scaled markers) + StromalScore + ImmuneScore + MCP_Fibroblasts + MCP_Endothelial + MCP_Tcells\n")
cat("Note: ssGSEA scores are not used as composition covariates.\n\n")

R <- rbindlist(lapply(c(SIG4, "SCORE"), function(g) {
  y <- if (g=="SCORE") COMP$risk else as.numeric(E[g, ])
  d <- data.table(y=y, grp=grp, THY=COMP$THY, STR=COMP$STR, IMM=COMP$IMM,
                  FIB=COMP$FIB, END=COMP$END, TCL=COMP$TCL)
  d <- d[complete.cases(d)]
  f0 <- lm(y ~ grp, d)
  f1 <- lm(y ~ grp + THY + STR + IMM + FIB + END + TCL, d)
  pid <- substr(COMP$sample[complete.cases(data.table(y=y, grp=grp, THY=COMP$THY, STR=COMP$STR, IMM=COMP$IMM,FIB=COMP$FIB,END=COMP$END,TCL=COMP$TCL))],1,12)
  a0 <- patient_cluster_coefficient(f0,pid);a1 <- patient_cluster_coefficient(f1,pid)
  c0 <- a0[c("lo","hi")]; c1 <- a1[c("lo","hi")]
  b0 <- coef(f0)["grpTumor"]; b1 <- coef(f1)["grpTumor"]
  # Published coefficients are written at full precision; rounding is applied
  # only in the display tables. n_patients and the patient-clustered P value
  # are part of the published M1 schema.
  data.table(unit=g, n=nrow(d), n_patients=unname(a1["patients"]),
             beta_raw=b0, lo_raw=c0[1], hi_raw=c0[2],
             beta_adj=b1, lo_adj=c1[1], hi_adj=c1[2],
             pct_beta_retained=100*abs(b1)/abs(b0),
             p_adj=unname(a1["p"]))
}))
print(R)
stopifnot(identical(names(R), M1_COLS))
fwrite(R[, ..M1_COLS], file.path(OUT, "M1_coefficients_with_CI.csv"))


ci2 <- ci[, .(sample, BRAF_V600E, RAS_MUT, SEQ)]
D3 <- data.table(sample = colnames(E), grp = as.character(grp))
D3 <- merge(D3, ci2, by="sample", all.x=TRUE)
D3[, sub := ifelse(grp=="Normal", "Normal",
            ifelse(!is.na(SEQ) & SEQ & BRAF_V600E==1, "BRAFmut",
            ifelse(!is.na(SEQ) & SEQ & RAS_MUT==1, "RASmut",
            ifelse(!is.na(SEQ) & SEQ, "driverneg", "unsequenced"))))]
print(table(D3$sub, useNA="ifany"))

M3 <- rbindlist(lapply(c(SIG4, "SCORE"), function(g) {
  y <- if (g=="SCORE") COMP$risk[match(D3$sample, COMP$sample)] else as.numeric(E[g, match(D3$sample, colnames(E))])
  d <- data.table(y=y, sub=D3$sub,patient=substr(D3$sample,1,12)); d <- d[!is.na(y) & !is.na(sub)]
  nrm <- d$y[d$sub=="Normal"]
  rbindlist(lapply(c("BRAFmut","RASmut","driverneg"), function(s) {
    t_ <- d$y[d$sub==s]
    if (length(t_) < 10) return(NULL)
    dz <- copy(d[sub %in% c("Normal",s)]);dz[,grp:=factor(ifelse(sub=="Normal","Normal","Tumor"),levels=c("Normal","Tumor"))]
    test <- patient_cluster_coefficient(lm(y~grp,dz),dz$patient)
    # Cohen's d is a descriptive effect size and is rounded for display; the
    # coefficient, interval and P values are written at full precision. The
    # published field is named beta (not mean_difference) so that the code
    # output and the frozen table share one schema. The asymptotic Wilcoxon
    # P is retained under its historical name for traceability; the stable
    # tail probability is supplied in S10_P_value_trace.csv.
    w <- suppressWarnings(stats::wilcox.test(t_, nrm))
    data.table(unit=g, subgroup=s, n_tumor=length(t_), n_normal=length(nrm),
               d_cohen = round((mean(t_)-mean(nrm))/sqrt((var(t_)+var(nrm))/2), 2),
               p = unname(test["p"]), beta=unname(test["beta"]),
               ci_lo=unname(test["lo"]), ci_hi=unname(test["hi"]),
               n_patients=unname(test["patients"]),
               wilcoxon_p_recomputed=w$p.value, statistic=unname(w$statistic))
  }))
}))
M3[, p_BH := p.adjust(p, "BH"), by=unit]
print(dcast(M3, unit ~ subgroup, value.var="d_cohen"))
cat("\nBH-adjusted P =\n"); print(dcast(M3, unit ~ subgroup, value.var="p_BH")[, lapply(.SD, function(x) signif(x,2)), by=unit])
stopifnot(identical(names(M3), M3_COLS))
fwrite(M3[, ..M3_COLS], file.path(OUT, "M3_by_driver.csv"))

# Retain the actual proxy matrix, collinearity checks and two bounded alternatives.
DD<-data.table(sample=COMP$sample,pid=substr(COMP$sample,1,12),grp=grp,
 THY=COMP$THY,STR=COMP$STR,IMM=COMP$IMM,FIB=COMP$FIB,END=COMP$END,TCL=COMP$TCL)
fwrite(DD,file.path(OUT,"composition_inputs.csv"))
zz<-as.matrix(scale(DD[,.(THY,STR,IMM,FIB,END,TCL)]))
cc<-cor(zz,use="complete.obs")
fwrite(data.table(variable=rownames(cc),as.data.table(cc)),file.path(OUT,"composition_proxy_correlations.csv"))
vif<-sapply(seq_len(ncol(zz)),function(j)1/(1-summary(lm(zz[,j]~zz[,-j]))$r.squared))
dm<-model.matrix(~grp+THY+STR+IMM+FIB+END+TCL,DD)
jsonlite::write_json(list(n=nrow(DD),patients=uniqueN(DD$pid),rank=qr(dm)$rank,
 columns=ncol(dm),vif=setNames(vif,colnames(zz)),scaled_condition_number=kappa(cbind(model.matrix(~grp,DD),zz)),
 max_abs_correlation=max(abs(cc[upper.tri(cc)]))),file.path(OUT,"composition_diagnostics.json"),auto_unbox=TRUE,pretty=TRUE,digits=16)
fms<-list(thyroid_ESTIMATE=y~grp+THY+STR+IMM,thyroid_lineage_markers=y~grp+THY+FIB+END+TCL)
sensitivity<-rbindlist(lapply(c(SIG4,"SCORE"),function(g)rbindlist(lapply(names(fms),function(nm){
 d<-copy(DD);d[,y:=if(g=="SCORE")COMP$risk else as.numeric(E[g,])]
 d<-d[complete.cases(d)];r<-patient_cluster_coefficient(lm(fms[[nm]],d),d$pid)
 data.table(outcome=g,model=nm,beta=unname(r["beta"]),ci_lo=unname(r["lo"]),ci_hi=unname(r["hi"]),p=unname(r["p"]))
}))))
fwrite(sensitivity,file.path(OUT,"composition_reduced_models.csv"))


# Descriptive residual diagnostics; HC1 inference does not assume constant variance.
residual_checks<-rbindlist(lapply(c(SIG4,"SCORE"),function(g){
 d<-copy(DD);d[,y:=if(g=="SCORE")COMP$risk else as.numeric(E[g,])]
 d<-d[complete.cases(d)];f<-lm(y~grp+THY+STR+IMM+FIB+END+TCL,d)
 r<-resid(f);h<-hatvalues(f);cook<-cooks.distance(f)
 data.table(outcome=g,n=nrow(d),residual_mean=mean(r),residual_sd=sd(r),
  max_abs_studentized_residual=max(abs(rstudent(f))),max_leverage=max(h),
  max_cooks_distance=max(cook),spearman_abs_residual_fitted=cor(abs(r),fitted(f),method="spearman"))
}))
fwrite(residual_checks,file.path(OUT,"composition_residual_diagnostics.csv"))
