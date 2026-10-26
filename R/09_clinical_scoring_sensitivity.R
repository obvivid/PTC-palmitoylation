# Clinical effects and scoring sensitivity
# Run through run_analysis.R; see docs/METHODS.md.

suppressMessages({ library(data.table); library(pROC) })
ROOT <- PROJECT_ROOT
DIR  <- RESULT_DIR
set.seed(11)
sec <- function(t){ cat("\n"); cat(strrep("=",78),"\n",sep=""); cat(t,"\n"); cat(strrep("=",78),"\n",sep="") }

cliff <- function(x, y) {
  n1 <- length(x); n2 <- length(y)
  (sum(outer(x, y, ">")) - sum(outer(x, y, "<"))) / (n1 * n2)
}
med_ci <- function(x, y, n=2000) {
  a <- replicate(n, median(sample(x, replace=TRUE)) - median(sample(y, replace=TRUE)))
  quantile(a, c(.025,.975))
}

sec("[a] Clinical associations: NX excluded / histology sensitivity / BH recomputed")
d <- fread(file.path(DIR, "clinical_immune_fixed.csv"), na.strings=c("","NA"))
d <- d[!is.na(risk)]
cat(sprintf("Tumor samples %d; n_stage distribution: %s\n", nrow(d),
            paste(names(table(d$n_stage, useNA="ifany")), table(d$n_stage, useNA="ifany"), sep="=", collapse=" ")))

out <- list()
addrow <- function(lbl, n, type, eff_name, eff, lo, hi, p, note) {
  out[[length(out)+1]] <<- data.table(var=lbl, n=n, type=type, effect=eff_name,
                                      eff=eff, ci_lo=lo, ci_hi=hi, p=p, note=note)
}

# Effect direction is stated explicitly and follows the manuscript table:
# mutant minus wild type for the two mutation flags (coded 1 and 0), female minus
# male for sex. Relying on factor() defaults orders levels alphabetically and would
# silently reverse the two mutation comparisons.
BINVAR <- list(
  list(lbl="BRAF V600E", f="BRAF_V600E", pos="1", neg="0", pos_lbl="mutant", neg_lbl="wild type"),
  list(lbl="RAS mutation", f="RAS_MUT", pos="1", neg="0", pos_lbl="mutant", neg_lbl="wild type"),
  list(lbl="Sex", f="sex", pos="Female", neg="Male", pos_lbl="female", neg_lbl="male"))
for (v in BINVAR) {
  x <- d[[v$f]]; y <- d$risk; k <- !is.na(x)
  g <- factor(x[k]); stopifnot(v$pos %in% levels(g), v$neg %in% levels(g))
  a <- y[k][g==v$pos]; b <- y[k][g==v$neg]
  w <- stable_wilcox_test(a, b); ci <- med_ci(a, b)
  cat(sprintf("  %-12s n=%3d | median difference(%s-%s) %+7.3f (%.3f~%.3f) | delta %+.3f | P=%.3g\n",
              v$lbl, sum(k), v$pos_lbl, v$neg_lbl, median(a)-median(b), ci[1], ci[2],
              cliff(a,b), w$p.value))
  addrow(v$lbl, sum(k), "two-group", "Cliff\'s δ", cliff(a,b), ci[1], ci[2], w$p.value,
         paste(v$pos_lbl, "minus", v$neg_lbl))
}

x <- d$n_stage; y <- d$risk
k <- !is.na(x) & x != "NX"
g <- factor(x[k]); kk <- kruskal.test(y[k] ~ g)
eps2 <- as.numeric(kk$statistic) / (sum(k) - 1)
cat(sprintf("  N status (NX excluded) n=%3d | groups: %s | KW P=%.3g | eps2=%.4f\n",
            sum(k), paste(levels(g), collapse="/"), kk$p.value, eps2))
addrow("N status (NX excluded)", sum(k), "multi-group", "ε²", eps2, NA, NA, kk$p.value,
       paste(nlevels(g), "groups: N0/N1/N1A/N1B"))

n0 <- y[k][g=="N0"]; np <- y[k][g!="N0"]
w2 <- stable_wilcox_test(n0, np)
cat(sprintf("    [sensitivity] N0 vs N+ : %d vs %d, Wilcoxon P=%.3g, delta(N0 vs N+)=%+.3f\n",
            length(n0), length(np), w2$p.value, cliff(n0, np)))

for (v in list(list(lbl="T stage", f="t4"), list(lbl="AJCC stage", f="stage4"),
               list(lbl="Histology", f="histgrp"))) {
  x <- d[[v$f]]; k <- !is.na(x); g <- factor(x[k])
  kk <- kruskal.test(d$risk[k] ~ g); eps2 <- as.numeric(kk$statistic)/(sum(k)-1)
  cat(sprintf("  %-12s n=%3d | KW P=%.3g | eps2=%.4f (%d groups: %s)\n",
              v$lbl, sum(k), kk$p.value, eps2, nlevels(g), paste(levels(g), collapse="/")))
  addrow(v$lbl, sum(k), "multi-group", "ε²", eps2, NA, NA, kk$p.value, paste0(nlevels(g)," groups"))
}

kh <- !is.na(d$histgrp) & !d$histgrp %in% c("Hurthle","Other")
gh <- factor(d$histgrp[kh]); kkh <- kruskal.test(d$risk[kh] ~ gh)
cat(sprintf("    [sensitivity] histology, three main types only n=%3d | KW P=%.3g | eps2=%.4f\n",
            sum(kh), kkh$p.value, as.numeric(kkh$statistic)/(sum(kh)-1)))

k <- !is.na(d$age)
r <- cor.test(d$risk[k], d$age[k], method="spearman")
bs <- replicate(1000, { s <- sample(sum(k), replace=TRUE)
  suppressWarnings(cor(d$risk[k][s], d$age[k][s], method="spearman")) })
ci <- quantile(bs, c(.025,.975), na.rm=TRUE)
cat(sprintf("  %-12s n=%3d | Spearman ρ %+.4f (%.3f~%.3f) | P=%.3g\n",
            "Age", sum(k), r$estimate, ci[1], ci[2], r$p.value))
addrow("Age", sum(k), "continuous", "Spearman ρ", as.numeric(r$estimate), ci[1], ci[2], r$p.value, "")

res <- rbindlist(out)
res[,median_difference:=NA_real_]
for(v in BINVAR){
 a<-d[!is.na(get(v$f)) & as.character(get(v$f))==v$pos]$risk
 b<-d[!is.na(get(v$f)) & as.character(get(v$f))==v$neg]$risk
 res[var==v$lbl,median_difference:=median(a)-median(b)]
}
res[, p_BH := p.adjust(p, "BH")]
res[, sigBH := p_BH < 0.05]
res[, ci_quantity := fifelse(type=="two-group","median_difference",fifelse(type=="continuous","spearman_rho","not_estimated"))]
sec("[BH-corrected (8 variables, N with NX excluded)]")
print(res[, .(var, n, eff=round(eff,4), p=signif(p,3), p_BH=signif(p_BH,3), sigBH)])
cat(sprintf("  Significant after BH: %d / %d\n", sum(res$sigBH), nrow(res)))

sec("[b] Scoring sensitivity check (final definition: standardisation within the target samples)")
SIG   <- fread(file.path(DIR, "signature_genes.csv"))$gene
model <- readRDS(file.path(DIR, "signature_model.rds"))
COEF  <- coef(model$model); INT <- COEF[1]; COEF <- COEF[-1]; names(COEF) <- SIG
cat(sprintf("  Full model definition: intercept %.6f ; %s\n", INT,
            paste(sprintf("%s=%+.6f", SIG, COEF), collapse="  ")))

ann <- fread(file.path(ROOT,"data/geo/GPL570.annot.gz"), skip=27, header=TRUE, sep="\t",
             quote="", na.strings=c("","NA"), fill=TRUE); setnames(ann,1,"probe")
A <- ann[!is.na(`Gene symbol`) & `Gene symbol` %in% SIG,
         .(probe=as.character(probe), symbol=as.character(`Gene symbol`))]

read_gpl570 <- function(fp) {
  hl <- grep("!series_matrix_table_begin", readLines(gzfile(fp), n=200))[1]
  M  <- fread(fp, skip=hl, header=TRUE, na.strings=c("","NA")); setnames(M,1,"probe")
  M[, probe := gsub('"', '', as.character(probe))]; M
}
getchar <- function(fp) {
  tl <- readLines(gzfile(fp), n=200)
  ln <- grep("^!Sample_characteristics_ch1", tl, value=TRUE)[1]
  strsplit(gsub('"','',sub("^!Sample_characteristics_ch1\\t?","",ln)), "\t")[[1]]
}
extract <- function(fp, log2_first=FALSE) {
  M <- read_gpl570(fp); sc <- setdiff(names(M), "probe")
  if (log2_first) for (cc in sc) set(M, j=cc, value=as.numeric(suppressWarnings(log2(as.numeric(M[[cc]])))))
  sub <- merge(A, M, by="probe")
  G <- sub[, lapply(.SD, function(x) mean(as.numeric(x), na.rm=TRUE)), by=symbol, .SDcols=sc]
  Gm <- as.matrix(G[, ..sc]); rownames(Gm) <- G$symbol; Gm
}
score_target <- function(Gm, keep, zscore=TRUE) {
  gs <- intersect(SIG, rownames(Gm))
  X  <- Gm[gs, keep, drop=FALSE]
  Z  <- if (zscore) t(scale(t(X))) else X
  as.numeric(t(Z) %*% COEF[gs])
}
auc_of <- function(lab, lp) {
  r <- roc(lab, lp, levels=c("Normal","Tumor"), direction="<", quiet=TRUE)
  as.numeric(r$auc)
}

fp30 <- file.path(ROOT,"data/geo/GSE33630_series_matrix.txt.gz")
G30  <- extract(fp30)
dx   <- gsub("pathological diagnostic: ", "", getchar(fp30))
lab  <- ifelse(grepl("non-tumor|normal", dx, ignore.case=TRUE), "Normal",
        ifelse(grepl("papillary", dx, ignore.case=TRUE), "PTC",
        ifelse(grepl("anaplastic", dx, ignore.case=TRUE), "ATC", NA)))
keep <- which(lab %in% c("PTC","Normal"))
cat(sprintf("  GSE33630 PTC/normal = %d/%d\n", sum(lab[keep]=="PTC"), sum(lab[keep]=="Normal")))
a1 <- auc_of(ifelse(lab[keep]=="PTC","Tumor","Normal"), score_target(G30, keep, TRUE))
a0 <- auc_of(ifelse(lab[keep]=="PTC","Tumor","Normal"), score_target(G30, keep, FALSE))
cat(sprintf("  GSE33630 PTC vs normal: within-target z = %.4f ; no z = %.4f\n", a1, a0))

fp65 <- file.path(ROOT,"data/geo/GSE65144_series_matrix.txt.gz")
ch   <- gsub("tissue type: ", "", getchar(fp65))
lab2 <- ifelse(grepl("normal|non-tumor|adjacent", ch, ignore.case=TRUE), "Normal",
        ifelse(grepl("anaplastic|carcinoma|tumor", ch, ignore.case=TRUE), "ATC", NA))
keep2 <- which(lab2 %in% c("ATC","Normal"))
G65l <- extract(fp65, log2_first=TRUE)
G65r <- extract(fp65, log2_first=FALSE)
b1 <- auc_of(ifelse(lab2[keep2]=="ATC","Tumor","Normal"), score_target(G65l, keep2, TRUE))
b0 <- auc_of(ifelse(lab2[keep2]=="ATC","Tumor","Normal"), score_target(G65r, keep2, TRUE))
cat(sprintf("  GSE65144 ATC/normal = %d/%d\n", sum(lab2[keep2]=="ATC"), sum(lab2[keep2]=="Normal")))
cat(sprintf("  GSE65144 ATC vs normal: z after log2 = %.4f ; z on linear scale = %.4f\n", b1, b0))

atc_keep <- which(lab %in% c("ATC","Normal"))
atc30 <- roc(ifelse(lab[atc_keep]=="ATC","Tumor","Normal"), score_target(G30,atc_keep,TRUE),
             levels=c("Normal","Tumor"),direction="<",quiet=TRUE)
atc65 <- roc(ifelse(lab2[keep2]=="ATC","Tumor","Normal"), score_target(G65l,keep2,TRUE),
             levels=c("Normal","Tumor"),direction="<",quiet=TRUE)
atc_out <- rbindlist(lapply(list(list(tag="GSE33630",ro=atc30,labels=lab[atc_keep]),
                               list(tag="GSE65144",ro=atc65,labels=lab2[keep2])),function(d) {
  ci <- ci.auc(d$ro,method="delong")
  data.table(dataset=d$tag,n_atc=sum(d$labels=="ATC"),n_normal=sum(d$labels=="Normal"),
             auc=as.numeric(d$ro$auc),delong_lo=ci[1],delong_hi=ci[3])
}))
fwrite(atc_out,file.path(DIR,"atc_validation_summary.csv"))
sens <- data.table(
  item=c("GSE33630_PTC_z_within_target","GSE33630_PTC_no_z","GSE65144_ATC_log2_z","GSE65144_ATC_linear_z"),
  auc=c(a1,a0,b1,b0))
fwrite(res,  file.path(DIR, "revision_clinical_exclNX_2026-10-02.csv"))
fwrite(sens, file.path(DIR, "revision_scoring_sensitivity_2026-10-02.csv"))
cat(sprintf("\nOutput:\n  %s\n  %s\n",
            file.path(DIR, "revision_clinical_exclNX_2026-10-02.csv"),
            file.path(DIR, "revision_scoring_sensitivity_2026-10-02.csv")))
