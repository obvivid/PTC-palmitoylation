# Clinical annotation and custom ssGSEA composition
# Run through run_analysis.R; see docs/METHODS.md.

suppressMessages({ library(data.table); library(GSVA); library(jsonlite) })
ROOT <- PROJECT_ROOT
FM   <- MATRIX_DIR
DIR  <- RESULT_DIR
sec  <- function(t){ cat("\n"); cat(strrep("=",78),"\n",sep=""); cat(t,"\n"); cat(strrep("=",78),"\n",sep="") }
sid15 <- function(x) substr(as.character(x), 1, 15)

E    <- readRDS(file.path(FM, "expr_full_logcpm.rds"))
meta <- as.data.table(readRDS(file.path(FM, "meta_full.rds")))
mi   <- match(colnames(E), meta$sample); E <- E[, !is.na(mi), drop=FALSE]; meta <- meta[mi[!is.na(mi)]]
tu   <- !is.na(meta$group) & meta$group == "Tumor"
E    <- E[, tu, drop=FALSE]; meta <- meta[tu]
stopifnot(all(meta$sample == colnames(E)))
meta[, sid := sid15(sample)]
cat(sprintf("Tumor samples %d\n", ncol(E)))

SIG   <- fread(file.path(DIR, "signature_genes.csv"))$gene
model <- readRDS(file.path(DIR, "signature_model.rds"))
COEF  <- coef(model$model)[-1]; names(COEF) <- SIG
meta[, risk := as.numeric(t(E[SIG, , drop=FALSE]) %*% COEF[SIG])]

sec("[1] Tissue composition scores (custom marker sets, ssGSEA, not ESTIMATE; recomputed on the full matrix)")
STROMA <- c("COL1A1","COL1A2","COL3A1","COL5A1","COL5A2","COL6A1","COL6A2","COL6A3",
            "FAP","PDGFRB","THY1","DCN","LUM","FN1","ACTA2","TAGLN","POSTN","MMP2","VIM","SPARC")
IMMUNE <- c("PTPRC","CD3D","CD3E","CD2","CD8A","CD19","MS4A1","CD79A","CD68","CD14",
            "FCGR3A","NKG7","KLRD1","IL2RB","LYZ","AIF1","TYROBP","CCL5","CXCL9","GZMA")
THYRO  <- c("TG","TSHR","TPO","SLC5A5","PAX8","FOXE1","NKX2-1")
for (nm in c("STROMA","IMMUNE","THYRO")) {
  orig <- get(nm)
  g    <- intersect(orig, rownames(E))
  assign(nm, g)
  cat(sprintf("  %-8s markers %2d/%2d present; missing: %s\n", nm, length(g), length(orig),
              paste(setdiff(orig, g), collapse=", ")))
}
ss <- gsva(ssgseaParam(E, list(STROMA=STROMA, IMMUNE=IMMUNE, THYRO=THYRO)), verbose=FALSE)
meta[, STROMA := as.numeric(ss["STROMA",])]
meta[, IMMUNE := as.numeric(ss["IMMUNE",])]
meta[, THYRO  := as.numeric(ss["THYRO", ])]
meta[, COMPOS := scale(STROMA)[,1] + scale(IMMUNE)[,1] - scale(THYRO)[,1]]

sec("[2] Mutation annotation: unsequenced samples have missing mutation status")
sq   <- unique(sid15(fromJSON(file.path(ROOT,"data/tcga/mutations/sequenced_samples.json"))))
mut  <- as.data.table(fromJSON(file.path(ROOT,"data/tcga/mutations/thca_mutations.json")))
b600 <- unique(sid15(mut[grepl("^BRAF V600", keyword)]$sampleId))
rasu <- unique(sid15(mut[grepl("^(HRAS|KRAS|NRAS) ", keyword)]$sampleId))
meta[, SEQ := sid %in% sq]
meta[, BRAF_V600E := fifelse(!SEQ, NA_integer_, as.integer(sid %in% b600))]
meta[, RAS_MUT    := fifelse(!SEQ, NA_integer_, as.integer(sid %in% rasu))]
cat(sprintf("  %d samples without sequencing information; mutation status recorded as missing\n", sum(!meta$SEQ)))
cat(sprintf("  Among sequenced samples: BRAF V600E %d, RAS mutation %d, both negative %d\n",
            sum(meta$BRAF_V600E == 1, na.rm=TRUE),
            sum(meta$RAS_MUT    == 1, na.rm=TRUE),
            sum(meta$BRAF_V600E == 0 & meta$RAS_MUT == 0, na.rm=TRUE)))

sec("[3] Signature vs tissue composition (effective N reported)")
for (f in c("STROMA","IMMUNE","THYRO","COMPOS")) {
  k <- !is.na(meta[[f]]) & !is.na(meta$risk)
  r <- cor.test(meta$risk[k], meta[[f]][k])
  cat(sprintf("  %-8s n=%d  r=%+.3f  P=%.3g\n", f, sum(k), r$estimate, r$p.value))
}

sec("[4] Signature vs clinical/mutation variables (unsequenced samples excluded)")
clu <- read_clinical_annotation(file.path(ROOT, "data/tcga/data_clinical_patient.txt"), meta$sample)
im <- match(meta$sample, clu$sample)
for (f in c("age","sex","stage4","t4","n_pos","n_stage","histgrp")) meta[[f]] <- clu[[f]][im]

assoc2 <- function(lbl, v) {
  k <- !is.na(v) & !is.na(meta$risk)
  if (sum(k) < 20) { cat(sprintf("  %-12s insufficient samples (%d)\n", lbl, sum(k))); return(invisible()) }
  v2 <- v[k]; y <- meta$risk[k]
  if (length(unique(v2)) == 2) {
    g <- factor(v2); w <- stable_wilcox_test(y ~ g)
    cat(sprintf("  %-12s n=%d  %s %.3f vs %.3f  P=%.3g\n", lbl, sum(k),
                paste(levels(g), collapse="/"),
                median(y[g == levels(g)[1]]), median(y[g == levels(g)[2]]), w$p.value))
  } else {
    g <- factor(v2); kk <- kruskal.test(y ~ g)
    cat(sprintf("  %-12s n=%d  Kruskal-Wallis P=%.3g\n", lbl, sum(k), kk$p.value))
  }
}
assoc2("BRAF V600E", meta$BRAF_V600E)
assoc2("RAS mutation", meta$RAS_MUT)
assoc2("N status",     meta$n_stage)
assoc2("T stage",      meta$t4)
assoc2("AJCC stage",   meta$stage4)
assoc2("Histology",    meta$histgrp)
assoc2("Sex",          meta$sex)
k <- !is.na(meta$age) & !is.na(meta$risk); r <- cor.test(meta$risk[k], meta$age[k])
cat(sprintf("  %-12s n=%d  r=%+.3f  P=%.3g\n", "Age", sum(k), r$estimate, r$p.value))

fwrite(meta[, .(sample, sid, risk, STROMA, IMMUNE, THYRO, COMPOS,
                SEQ, BRAF_V600E, RAS_MUT, age, sex, stage4, t4, n_pos, n_stage, histgrp)],
       file.path(DIR, "clinical_immune_fixed.csv"), na="NA")
cat(sprintf("  output: %s\n", file.path(DIR, "clinical_immune_fixed.csv")))
