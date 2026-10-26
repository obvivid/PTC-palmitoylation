# Epithelial-lineage identity review and raw-count aggregation.
# Run after module 18; input metadata and optional recorded doublet calls are
# supplied in Supplementary Data 1 / singlecell_lineage.
set.seed(20261022)
suppressPackageStartupMessages({library(data.table); library(Matrix); library(ggplot2)})

GATE <- file.path(PROJECT_ROOT, 'data/singlecell_lineage')
OUT <- file.path(RESULT_DIR, 'extended/singlecell_lineage/01_annotation')
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
fread <- function(..., tmpdir=OUT) data.table::fread(..., tmpdir=tmpdir)
gate_out <- file.path(RESULT_DIR, 'extended/singlecell_lineage/00_gate')
dir.create(gate_out, recursive=TRUE, showWarnings=FALSE)
file.copy(file.path(GATE, 'sample_donor_ledger.csv'), file.path(gate_out, 'sample_donor_ledger.csv'), overwrite=TRUE)
DATA <- file.path(PROJECT_ROOT, 'data/geo/GSE193581')
RDS_V7 <- if (exists('SINGLECELL_COUNT_OBJECT')) SINGLECELL_COUNT_OBJECT else file.path(RESULT_DIR, 'extended/singlecell/counts_annotated_193581.rds')
SIG4 <- c("ABHD17C","ZDHHC6","ABHD17B","ZDHHC20")

## ---------- 1. 载入与映射 ----------
L <- readRDS(RDS_V7)
C <- L$counts; CT <- as.character(L$celltype); TT <- as.character(L$tissue)
v7id <- colnames(C)
tmp  <- sub("^GSM[0-9]+_", "", v7id)
cand <- ifelse(grepl("^([^_]+)_\\1_", tmp), sub("^([^_]+)_\\1_", "\\1_", tmp), tmp)
SAMPLE <- sub("^GSM[0-9]+_([^_]+)_.*$", "\\1", v7id)

auth <- fread(file.path(DATA, "GSE193581_celltype_annotation.txt.gz"),
              header = FALSE, skip = 1, col.names = c("cell", "sample", "celltype"))
m <- match(cand, auth$cell)
AUTHOR <- auth$celltype[m]

sel <- TT %in% c("NORM", "PTC")
idx <- which(sel)
cat("NORM/PTC 细胞:", length(idx), "\n")

## ---------- 2. 模块得分（稀疏逐基因累计，避免稠密矩阵）----------
MODULES <- list(
  epithelial  = c("EPCAM","KRT8","KRT18","KRT19","CDH1"),
  thyroid_diff= c("TG","TPO","TFF3","SLC26A4","SLC5A5","TSHR","PAX8","FOXE1","DIO1","DIO2"),
  tfc         = c("TFF3","TPO","SLC26A4"),
  fibro       = c("COL1A1","COL1A2","DCN","LUM"),
  bcell       = c("CD79A","MS4A1","MZB1","IGKC","JCHAIN"),
  tcell       = c("CD3D","CD3E","TRAC","IL7R"),
  myeloid     = c("LYZ","CD68","AIF1","C1QA","FCN1"),
  endo        = c("PECAM1","VWF","EGFL7","CLDN5")
)
stopifnot(!any(SIG4 %in% unlist(MODULES)))  # 身份标志禁用四基因（既定规则）

lib <- Matrix::colSums(C)
lib[lib == 0] <- 1
nn <- length(idx)
SC <- matrix(0, nrow = nn, ncol = length(MODULES),
             dimnames = list(v7id[idx], names(MODULES)))
Csub <- C[, idx, drop = FALSE]
libsub <- lib[idx]
for (md in names(MODULES)) {
  gset <- intersect(MODULES[[md]], rownames(Csub))
  acc <- numeric(nn)
  for (g in gset) {
    vv <- as.numeric(Csub[g, , drop = TRUE])
    nz <- which(vv > 0)
    if (!length(nz)) next
    acc[nz] <- acc[nz] + log1p(vv[nz] / libsub[nz] * 1e4)
  }
  SC[, md] <- acc / length(MODULES[[md]])   # 均值按完整集合大小（未检出=0）
}
sc <- as.data.table(SC); sc$cell <- rownames(SC)

## ---------- 3. scDblFinder 按样本（缓存 dbl_calls_cache.csv.gz，避免重跑）----------
suppressPackageStartupMessages({library(SingleCellExperiment); library(scDblFinder)})
CACHE <- file.path(OUT, "dbl_calls_cache.csv.gz")
# Reuse recorded calls when supplied. To recompute them, omit the input copy
# in GATE and the generated copy in OUT before running this step.
if (!file.exists(CACHE) && file.exists(file.path(GATE, 'dbl_calls_cache.csv.gz')))
    file.copy(file.path(GATE, 'dbl_calls_cache.csv.gz'), CACHE)
if (file.exists(CACHE)) {
  dbl <- fread(CACHE)
  stopifnot(all(dbl$cell == v7id[idx]))
  dbl_call <- dbl$call; names(dbl_call) <- dbl$cell
  dbl_score <- dbl$score; names(dbl_score) <- dbl$cell
  qc <- data.table(sample=SAMPLE[idx], call=dbl$call)[,
       .(n_cells=.N, n_doublet=sum(call=='doublet'),
         pct_doublet=round(100*mean(call=='doublet'),2)),by=sample]
  fwrite(qc,file.path(OUT,'qc_summary.csv'))
  cat("载入 scDblFinder 缓存:", CACHE, "\n")
} else {
dbl_call <- character(nn); names(dbl_call) <- v7id[idx]
dbl_score <- numeric(nn); names(dbl_score) <- v7id[idx]
qc_rows <- list()
for (s in unique(SAMPLE[idx])) {
  j <- which(SAMPLE[idx] == s)
  sce <- SingleCellExperiment(list(counts = Csub[, j, drop = FALSE]))
  set.seed(20261022)  # 每样本同一种子起点；记录于脚本
  sce <- scDblFinder(sce, dbr = NULL, verbose = FALSE)
  dbl_call[j] <- as.character(sce$scDblFinder.class)
  dbl_score[j] <- as.numeric(sce$scDblFinder.score)
  qc_rows[[s]] <- data.table(sample = s, n_cells = length(j),
                             n_doublet = sum(dbl_call[j] == "doublet"),
                             pct_doublet = round(100 * mean(dbl_call[j] == "doublet"), 2))
  cat(sprintf("  %s: %d cells, %.2f%% doublets\n", s, length(j), qc_rows[[s]]$pct_doublet))
}
qc <- rbindlist(qc_rows)
fwrite(qc, file.path(OUT, "qc_summary.csv"))
fwrite(data.table(cell = v7id[idx], call = dbl_call, score = dbl_score), CACHE)
}

## ---------- 4. 身份分类与冻结 ----------
epi <- sc$epithelial; thy <- sc$thyroid_diff; tfc <- sc$tfc
fib <- sc$fibro; bc <- sc$bcell; tc <- sc$tcell; mye <- sc$myeloid; en <- sc$endo
imm <- pmax(tc, bc, mye)

author_sel <- AUTHOR[idx]
in_pool <- author_sel %in% c("Malignant cell", "Epithelial cell")
is_doublet <- dbl_call == "doublet"

immune_conflict <- imm > epi & imm > 1
fibro_conflict  <- fib > 1.5 * epi & fib > 1 & epi < 1
endo_conflict   <- en  > 1.5 * epi & en  > 1 & epi < 1
conflict_kind <- ifelse(immune_conflict, "immune",
                        ifelse(fibro_conflict, "fibro",
                               ifelse(endo_conflict, "endo", "")))

final_label <- rep("non_lineage", nn)
final_label[in_pool & is_doublet] <- "lineage_removed_doublet"
final_label[in_pool & !is_doublet & (immune_conflict | fibro_conflict | endo_conflict)] <-
  paste0("lineage_uncertain_", conflict_kind[in_pool & !is_doublet & (immune_conflict | fibro_conflict | endo_conflict)])
final_label[in_pool & !is_doublet & !(immune_conflict | fibro_conflict | endo_conflict)] <- "lineage_confident"

q2_class <- rep(NA_character_, nn)
q2_class[author_sel == "Malignant cell" & !is_doublet & !immune_conflict] <- "malignant"
q2_class[author_sel == "Malignant cell" & !is_doublet & immune_conflict] <- "malignant_uncertain"
q2_class[author_sel == "Malignant cell" & is_doublet] <- "malignant_uncertain"
q2_class[author_sel == "Epithelial cell" & !is_doublet & !immune_conflict & tfc > 0] <- "normal_follicular"
q2_class[author_sel == "Epithelial cell" & !is_doublet & (immune_conflict | tfc == 0)] <- "normal_follicular_uncertain"
q2_class[author_sel == "Epithelial cell" & is_doublet] <- "normal_follicular_uncertain"

ann <- data.table(
  cell = v7id[idx], sample = SAMPLE[idx], tissue = TT[idx],
  v7_celltype = CT[idx], author_celltype = author_sel,
  epi = epi, thyroid_diff = thy, tfc = tfc, fibro = fib, bcell = bc,
  tcell = tc, myeloid = mye, endo = en,
  scDblFinder_class = dbl_call, scDblFinder_score = dbl_score,
  final_label = final_label, q2_class = q2_class,
  in_candidate_pool = in_pool, conflict_kind = conflict_kind)
fwrite(ann, file.path(OUT, "cell_annotations.tsv.gz"), sep = "\t")

## 遗漏细胞扫描（只报告，不纳入）：池外但上皮得分高的细胞
missed <- ann[!in_candidate_pool & epi >= 2]
ev1 <- ann[, .(n = .N), by = .(sample, class = final_label)]
ev <- if (nrow(missed)) {
  rbind(ev1, missed[, .(n = .N), by = .(sample, class = paste0("missed_scan_epiGE2_", v7_celltype))], fill = TRUE)
} else ev1
cat("遗漏扫描（池外 epi>=2 细胞）:", nrow(missed), "个\n")
fwrite(ev, file.path(OUT, "annotation_evidence.csv"))

cat("\nfinal_label 分布:\n"); print(table(ann$final_label))
cat("\nq2_class 分布（谱系池内）:\n"); print(table(ann$q2_class[ann$final_label == "lineage_confident"]))
cat("\n各样本 lineage_confident 计数:\n")
ps <- ann[final_label == "lineage_confident", .(
  n_confident = .N,
  n_malignant = sum(q2_class == "malignant"),
  n_normal_follicular = sum(q2_class == "normal_follicular"),
  n_uncertain_q2 = sum(q2_class %in% c("malignant_uncertain", "normal_follicular_uncertain")),
  median_umi = median(libsub[match(cell, v7id[idx])])
), by = .(sample, tissue)]
setorder(ps, tissue, sample)
fwrite(ps, file.path(OUT, "frozen_per_sample_counts.csv"))
print(ps)

# Freeze the reviewed cell set and specimen-level raw-count sums.
frozen <- ann[final_label == "lineage_confident", .(cell, sample, tissue, q2_class)]
fwrite(frozen, file.path(OUT, "frozen_cell_list.csv"))

fr <- match(frozen$cell, v7id)
fr_samples <- sort(unique(frozen$sample))
pb_list <- split(seq_len(nrow(frozen)), factor(frozen$sample, levels = fr_samples))
pb <- vapply(pb_list, function(ii)
  as.numeric(Matrix::rowSums(C[, fr[ii], drop = FALSE])), numeric(nrow(C)))
rownames(pb) <- rownames(C)
storage.mode(pb) <- "integer"
meta <- fread(file.path(GATE, "sample_donor_ledger.csv"))
pb_meta <- meta[match(fr_samples, sample)]
fwrite(data.table(gene = rownames(pb), as.data.table(pb)),
       file.path(OUT, "pseudobulk_counts_frozen.csv.gz"))
fwrite(data.table(sample = fr_samples,
                  tissue = ifelse(grepl("^PTC", fr_samples), "PTC", "NORM"),
                  donor_id = pb_meta$donor_id,
                  treatment = pb_meta$treatment, n_cells = as.integer(table(frozen$sample)[fr_samples])),
       file.path(OUT, "pseudobulk_metadata_frozen.csv"))

