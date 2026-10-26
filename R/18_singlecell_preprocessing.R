# GSE193581 preprocessing and broad marker-based cell annotation.
DATA_DIR <- file.path(PROJECT_ROOT, "data")

suppressPackageStartupMessages({
  library(data.table); library(Matrix); library(Seurat); library(harmony)
})
SC <- file.path(DATA_DIR, "geo", "GSE193581", "umi")
OUT <- file.path(RESULT_DIR, "extended", "singlecell")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
set.seed(20261007)
SIG4 <- c("ABHD17C","ZDHHC6","ABHD17B","ZDHHC20")

files <- list.files(SC, pattern = "_UMI\\.txt\\.gz$")
cat("Samples:", length(files), "\n")

CACHE <- file.path(OUT, "seurat_GSE193581.rds")
if (file.exists(CACHE)) {
  cat("Loading cache:", CACHE, "\n")
  obj <- readRDS(CACHE); cat("  Cells:", ncol(obj), "\n")
} else {
obj_list <- list(); qc <- list()
for (f in files) {
  nm <- sub("_UMI\\.txt\\.gz$", "", f)
  gsm <- sub("_.*$", "", nm); lab <- sub("^[^_]+_", "", nm); grp <- sub("[0-9]+$", "", lab)
  ln1 <- readLines(gzfile(file.path(SC, f)), n = 1, warn = FALSE)
  bc  <- strsplit(ln1, "\t")[[1]]
  d <- fread(file.path(SC, f), sep = "\t", header = FALSE, skip = 1, showProgress = FALSE)
  gn <- as.character(d[[1]]); m <- as.matrix(d[, -1])
  if (ncol(m) != length(bc)) stop(sprintf("Column mismatch: m=%d, barcodes=%d", ncol(m), length(bc)))
  colnames(m) <- bc
  storage.mode(m) <- "double"; rownames(m) <- gn
  m <- as(m, "CsparseMatrix")
  ok <- !is.na(gn) & nzchar(gn); m <- m[ok, , drop = FALSE]; gn <- gn[ok]
  if (anyDuplicated(gn)) {
    ug <- unique(gn)
    A <- sparseMatrix(i = match(gn, ug), j = seq_along(gn), x = 1,
                      dims = c(length(ug), length(gn)))
    m <- as(A %*% m, "CsparseMatrix"); rownames(m) <- ug
  } else rownames(m) <- gn
  colnames(m) <- paste0(nm, "_", colnames(m))
  o <- CreateSeuratObject(counts = m, project = nm, min.cells = 3, min.features = 200)
  o$sample <- nm; o$gsm <- gsm; o$tissue <- grp; o$label <- lab
  o[["percent.mt"]] <- PercentageFeatureSet(o, pattern = "^MT-")
  n0 <- ncol(o); o <- subset(o, subset = nFeature_RNA >= 200 & nFeature_RNA <= 6000 & percent.mt < 20)
  qc[[nm]] <- data.table(sample = nm, tissue = grp, n_before = n0, n_after = ncol(o))
  obj_list[[nm]] <- o
  cat(sprintf("  %-24s %-6s %6d -> %6d\n", nm, grp, n0, ncol(o)))
  rm(d, m); gc(verbose = FALSE)
}
Q <- rbindlist(qc); fwrite(Q, file.path(OUT, "sc193581_qc.csv"))
cat(sprintf("\nTotal %d -> After QC %d (%.0f%%)\n", sum(Q$n_before), sum(Q$n_after),
            100*sum(Q$n_after)/sum(Q$n_before)))
print(Q[, .(n=.N, before=sum(n_before), after=sum(n_after)), by = tissue])

obj <- merge(obj_list[[1]], obj_list[-1]); rm(obj_list); gc(verbose = FALSE)
cat("Merged:", ncol(obj), "cells\n")

obj <- NormalizeData(obj, verbose = FALSE) |> FindVariableFeatures(nfeatures = 2000, verbose = FALSE) |>
  ScaleData(verbose = FALSE) |> RunPCA(npcs = 30, verbose = FALSE) |>
  RunHarmony("sample", verbose = FALSE) |> RunUMAP(reduction = "harmony", dims = 1:30, verbose = FALSE) |>
  FindNeighbors(reduction = "harmony", dims = 1:30, verbose = FALSE) |>
  FindClusters(resolution = 0.5, verbose = FALSE)
cat("Clusters:", nlevels(obj$seurat_clusters), "\n")
tmp <- paste0(CACHE, ".tmp")
saveRDS(obj, tmp); file.rename(tmp, CACHE)
cat("[cached]", CACHE, "\n")
}

MARK <- list(
  `T/NK`      = c("CD3D","CD3E","CD2","TRAC","IL7R","NKG7","KLRD1"),
  `B/Plasma`  = c("MS4A1","CD79A","CD79B","IGKC","MZB1","JCHAIN","IGHG1"),
  Myeloid     = c("LYZ","CD68","CD14","AIF1","ITGAX","C1QA","FCN1"),
  Mast        = c("TPSAB1","TPSB2","MS4A2","CPA3","KIT"),
  Endothelial = c("PECAM1","VWF","EGFL7","PLVAP","CLDN5"),
  Fibroblast  = c("COL1A1","COL1A2","DCN","FAP","PDGFRA","LUM"),
  Thyrocyte   = c("TG","TPO","SLC5A5","TSHR","PAX8","FOXE1","SLC26A4","DIO1","DIO2"),
  Prolif      = c("MKI67","TOP2A","PCNA","CCNB1","BIRC5"),
  Squamous    = c("KRT5","KRT14","KRT6A","KRT17","S100A2")
)
mk_use <- intersect(unique(c(unlist(MARK), SIG4)), rownames(obj))
Idents(obj) <- "seurat_clusters"
av <- AverageExpression(obj, features = mk_use, assays = "RNA", layer = "data")$RNA
sc <- sapply(names(MARK), function(ct) {
  g <- intersect(MARK[[ct]], rownames(av)); if (!length(g)) return(rep(-Inf, ncol(av)))
  colMeans(av[g, , drop = FALSE]) })
rownames(sc) <- colnames(av)
ctv <- colnames(sc)[max.col(sc, ties.method = "first")]
names(ctv) <- sub("^g([0-9]+)$", "\\1", rownames(sc))   # Use a capturing group supported by the default R regular expressions.
lab_vec <- ctv[as.character(obj$seurat_clusters)]
if (anyNA(lab_vec)) {   # Fallback: try the cluster identifier with a g prefix.
  lab_vec <- ctv[paste0("g", obj$seurat_clusters)]
}
stopifnot(!anyNA(lab_vec))
obj$celltype <- unname(lab_vec)
cat("\nCell-type composition:\n"); print(table(obj$celltype))
fwrite(data.table(cluster = names(ctv), celltype = ctv), file.path(OUT, "sc193581_annotation.csv"))
saveRDS(obj, file.path(OUT, "seurat_GSE193581.rds"))
cat("[saved] seurat_GSE193581.rds\n")


obj <- tryCatch(JoinLayers(obj), error=function(e) obj)
cnt <- GetAssayData(obj, assay="RNA", layer="counts")
stopifnot(!anyNA(obj$celltype), ncol(cnt)==length(obj$celltype))
checkpoint <- file.path(OUT, "counts_annotated_193581.rds")
meta_index <- match(colnames(cnt),rownames(obj[[]]))
stopifnot(!anyNA(meta_index))
saveRDS(list(counts=cnt,celltype=obj$celltype[meta_index],tissue=obj$tissue[meta_index],sample=obj$sample[meta_index],cluster=as.character(obj$seurat_clusters[meta_index])),paste0(checkpoint,".tmp"))
if (!file.rename(paste0(checkpoint,".tmp"),checkpoint)) stop("Could not save annotated checkpoint.")

# Export the same coordinates and annotation used for the UMAP figure.
um <- Embeddings(obj,"umap")
um <- data.table(cell=rownames(um),UMAP1=um[,1],UMAP2=um[,2],celltype=obj$celltype,tissue=obj$tissue,sample=obj$sample)
fwrite(um,file.path(OUT,"umap_coords_current_annotations.csv"))
