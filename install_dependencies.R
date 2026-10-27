# Run explicitly to install missing dependencies into your normal R library.
# This uses currently available releases; recorded versions are not automatically pinned.
cran <- c("data.table","glmnet","randomForest","e1071","pROC","httr","jsonlite",
          "msigdbr","ggplot2","patchwork","ggrepel","ggtext","scales","readxl","statmod","Seurat","harmony","Matrix","survival","sandwich","png")
bioc <- c("limma","GSVA","clusterProfiler","org.Hs.eg.db","edgeR","scDblFinder","SingleCellExperiment")
missing <- cran[!vapply(cran,requireNamespace,logical(1),quietly=TRUE)]
if (length(missing)) install.packages(missing,repos="https://cloud.r-project.org")
missing <- bioc[!vapply(bioc,requireNamespace,logical(1),quietly=TRUE)]
if (length(missing)) {
  if (!requireNamespace("BiocManager",quietly=TRUE)) install.packages("BiocManager",repos="https://cloud.r-project.org")
  BiocManager::install(missing,update=FALSE,ask=FALSE)
}

# ESTIMATE is distributed through R-Forge; this is an explicit installation step.
if (!requireNamespace("estimate",quietly=TRUE)) install.packages("estimate",repos="https://R-Forge.R-project.org")
