# External PTC validation with deduplication
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
read_gpl570 <- function(fp) {
  hl <- grep("!series_matrix_table_begin", readLines(gzfile(fp), n=200))[1]
  M  <- fread(fp, skip=hl, header=TRUE, na.strings=c("","NA")); setnames(M,1,"probe")
  M[, probe := gsub('"', '', as.character(probe))]; M
}
getmeta <- function(fp) {
  tl <- readLines(gzfile(fp), n=200)
  g <- function(tag){ ln <- grep(paste0("^!",tag), tl, value=TRUE); if (!length(ln)) return(NULL)
    lapply(ln, function(l) strsplit(gsub('"','',sub(paste0("^!",tag,"\\t?"),"",l)), "\t")[[1]]) }
  list(gsm=g("Sample_geo_accession")[[1]], title=g("Sample_title")[[1]], char_lines=g("Sample_characteristics_ch1"))
}
pick <- function(lines, pat){ if (is.null(lines)) return(NULL)
  for (l in lines) if (any(grepl(pat, l, ignore.case=TRUE))) return(l); NULL }

load_ds <- function(fp, f_lab, f_pid) {
  M <- read_gpl570(fp); md <- getmeta(fp); sc <- setdiff(names(M), "probe")
  lab <- f_lab(md); names(lab) <- md$gsm; lab <- lab[sc]
  pid <- f_pid(md$title); names(pid) <- md$gsm; pid <- pid[sc]
  sub <- merge(A, M, by="probe")
  G <- sub[, lapply(.SD, function(x) mean(as.numeric(x), na.rm=TRUE)), by=symbol, .SDcols=sc]
  Gm <- as.matrix(G[, ..sc]); rownames(Gm) <- G$symbol
  list(G=Gm, gsm=sc, title=setNames(md$title, md$gsm)[sc], lab=lab, pid=pid)
}
pid_lead  <- function(t) ifelse(grepl("^ATC", t), t, sub("^([0-9]+).*$", "\\1", t))
pid_last  <- function(t) sapply(regmatches(t, gregexpr("[0-9]+", t)), function(v) if(length(v)) tail(v,1) else NA_character_)
pid_pair  <- function(t) { tk <- sub(".* ", "", t); tk <- sub("^N(PTC)", "\\1", tk); sub("[NT]$", "", tk) }

sec("[0] Load the four datasets")
g30 <- load_ds(file.path(GEO,"GSE33630_series_matrix.txt.gz"),
  function(md){ x <- gsub("^[^:]*: *","",pick(md$char_lines,"diagnostic"))
    ifelse(grepl("non-tumor|normal",x,ignore.case=TRUE),"Normal",
    ifelse(grepl("papillary",x,ignore.case=TRUE),"PTC",
    ifelse(grepl("anaplastic",x,ignore.case=TRUE),"ATC",NA))) }, pid_lead)
g55 <- load_ds(file.path(GEO,"GSE35570_matrix.txt.gz"),
  function(md){ x <- gsub("^[^:]*: *","",pick(md$char_lines,"type of thyroid"))
    ifelse(grepl("normal",x,ignore.case=TRUE),"Normal", ifelse(grepl("ptc",x,ignore.case=TRUE),"PTC",NA)) }, pid_last)
g29 <- load_ds(file.path(GEO,"GSE29265_matrix.txt.gz"),
  function(md){ x <- gsub("^[^:]*: *","",pick(md$char_lines,"^tissue"))
    ifelse(grepl("nontumor|non-tumor|normal",x,ignore.case=TRUE),"Normal",
    ifelse(grepl("papillary",x,ignore.case=TRUE),"PTC",
    ifelse(grepl("anaplastic",x,ignore.case=TRUE),"ATC",NA))) }, pid_pair)
g60 <- load_ds(file.path(GEO,"GSE60542_matrix.txt.gz"),
  function(md){ x <- md$title
    ifelse(grepl("Normal thyroid",x,ignore.case=TRUE),"Normal",
    ifelse(grepl("Papillary thyroid carcinoma",x,ignore.case=TRUE),"PTC",NA)) }, pid_lead)
for (nm in c("g30","g55","g29","g60")) { d <- get(nm)
  cat(sprintf("  %-5s %3d samples | %s\n", nm, length(d$gsm),
              paste(names(table(d$lab,useNA="ifany")), table(d$lab,useNA="ifany"), sep="=", collapse=" "))) }

sec("[1] Deduplication")
shared <- intersect(unique(na.omit(g30$pid)), unique(na.omit(g55$pid)))
cat(sprintf("  GSE35570 x GSE33630 shared patient IDs %d; %d records removed\n",
            length(shared), sum(g55$pid %in% shared)))
keep55 <- which(!(g55$pid %in% shared))
sel60 <- which(!is.na(g60$lab))
k60 <- paste(g60$pid[sel60], g60$lab[sel60], sep="_")
agg60 <- sapply(split(seq_along(k60), k60), function(ii) rowMeans(g60$G[, sel60[ii], drop=FALSE]))
lab60 <- sub("^[^_]*_", "", colnames(agg60)); pid60 <- sub("_.*$", "", colnames(agg60))
cat(sprintf("  GSE60542 after merging %d records (%d patients)\n", ncol(agg60), length(unique(pid60))))

mk <- function(G, lab, pid, tag) {
  keep <- !is.na(lab) & lab %in% c("PTC","Normal")
  list(G=G[, keep, drop=FALSE], lab=lab[keep], pid=pid[keep], tag=tag)
}
D30 <- mk(g30$G, g30$lab, g30$pid, "GSE33630")
D55 <- mk(g55$G[, keep55, drop=FALSE], g55$lab[keep55], g55$pid[keep55], "GSE35570 (deduplicated)")
D29 <- mk(g29$G, g29$lab, g29$pid, "GSE29265")
D60 <- mk(agg60, lab60, pid60, "GSE60542 (merged)")
for (d in list(D30,D55,D29,D60))
  cat(sprintf("  %-14s target samples %3d (PTC %d / normal %d, patients %d)\n", d$tag, length(d$lab),
              sum(d$lab=="PTC"), sum(d$lab=="Normal"), length(unique(d$pid))))

sec("[2] Score after z-standardisation within the target samples")
score_within <- function(d, w) {
  gs <- intersect(names(w), rownames(d$G))
  Z  <- t(scale(t(d$G[gs, , drop=FALSE])))
  as.numeric(t(Z) %*% w[gs])
}
boot_pid <- function(d, lp, n=2000) {
  set.seed(11); u <- unique(d$pid)
  a <- replicate(n, { ps <- sample(u, replace=TRUE)
    ii <- unlist(lapply(ps, function(p) which(d$pid == p)))
    if (length(unique(d$lab[ii])) < 2) return(NA_real_)
    as.numeric(auc(roc(d$lab[ii], lp[ii], levels=c("Normal","PTC"), direction="<", quiet=TRUE))) })
  a <- a[!is.na(a)]; c(quantile(a,.025), quantile(a,.975))
}

RES <- list()
for (d in list(D30,D55,D29,D60)) {
  lp <- score_within(d, COEF)
  ro <- roc(d$lab, lp, levels=c("Normal","PTC"), direction="<", quiet=TRUE); ci <- ci.auc(ro, method="delong")
  bi <- boot_pid(d, lp)
  cat(sprintf("  %-14s n=%3d  AUC=%.4f (DeLong %.3f-%.3f | patient bootstrap %.3f-%.3f)\n",
              d$tag, length(d$lab), as.numeric(ro$auc), ci[1], ci[3], bi[1], bi[2]))
  RES[[d$tag]] <- data.table(dataset=d$tag, n=length(d$lab), n_ptc=sum(d$lab=="PTC"),
                             n_normal=sum(d$lab=="Normal"), n_pid=length(unique(d$pid)),
                             auc=as.numeric(ro$auc), delong_lo=ci[1], delong_hi=ci[3],
                             boot_lo=bi[1], boot_hi=bi[2])
}

res <- rbindlist(RES)
print(res)
fwrite(res, file.path(DIR, "external_final_summary.csv"))
cat(sprintf("\n  output: %s\n", DIR))

saveRDS(list(D30,D55,D29,D60),file.path(DIR,"external_four_gene_inputs.rds"))
