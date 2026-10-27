# Figures 1 to 5
# Run through run_analysis.R; see docs/METHODS.md.

suppressMessages({
  library(data.table); library(ggplot2); library(pROC); library(patchwork)
  library(ggrepel); library(grid)
})
ROOT <- PROJECT_ROOT
GEO  <- file.path(ROOT, "data/geo")
DIR  <- RESULT_DIR
FIG  <- if (exists("FIGURE_DIR")) FIGURE_DIR else file.path(DIR, "figures"); dir.create(FIG, showWarnings=FALSE, recursive=TRUE)
sec  <- function(t){ cat("\n==",t,"==\n") }

SZ <- list(tag = 9, axis_title = 9, axis_text = 8, legend = 7, annot = 7, gene = 7)
mp <- function(pt) pt / 2.845            # ggplot geom_text size is in mm
LW  <- 0.6                                # data lines -> measured ~1.28 pt
LWG <- 0.5                                # grid/panel/ticks -> measured ~1.07 pt
DIM <- list(Fig1 = c(6.8, 3.2), Fig2 = c(6.8, 3.0), Fig3 = c(6.8, 5.0),
            Fig4 = c(6.8, 3.4), Fig5 = c(6.8, 3.1))

theme_set(
  theme_bw(base_size = SZ$axis_text) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = LWG, colour = "grey90"),
    panel.border    = element_rect(linewidth = LWG, colour = "grey30"),
    axis.ticks      = element_line(linewidth = LWG, colour = "grey30"),
    axis.title      = element_text(size = SZ$axis_title),
    axis.text       = element_text(size = SZ$axis_text, colour = "grey20"),
    legend.text     = element_text(size = SZ$legend),
    legend.title    = element_text(size = SZ$legend),
    legend.key.size = unit(0.32, "cm"),
    plot.tag        = element_text(size = SZ$tag, face = "bold"),
    plot.subtitle   = element_text(size = SZ$axis_text, hjust = 0, colour = "grey20"),
    strip.text      = element_text(size = SZ$axis_text)
  )
)
TAG <- theme(plot.tag = element_text(size = SZ$tag, face = "bold"))

save_fig <- function(p, i) {
  if (exists("BUILD_PLOTS_ONLY") && isTRUE(BUILD_PLOTS_ONLY)) return(invisible(NULL))
  dim <- DIM[[paste0("Fig", i)]]
  ggsave(file.path(FIG, sprintf("Fig%d.png", i)), p, width=dim[1], height=dim[2], dpi=300)
  ggsave(file.path(FIG, sprintf("Fig%d.tif", i)), p, width=dim[1], height=dim[2], dpi=300,
         device="tiff", compression="lzw")
  ggsave(file.path(FIG, sprintf("Fig%d.pdf", i)), p, width=dim[1], height=dim[2])
  cat(sprintf("  Fig%d written %.1f x %.1f in\n", i, dim[1], dim[2]))
}

set.seed(7)

SIG   <- fread(file.path(DIR,"signature_genes.csv"))$gene
model <- readRDS(file.path(DIR,"signature_model.rds"))
CALL  <- coef(summary(model$model))
INT   <- CALL[1,1]; COEF <- CALL[-1,1]; names(COEF) <- SIG

sec("Fig 1")
L <- fread(file.path(DIR,"sel_lasso.csv"))$gene
R <- fread(file.path(DIR,"sel_rf.csv"))$gene
S <- fread(file.path(DIR,"sel_svm.csv"))$gene
reg <- c(
  onlyL=length(setdiff(L, union(R,S))), onlyR=length(setdiff(R, union(L,S))),
  onlyS=length(setdiff(S, union(L,R))),
  LR=length(setdiff(intersect(L,R), S)), LS=length(setdiff(intersect(L,S), R)),
  RS=length(setdiff(intersect(R,S), L)), LRS=length(Reduce(intersect, list(L,R,S))))

circ <- function(x,y,r,n=200) data.frame(x=x+r*cos(seq(0,2*pi,length.out=n)),
                                         y=y+r*sin(seq(0,2*pi,length.out=n)))
c1 <- circ(0.37,0.60,0.27); c2 <- circ(0.63,0.60,0.27); c3 <- circ(0.50,0.38,0.27)
lab <- data.frame(x=c(0.19,0.81,0.50,0.50,0.29,0.71,0.50),
                  y=c(0.70,0.70,0.26,0.735,0.44,0.44,0.585),
                  n=c(reg["onlyL"],reg["onlyR"],reg["onlyS"],reg["LR"],reg["LS"],reg["RS"],reg["LRS"]))
f1a <- ggplot() +
  geom_path(data=c1,aes(x,y),colour="#2C7FB8",linewidth=LW) +
  geom_path(data=c2,aes(x,y),colour="#D95F02",linewidth=LW) +
  geom_path(data=c3,aes(x,y),colour="#1B9E77",linewidth=LW) +
  geom_text(data=lab,aes(x,y,label=n),size=mp(SZ$annot)) +
  annotate("text",x=0.17,y=0.90,label="LASSO (12)",colour="#2C7FB8",size=mp(SZ$annot),hjust=0) +
  annotate("text",x=0.83,y=0.90,label="Random forest (10)",colour="#D95F02",size=mp(SZ$annot),hjust=1) +
  annotate("text",x=0.50,y=0.055,label="SVM-RFE (5)",colour="#1B9E77",size=mp(SZ$annot)) +
  annotate("text",x=0.50,y=-0.035,label="ABHD17C  ZDHHC6  ABHD17B  ZDHHC20",
           size=mp(SZ$gene), fontface="italic") +
  coord_equal(xlim=c(0.06,0.94),ylim=c(-0.09,0.95),clip="off") + theme_void()

coefficient_intervals <- fread(file.path(DIR,"signature_coefficients.csv"))
cf <- data.frame(gene=factor(SIG,levels=rev(SIG)),est=COEF[SIG],
 lo=coefficient_intervals$ci_lo[match(SIG,coefficient_intervals$term)],
 hi=coefficient_intervals$ci_hi[match(SIG,coefficient_intervals$term)])

f1b <- ggplot(cf,aes(est,gene)) +
  geom_vline(xintercept=0,linetype=2,colour="grey50",linewidth=LWG) +
  geom_errorbarh(aes(xmin=lo,xmax=hi),height=0.15,colour="grey30",linewidth=LW) +
  geom_point(size=2.6,colour="#2C7FB8") +
  labs(x="Logistic regression coefficient", y=NULL) +
  theme(axis.text.y=element_text(face="bold.italic"))
p1 <- f1a + f1b + plot_layout(widths=c(1,1)) +
  plot_annotation(tag_levels="a") & TAG
save_fig(p1, 1)

sec("Fig 2")
oof <- fread(file.path(DIR,"nested_cv_oof_predictions.csv"))
rep_auc <- sapply(paste0("pred_r",1:5), function(c)
  as.numeric(auc(roc(oof$group, oof[[c]], levels=c("Normal","Tumor"), direction="<", quiet=TRUE))))
f2a <- ggplot(data.frame(rep=factor(1:5), auc=rep_auc), aes(rep,auc)) +
  geom_hline(yintercept=mean(rep_auc),linetype=2,colour="red",linewidth=LWG) +
  geom_point(size=3.2,colour="#2C7FB8") +
  annotate("text",x=5,y=mean(rep_auc),label=sprintf("mean %.3f",mean(rep_auc)),
           hjust=1.05,vjust=-0.8,size=mp(SZ$annot),colour="red") +
  coord_cartesian(ylim=c(0.945,0.975)) +
  labs(x="Repeat", y="AUC (pooled out-of-fold)")

fr <- fread(file.path(DIR,"nested_cv_gene_selection_freq.csv"))
fr[, sig := gene %in% SIG]; setorder(fr, -n_folds)
fr[, gene := factor(gene, levels=rev(as.character(gene)))]
f2b <- ggplot(fr, aes(n_folds, gene, fill=sig)) +
  geom_col(width=0.7) + scale_fill_manual(values=c(`FALSE`="grey75",`TRUE`="#C0392B"), guide="none") +
  labs(x="Selected in N of 50 training folds", y=NULL) +
  geom_text(aes(label=n_folds), hjust=-0.25, size=mp(SZ$annot)) +
  theme(axis.text.y=element_text(face="italic"))
p2 <- f2a + f2b + plot_layout(widths=c(1,1.15)) +
  plot_annotation(tag_levels="a") & TAG
save_fig(p2, 2)

sec("Fig 3")
E  <- readRDS(file.path(MATRIX_DIR,"expr_full_logcpm.rds"))
dg <- fread(file.path(DIR,"deg_full_palmitoylation_genes.csv"))
gn <- dg[order(logFC)]$gene
E29 <- E[gn, colnames(E) %in% oof$sample, drop=FALSE]
grp <- oof$group[match(colnames(E29), oof$sample)]
stopifnot(!any(is.na(grp)))
ord <- order(grp, E29[1,])
stopifnot(all(grp[ord] == sort(grp)))
Eo <- E29[, ord]; grpo <- grp[ord]
Z <- t(scale(t(Eo)))
hm <- data.table(gene=rep(gn, ncol(Z)), x=rep(seq_len(ncol(Z)), each=nrow(Z)),
                 z=as.numeric(Z), group=rep(grpo, each=nrow(Z)))
hm[, gene := factor(gene, levels=gn)]
nn <- sum(grpo=="Normal"); nT <- ncol(Z) - nn
ann_bar <- ggplot(unique(hm[, .(x, group)]), aes(x, 1, fill=group)) + geom_raster() +
  geom_vline(xintercept=nn+0.5, colour="white", linewidth=LWG) +
  scale_fill_manual(values=c(Normal="#1B9E77", Tumor="#C0392B"), name=NULL) +
  scale_x_continuous(breaks=c((nn+1)/2, nn+(nT+1)/2), labels=c("Normal","Tumour")) +
  theme_void() + theme(legend.position="none",
                       axis.text.x=element_text(size=SZ$gene, colour="grey20"))
f3a_hm <- ggplot(hm, aes(x, gene, fill=z)) + geom_raster() +
  geom_vline(xintercept=nn+0.5, colour="black", linewidth=LWG) +
  scale_fill_gradient2(low="#2166AC",mid="white",high="#B2182B",midpoint=0,
                       name="row z", limits=c(-3,3), oob=scales::squish) +
labs(x=sprintf("normal (n = %d) | tumour (n = %d)", nn, nT), y=NULL) +   # panel label showing n
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text.y=element_text(size=SZ$gene, face="italic"))
f3a <- (ann_bar / f3a_hm + plot_layout(heights=c(0.06, 1)))
vol <- copy(dg); vol[, sig := gene %in% SIG]
vol[, lab := fifelse(gene %in% SIG, gene, fifelse(abs(logFC)>0.5, gene, NA_character_))]
f3b <- ggplot(vol, aes(logFC, -log10(adj.P.Val))) +
  geom_point(aes(colour=sig),size=1.8) +
  geom_vline(xintercept=c(-0.5,0.5),linetype=3,colour="grey60",linewidth=LW) +
  scale_colour_manual(values=c(`FALSE`="grey55",`TRUE`="#C0392B"), guide="none") +
  geom_text_repel(data=vol[!is.na(lab)], aes(label=lab), size=mp(SZ$annot),
                  max.overlaps=30, segment.size=0.3, fontface="italic") +
  labs(x="log2 fold change (tumour vs normal)", y="-log10 FDR")
p3 <- (f3a | f3b) + plot_layout(widths=c(1.25,1)) +
  plot_annotation(tag_levels=list(c("a","","b"))) & TAG
save_fig(p3, 3)

sec("Fig 4 load external cohorts")
ann <- fread(file.path(ROOT,"data/geo/GPL570.annot.gz"), skip=27, header=TRUE, sep="\t",
             quote="", na.strings=c("","NA"), fill=TRUE); setnames(ann,1,"probe")
A <- ann[!is.na(`Gene symbol`) & `Gene symbol` %in% SIG,
         .(probe=as.character(probe), symbol=as.character(`Gene symbol`))]
read_gpl570 <- function(fp) {
  hl <- grep("!series_matrix_table_begin", readLines(gzfile(fp), n=200))[1]
  M  <- fread(fp, skip=hl, header=TRUE, na.strings=c("","NA")); setnames(M,1,"probe")
  M[, probe := gsub('"','',as.character(probe))]; M }
getmeta <- function(fp) {
  tl <- readLines(gzfile(fp), n=200)
  g <- function(tag){ ln <- grep(paste0("^!",tag), tl, value=TRUE); if(!length(ln)) return(NULL)
    lapply(ln, function(l) strsplit(gsub('"','',sub(paste0("^!",tag,"\\t?"),"",l)),"\t")[[1]]) }
  list(gsm=g("Sample_geo_accession")[[1]], title=g("Sample_title")[[1]], char_lines=g("Sample_characteristics_ch1")) }
pick <- function(lines, pat){ if(is.null(lines)) return(NULL)
  for (l in lines) if (any(grepl(pat,l,ignore.case=TRUE))) return(l); NULL }
load_ds <- function(fp, f_lab, f_pid) {
  M <- read_gpl570(fp); md <- getmeta(fp); sc <- setdiff(names(M),"probe")
  lab <- f_lab(md); names(lab) <- md$gsm; lab <- lab[sc]
  pid <- f_pid(md$title); names(pid) <- md$gsm; pid <- pid[sc]
  sub <- merge(A, M, by="probe")
  G <- sub[, lapply(.SD, function(x) mean(as.numeric(x),na.rm=TRUE)), by=symbol, .SDcols=sc]
  Gm <- as.matrix(G[, ..sc]); rownames(Gm) <- G$symbol
  list(G=Gm, gsm=sc, lab=lab, pid=pid) }
pid_lead <- function(t) ifelse(grepl("^ATC",t), t, sub("^([0-9]+).*$","\\1",t))
pid_last <- function(t) sapply(regmatches(t,gregexpr("[0-9]+",t)), function(v) if(length(v)) tail(v,1) else NA_character_)
pid_pair <- function(t){ tk <- sub(".* ","",t); tk <- sub("^N(PTC)","\\1",tk); sub("[NT]$","",tk) }

g30 <- load_ds(file.path(GEO,"GSE33630_series_matrix.txt.gz"),
  function(md){ x <- gsub("^[^:]*: *","",pick(md$char_lines,"diagnostic"))
    ifelse(grepl("non-tumor|normal",x,ignore.case=TRUE),"Normal",
    ifelse(grepl("papillary",x,ignore.case=TRUE),"PTC", ifelse(grepl("anaplastic",x,ignore.case=TRUE),"ATC",NA))) }, pid_lead)
g55 <- load_ds(file.path(GEO,"GSE35570_matrix.txt.gz"),
  function(md){ x <- gsub("^[^:]*: *","",pick(md$char_lines,"type of thyroid"))
    ifelse(grepl("normal",x,ignore.case=TRUE),"Normal", ifelse(grepl("ptc",x,ignore.case=TRUE),"PTC",NA)) }, pid_last)
g29 <- load_ds(file.path(GEO,"GSE29265_matrix.txt.gz"),
  function(md){ x <- gsub("^[^:]*: *","",pick(md$char_lines,"^tissue"))
    ifelse(grepl("nontumor|non-tumor|normal",x,ignore.case=TRUE),"Normal",
    ifelse(grepl("papillary",x,ignore.case=TRUE),"PTC", ifelse(grepl("anaplastic",x,ignore.case=TRUE),"ATC",NA))) }, pid_pair)
g60 <- load_ds(file.path(GEO,"GSE60542_matrix.txt.gz"),
  function(md){ x <- md$title
    ifelse(grepl("Normal thyroid",x,ignore.case=TRUE),"Normal",
    ifelse(grepl("Papillary thyroid carcinoma",x,ignore.case=TRUE),"PTC",NA)) }, pid_lead)

shared <- intersect(unique(na.omit(g30$pid)), unique(na.omit(g55$pid)))
keep55 <- which(!(g55$pid %in% shared))
sel60  <- which(!is.na(g60$lab)); k60 <- paste(g60$pid[sel60], g60$lab[sel60], sep="_")
agg60  <- sapply(split(seq_along(k60), k60), function(ii) rowMeans(g60$G[, sel60[ii], drop=FALSE]))
lab60  <- sub("^[^_]*_","",colnames(agg60)); pid60 <- sub("_.*$","",colnames(agg60))
mk <- function(G,lab,pid,tag){ keep <- !is.na(lab) & lab %in% c("PTC","Normal")
  list(G=G[,keep,drop=FALSE], lab=lab[keep], pid=pid[keep], tag=tag) }
DS <- list(mk(g30$G,g30$lab,g30$pid,"GSE33630"),
           mk(g55$G[,keep55,drop=FALSE],g55$lab[keep55],g55$pid[keep55],"GSE35570 (deduplicated)"),
           mk(g29$G,g29$lab,g29$pid,"GSE29265"),
           mk(agg60,lab60,pid60,"GSE60542 (merged)"))
sc_within <- function(d,w){ gs <- intersect(names(w), rownames(d$G))
  Z <- t(scale(t(d$G[gs,,drop=FALSE]))); as.numeric(t(Z) %*% w[gs]) }

rocdf <- list(); aucs <- c()
for (d in DS) {
  lp <- sc_within(d, COEF)
  ro <- roc(d$lab, lp, levels=c("Normal","PTC"), direction="<", quiet=TRUE)
  aucs[d$tag] <- as.numeric(ro$auc)
  rocdf[[d$tag]] <- data.table(tag=d$tag, sens=rev(ro$sensitivities), spec=rev(ro$specificities))
}
rd <- rbindlist(rocdf); rd[, tag := factor(tag, levels=names(aucs))]
leg <- sprintf("%s: AUC %.3f", names(aucs), aucs)
f4a <- ggplot(rd, aes(1-spec, sens, colour=tag)) + geom_line(linewidth=LW) +
  geom_abline(slope=1,intercept=0,linetype=2,colour="grey60",linewidth=LWG) +
  scale_colour_brewer(palette="Dark2", labels=leg, name=NULL) +
  labs(x="1 - specificity", y="Sensitivity") +
  coord_equal() + theme(legend.position=c(0.98,0.02), legend.justification=c(1,0),
                        legend.background=element_rect(fill=alpha("white",0.7)),
                        legend.text=element_text(size=SZ$legend),
                        legend.key.size=unit(0.3,"cm"))

sec("Fig 4b GSE82208")
g82 <- load_ds(file.path(GEO,"GSE82208_series_matrix.txt.gz"),
  function(md){ x <- gsub("^[^:]*: *","",pick(md$char_lines,"^class:"))
    ifelse(grepl("FTA|adenoma",x,ignore.case=TRUE),"FTA",
    ifelse(grepl("FTC|carcinoma|cancer",x,ignore.case=TRUE),"FTC",NA)) }, pid_lead)
k82 <- !is.na(g82$lab)
G82 <- g82$G[,k82,drop=FALSE]; lab82 <- g82$lab[k82]
gs82 <- intersect(SIG, rownames(G82))
Z82 <- t(scale(t(G82[gs82,,drop=FALSE]))); lp82 <- as.numeric(t(Z82) %*% COEF[gs82])
ro82 <- roc(lab82, lp82, levels=c("FTA","FTC"), direction="<", quiet=TRUE)
cat(sprintf("  GSE82208 n=%d (%d/%d) AUC=%.4f\n", length(lab82), sum(lab82=="FTC"), sum(lab82=="FTA"), as.numeric(ro82$auc)))
f4b <- ggplot(data.frame(lp=lp82, grp=lab82), aes(grp, lp, fill=grp)) +
  geom_boxplot(width=0.55, outlier.shape=NA, alpha=0.85, linewidth=LW) +
  geom_jitter(width=0.12, size=1.2, alpha=0.75) +
  scale_fill_manual(values=c(FTA="#7FB3D5",FTC="#C0392B"), guide="none") +
  scale_x_discrete(labels=c(FTA="Adenoma", FTC="Carcinoma")) +
  labs(x=NULL, y="Signature score", subtitle="Adenoma vs carcinoma")
f4c <- ggplot(data.frame(sens=rev(ro82$sensitivities), spec=rev(ro82$specificities)), aes(1-spec,sens)) +
  geom_line(linewidth=LW,colour="#C0392B") + geom_abline(slope=1,intercept=0,linetype=2,colour="grey60",linewidth=LWG) +
  labs(x="1 - specificity", y="Sensitivity", subtitle=sprintf("AUC %.3f",as.numeric(ro82$auc))) +
  coord_equal()
p4 <- f4a + (f4b | f4c) + plot_layout(widths=c(1,1.1)) +
  plot_annotation(tag_levels="a") & TAG
save_fig(p4, 4)

sec("Fig 5")
cl <- fread(file.path(DIR,"clinical_immune_fixed.csv"), na.strings=c("","NA"))
cl <- cl[!is.na(risk)]
mkbox <- function(var, lab0, lab1) {
  k <- !is.na(cl[[var]]); d <- cl[k]
  data.table(g=ifelse(d[[var]]==1, lab1, lab0), v=d$risk, var=lab1)
}
b <- rbind(mkbox("BRAF_V600E","BRAF WT","BRAF V600E"),
           mkbox("RAS_MUT","RAS WT","RAS mutant"))
b[, var := factor(var, levels=c("BRAF V600E","RAS mutant"))]
b[, g := factor(g, levels=c("BRAF WT","BRAF V600E","RAS WT","RAS mutant"))]
nn_lab <- b[, .(n=.N), by=.(var,g)][, lbl := sprintf("%s\n(n = %d)", g, n)]   # panel label showing n
b <- merge(b, nn_lab[, .(var,g,lbl)], by=c("var","g"))
b[, lbl := factor(lbl, levels=unique(nn_lab$lbl))]
f5a <- ggplot(b, aes(lbl, v, fill=g)) +
  facet_wrap(~var, scales="free_x") +
  geom_boxplot(width=0.55, outlier.shape=NA, alpha=0.85, linewidth=LW) +
  geom_jitter(width=0.12, size=0.6, alpha=0.5) +
  scale_fill_manual(values=c("BRAF WT"="grey70","BRAF V600E"="#C0392B",
                             "RAS WT"="grey70","RAS mutant"="#2C7FB8"), guide="none") +
  labs(x=NULL, y="Signature score") +
  theme(axis.text.x=element_text(size=SZ$gene),
        strip.background=element_rect(fill="grey95", colour="grey30", linewidth=LWG))
cmp <- data.table(score=c("STROMA","IMMUNE","THYRO","COMPOS"),
                  r=sapply(c("STROMA","IMMUNE","THYRO","COMPOS"), function(v) cor(cl$risk, cl[[v]], use="complete.obs")))
cmp[, score := factor(score, levels=c("STROMA","IMMUNE","THYRO","COMPOS"))]
f5b <- ggplot(cmp, aes(score, r, fill=r>0)) + geom_col(width=0.62) +
  geom_hline(yintercept=0, linewidth=LWG) +
  geom_text(aes(label=sprintf("%+.3f",r)), vjust=ifelse(cmp$r>0,-0.5,1.3), size=mp(SZ$annot)) +
  scale_fill_manual(values=c(`TRUE`="#4C9AC4",`FALSE`="#E08214"), guide="none") +
  labs(x=NULL, y="Pearson r with signature score") + coord_cartesian(ylim=c(-0.55,0.55))
p5 <- f5a + f5b + plot_layout(widths=c(1.25,1)) +
  plot_annotation(tag_levels="a") & TAG
save_fig(p5, 5)

cat("\n== done, output directory:", FIG, "==\n")
print(list.files(FIG))
