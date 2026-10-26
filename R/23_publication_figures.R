# Final publication figures. Run after steps 16 and 22, or source with the
# documented PROJECT_ROOT, MATRIX_DIR, RESULT_DIR and OUTPUT_DIR variables.
# Frozen supplementary coefficients reproduce the reported adjustment results.
BUILD_PLOTS_ONLY <- TRUE
FIGURE_DIR <- file.path(OUTPUT_DIR, "figures")
source(file.path(SOURCE_ROOT, "R", "16_figures.R"), local=TRUE, encoding="UTF-8")
FIG <- FIGURE_DIR
TAB <- file.path(SOURCE_ROOT, "data", "figure_tables")
theme_set(theme_bw(base_size=9) + theme(panel.grid.minor=element_blank(),
  panel.grid.major=element_line(colour="grey93", linewidth=.3),
  axis.text=element_text(colour="grey20"), plot.tag=element_text(face="bold",size=11),
  legend.position="bottom", strip.background=element_rect(fill="grey96")))
cols <- c(Normal="#0072B2", Tumor="#D55E00", PTC="#D55E00")
export <- function(p,name,h) {
  ggsave(file.path(FIG,paste0(name,".png")),p,width=6.8,height=h,dpi=300,bg="white")
  ggsave(file.path(FIG,paste0(name,".tif")),p,width=6.8,height=h,dpi=600,compression="lzw",bg="white")
  ggsave(file.path(FIG,paste0(name,".pdf")),p,width=6.8,height=h,device=cairo_pdf,bg="white")
}
box <- function(dat,ylab) ggplot(dat,aes(group,value,fill=group)) +
  geom_boxplot(width=.55,outlier.shape=NA,linewidth=.4) +
  geom_point(position=position_jitter(width=.13,seed=11),size=.45,alpha=.28) +
  scale_fill_manual(values=cols,guide="none") + labs(x=NULL,y=ylab)
expr <- rbindlist(lapply(SIG,function(g) data.table(gene=g,group=grp,value=as.numeric(E[g,colnames(E29)]))))
expr[,group:=factor(group,levels=c("Normal","Tumor"),labels=c("Normal","PTC"))]
px <- box(expr,"log2(CPM + 1)") + facet_wrap(~gene,nrow=1)
score <- fread(file.path(DIR,"signature_scores_tcga.csv"))
score <- score[sample %in% oof$sample];score[,group:=factor(group,levels=c("Normal","Tumor"),labels=c("Normal","PTC"))]
ps <- box(score[,.(group,value=risk)],"Signature score")
# Candidate-level results precede introduction of the selected signature.
vol[,lab:=fifelse(abs(logFC)>0.5,gene,NA_character_)]
f3b <- ggplot(vol,aes(logFC,-log10(adj.P.Val))) +
 geom_point(colour="#0072B2",size=1.8) +
 geom_vline(xintercept=c(-.5,.5),linetype=3,colour="grey60",linewidth=.6) +
 geom_text_repel(data=vol[!is.na(lab)],aes(label=lab),size=7/2.845,
                 max.overlaps=30,segment.size=.3,seed=7) +
 labs(x="log2 fold change (tumour vs normal)",y="-log10 FDR")
fig1 <- (wrap_elements(full=f3a) | f3b) + plot_annotation(tag_levels="a")
export(fig1,"Fig1",5.0)

enrich <- function(file,n,title) {
 d <- fread(file.path(DIR,file));d <- head(d[order(p.adjust)],n)
 d[,ratio:=sapply(strsplit(GeneRatio,"/"),function(v)as.numeric(v[1])/as.numeric(v[2]))]
 d[,Description:=factor(Description,levels=rev(Description))]
 ggplot(d,aes(ratio,Description,size=Count,colour=-log10(p.adjust))) + geom_point() +
 scale_colour_viridis_c(name="-log10(FDR)",breaks=scales::breaks_pretty(n=3),guide=guide_colourbar(barwidth=unit(2.2,"cm"))) + scale_size_continuous(range=c(2,5),name="Genes") +
 scale_y_discrete(labels=function(x)vapply(x,function(y)paste(strwrap(y,44),collapse="\n"),"")) + labs(x="Input gene ratio",y=NULL,title=title)+theme(axis.text.y=element_text(size=7),legend.text=element_text(size=7),legend.title=element_text(size=7),plot.title=element_text(size=10),legend.box="vertical")
}
go <- enrich("enrich_GO_BP_clusterProfiler.csv",8,"GO biological processes")
kg <- enrich("enrich_KEGG_clusterProfiler.csv",5,"KEGG pathways")
ed <- fread(file.path(DIR,"ppi_edges_string.csv"));nd <- sort(unique(c(ed$a,ed$b)))
# Equal spacing around the complete ring, with one node at the top.
angle <- pi/2 + (seq_along(nd)-6)*2*pi/length(nd)
node <- data.table(gene=nd,x=cos(angle),y=sin(angle))
node[,degree:=sapply(gene,function(g)sum(ed$a==g|ed$b==g))];node[,signature:=gene %in% SIG]
node[,cap:=y>0 & abs(x)<.16]
node[,label_x:=ifelse(cap,x,x+sign(x)*.085)]
node[,label_y:=ifelse(cap,y+sign(y)*.085,y)]
node[gene=="PPT2",label_y:=label_y+.035]
node[,label_hjust:=ifelse(cap,.5,ifelse(x>=0,0,1))]
node[,label_vjust:=ifelse(cap,ifelse(y>0,0,1),.5)]
node[y < -.97, `:=`(label_y=y-.07,label_vjust=1)]
# All names remain horizontal and adjacent to their nodes, without leaders.
edge <- merge(ed,node[,.(a=gene,x,y)],by="a");edge <- merge(edge,node[,.(b=gene,xend=x,yend=y)],by="b")
net <- ggplot()+geom_curve(data=edge,aes(x,y,xend=xend,yend=yend),curvature=.045,linewidth=.25,alpha=.15)+
 geom_point(data=node,aes(x,y,size=degree),colour="#0072B2")+
 geom_text(data=node,aes(label_x,label_y,label=gene,hjust=label_hjust,vjust=label_vjust),size=2.4)+
 scale_size_continuous(range=c(2,4.2),guide="none")+coord_equal(xlim=c(-1.7,1.7),ylim=c(-1.22,1.22))+theme_void(base_size=9)+theme(legend.position="bottom",legend.text=element_text(size=7),plot.title=element_text(size=10))+labs(title="STRING associations\n25 nodes; 200 edges")
export((wrap_elements(full=go) / (wrap_elements(full=kg) | wrap_elements(full=net)))+plot_layout(heights=c(1.3,1.2))+plot_annotation(tag_levels="a"),"Fig2",7.4)
f2b <- f2b + scale_x_continuous(limits=c(0,55),breaks=seq(0,50,10))
# Panel order follows selection, model construction, expression, then evaluation.
f1a <- ggplot() +
 geom_path(data=c1,aes(x,y),colour="#2C7FB8",linewidth=.6) +
 geom_path(data=c2,aes(x,y),colour="#D95F02",linewidth=.6) +
 geom_path(data=c3,aes(x,y),colour="#1B9E77",linewidth=.6) +
 geom_text(data=lab,aes(x,y,label=n),size=7/2.845) +
 annotate("text",x=.26,y=.96,label="LASSO\n(12)",colour="#2C7FB8",size=7/2.845) +
 annotate("text",x=.75,y=.96,label="Random forest\n(10)",colour="#D95F02",size=7/2.845) +
 annotate("text",x=.5,y=.045,label="SVM-RFE (5)",colour="#1B9E77",size=7/2.845) +
 annotate("text",x=.5,y=-.075,label="ABHD17C  ZDHHC6\nABHD17B  ZDHHC20",size=7/2.845) +
 coord_equal(xlim=c(.06,.94),ylim=c(-.16,1.07),clip="off") + theme_void()
fig3 <- (f1a | f1b) / ((px | ps)+plot_layout(widths=c(3,1))) /
 (f2a | f2b) + plot_layout(heights=c(1,.8,1.1)) + plot_annotation(tag_levels="a")
export(fig3,"Fig3",8.0)

ex <- fread(file.path(DIR,"external_final_summary.csv"));ex[,dataset:=c("GSE33630","GSE35570","GSE29265","GSE60542")]
stopifnot(max(abs(as.numeric(aucs)-ex$auc))<1e-7)
forest <- ggplot(ex,aes(auc,factor(dataset,levels=rev(dataset))))+
 geom_segment(aes(x=boot_lo,xend=boot_hi,yend=factor(dataset,levels=rev(dataset))),linewidth=2,colour="grey75")+
 geom_segment(aes(x=delong_lo,xend=delong_hi,yend=factor(dataset,levels=rev(dataset))),linewidth=.65,colour="#0072B2")+
 geom_point(size=2.5,colour="#0072B2")+geom_text(aes(label=sprintf("%.4f",auc)),nudge_y=.25,size=2.8)+
 coord_cartesian(xlim=c(.80,1.015))+labs(x="AUC with 95% confidence intervals",y=NULL)
# Expression distributions use the same deduplicated/aggregated data as ROC analysis.
external_expression <- rbindlist(lapply(DS,function(d) {
  rbindlist(lapply(SIG,function(g) data.table(dataset=sub(" .*","",d$tag),
    gene=g,group=d$lab,value=as.numeric(d$G[g,]))))
}))
external_expression[,dataset:=factor(dataset,levels=c("GSE33630","GSE35570","GSE29265","GSE60542"))]
external_expression[,gene:=factor(gene,levels=SIG)]
external_expression[,group:=factor(group,levels=c("Normal","PTC"))]
stopifnot(nrow(external_expression)==948,all(is.finite(external_expression$value)))
# Check all 32 group medians against the reported supplementary table.
gd <- fread(file.path(DIR,"per_gene_external_direction.csv"))
gd[,dataset:=sub(" .*","",dataset)]
observed <- dcast(external_expression[,.(median=median(value)),by=.(dataset,gene,group)],
                  dataset+gene~group,value.var="median")
check <- merge(gd,observed,by=c("dataset","gene"))
stopifnot(nrow(check)==16,max(abs(check$Normal-check$median_Normal))<1e-10,
          max(abs(check$PTC-check$median_PTC))<1e-10)
# Narrow outlines and paired offsets leave clear gaps between cohorts.
external_expression[,cohort_x:=as.numeric(dataset)]
external_expression[,plot_x:=cohort_x+ifelse(group=="Normal",-.19,.19)]
summaries <- external_expression[,.(median=median(value),
  q1=quantile(value,.25),q3=quantile(value,.75)),by=.(dataset,gene,group,plot_x)]
expression_plot <- ggplot(external_expression,
    aes(plot_x,value,group=interaction(dataset,group),colour=group)) +
 geom_point(position=position_jitter(width=.095,height=0,seed=11),
            shape=16,size=.65,alpha=.45,show.legend=FALSE) +
 geom_errorbar(data=summaries,aes(x=plot_x,y=NULL,ymin=q1,ymax=q3),
               width=.13,linewidth=.4) +
 geom_segment(data=summaries,aes(x=plot_x-.13,xend=plot_x+.13,
              y=median,yend=median),linewidth=.65) +
 facet_wrap(~gene,ncol=2,scales="free_y") +
 scale_colour_manual(values=c(Normal="#0072B2",PTC="#D55E00"),name=NULL) +
 scale_x_continuous(breaks=1:4,labels=levels(external_expression$dataset),
                    expand=expansion(add=.5)) +
 scale_y_continuous(breaks=scales::breaks_pretty(n=4),expand=expansion(mult=c(.08,.12))) +
 labs(x=NULL,y="log2 expression") +
 theme(panel.grid=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.border=element_blank(),
       axis.line=element_line(colour="grey25",linewidth=.4),
       axis.ticks=element_line(colour="grey25",linewidth=.35),
       strip.background=element_blank(),strip.text=element_text(face="bold.italic",size=10),
       axis.text.x=element_text(size=7),axis.text.y=element_text(size=8),
       axis.title.y=element_text(size=9),legend.position="bottom",
       legend.key.width=unit(.65,"cm"),legend.key.height=unit(.32,"cm"),
       panel.spacing.x=unit(.65,"cm"),panel.spacing.y=unit(.6,"cm"))

export((f4a | forest) / expression_plot + plot_layout(heights=c(1,1.65)) +
       plot_annotation(tag_levels="a"),"Fig4",8.0)
fwrite(external_expression[,.(dataset,gene,group,value)],file.path(FIG,"external_gene_distribution_plot_data.csv"))

pub <- fread(file.path(TAB,"tcga_published_scores.csv"));pub[,key:=substr(sample,1,15)]
ss <- copy(score);ss[,key:=substr(sample,1,15)];pub <- merge(pub,ss[,.(key,risk)],by="key")
scatter <- function(v,label) {
 d <- pub[is.finite(get(v)) & is.finite(risk)];r <- cor(d[[v]],d$risk)
 ggplot(d,aes(x=.data[[v]],y=risk))+geom_point(size=.6,alpha=.35,colour="#0072B2")+
 geom_smooth(method="lm",se=FALSE,colour="#D55E00",linewidth=.5)+
 labs(x=label,y="Signature score",subtitle=sprintf("Pearson r = %+.3f; n = %d",r,nrow(d)))
}
b[,v:=v+INT];f5a$data <- b
f5b <- f5b + labs(y="Pearson r")
export(f5a / ((scatter("BRAF_RAF_score","BRAF-RAS score") | scatter("differentiation_score","Thyroid differentiation score"))) /
 (scatter("ERK_score","ERK-output score") | scatter("purity","Tumour purity")) / f5b +
 plot_layout(heights=c(1,1.2,1.2,.7))+plot_annotation(tag_levels="a"),"Fig5",8.0)

ad <- fread(file.path(TAB,"Supplementary_Table_S11_composition_coefficients.csv"))
al <- rbind(ad[,.(unit,n,beta=beta_raw,lo=lo_raw,hi=hi_raw,model="Unadjusted")],ad[,.(unit,n,beta=beta_adj,lo=lo_adj,hi=hi_adj,model="Adjusted")])
adjplot <- function(d,title)ggplot(d,aes(beta,factor(unit,levels=rev(c(SIG,"SCORE"))),colour=model))+
 geom_vline(xintercept=0,linetype=2,colour="grey55")+
 geom_errorbar(aes(xmin=lo,xmax=hi),orientation="y",width=.2,position=position_dodge(.5))+
 geom_point(position=position_dodge(.5),size=2)+scale_colour_manual(values=c(Unadjusted="#0072B2",Adjusted="#D55E00"))+
 labs(x="Tumour-normal coefficient (95% CI)",y=NULL,title=title,colour=NULL)
md <- fread(file.path(TAB,"Supplementary_Table_S10_mutation_stratification.csv"))
md[,unit:=factor(unit,levels=rev(c(SIG,"SCORE")))]
mp <- ggplot(md,aes(subgroup,unit,fill=d_cohen))+geom_tile(colour="white",linewidth=.7)+
 geom_text(aes(label=sprintf("%+.2f",d_cohen)),size=3)+scale_fill_gradient2(low="#0072B2",mid="white",high="#D55E00",midpoint=0,name="SMD")+
 labs(x=NULL,y=NULL,title="Mutation-stratified tumour-normal differences")+
 scale_x_discrete(labels=c(BRAFmut="BRAF V600E\n(n = 283)",RASmut="RAS mutant\n(n = 58)",driverneg="No recorded BRAF/RAS\n(n = 146)"))+
 theme(panel.grid=element_blank())
export((adjplot(al[unit!="SCORE"],"Signature genes") | adjplot(al[unit=="SCORE"],"Signature score")) / mp+
 plot_layout(heights=c(1,1))+plot_annotation(tag_levels="a"),"Fig7",5.7)

k <- which(g30$lab %in% c("Normal","ATC"));a30 <- g30;a30$G <- g30$G[,k];a30$lab <- g30$lab[k]
r30 <- roc(a30$lab,sc_within(a30,COEF),levels=c("Normal","ATC"),direction="<",quiet=TRUE)
g65 <- load_ds(file.path(GEO,"GSE65144_series_matrix.txt.gz"),
 function(md)ifelse(grepl("^N",md$title),"Normal","ATC"),function(t)t)
M65 <- read_gpl570(file.path(GEO,"GSE65144_series_matrix.txt.gz"))
sc65 <- setdiff(names(M65),"probe")
M65[, (sc65):=lapply(.SD,function(x)log2(as.numeric(x)+1)),.SDcols=sc65]
sub65 <- merge(A,M65,by="probe")
avg65 <- sub65[,lapply(.SD,mean,na.rm=TRUE),by=symbol,.SDcols=sc65]
g65$G <- as.matrix(avg65[,..sc65]);rownames(g65$G) <- avg65$symbol
r65 <- roc(g65$lab,sc_within(g65,COEF),levels=c("Normal","ATC"),direction="<",quiet=TRUE)
cat("ATC AUCs:",as.numeric(auc(r30)),as.numeric(auc(r65)),"\n")
stopifnot(abs(as.numeric(auc(r30))-.8707070707)<.001, abs(as.numeric(auc(r65))-.8782051282)<.001)
ar <- rbindlist(lapply(list(GSE33630=r30,GSE65144=r65),function(r)data.table(fpr=1-r$specificities,tpr=r$sensitivities)),idcol="dataset")
atc <- ggplot(ar,aes(fpr,tpr,colour=dataset))+geom_path(linewidth=.7)+geom_abline(slope=1,intercept=0,linetype=2,colour="grey60")+
 coord_equal()+scale_colour_manual(values=c("#0072B2","#D55E00"),labels=c("GSE33630: AUC 0.871","GSE65144: AUC 0.878"))+
 labs(x="1 - specificity",y="Sensitivity",colour=NULL,title="ATC versus normal")
export(atc / (f4b | f4c)+plot_layout(heights=c(1,1))+plot_annotation(tag_levels="a"),"Supplementary_Fig_S8",5.8)

# Patient-paired differences: all available matched tissue pairs.
ee <- as.data.table(t(E[SIG,colnames(E29)]));ee[,sample:=colnames(E29)];ee[,pid:=substr(sample,1,12)];ee[,group:=grp]
dd <- merge(ee[group=="Tumor"],ee[group=="Normal"],by="pid",suffixes=c("_T","_N"))
pd <- rbindlist(lapply(SIG,function(g)data.table(gene=g,value=dd[[paste0(g,"_T")]]-dd[[paste0(g,"_N")]])))
pp <- ggplot(pd,aes(gene,value))+geom_hline(yintercept=0,linetype=2,colour="grey55")+
 geom_boxplot(fill="#0072B2",alpha=.3,outlier.shape=NA,width=.5)+geom_point(position=position_jitter(width=.12,seed=11),size=.8,alpha=.5)+labs(x=NULL,y="Paired difference in log2(CPM + 1)")
sp <- copy(score);sp[,pid:=substr(sample,1,12)];sd <- merge(sp[group=="PTC"],sp[group=="Normal"],by="pid",suffixes=c("_T","_N"))
stopifnot(nrow(sd)==59,nrow(dd)==59)
pp2 <- ggplot(sd,aes(x="59 matched patients",y=risk_T-risk_N))+geom_hline(yintercept=0,linetype=2,colour="grey55")+
 geom_boxplot(fill="#D55E00",alpha=.3,outlier.shape=NA,width=.5)+geom_point(position=position_jitter(width=.12,seed=11),size=.9,alpha=.5)+labs(x=NULL,y="Paired difference in signature score")
export((pp | pp2)+plot_layout(widths=c(2,1))+plot_annotation(tag_levels="a"),"Supplementary_Fig_S1",3.5)
clinicalbox <- function(v,lev,title) {
 d <- cl[!is.na(get(v)) & get(v) %in% lev];d[,g:=factor(get(v),levels=lev)]
 ggplot(d,aes(g,risk+INT))+geom_boxplot(fill="#0072B2",alpha=.3,outlier.shape=NA,width=.55)+
 geom_point(position=position_jitter(width=.12,seed=11),size=.65,alpha=.35)+
 labs(x=NULL,y="Signature score",title=title)+scale_x_discrete(labels=c(PTC_classic="Classical",FVPTC="Follicular\nvariant",Tall_columnar="Tall cell",Hurthle="Hurthle",Other="Other",T1="T1",T2="T2",T3="T3",T4="T4",I="I",II="II",III="III",IV="IV",N0="N0",N1="N1",N1A="N1A",N1B="N1B"))+theme(axis.text.x=element_text(size=7),plot.title=element_text(size=10))
}
export((clinicalbox("histgrp",c("PTC_classic","FVPTC","Tall_columnar","Hurthle","Other"),"Histological category") |
 clinicalbox("t4",paste0("T",1:4),"T stage")) /
 (clinicalbox("stage4",c("I","II","III","IV"),"AJCC stage") | clinicalbox("n_stage",c("N0","N1","N1A","N1B"),"Nodal status")) +
 plot_annotation(tag_levels="a"),"Supplementary_Fig_S4",5.0)
cc <- fread(file.path(DIR,"drug_prediction_control.csv"));cv <- setNames(cc$value,cc$metric)
ch <- data.table(metric=rep(c("All associated chemicals","Chemicals linked to >=3 genes"),each=2),
 type=rep(c("Observed","Random-set median"),2),value=c(cv["n_chem_hits"],cv["null_n_chem_median"],cv["n_chem_ge3"],cv["null_ge3_median"]))
cp <- ggplot(ch,aes(type,value,fill=type))+geom_col(width=.55)+geom_text(aes(label=value),vjust=-.3,size=3)+
 facet_wrap(~metric,scales="free_y")+scale_fill_manual(values=c(Observed="#D55E00",`Random-set median`="#0072B2"),guide="none")+
 scale_y_continuous(expand=expansion(mult=c(0,.12)))+labs(x=NULL,y="Chemical count")+
 theme(axis.text.x=element_text(size=8))
export(cp,"Supplementary_Fig_S3",3.1)
fwrite(ex,file.path(FIG,"external_auc_plot_data.csv"))
fwrite(pub[,.(BRAF_RAF_score,differentiation_score,ERK_score,purity,risk)],file.path(FIG,"molecular_axis_plot_data.csv"))
cat("Final publication figures completed. Fig6 is produced by step 22.\n")
