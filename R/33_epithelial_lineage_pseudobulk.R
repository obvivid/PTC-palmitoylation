# Independent-donor pseudobulk inference and bounded sensitivities.
options(warn=1)
Sys.setenv(OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2',MKL_NUM_THREADS='2')
suppressPackageStartupMessages({library(data.table);library(Matrix);library(edgeR);library(limma)})
SRC <- file.path(RESULT_DIR, 'extended/singlecell_lineage')
O <- file.path(SRC, '02_pseudobulk');dir.create(O,recursive=TRUE,showWarnings=FALSE)
SIG <- c('ABHD17C','ZDHHC6','ABHD17B','ZDHHC20')
rd <- function(p) fread(p,tmpdir=O)
frozen <- rd(file.path(SRC,'01_annotation/frozen_cell_list.csv'))
md <- rd(file.path(SRC,'01_annotation/pseudobulk_metadata_frozen.csv'))
ledger <- rd(file.path(SRC,'00_gate/sample_donor_ledger.csv'))
raw <- rd(file.path(SRC,'01_annotation/pseudobulk_counts_frozen.csv.gz'))
pb <- as.matrix(raw[,-1]);rownames(pb)<-raw$gene
stopifnot(identical(colnames(pb),md$sample))
eligible <- md[treatment=='untreated' & n_cells>=50]
stopifnot(nrow(eligible)==8,uniqueN(eligible$donor_id)==7)
# Donor 7 supplies both tissues. Preserve all three eligible normal donors.
# This is a design correction after review, not a preregistered rule.
primary <- eligible[sample!='PTC07']
stopifnot(nrow(primary)==7,uniqueN(primary$donor_id)==7,
          sum(primary$tissue=='NORM')==3,sum(primary$tissue=='PTC')==4)

diags <- list(); logs <- list(); units <- list()
fit_pb <- function(counts,met,label,groupcol='tissue',levels=c('NORM','PTC'),blocked=FALSE) {
    met<-copy(met)
    counts<-counts[,met$unit,drop=FALSE]
    grp<-factor(met[[groupcol]],levels=levels);X<-model.matrix(~grp)
    stopifnot(!anyNA(grp),qr(X)$rank==ncol(X),all(colnames(counts)==met$unit))
    y<-DGEList(counts);keep<-filterByExpr(y,design=X)
    stopifnot(all(keep[SIG]))
    y<-calcNormFactors(y[keep,,keep.lib.sizes=FALSE])
    v<-voom(y,X,plot=FALSE)
    cor1<-cor2<-NA_real_
    if(blocked){
        cor1<-duplicateCorrelation(v,X,block=met$donor_id)$consensus
        stopifnot(is.finite(cor1),abs(cor1)<1)
        v<-voom(y,X,block=met$donor_id,correlation=cor1,plot=FALSE)
        cor2<-duplicateCorrelation(v,X,block=met$donor_id)$consensus
        stopifnot(is.finite(cor2),abs(cor2)<1)
        fit<-eBayes(lmFit(v,X,block=met$donor_id,correlation=cor2))
    } else fit<-eBayes(lmFit(v,X))
    tt<-topTable(fit,coef=2,number=Inf,sort.by='none',confint=TRUE)
    z<-as.data.table(tt[SIG,c('logFC','CI.L','CI.R','P.Value')],keep.rownames='gene')
    z[,`:=`(BH4=p.adjust(P.Value,'BH',n=4),analysis=label,n_units=nrow(met),
             n_donors=uniqueN(met$donor_id),n_control=sum(grp==levels[1]),n_case=sum(grp==levels[2]))]
    ii<-match(SIG,rownames(fit$coefficients))
    lo<-fit$coefficients[ii,2]-qt(.975,fit$df.total[ii])*fit$stdev.unscaled[ii,2]*sqrt(fit$s2.post[ii])
    stopifnot(isTRUE(all.equal(as.numeric(lo),z$CI.L,tolerance=1e-10)))
    diags[[label]]<<-data.table(analysis=label,n_units=nrow(met),n_donors=uniqueN(met$donor_id),
                              n_genes=sum(keep),n_control=sum(grp==levels[1]),n_case=sum(grp==levels[2]),
                              residual_df=fit$df.residual[1],moderated_df=fit$df.total[1],
                              block=blocked,correlation_round1=cor1,correlation_round2=cor2)
    lc<-as.data.table(t(cpm(y,log=TRUE,prior.count=1)[SIG,,drop=FALSE]),keep.rownames='unit')
    logs[[label]]<<-merge(met,lc,by='unit')[,analysis:=label]
    units[[label]]<<-copy(met)[,analysis:=label]
    z
}
md[,unit:=sample];eligible[,unit:=sample];primary[,unit:=sample]
q1cases<-list(Q1_primary_independent=primary,
             Q1_all_eight_unblocked=eligible,
             Q1_retain_PTC07=eligible[sample!='NORM07'],
             Q1_remove_donor7=eligible[donor_id!=7],
             Q1_min100=primary[n_cells>=100])
q1<-rbindlist(lapply(names(q1cases),function(k)fit_pb(pb,q1cases[[k]],k)))
loo<-rbindlist(lapply(primary$donor_id,function(id)fit_pb(pb,primary[donor_id!=id],paste0('Q1_leave_out_',id))))
treated<-md[n_cells>=50]
stopifnot(nrow(treated)==10)
tr<-fit_pb(pb,treated,'Q1_include_treated_blocked',blocked=TRUE)

# Aggregate identity-defined cells once; retain both the complete and primary units.
counts_path <- if (exists('SINGLECELL_COUNT_OBJECT')) SINGLECELL_COUNT_OBJECT else file.path(RESULT_DIR, 'extended/singlecell/counts_annotated_193581.rds')
L<-readRDS(counts_path);C<-L$counts
cells<-frozen[q2_class %in% c('malignant','normal_follicular')]
mm<-match(cells$cell,colnames(C));stopifnot(!anyNA(mm))
ukey<-paste(cells$sample,cells$q2_class,sep='|')
spl<-split(seq_len(nrow(cells)),factor(ukey,levels=unique(ukey)))
pb2<-vapply(spl,function(ii)as.numeric(Matrix::rowSums(C[,mm[ii],drop=FALSE])),numeric(nrow(C)))
rownames(pb2)<-rownames(C);storage.mode(pb2)<-'integer'
m2<-data.table(unit=names(spl),n_cells=lengths(spl))
m2[,`:=`(sample=sub('\\|.*$','',unit),identity=sub('^.*\\|','',unit))]
m2[,`:=`(donor_id=md$donor_id[match(sample,md$sample)],
          treatment=md$treatment[match(sample,md$sample)],tissue=md$tissue[match(sample,md$sample)])]
fwrite(data.table(gene=rownames(pb2),pb2),file.path(O,'identity_pseudobulk_counts.csv.gz'))
fwrite(m2,file.path(O,'identity_pseudobulk_metadata.csv'))
rm(L,C);gc()
q2main<-m2[sample %in% primary$sample & n_cells>=50]
stopifnot(nrow(q2main)==7,uniqueN(q2main$donor_id)==7)
q2<-rbindlist(list(
    fit_pb(pb2,q2main,'Q2_primary_independent','identity',c('normal_follicular','malignant')),
    fit_pb(pb2,m2[sample %in% eligible$sample & n_cells>=50],'Q2_all_eight_unblocked','identity',c('normal_follicular','malignant')),
    fit_pb(pb2,m2[sample %in% eligible$sample & sample!='NORM07' & n_cells>=50],'Q2_retain_PTC07','identity',c('normal_follicular','malignant')),
    fit_pb(pb2,m2[sample %in% eligible$sample & donor_id!=7 & n_cells>=50],'Q2_remove_donor7','identity',c('normal_follicular','malignant')),
    fit_pb(pb2,m2[sample %in% eligible$sample],'Q2_lowcell_unit_blocked','identity',c('normal_follicular','malignant'),blocked=TRUE)))
fwrite(q1,file.path(O,'Q1_results.csv'));fwrite(q2,file.path(O,'Q2_results.csv'))
fwrite(loo,file.path(O,'Q1_leave_one_donor_out.csv'));fwrite(tr,file.path(O,'Q1_treated_sensitivity.csv'))
fwrite(rbindlist(diags),file.path(O,'model_diagnostics.csv'))
fwrite(rbindlist(logs,fill=TRUE),file.path(O,'per_unit_expression.csv'))
fwrite(rbindlist(units,fill=TRUE),file.path(O,'analysis_unit_ledger.csv'))
capture.output(sessionInfo(),file=file.path(O,'sessionInfo.txt'))
cat('\nPRIMARY Q1\n');print(q1[analysis=='Q1_primary_independent'])
cat('\nPRIMARY Q2\n');print(q2[analysis=='Q2_primary_independent'])
cat('\nModel diagnostics\n');print(rbindlist(diags))
