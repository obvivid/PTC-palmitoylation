suppressPackageStartupMessages({library(data.table);library(pROC)})
ROOT <- file.path(RESULT_DIR,'single_sample')
GEO <- file.path(PROJECT_ROOT,'data','geo')
dir.create(file.path(ROOT,'tables'),recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(ROOT,'workflow'),recursive=TRUE,showWarnings=FALSE)
source(file.path(SOURCE_ROOT,'single_sample','single_sample_score.R'))
PUBLIC <- file.path(SOURCE_ROOT,'R','05_external_ptc.R')
# Load the published parser definitions, without running the original analysis.
env<-new.env(parent=.GlobalEnv)
for(e in parse(PUBLIC,encoding='UTF-8')){
 if(is.call(e)&&identical(e[[1]],as.name('<-'))&&is.call(e[[3]])&&identical(e[[3]][[1]],as.name('function')))eval(e,envir=env)
}
ann<-fread(file.path(GEO,'GPL570.annot.gz'),skip=27,header=TRUE,sep='\t',quote='',na.strings=c('','NA'),fill=TRUE)
setnames(ann,1,'probe')
mapping<-ann[!is.na(`Gene symbol`)&grepl('^[A-Za-z0-9][A-Za-z0-9_.-]*$',`Gene symbol`),
             .(probe=as.character(probe),symbol=as.character(`Gene symbol`))]
stopifnot(!anyDuplicated(mapping$probe))
read_all<-function(filename,linear=FALSE){
 M<-env$read_gpl570(file.path(GEO,filename));cols<-setdiff(names(M),'probe')
 if(linear)M[,(cols):=lapply(.SD,function(v)log2(as.numeric(v)+1)),.SDcols=cols]
 sub<-merge(mapping,M,by='probe',allow.cartesian=FALSE)
 G<-sub[,lapply(.SD,function(v)mean(as.numeric(v),na.rm=TRUE)),by=symbol,.SDcols=cols]
 mat<-as.matrix(G[,-1]);rownames(mat)<-G$symbol
 stopifnot(all(is.finite(mat)));mat
}
message('Reading full expression profiles for the four PTC cohorts')
g30<-read_all('GSE33630_series_matrix.txt.gz')
g55<-read_all('GSE35570_matrix.txt.gz')
g29<-read_all('GSE29265_matrix.txt.gz')
g60<-read_all('GSE60542_matrix.txt.gz')
make_profile <- function(G,md,lab,pid) {
  k <- !is.na(lab)&lab%in%c('PTC','Normal')
  list(G=G[,md$gsm[k],drop=FALSE],lab=lab[k],pid=pid[k])
}
m30 <- env$getmeta(file.path(GEO,'GSE33630_series_matrix.txt.gz'))
l30 <- gsub('^[^:]*: *','',env$pick(m30$char_lines,'diagnostic'))
l30 <- ifelse(grepl('non-tumor|normal',l30,ignore.case=TRUE),'Normal',ifelse(grepl('papillary',l30,ignore.case=TRUE),'PTC',NA))
p30 <- env$pid_lead(m30$title)
m55 <- env$getmeta(file.path(GEO,'GSE35570_matrix.txt.gz'))
l55 <- gsub('^[^:]*: *','',env$pick(m55$char_lines,'type of thyroid'))
l55 <- ifelse(grepl('normal',l55,ignore.case=TRUE),'Normal',ifelse(grepl('ptc',l55,ignore.case=TRUE),'PTC',NA))
p55 <- env$pid_last(m55$title); l55[p55%in%intersect(unique(p30),unique(p55))] <- NA
m29 <- env$getmeta(file.path(GEO,'GSE29265_matrix.txt.gz'))
l29 <- gsub('^[^:]*: *','',env$pick(m29$char_lines,'^tissue'))
l29 <- ifelse(grepl('nontumor|non-tumor|normal',l29,ignore.case=TRUE),'Normal',ifelse(grepl('papillary',l29,ignore.case=TRUE),'PTC',NA))
datasets <- list(GSE33630=make_profile(g30,m30,l30,p30),
                 GSE35570=make_profile(g55,m55,l55,p55),
                 GSE29265=make_profile(g29,m29,l29,env$pid_pair(m29$title)))
md60<-env$getmeta(file.path(GEO,'GSE60542_matrix.txt.gz'))
lab60<-ifelse(grepl('Normal thyroid',md60$title,ignore.case=TRUE),'Normal',ifelse(grepl('Papillary thyroid carcinoma',md60$title,ignore.case=TRUE),'PTC',NA))
pid60<-sub('^([0-9]+).*$','\\1',md60$title);keep60<-which(!is.na(lab60))
agg60<-sapply(split(keep60,paste(pid60[keep60],lab60[keep60],sep='_')),function(ii)rowMeans(g60[,md60$gsm[ii],drop=FALSE]))
datasets$GSE60542 <- list(G=agg60,lab=sub('^[^_]*_','',colnames(agg60)),pid=sub('_.*$','',colnames(agg60)))
for(nm in names(datasets)){
 d<-datasets[[nm]]
 d$role<-if(nm=='GSE33630')'development'else if(nm=='GSE29265')'sensitivity_unresolved_independence'else'primary_transfer'
 d$positive_label<-'PTC';d$negative_label<-'Normal';d$bootstrap_unit<-'patient';datasets[[nm]]<-d
}
universe<-sort(rownames(datasets$GSE33630$G))
model<-list(version='1.0.0',date='2026-10-13',input_scale='GPL570_gene_log2',
 up=c('ABHD17C','ZDHHC20'),down=c('ABHD17B','ZDHHC6'),universe=universe,
 development='GSE33630_PTC_Normal',thresholds=list(rank=list(youden=0,sensitivity95=0),contrast=list(youden=0,sensitivity95=0)),
description='Exploratory frozen single-sample scoring; not a clinically validated assay')
mapping<-mapping[symbol%in%universe]
writeLines(universe,file.path(ROOT,'workflow/frozen_gene_universe.txt'))
fwrite(mapping,file.path(ROOT,'workflow/frozen_probe_map.csv'))

calc<-function(y,s,t){
 p<-s>=t;tp<-sum(p&y);fn<-sum(!p&y);tn<-sum(!p&!y);fp<-sum(p&!y)
 se<-tp/(tp+fn);sp<-tn/(tn+fp)
 c(tp=tp,fn=fn,tn=tn,fp=fp,sensitivity=se,specificity=sp,accuracy=(tp+tn)/length(y),balanced_accuracy=(se+sp)/2)
}
choose<-function(y,s,rule){
 u<-sort(unique(s));ts<-c(-Inf,head(u,-1)+(tail(u,-1)-head(u,-1))/2,Inf)
 x<-rbindlist(lapply(ts,function(t)as.list(c(threshold=t,calc(y,s,t)))))
 if(rule=='youden'){
  x[,j:=sensitivity+specificity-1];x<-x[abs(j-max(j))<1e-12];setorder(x,-sensitivity,-threshold)
 }else{
  x<-x[sensitivity>=.95-1e-12];x<-x[abs(specificity-max(specificity))<1e-12];setorder(x,-threshold)
 }
 t<-x$threshold[1];if(!is.finite(t))stop('No finite development threshold; branch cannot be frozen.')
 t
}
wilson<-function(k,n){z<-qnorm(.975);p<-k/n;v<-1+z*z/n;c((p+z*z/(2*n)-z*sqrt(p*(1-p)/n+z*z/(4*n*n)))/v,(p+z*z/(2*n)+z*sqrt(p*(1-p)/n+z*z/(4*n*n)))/v)}
bootstrap<-function(y,s,t,pid,B=2000){
 set.seed(20261013);p<-s>=t
 m<-rowsum(cbind(tp=as.integer(p&y),fn=as.integer(!p&y),tn=as.integer(!p&!y),fp=as.integer(p&!y)),pid,reorder=FALSE)
 z<-replicate(B,colSums(m[sample.int(nrow(m),nrow(m),replace=TRUE),,drop=FALSE]))
 a<-(z[1,]+z[3,])/colSums(z);ba<-(z[1,]/(z[1,]+z[2,])+z[3,]/(z[3,]+z[4,]))/2
 q<-function(v)as.numeric(quantile(v[is.finite(v)],c(.025,.975)))
 c(accuracy_lo=q(a)[1],accuracy_hi=q(a)[2],balanced_accuracy_lo=q(ba)[1],balanced_accuracy_hi=q(ba)[2])
}
metrics<-function(d,s,t,method,op,tag=d$tag,role=d$role){
 y<-d$lab==d$positive_label;c<-calc(y,s,t);se<-wilson(c['tp'],c['tp']+c['fn']);sp<-wilson(c['tn'],c['tn']+c['fp'])
 ro<-roc(y,s,levels=c(FALSE,TRUE),direction='<',quiet=TRUE);ac<-as.numeric(ci.auc(ro,method='delong'))
 as.data.table(c(list(dataset=tag,role=role,method=method,operating_point=op,n=length(y),n_positive=sum(y),n_negative=sum(!y),
 positive_label=d$positive_label,negative_label=d$negative_label,bootstrap_unit=d$bootstrap_unit,threshold=t,
 auc=as.numeric(auc(ro)),auc_lo=ac[1],auc_hi=ac[3]),as.list(c),
 list(sensitivity_lo=se[1],sensitivity_hi=se[2],specificity_lo=sp[1],specificity_hi=sp[2]),as.list(bootstrap(y,s,t,d$pid))))
}
scores<-function(d,method)score_sample_matrix(d$G,model,method)$score
dev<-datasets$GSE33630;dev$tag<-'GSE33630'
for(method in c('rank','contrast')){
 s<-scores(dev,method);y<-dev$lab=='PTC'
 for(op in c('youden','sensitivity95'))model$thresholds[[method]][[op]]<-choose(y,s,op)
}
saveRDS(model,file.path(ROOT,'workflow/model.rds'))
fwrite(rbindlist(lapply(names(model$thresholds),function(m) data.table(method=m,operating_point=names(model$thresholds[[m]]),threshold=unlist(model$thresholds[[m]])))),file.path(ROOT,'tables/frozen_thresholds.csv'))
fwrite(data.table(gene=universe,log2_expression=dev$G[universe,1]),file.path(ROOT,'workflow/example_single_sample.csv'))

message('Reading the predefined ATC and follicular extension profiles')
meta30<-env$getmeta(file.path(GEO,'GSE33630_series_matrix.txt.gz'))
diag30<-gsub('^[^:]*: *','',env$pick(meta30$char_lines,'diagnostic'))
labs30<-ifelse(grepl('non-tumor|normal',diag30,ignore.case=TRUE),'Normal',ifelse(grepl('anaplastic',diag30,ignore.case=TRUE),'ATC',NA))
keep<-which(!is.na(labs30))
datasets$ATC_GSE33630<-list(G=g30[,meta30$gsm[keep],drop=FALSE],lab=labs30[keep],pid=env$pid_lead(meta30$title[keep]),
 tag='ATC_GSE33630',role='histology_extension_shared_controls',positive_label='ATC',negative_label='Normal',bootstrap_unit='patient_or_ATC_record')
g65<-read_all('GSE65144_series_matrix.txt.gz',linear=TRUE);m65<-env$getmeta(file.path(GEO,'GSE65144_series_matrix.txt.gz'))
l65<-ifelse(grepl('^N',m65$title),'Normal','ATC')
datasets$ATC_GSE65144<-list(G=g65[,m65$gsm,drop=FALSE],lab=l65,pid=m65$gsm,tag='ATC_GSE65144',role='histology_extension',
 positive_label='ATC',negative_label='Normal',bootstrap_unit='record')
g82<-read_all('GSE82208_series_matrix.txt.gz');m82<-env$getmeta(file.path(GEO,'GSE82208_series_matrix.txt.gz'))
l82<-ifelse(grepl('adenoma',m82$title,ignore.case=TRUE),'FTA',ifelse(grepl('cancer',m82$title,ignore.case=TRUE),'FTC',NA))
k82<-which(!is.na(l82));stopifnot(sum(l82=='FTA',na.rm=TRUE)==25,sum(l82=='FTC',na.rm=TRUE)==27)
datasets$FTC_GSE82208<-list(G=g82[,m82$gsm[k82],drop=FALSE],lab=l82[k82],pid=m82$gsm[k82],tag='FTC_GSE82208',role='benign_malignant_extension',
 positive_label='FTC',negative_label='FTA',bootstrap_unit='record')
stopifnot(sum(datasets$ATC_GSE33630$lab=='ATC')==11,sum(l65=='ATC')==12,sum(l65=='Normal')==13)

results<-list();record_scores<-list();validation<-list()
for(nm in names(datasets)){
 d<-datasets[[nm]];d$tag<-nm;datasets[[nm]]<-d
 for(method in c('rank','contrast')){
  s<-scores(d,method)
  for(op in c('youden','sensitivity95')){
   t<-model$thresholds[[method]][[op]]
   results[[paste(nm,method,op)]]<-metrics(d,s,t,method,op)
  }
  record_scores[[paste(nm,method)]]<-data.table(dataset=nm,record=colnames(d$G),pid=d$pid,truth=d$lab,method=method,score=s)
 }
 validation[[nm]]<-data.table(dataset=nm,background_genes=length(universe),required_present=all(universe%in%rownames(d$G)),n=ncol(d$G))
}
result<-rbindlist(results);fwrite(result,file.path(ROOT,'tables/fixed_pipeline_results.csv'))
fwrite(rbindlist(record_scores),file.path(ROOT,'tables/record_scores.csv'))
fwrite(rbindlist(validation),file.path(ROOT,'tables/input_checks.csv'))

# Patient-grouped internal threshold validation, five folds; score definition stays frozen.
set.seed(20261013);u<-unique(dev$pid);hasnormal<-u%in%dev$pid[dev$lab=='Normal'];f<-integer(length(u))
for(h in c(TRUE,FALSE)){ii<-which(hasnormal==h);f[ii]<-rep(1:5,length.out=length(ii))[order(runif(length(ii)))]}
fid<-f[match(dev$pid,u)];cv<-list();folds<-list()
for(method in c('rank','contrast'))for(op in c('youden','sensitivity95')){
 s<-scores(dev,method);pred<-logical(length(s));th<-numeric(length(s));y<-dev$lab=='PTC'
 for(k in 1:5){tr<-fid!=k;te<-fid==k;t<-choose(y[tr],s[tr],op);pred[te]<-s[te]>=t;th[te]<-t
  folds[[paste(method,op,k)]]<-data.table(method=method,operating_point=op,fold=k,threshold=t,n_test=sum(te))
 }
 c<-calc(y,as.numeric(pred),.5)
 cv[[paste(method,op)]]<-as.data.table(c(list(method=method,operating_point=op,n=length(y)),as.list(c)))
}
fwrite(rbindlist(cv),file.path(ROOT,'tables/development_patient_grouped_cv.csv'))
fwrite(rbindlist(folds),file.path(ROOT,'tables/development_fold_thresholds.csv'))

# Same sample scored alone or in different batch compositions.
dd<-datasets$GSE60542
invariance<-list()
for(method in c('rank','contrast')){
 s<-scores(dd,method);single<-score_single_sample(setNames(dd$G[,1],rownames(dd$G)),model,method)$score
 idx<-which(dd$lab=='PTC');only_positive<-score_sample_matrix(dd$G[,idx,drop=FALSE],model,method)$score
 idxn<-which(dd$lab=='Normal');only_negative<-score_sample_matrix(dd$G[,idxn,drop=FALSE],model,method)$score
 invariance[[method]]<-data.table(method=method,single_difference=abs(single-s[1]),positive_batch_max_difference=max(abs(only_positive-s[idx])),normal_batch_max_difference=max(abs(only_negative-s[idxn])))
}
inv<-rbindlist(invariance);stopifnot(all(as.matrix(inv[,-1])<1e-12));fwrite(inv,file.path(ROOT,'tables/batch_invariance.csv'))
errors<-list()
expect_rejection<-function(name,fun){e<-tryCatch({fun();NULL},error=function(e)conditionMessage(e));if(is.null(e))stop('Expected rejection did not occur: ',name);data.table(check=name,rejected=TRUE,message=e)}
x<-setNames(dd$G[,1],rownames(dd$G))
errors[[1]]<-expect_rejection('missing_signature_gene',function()score_single_sample(x[names(x)!='ABHD17C'],model,'contrast'))
errors[[2]]<-expect_rejection('missing_rank_background_gene',function()score_single_sample(x[-which(names(x)==setdiff(universe,c(model$up,model$down))[1])],model,'rank'))
errors[[3]]<-expect_rejection('nonfinite_gene',function(){v<-x;v['ZDHHC20']<-NA;score_single_sample(v,model,'contrast')})
errors[[4]]<-expect_rejection('unsupported_RNAseq_scale',function()score_single_sample(x,model,'contrast',input_scale='RNAseq_logCPM'))
fwrite(rbindlist(errors),file.path(ROOT,'tables/input_rejections.csv'))

# Prespecified common and gene-specific measurement perturbations.
stress<-list()
for(method in c('rank','contrast')){
 s0<-scores(dd,method)
 for(case in c('common_plus5','common_times2_plus3','ZDHHC20_minus0.5','ZDHHC20_minus0.25','ZDHHC20_plus0.25','ZDHHC20_plus0.5')){
  dx<-dd
  if(case=='common_plus5')dx$G<-dd$G+5
  else if(case=='common_times2_plus3')dx$G<-dd$G*2+3
  else {v<-switch(case,ZDHHC20_minus0.5=-.5,ZDHHC20_minus0.25=-.25,ZDHHC20_plus0.25=.25,ZDHHC20_plus0.5=.5);dx$G['ZDHHC20',]<-dx$G['ZDHHC20',]+v}
  sx<-scores(dx,method);r<-metrics(dd,sx,model$thresholds[[method]]$youden,method,'youden',tag=case,role='synthetic_measurement_stress')
  r[,max_score_difference:=max(abs(sx-s0))];r[,changed_calls:=sum((sx>=model$thresholds[[method]]$youden)!=(s0>=model$thresholds[[method]]$youden))]
  stress[[paste(method,case)]]<-r
  if(case=='common_plus5')stopifnot(max(abs(sx-s0))<1e-12)
  if(case=='common_times2_plus3'&&method=='rank')stopifnot(max(abs(sx-s0))<1e-12)
 }
}
fwrite(rbindlist(stress),file.path(ROOT,'tables/measurement_stress.csv'))

# Synthetic matched-profile signal mixing; not a physical cell-purity assay.
matched<-intersect(dd$pid[dd$lab=='PTC'],dd$pid[dd$lab=='Normal']);stopifnot(length(matched)>0)
mix<-list();mixrecords<-list()
for(method in c('rank','contrast'))for(a in c(0,.25,.5,.75,1)){
 vals<-sapply(matched,function(pid){
  ti<-which(dd$pid==pid&dd$lab=='PTC');ni<-which(dd$pid==pid&dd$lab=='Normal');stopifnot(length(ti)==1,length(ni)==1)
  xx<-log2(a*2^dd$G[,ti]+(1-a)*2^dd$G[,ni]);names(xx)<-rownames(dd$G)
  score_single_sample(xx,model,method)$score
 })
 t<-model$thresholds[[method]]$youden
 mix[[paste(method,a)]]<-data.table(method=method,tumor_signal_weight=a,n_profiles=length(vals),positive=sum(vals>=t),positive_fraction=mean(vals>=t),threshold=t,median_score=median(vals))
 mixrecords[[paste(method,a)]]<-data.table(method=method,tumor_signal_weight=a,pid=matched,score=vals,threshold=t,positive=vals>=t)
}
fwrite(rbindlist(mix),file.path(ROOT,'tables/synthetic_signal_mixing.csv'))
fwrite(rbindlist(mixrecords),file.path(ROOT,'tables/synthetic_signal_mixing_records.csv'))

saveRDS(list(model=model,datasets=datasets,results=result),file.path(ROOT,'analysis_cache.rds'))
writeLines(capture.output(sessionInfo()),file.path(ROOT,'session_info.txt'))
cat('Fixed single-sample analysis completed. Tables and excluded checkpoints are under outputs.\n')
