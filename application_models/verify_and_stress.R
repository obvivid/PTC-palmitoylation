ROOT<-file.path(RESULT_DIR,'weighted_single_sample')
source(file.path(SOURCE_ROOT,'application_models','models.R'))
suppressPackageStartupMessages(library(jsonlite))
cache<-readRDS(file.path(RESULT_DIR,'single_sample','analysis_cache.rds'))
state<-readRDS(file.path(ROOT,'workflow/development_locked_models.rds'));locked<-state$models
D<-fread(file.path(ROOT,'tables/input_features.csv'))
# Export fixed parameters for a second, non-R implementation of prediction.
exports<-list()
for(id in names(locked)){
 m<-locked[[id]];e<-list(id=id,candidate=m$candidate,parameter=m$parameter,thresholds=as.list(m$thresholds))
 if(m$candidate$family!='baseline'){
  e$type<-m$type;e$feature_names<-m$feature_names;e$keep<-m$keep;e$mean<-unname(m$mean);e$sd<-unname(m$sd)
  if(m$type=='linear'){e$intercept<-m$intercept;e$weights<-m$weights}
  if(m$type=='constant')e$constant<-m$constant
  if(m$type=='svm'){e$SV<-m$fitted$SV;e$coefs<-as.numeric(m$fitted$coefs);e$rho<-m$fitted$rho;e$gamma<-m$fitted$gamma;e$direction<-m$svm_direction}
  if(m$type=='forest'){
   trees<-rbindlist(lapply(seq_len(m$fitted$ntree),function(i){a<-as.data.table(getTree(m$fitted,i));setnames(a,c('left','right','variable','cut','status','prediction'));a[,`:=`(tree=i,node=seq_len(.N))];a}))
   # Preserve a split exactly when a sample lies on the numerical boundary.
   # Default CSV precision can otherwise reverse a random-forest branch.
   trees[,cut:=sprintf('%.17g',cut)]
   e$ntree<-m$fitted$ntree;e$trees_file<-paste0('forest_',id,'.csv');fwrite(trees,file.path(ROOT,'tables',e$trees_file))
  }
 }
 exports[[id]]<-e
}
write_json(exports,file.path(ROOT,'tables/model_parameters_for_verification.json'),auto_unbox=TRUE,digits=17,pretty=TRUE,na='null')
genes<-GENES;universe<-cache$model$universe
for(e in parse(file.path(SOURCE_ROOT,'application_models','prepare_data.R'),encoding='UTF-8'))if(is.call(e)&&identical(e[[1]],as.name('<-'))&&identical(e[[2]],as.name('extract_features')))eval(e)
targets<-unique(c('equal_contrast','equal_rank',state$selections$candidate,'raw__radial_svm','difference__linear_svm'))
g<-cache$datasets$GSE60542;patients<-intersect(g$pid[g$lab=='PTC'],g$pid[g$lab=='Normal'])
stopifnot(length(patients)==24)
ix<-sapply(patients,function(pid)which(g$pid==pid&g$lab=='PTC'));inx<-sapply(patients,function(pid)which(g$pid==pid&g$lab=='Normal'))
mix<-list();mixrecords<-list()
for(a in c(0,.25,.5,.75,1)){
 G<-log2(a*2^g$G[,ix,drop=FALSE]+(1-a)*2^g$G[,inx,drop=FALSE]);d<-extract_features(G)
 for(id in targets){
  m<-locked[[id]];s<-predict_model(m,d);t<-m$thresholds['youden'];pred<-s>=t
  mix[[length(mix)+1]]<-data.table(candidate=id,ptc_signal_weight=a,n=24,positive=sum(pred),positive_fraction=mean(pred),threshold=t,median_score=median(s))
  mixrecords[[length(mixrecords)+1]]<-data.table(candidate=id,ptc_signal_weight=a,patient=patients,score=s,threshold=t,positive=as.integer(pred))
 }
}
fwrite(rbindlist(mix),file.path(ROOT,'tables/weighted_signal_mixing.csv'));fwrite(rbindlist(mixrecords),file.path(ROOT,'tables/weighted_signal_mixing_records.csv'))
stress<-list();invariant<-list()
for(id in targets){
 m<-locked[[id]];d<-D[dataset=='GSE60542'];base<-predict_model(m,d)
 single<-sapply(seq_len(nrow(d)),function(i)predict_model(m,d[i]));only_cases<-predict_model(m,d[y==1]);only_controls<-predict_model(m,d[y==0])
 invariant[[id]]<-data.table(candidate=id,single_max_difference=max(abs(single-base)),case_only_max_difference=max(abs(only_cases-base[d$y==1])),control_only_max_difference=max(abs(only_controls-base[d$y==0])))
 for(change in c('unchanged','common_plus5','common_times2_plus3','ZDHHC20_minus0.5','ZDHHC20_plus0.5')){
  G<-g$G
  if(change=='common_plus5')G<-G+5
  if(change=='common_times2_plus3')G<-2*G+3
  if(change=='ZDHHC20_minus0.5')G['ZDHHC20',]<-G['ZDHHC20',]-.5
  if(change=='ZDHHC20_plus0.5')G['ZDHHC20',]<-G['ZDHHC20',]+.5
  dx<-extract_features(G);s<-predict_model(m,dx);cc<-counts(d$y,s,m$thresholds['youden'])
  stress[[length(stress)+1]]<-as.data.table(c(list(candidate=id,perturbation=change,score_max_difference=max(abs(s-base))),as.list(cc)))
 }
}
fwrite(rbindlist(stress),file.path(ROOT,'tables/weighted_measurement_stress.csv'));fwrite(rbindlist(invariant),file.path(ROOT,'tables/weighted_batch_invariance.csv'))
stopifnot(all(as.matrix(rbindlist(invariant)[,-1])<1e-12))
# Export a compact fully specified learned contrast on the original input scale.
id<-'difference__linear_svm';m<-locked[[id]];cf<-formula_coefficients(m)
total_up<-sum(cf$coefficients[c('ABHD17C','ZDHHC20')])
coef<-cf$coefficients/total_up;cut<-(m$thresholds['youden']-cf$intercept)/total_up
stopifnot(abs(sum(coef))<1e-12)
fwrite(data.table(gene=names(coef),coefficient=unname(coef),threshold=as.numeric(cut),input_scale='GPL570_gene_log2',model_origin='GSE33630 patient-grouped CV, linear SVM'),file.path(ROOT,'workflow/learned_contrast_frozen.csv'))
raw<-as.matrix(D[,..GENES]);equiv<-as.numeric(raw%*%coef)
stopifnot(identical(equiv>=cut,predict_model(m,D)>=m$thresholds['youden']))
cf<-formula_coefficients(locked[['raw__enet']]);fwrite(data.table(term=c('intercept',names(cf$coefficients)),coefficient=c(cf$intercept,cf$coefficients),threshold=locked[['raw__enet']]$thresholds['youden']),file.path(ROOT,'workflow/learned_logistic_frozen.csv'))
fwrite(data.table(gene=GENES,log2_expression=as.numeric(g$G[GENES,1])),file.path(ROOT,'workflow/example_four_gene_input.csv'))
cat('Frozen predictors exported, single-sample checks passed, and bounded stress/mixing completed.\n')
