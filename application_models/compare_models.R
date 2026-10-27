ROOT<-file.path(RESULT_DIR,'weighted_single_sample')
source(file.path(SOURCE_ROOT,'application_models','models.R'))
D<-fread(file.path(ROOT,'tables/input_features.csv'))
candidates<-make_candidates();stopifnot(length(candidates)==30)
fwrite(rbindlist(lapply(candidates,function(x)as.data.table(x))),file.path(ROOT,'tables/candidate_registry.csv'))
dev<-D[dataset=='GSE33630']
repeats<-list();outer_records<-list();outer_folds<-list();choices<-list();failures<-character()
for(rep in 1:5){
 fid<-patient_folds(dev,5,20261014+rep)
 for(k in 1:5){
  train<-dev[fid!=k];test<-dev[fid==k]
  inner<-patient_folds(train,4,20261014+100*rep+k)
  tuned<-lapply(candidates,function(c)tune_candidate(train,c,inner,20261014+rep*10000+k*100));names(tuned)<-names(candidates)
  winner<-select_candidate(tuned);choices[[paste(rep,k)]]<-data.table(repeat_id=rep,fold=k,selected_candidate=winner$candidate$id,inner_ba=winner$selection$balanced_accuracy,inner_auc=winner$selection$auc)
  for(id in names(tuned)){
   t<-tuned[[id]];f<-fit_model(train,t$candidate,t$parameter,20261014+rep*10000+k*100+99)
   s<-predict_model(f,test);failures<-c(failures,t$failures)
   for(op in c('youden','sensitivity95')){
    threshold_value<-t$thresholds[[op]]
    row<-data.table(repeat_id=rep,fold=k,candidate=id,operating_point=op,dataset=test$dataset,record=test$record,patient_key=test$patient_key,y=test$y,score=s,threshold=threshold_value,positive=as.integer(s>=threshold_value))
    outer_records[[length(outer_records)+1]]<-row
    cc<-counts(test$y,s,threshold_value)
    outer_folds[[length(outer_folds)+1]]<-as.data.table(c(list(repeat_id=rep,fold=k,candidate=id,operating_point=op,auc=auc_fixed(test$y,s)),as.list(cc)))
    if(id==winner$candidate$id){rr<-copy(row);rr[,`:=`(candidate='all_model_selection',selected_candidate=id)];outer_records[[length(outer_records)+1]]<-rr;ff<-copy(outer_folds[[length(outer_folds)]]);ff[,candidate:='all_model_selection'];outer_folds[[length(outer_folds)+1]]<-ff}
   }
  }
  fwrite(rbindlist(choices),file.path(ROOT,'tables/outer_model_choices.csv'))
  cat(sprintf('Nested development: repeat %d/5, fold %d/5 completed; selected %s\n',rep,k,winner$candidate$id));flush.console()
 }
}
records<-rbindlist(outer_records,fill=TRUE);folds<-rbindlist(outer_folds)
fwrite(records,file.path(ROOT,'tables/nested_oof_records.csv'));fwrite(folds,file.path(ROOT,'tables/nested_fold_results.csv'))
rr<-records[,{
 c<-counts(y,positive,.5)
 a<-if(candidate[1]=='all_model_selection')NA_real_ else auc_fixed(y,score)
 c(as.list(c),list(auc_pooled_oof=a))
},by=.(repeat_id,candidate,operating_point)]
aa<-folds[,.(auc_mean_folds=mean(auc)),by=.(repeat_id,candidate,operating_point)]
rr<-merge(rr,aa,by=c('repeat_id','candidate','operating_point'))
fwrite(rr,file.path(ROOT,'tables/nested_repeat_results.csv'))
cv<-rr[,lapply(.SD,function(x)mean(x,na.rm=TRUE)),by=.(candidate,operating_point),.SDcols=c('sensitivity','specificity','accuracy','balanced_accuracy','auc_pooled_oof','auc_mean_folds')]
cv[,accuracy_sd:=rr[,sd(accuracy),by=.(candidate,operating_point)]$V1]
cv[,ba_sd:=rr[,sd(balanced_accuracy),by=.(candidate,operating_point)]$V1]
fwrite(cv,file.path(ROOT,'tables/nested_summary.csv'))
cat('Development-only tuning and locking all candidates\n');flush.console()
final_fold<-patient_folds(dev,4,20261014+999)
final_tuned<-lapply(candidates,function(c)tune_candidate(dev,c,final_fold,20261014+99900));names(final_tuned)<-names(candidates)
locked<-lapply(final_tuned,function(t){m<-fit_model(dev,t$candidate,t$parameter,20261014+888);m$thresholds<-t$thresholds;m$training_datasets<-'GSE33630';m})
names(locked)<-names(candidates)
winner<-select_candidate(final_tuned)
linear<-select_candidate(final_tuned,linear_only=TRUE,learning_only=TRUE)
learning<-select_candidate(final_tuned,learning_only=TRUE)
selected<-data.table(role=c('development_selected','best_learned_linear','best_learned_overall'),candidate=c(winner$candidate$id,linear$candidate$id,learning$candidate$id))
fwrite(selected,file.path(ROOT,'tables/locked_selections.csv'))
params<-rbindlist(lapply(final_tuned,function(t)data.table(candidate=t$candidate$id,parameter=paste(names(t$parameter),unlist(t$parameter),sep='=',collapse=';'),inner_ba=t$selection$balanced_accuracy,inner_auc=t$selection$auc,inner_sensitivity=t$selection$sensitivity,threshold_youden=t$thresholds['youden'],threshold_sensitivity95=t$thresholds['sensitivity95'])))
fwrite(params,file.path(ROOT,'tables/locked_parameters.csv'))
saveRDS(list(models=locked,selections=selected,training=dev,candidates=candidates),file.path(ROOT,'workflow/development_locked_models.rds'))
for(id in names(locked)){
 coef<-formula_coefficients(locked[[id]])
 if(!is.null(coef))fwrite(data.table(term=c('intercept',names(coef$coefficients)),coefficient=c(coef$intercept,coef$coefficients),input=coef$input),file.path(ROOT,'workflow',paste0(id,'_coefficients.csv')))
}
wilson<-function(k,n){z<-qnorm(.975);p<-k/n;v<-1+z*z/n;c((p+z*z/(2*n)-z*sqrt(p*(1-p)/n+z*z/(4*n*n)))/v,(p+z*z/(2*n)+z*sqrt(p*(1-p)/n+z*z/(4*n*n)))/v)}
evaluation<-function(d,s,t,id,op,stage){
 cc<-counts(d$y,s,t);se<-wilson(cc['tp'],cc['tp']+cc['fn']);sp<-wilson(cc['tn'],cc['tn']+cc['fp'])
 as.data.table(c(list(stage=stage,dataset=d$dataset[1],candidate=id,operating_point=op,n=nrow(d),n_positive=sum(d$y),n_negative=sum(!d$y),threshold=t,auc=auc_fixed(d$y,s)),as.list(cc),list(sensitivity_lo=se[1],sensitivity_hi=se[2],specificity_lo=sp[1],specificity_hi=sp[2])))
}
transfer<-list();transfer_records<-list()
for(nm in unique(D$dataset))for(id in names(locked)){
 d<-D[dataset==nm];m<-locked[[id]];s<-predict_model(m,d)
 for(op in names(m$thresholds)){
  t<-m$thresholds[[op]]
  transfer[[length(transfer)+1]]<-evaluation(d,s,t,id,op,if(nm=='GSE33630')'development_apparent'else'fixed_transfer')
  transfer_records[[length(transfer_records)+1]]<-data.table(dataset=nm,candidate=id,operating_point=op,record=d$record,patient_key=d$patient_key,y=d$y,score=s,threshold=t,positive=as.integer(s>=t))
 }
}
fwrite(rbindlist(transfer),file.path(ROOT,'tables/fixed_transfer_results.csv'));fwrite(rbindlist(transfer_records),file.path(ROOT,'tables/fixed_transfer_records.csv'))
cat('Frozen single-development comparison complete. Selected ',winner$candidate$id,'; learned linear ',linear$candidate$id,'\n');flush.console()
primary<-c('GSE33630','GSE35570','GSE60542');loco<-list();loco_records<-list();loco_choices<-list()
for(held in primary){
 train<-D[dataset%in%setdiff(primary,held)];test<-D[dataset==held]
 inner<-patient_folds(train,4,20261014+match(held,primary)*77)
 tt<-lapply(candidates,function(c)tune_candidate(train,c,inner,20261014+match(held,primary)*777));names(tt)<-names(candidates)
 best<-select_candidate(tt)
 loco_choices[[held]]<-data.table(held_dataset=held,selected_candidate=best$candidate$id,inner_ba=best$selection$balanced_accuracy,inner_auc=best$selection$auc)
 for(id in names(tt)){
  t<-tt[[id]];m<-fit_model(train,t$candidate,t$parameter,20261014+666);s<-predict_model(m,test)
  for(op in names(t$thresholds)){
   th<-t$thresholds[[op]]
   row<-evaluation(test,s,th,id,op,'leave_cohort_out');loco[[length(loco)+1]]<-row
   record<-data.table(dataset=held,candidate=id,operating_point=op,record=test$record,patient_key=test$patient_key,y=test$y,score=s,threshold=th,positive=as.integer(s>=th));loco_records[[length(loco_records)+1]]<-record
   if(id==best$candidate$id){r<-copy(row);r[,candidate:='all_model_selection'];loco[[length(loco)+1]]<-r;pr<-copy(record);pr[,`:=`(candidate='all_model_selection',selected_candidate=id)];loco_records[[length(loco_records)+1]]<-pr}
  }
 }
 cat('Leave-cohort-out evaluation complete: ',held,' selected ',best$candidate$id,'\n');flush.console()
}
fwrite(rbindlist(loco),file.path(ROOT,'tables/leave_cohort_out_results.csv'));fwrite(rbindlist(loco_records,fill=TRUE),file.path(ROOT,'tables/leave_cohort_out_records.csv'));fwrite(rbindlist(loco_choices),file.path(ROOT,'tables/leave_cohort_out_choices.csv'))
# An optional consolidated development model. Its training cohorts are not validation.
pool<-D[dataset%in%primary];inner<-patient_folds(pool,4,20261014+2026)
tt<-lapply(candidates,function(c)tune_candidate(pool,c,inner,20261014+22026));names(tt)<-names(candidates)
best<-select_candidate(tt);bestlin<-select_candidate(tt,linear_only=TRUE,learning_only=TRUE)
poolids<-unique(c(best$candidate$id,bestlin$candidate$id,'equal_contrast'))
poolmodels<-lapply(poolids,function(id){t<-tt[[id]];m<-fit_model(pool,t$candidate,t$parameter,20261014+77);m$thresholds<-t$thresholds;m$training_datasets<-primary;m});names(poolmodels)<-poolids
saveRDS(list(models=poolmodels,selected=best$candidate$id,best_linear=bestlin$candidate$id,training_datasets=primary),file.path(ROOT,'workflow/pooled_development_models.rds'))
fwrite(data.table(role=c('pooled_selected','pooled_best_linear'),candidate=c(best$candidate$id,bestlin$candidate$id)),file.path(ROOT,'tables/pooled_selections.csv'))
for(id in names(poolmodels)){
 cf<-formula_coefficients(poolmodels[[id]])
 if(!is.null(cf))fwrite(data.table(term=c('intercept',names(cf$coefficients)),coefficient=c(cf$intercept,cf$coefficients),input=cf$input),file.path(ROOT,'workflow',paste0('pooled_',id,'_coefficients.csv')))
}
writeLines(unique(failures),file.path(ROOT,'tables/fit_failures.txt'))
writeLines(capture.output(sessionInfo()),file.path(ROOT,'session_info_models.txt'))
cat('All 30 candidate comparisons completed.\n')
