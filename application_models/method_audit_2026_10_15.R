# Bounded audit of frozen models and previously generated mixtures; no model search.
ROOT<-file.path(RESULT_DIR,'weighted_single_sample')
source(file.path(SOURCE_ROOT,'application_models','models.R'))
suppressPackageStartupMessages(library(jsonlite))
cache<-readRDS(file.path(RESULT_DIR,'single_sample','analysis_cache.rds'))
state<-readRDS(file.path(ROOT,'workflow/development_locked_models.rds'))
D<-fread(file.path(ROOT,'tables/input_features.csv'))
g<-cache$datasets$GSE60542
ids<-intersect(g$pid[g$lab=='PTC'],g$pid[g$lab=='Normal'])
stopifnot(length(ids)==24L)
ip<-vapply(ids,function(z){ix<-which(g$pid==z&g$lab=='PTC');stopifnot(length(ix)==1);ix},integer(1))
inorm<-vapply(ids,function(z){ix<-which(g$pid==z&g$lab=='Normal');stopifnot(length(ix)==1);ix},integer(1))
stopifnot(!length(intersect(ip,inorm)))
m<-state$models[['difference__linear_svm']]
stopifnot(identical(m$training_datasets,'GSE33630'))
stopifnot(identical(unique(state$training$dataset),'GSE33630'))
X<-features(state$training,'difference')
stopifnot(max(abs(colMeans(X)-m$mean))<1e-12,max(abs(apply(X,2,sd)-m$sd))<1e-12)

records<-fread(file.path(ROOT,'tables/weighted_signal_mixing_records.csv'),colClasses=c(patient='character'))
new<-list();max_delta<-0
for(a in c(0,.25,.5,.75,1)){
  Xmix<-t(log2(a*2^g$G[GENES,ip,drop=FALSE]+(1-a)*2^g$G[GENES,inorm,drop=FALSE]))
  d<-as.data.table(Xmix);d[,baseline_contrast:=(ABHD17C+ZDHHC20-ABHD17B-ZDHHC6)/2]
  if(a==0)stopifnot(max(abs(Xmix-t(g$G[GENES,inorm,drop=FALSE])))<1e-12)
  if(a==1)stopifnot(max(abs(Xmix-t(g$G[GENES,ip,drop=FALSE])))<1e-12)
  for(id in c('equal_contrast','difference__linear_svm')){
    model<-state$models[[id]];s<-predict_model(model,d)
    saved<-records[candidate==id&abs(ptc_signal_weight-a)<1e-12]
    saved<-saved[match(ids,patient)]
    max_delta<-max(max_delta,max(abs(s-saved$score)))
    stopifnot(max(abs(s-saved$score))<1e-10)
    for(op in c('youden','sensitivity95')){
      t<-model$thresholds[[op]]
      new[[length(new)+1]]<-data.table(candidate=id,operating_point=op,ptc_signal_weight=a,n=24,threshold=t,positive=sum(s>=t))
    }
  }
}
fwrite(rbindlist(new),file.path(ROOT,'tables/audit_frozen_operating_points_mixing.csv'))

paired<-rbindlist(lapply(c(0,.25,.5,.75,1),function(a){
  old<-records[candidate=='equal_contrast'&abs(ptc_signal_weight-a)<1e-12]
  alt<-records[candidate=='difference__linear_svm'&abs(ptc_signal_weight-a)<1e-12]
  z<-merge(old[,.(patient,old=positive)],alt[,.(patient,new=positive)],by='patient')
  stopifnot(nrow(z)==24,uniqueN(z$patient)==24)
  gain<-sum(z$old==0&z$new==1);loss<-sum(z$old==1&z$new==0)
  p<-if(gain+loss)binom.test(gain,gain+loss,p=.5,alternative='two.sided')$p.value else 1
  data.table(ptc_signal_weight=a,n=24,old_positive=sum(z$old),new_positive=sum(z$new),gain=gain,loss=loss,paired_difference=(gain-loss)/24,exact_mcnemar_p=p)
}))
fwrite(paired,file.path(ROOT,'tables/audit_mixing_paired_comparison.csv'))

# Check exact whole-profile reuse in these cached inputs; correlation does not prove identity.
train<-cache$datasets$GSE33630$G;test<-g$G
genes<-intersect(rownames(train),rownames(test));A<-train[genes,,drop=FALSE];B<-test[genes,,drop=FALSE]
rho<-cor(A,B)
closest<-data.table(test_record=colnames(B),closest_training_record=colnames(A)[max.col(t(rho),ties.method='first')],max_pearson=apply(rho,2,max))
closest[,max_absolute_difference:=vapply(seq_len(.N),function(j)max(abs(A[,closest_training_record[j]]-B[,test_record[j]])),numeric(1))]
fwrite(closest,file.path(ROOT,'tables/audit_training_test_profile_similarity.csv'))
stopifnot(!any(closest$max_absolute_difference<1e-12))

write_json(list(passed=TRUE,n_mixture_patients=24,training_dataset='GSE33630',mixture_dataset='GSE60542',training_scaler_matches_development=TRUE,mixture_recalculation_max_difference=max_delta,no_exact_training_profile_copies_detected=TRUE,max_cross_dataset_pearson=max(rho),scope='Frozen-model audit; paired calls and exact McNemar tests are exploratory. Whole-profile comparison is not a definitive patient identity audit.'),file.path(ROOT,'method_audit_2026-10-15.json'),auto_unbox=TRUE,pretty=TRUE,digits=16)
print(paired)
cat('Maximum cross-dataset correlation:',max(rho),'\n')
cat('Frozen method audit completed; no model parameters changed.\n')
