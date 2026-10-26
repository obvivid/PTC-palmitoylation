suppressPackageStartupMessages({library(data.table);library(glmnet);library(e1071);library(randomForest);library(pROC)})
GENES<-c('ABHD17C','ZDHHC20','ABHD17B','ZDHHC6')
PAIRS<-combn(GENES,2,simplify=FALSE)
features<-function(d,encoding){
 X<-as.matrix(d[,..GENES])
 if(encoding=='raw')return(X)
 if(encoding=='difference'){
  out<-sweep(X[,1:3,drop=FALSE],1,X[,4],'-')
  colnames(out)<-paste0(GENES[1:3],'_minus_',GENES[4]);return(out)
 }
 if(encoding=='rank')return(as.matrix(d[,paste0('pct_',GENES),with=FALSE]))
 if(encoding=='order'){
  out<-sapply(PAIRS,function(pair)as.numeric(X[,pair[1]]>X[,pair[2]]))
  colnames(out)<-vapply(PAIRS,function(pair)paste(pair,collapse='_gt_'),character(1));return(out)
 }
 stop('Unknown encoding')
}
make_candidates<-function(){
 out<-list(list(id='equal_contrast',family='baseline',encoding='contrast',complexity=0),list(id='equal_rank',family='baseline',encoding='rank_baseline',complexity=0))
 for(enc in c('raw','difference','rank','order'))for(fam in c('ridge','lasso','enet','linear_svm','radial_svm','shrink_lda','random_forest')){
  out[[length(out)+1]]<-list(id=paste(enc,fam,sep='__'),family=fam,encoding=enc,complexity=if(fam%in%c('ridge','lasso','enet'))1 else if(fam=='shrink_lda')2 else if(fam=='linear_svm')3 else if(fam=='radial_svm')4 else 5)
 }
 names(out)<-vapply(out,`[[`,character(1),'id');out
}
grid<-function(cand,p){
 if(cand$family=='baseline')return(list(list()))
 if(cand$family%in%c('ridge','lasso','enet'))return(lapply(c(1,.3,.1,.03,.01,.003,.001),function(l)list(lambda=l,alpha=switch(cand$family,ridge=0,lasso=1,enet=.5))))
 if(cand$family=='linear_svm')return(lapply(c(.01,.1,1,10,100),function(C)list(cost=C)))
 if(cand$family=='radial_svm')return(unname(lapply(seq_len(9),function(i){g<-expand.grid(cost=c(.1,1,10),gamma=c(.1,.5,2));as.list(g[i,])})))
 if(cand$family=='shrink_lda')return(lapply(c(1,.5,.25,0),function(a)list(shrink=a)))
 if(cand$family=='random_forest'){
  g<-unique(expand.grid(nodesize=c(15,8,3),mtry=unique(c(1,ceiling(p/2),p))))
  return(lapply(seq_len(nrow(g)),function(i)as.list(g[i,])))
 }
 stop('Unknown family')
}
counts<-function(y,score,t){
 pred<-score>=t;tp<-sum(pred&y==1);fn<-sum(!pred&y==1);tn<-sum(!pred&y==0);fp<-sum(pred&y==0)
 se<-tp/(tp+fn);sp<-tn/(tn+fp)
 c(tp=tp,fn=fn,tn=tn,fp=fp,sensitivity=se,specificity=sp,accuracy=(tp+tn)/length(y),balanced_accuracy=(se+sp)/2)
}
auc_fixed<-function(y,score){if(length(unique(score))==1)return(.5);as.numeric(pROC::auc(pROC::roc(y,score,levels=c(0,1),direction='<',quiet=TRUE)))}
threshold<-function(y,score,point='youden'){
 u<-sort(unique(score));gap<-max(1,abs(range(u)),diff(range(u)))*1e-7
 ts<-c(min(u)-gap,if(length(u)>1)head(u,-1)+(tail(u,-1)-head(u,-1))/2,max(u)+gap)
 m<-rbindlist(lapply(ts,function(t)as.list(c(threshold=t,counts(y,score,t)))))
 if(point=='youden'){
  m<-m[abs(balanced_accuracy-max(balanced_accuracy))<1e-12];setorder(m,-sensitivity,-threshold)
 }else{
  m<-m[sensitivity>=.95-1e-12];m<-m[abs(specificity-max(specificity))<1e-12];setorder(m,-threshold)
 }
 m$threshold[1]
}
fit_model<-function(d,cand,parameter,seed){
 if(cand$family=='baseline')return(list(candidate=cand,parameter=parameter))
 X<-features(d,cand$encoding);y<-d$y
 mu<-colMeans(X);sdv<-apply(X,2,sd);keep<-which(is.finite(sdv)&sdv>1e-12)
 obj<-list(candidate=cand,parameter=parameter,feature_names=colnames(X),keep=keep,mean=mu,sd=sdv)
 if(!length(keep)){
  obj$type<-'constant';obj$constant<-qlogis((sum(y)+.5)/(length(y)+1));return(obj)
 }
 Z<-sweep(sweep(X[,keep,drop=FALSE],2,mu[keep],'-'),2,sdv[keep],'/')
 set.seed(seed)
 fam<-cand$family
 if(fam%in%c('ridge','lasso','enet')){
  # glmnet requires at least two columns, including for effectively one-feature pairs.
  zz<-if(ncol(Z)==1)cbind(Z,zero_padding=0)else Z
  f<-glmnet(zz,y,family='binomial',alpha=parameter$alpha,lambda=parameter$lambda,standardize=FALSE,thresh=1e-10,maxit=100000)
  b<-as.numeric(coef(f,s=parameter$lambda));obj$type<-'linear';obj$intercept<-b[1];obj$weights<-b[seq_len(ncol(Z))+1]
 }else if(fam=='shrink_lda'){
  m0<-colMeans(Z[y==0,,drop=FALSE]);m1<-colMeans(Z[y==1,,drop=FALSE])
  covmat<-function(a){if(ncol(a)==1)matrix(var(a[,1]),1,1)else cov(a)}
  pooled<-((sum(y==0)-1)*covmat(Z[y==0,,drop=FALSE])+(sum(y==1)-1)*covmat(Z[y==1,,drop=FALSE]))/(length(y)-2)
  Sigma<-(1-parameter$shrink)*pooled+parameter$shrink*diag(diag(pooled),nrow=ncol(Z))+diag(1e-6,ncol(Z))
  obj$type<-'linear';obj$weights<-as.numeric(solve(Sigma,m1-m0));obj$intercept<--.5*sum((m0+m1)*obj$weights)
 }else if(fam%in%c('linear_svm','radial_svm')){
  fy<-factor(y,levels=c(0,1));f<-svm(x=Z,y=fy,type='C-classification',kernel=if(fam=='linear_svm')'linear'else'radial',cost=parameter$cost,gamma=if(fam=='radial_svm')parameter$gamma else 1/ncol(Z),scale=FALSE,probability=FALSE)
  pr<-predict(f,Z,decision.values=TRUE);v<-as.numeric(attr(pr,'decision.values'))
  direction<-if(mean(v[y==1])>=mean(v[y==0]))1 else -1
  obj$svm_direction<-direction
  if(fam=='linear_svm'){
   obj$type<-'linear';obj$weights<-as.numeric(t(f$coefs)%*%f$SV)*direction;obj$intercept<--f$rho*direction
  }else{obj$type<-'svm';obj$fitted<-f}
 }else if(fam=='random_forest'){
  obj$type<-'forest';obj$fitted<-randomForest(x=Z,y=factor(y,levels=c(0,1)),ntree=300,mtry=min(parameter$mtry,ncol(Z)),nodesize=parameter$nodesize)
 }else stop('Unknown model')
 obj
}
predict_model<-function(obj,d){
 c<-obj$candidate
 if(c$family=='baseline')return(d[[if(c$encoding=='contrast')'baseline_contrast'else'baseline_rank']])
 if(obj$type=='constant')return(rep(obj$constant,nrow(d)))
 X<-features(d,c$encoding)[,obj$keep,drop=FALSE]
 Z<-sweep(sweep(X,2,obj$mean[obj$keep],'-'),2,obj$sd[obj$keep],'/')
 if(obj$type=='linear')return(as.numeric(obj$intercept+Z%*%obj$weights))
 if(obj$type=='svm')return(as.numeric(attr(predict(obj$fitted,Z,decision.values=TRUE),'decision.values'))*obj$svm_direction)
 if(obj$type=='forest')return(as.numeric(predict(obj$fitted,Z,type='prob')[,'1']))
 stop('Unknown stored model')
}
patient_folds<-function(d,K,seed){
 set.seed(seed);p<-d[,.(dataset=dataset[1],has_normal=any(y==0)),by=patient_key];p[,fold:=0L]
 for(ds in unique(p$dataset))for(h in c(TRUE,FALSE)){
  ids<-which(p$dataset==ds&p$has_normal==h)
  if(length(ids))p$fold[ids[sample.int(length(ids))]]<-rep(seq_len(K),length.out=length(ids))
 }
 folds<-p$fold[match(d$patient_key,p$patient_key)]
 stopifnot(all(folds%in%seq_len(K)),all(vapply(split(folds,d$patient_key),function(x)length(unique(x))==1,logical(1))))
 for(k in seq_len(K))stopifnot(length(unique(d$y[folds==k]))==2,length(unique(d$y[folds!=k]))==2)
 folds
}
tune_candidate<-function(d,cand,fold,seed){
 pp<-if(cand$family=='baseline')1 else ncol(features(d,cand$encoding));params<-grid(cand,pp)
 pred<-matrix(NA_real_,nrow(d),length(params));failures<-character()
 for(k in sort(unique(fold)))for(j in seq_along(params)){
  f<-tryCatch(fit_model(d[fold!=k],cand,params[[j]],seed+k*100+j),error=function(e)e)
  if(inherits(f,'error')){failures<-c(failures,paste(cand$id,k,j,conditionMessage(f),sep=':'));next}
  pred[fold==k,j]<-predict_model(f,d[fold==k])
 }
 measures<-rbindlist(lapply(seq_along(params),function(j){
  s<-pred[,j]
  if(any(!is.finite(s)))return(data.table(parameter=j,balanced_accuracy=-Inf,auc=-Inf,sensitivity=-Inf))
  t<-threshold(d$y,s);cc<-counts(d$y,s,t)
  data.table(parameter=j,balanced_accuracy=cc['balanced_accuracy'],auc=auc_fixed(d$y,s),sensitivity=cc['sensitivity'],threshold=t)
 }),fill=TRUE)
 setorder(measures,-balanced_accuracy,-auc,-sensitivity,parameter)
 best<-measures[1];if(!is.finite(best$balanced_accuracy))stop(paste('Candidate failed',cand$id,paste(failures,collapse='; ')))
 j<-best$parameter;score<-pred[,j]
 list(candidate=cand,parameter=params[[j]],parameter_index=j,selection=best,oof=score,thresholds=c(youden=best$threshold,sensitivity95=threshold(d$y,score,'sensitivity95')),failures=failures,grid_metrics=measures)
}
select_candidate<-function(tuned,linear_only=FALSE,learning_only=FALSE){
 allowed<-vapply(tuned,function(x){c<-x$candidate;(!linear_only||c$family%in%c('ridge','lasso','enet','linear_svm','shrink_lda'))&&(!learning_only||c$family!='baseline')},logical(1))
 a<-tuned[allowed]
 scores<-rbindlist(lapply(a,function(x)data.table(id=x$candidate$id,ba=x$selection$balanced_accuracy,auc=x$selection$auc,se=x$selection$sensitivity,complexity=x$candidate$complexity,parameter=x$parameter_index)))
 setorder(scores,-ba,-auc,-se,complexity,parameter,id)
 tuned[[scores$id[1]]]
}
formula_coefficients<-function(obj){
 if(obj$candidate$family=='baseline'){
  if(obj$candidate$encoding!='contrast')return(NULL)
  return(list(intercept=0,coefficients=setNames(c(.5,.5,-.5,-.5),GENES),input='log2_expression'))
 }
 if(obj$type!='linear')return(NULL)
 b<-rep(0,length(obj$feature_names));b[obj$keep]<-obj$weights/obj$sd[obj$keep]
 b0<-obj$intercept-sum(b[obj$keep]*obj$mean[obj$keep])
 names(b)<-obj$feature_names
 if(obj$candidate$encoding=='difference'){b<-setNames(c(b,-sum(b)),GENES);input<-'log2_expression'}
 else if(obj$candidate$encoding=='raw'){input<-'log2_expression'}
 else input<-obj$candidate$encoding
 list(intercept=b0,coefficients=b,input=input)
}
