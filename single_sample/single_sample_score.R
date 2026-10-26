# Single-sample scoring of a frozen four-gene thyroid expression panel.
# Input: named gene-level log2 expression from GPL570, with mapping defined by model.
# This scoring module never estimates parameters from a target cohort.
score_single_sample <- function(expression, model, method=c('rank','contrast'),
                                operating_point=c('youden','sensitivity95'),
                                input_scale='GPL570_gene_log2') {
  method<-match.arg(method); operating_point<-match.arg(operating_point)
  if(!identical(input_scale,model$input_scale))stop('Unsupported input scale or assay; calibration has not been established.')
  if(!is.numeric(expression)||is.null(names(expression)))stop('Supply a named numeric gene expression vector.')
  if(anyNA(names(expression))||any(names(expression)=='')||anyDuplicated(names(expression)))stop('Gene identifiers must be nonempty and unique.')
  required<-if(method=='rank')model$universe else c(model$up,model$down)
  missing<-setdiff(required,names(expression))
  if(length(missing))stop(sprintf('Missing %d required genes; background substitution or imputation is not allowed.',length(missing)))
  x<-expression[required]
  if(any(!is.finite(x)))stop('Required expression values must all be finite.')
  if(method=='rank') {
    ranks<-rank(x,ties.method='average');names(ranks)<-required
    score<-(mean(ranks[model$up])-mean(ranks[model$down]))/(length(required)-2)
  }else score<-(sum(x[model$up])-sum(x[model$down]))/2
  threshold<-model$thresholds[[method]][[operating_point]]
  if(length(threshold)!=1||!is.finite(threshold))stop('Frozen threshold is absent or invalid.')
  data.frame(method=method,operating_point=operating_point,score=as.numeric(score),threshold=threshold,
             positive=as.numeric(score)>=threshold,stringsAsFactors=FALSE)
}

score_sample_matrix <- function(expression,model,method='rank',operating_point='youden',
                                input_scale='GPL570_gene_log2'){
  if(!is.matrix(expression)||is.null(rownames(expression)))stop('Supply a numeric gene-by-sample matrix with gene names.')
  results<-lapply(seq_len(ncol(expression)),function(j){
    x<-expression[,j];names(x)<-rownames(expression)
    r<-score_single_sample(x,model,method,operating_point,input_scale)
    r$record<-if(is.null(colnames(expression)))as.character(j)else colnames(expression)[j]
    r
  })
  do.call(rbind,results)
}
