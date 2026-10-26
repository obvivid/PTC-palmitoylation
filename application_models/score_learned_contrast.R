# Input is a named numeric vector of deposited GPL570 gene-level log2 expression.
score_learned_contrast<-function(expression,parameter_csv,input_scale='GPL570_gene_log2'){
 p<-read.csv(parameter_csv,stringsAsFactors=FALSE)
 if(!identical(input_scale,unique(p$input_scale)))stop('Unsupported input scale.')
 if(!is.numeric(expression)||is.null(names(expression))||anyDuplicated(names(expression)))stop('Supply uniquely named numeric expression values.')
 if(!all(p$gene%in%names(expression)))stop('All four panel genes are required.')
 x<-expression[p$gene];if(any(!is.finite(x)))stop('Expression values must be finite.')
 score<-sum(x*p$coefficient);t<-unique(p$threshold);stopifnot(length(t)==1)
 data.frame(score=score,threshold=t,positive=score>=t)
}
