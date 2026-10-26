suppressPackageStartupMessages(library(data.table))
ROOT <- file.path(RESULT_DIR,'weighted_single_sample')
PREVIOUS <- file.path(RESULT_DIR,'single_sample')
cache <- readRDS(file.path(PREVIOUS,'analysis_cache.rds'))
genes <- c('ABHD17C','ZDHHC20','ABHD17B','ZDHHC6')
universe <- cache$model$universe
stopifnot(length(universe)==20845)
extract_features <- function(G){
 stopifnot(all(universe%in%rownames(G)),all(genes%in%rownames(G)))
 X <- t(G[genes,,drop=FALSE])
 ranks <- apply(G[universe,,drop=FALSE],2,rank,ties.method='average')
 rownames(ranks)<-universe
 R <- t((ranks[genes,,drop=FALSE]-1)/(length(universe)-1))
 colnames(R)<-paste0('pct_',genes)
 out <- as.data.table(cbind(X,R))
 out[,baseline_contrast:=(ABHD17C+ZDHHC20-ABHD17B-ZDHHC6)/2]
 out[,baseline_rank:=(rowMeans(t(ranks[c('ABHD17C','ZDHHC20'),,drop=FALSE]))-rowMeans(t(ranks[c('ABHD17B','ZDHHC6'),,drop=FALSE])))/(length(universe)-2)]
 out
}
rows <- lapply(names(cache$datasets),function(nm){
 d<-cache$datasets[[nm]]
 out<-extract_features(d$G)
 out[,`:=`(dataset=nm,record=colnames(d$G),pid=as.character(d$pid),patient_key=paste(nm,d$pid,sep='::'),truth=d$lab,y=as.integer(d$lab==d$positive_label))]
 setcolorder(out,c('dataset','record','pid','patient_key','truth','y',genes))
 out
})
data<-rbindlist(rows)
fwrite(data,file.path(ROOT,'tables/input_features.csv'))
manifest<-data[,.(n=.N,positive=sum(y),negative=sum(!y),patients=uniqueN(patient_key)),by=dataset]
fwrite(manifest,file.path(ROOT,'tables/input_manifest.csv'))
source(file.path(SOURCE_ROOT,'single_sample','single_sample_score.R'))
for(nm in c('GSE33630','GSE35570','GSE29265','GSE60542'))for(method in c('contrast','rank')){
 x<-data[dataset==nm][[paste0('baseline_',method)]]
 expected<-score_sample_matrix(cache$datasets[[nm]]$G,cache$model,method)$score
 stopifnot(max(abs(x-expected))<1e-12)
}
writeLines(universe,file.path(ROOT,'workflow/frozen_gene_universe.txt'))
writeLines(capture.output(sessionInfo()),file.path(ROOT,'session_info_prepare.txt'))
cat('Prepared verified features and patient groups for ',nrow(data),' records.\n')
