read_frozen_model <- function(directory) {
  t <- read.csv(file.path(directory,'frozen_thresholds.csv'),stringsAsFactors=FALSE)
  th <- lapply(split(t,t$method),function(x)as.list(setNames(x$threshold,x$operating_point)))
  list(version='1.0.0',date='2026-10-13',input_scale='GPL570_gene_log2',
       up=c('ABHD17C','ZDHHC20'),down=c('ABHD17B','ZDHHC6'),
       universe=readLines(file.path(directory,'frozen_gene_universe.txt')),thresholds=th)
}
