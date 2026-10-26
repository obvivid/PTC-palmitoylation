# External gene-expression contrasts with patient-clustered uncertainty.
suppressPackageStartupMessages({library(data.table);library(sandwich)})
DS <- readRDS(file.path(RESULT_DIR,"external_four_gene_inputs.rds"))
SIG <- fread(file.path(RESULT_DIR,"signature_genes.csv"))$gene
rows <- rbindlist(lapply(DS,function(d)rbindlist(lapply(SIG,function(g){
  y <- as.numeric(d$G[g,])
  dat <- data.frame(y=y,grp=factor(ifelse(d$lab=="PTC","Tumor","Normal"),levels=c("Normal","Tumor")),pid=d$pid)
  z <- patient_cluster_coefficient(lm(y~grp,dat),dat$pid)
  data.table(dataset=sub(" .*","",d$tag),gene=g,median_PTC=median(y[d$lab=="PTC"]),
    median_Normal=median(y[d$lab=="Normal"]),direction=ifelse(median(y[d$lab=="PTC"])>median(y[d$lab=="Normal"]),"up","down"),
    p=unname(z["p"]),mean_difference=unname(z["beta"]),ci_lo=unname(z["lo"]),ci_hi=unname(z["hi"]),n_patients=unname(z["patients"]))
}))))
rows[,p_BH:=p.adjust(p,"BH"),by=dataset]
fwrite(rows,file.path(RESULT_DIR,"per_gene_external_direction.csv"))
print(rows)
