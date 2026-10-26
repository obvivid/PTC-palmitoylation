# Exploratory survival associations of the frozen diagnostic expression score.
# Run through run_analysis.R after step 08. No threshold is optimised.
suppressPackageStartupMessages({library(data.table); library(survival)})
sc <- fread(file.path(RESULT_DIR,"clinical_immune_fixed.csv"))
sc[, patient_id := substr(sample,1,12)]
stopifnot(!anyDuplicated(sc$patient_id))
cl <- fread(file.path(PROJECT_ROOT,"data/tcga/data_clinical_patient.txt"),
            skip="PATIENT_ID",na.strings=c("","NA"))
stopifnot(!anyDuplicated(cl$PATIENT_ID))
d <- merge(sc[,.(patient_id,score=risk)],cl,by.x="patient_id",by.y="PATIENT_ID")
endpoints <- list(OS=c("OS_STATUS","OS_MONTHS"),DSS=c("DSS_STATUS","DSS_MONTHS"),
                  DFS=c("DFS_STATUS","DFS_MONTHS"),PFS=c("PFS_STATUS","PFS_MONTHS"))
out <- rbindlist(lapply(names(endpoints),function(nm){
 status <- as.character(d[[endpoints[[nm]][1]]])
 event <- suppressWarnings(as.integer(sub(":.*$","",status)))
 time <- suppressWarnings(as.numeric(d[[endpoints[[nm]][2]]]))
 keep <- is.finite(d$score)&is.finite(time)&time>0&!is.na(event)&event %in% c(0L,1L)
 z <- data.frame(time=time[keep],event=event[keep],score=d$score[keep])
 cut <- median(z$score)
 z$group <- factor(z$score>cut,levels=c(FALSE,TRUE),labels=c("Low","High"))
 stopifnot(nlevels(droplevels(z$group))==2L)
 fit <- coxph(Surv(time,event)~group,data=z,ties="efron")
 s <- summary(fit); lr <- survdiff(Surv(time,event)~group,data=z,rho=0)
 ph <- cox.zph(fit)$table["GLOBAL","p"]
 data.table(endpoint=nm,n=nrow(z),events=sum(z$event),n_low=sum(z$group=="Low"),
            n_high=sum(z$group=="High"),events_low=sum(z$event[z$group=="Low"]),
            events_high=sum(z$event[z$group=="High"]),score_median=cut,
            HR_high_vs_low=s$conf.int[1,"exp(coef)"],CI_low=s$conf.int[1,"lower .95"],
            CI_high=s$conf.int[1,"upper .95"],logrank_p=pchisq(lr$chisq,1,lower.tail=FALSE),
            cox_wald_p=s$coefficients[1,"Pr(>|z|)"],PH_global_p=ph)
}))
out[,logrank_BH_adjusted_P:=p.adjust(logrank_p,"BH")]
fwrite(out,file.path(RESULT_DIR,"survival_endpoints.csv"))
print(out)
# Scores from step 08 omit the fitted intercept. Adding a constant leaves the
# endpoint-specific median split, hazard ratios and tests unchanged.

flow <- rbindlist(lapply(names(endpoints),function(ep){
  status<-as.character(d[[endpoints[[ep]][1]]]);ev<-suppressWarnings(as.integer(sub(":.*$","",status)))
  tt<-suppressWarnings(as.numeric(d[[endpoints[[ep]][2]]]))
  vs<-is.finite(d$score);ve<-!is.na(ev)&ev %in% 0:1;vt<-is.finite(tt);vp<-tt>0
  keep<-vs&ve&vt&!is.na(vp)&vp
  data.table(endpoint=ep,discovery_tumours=nrow(sc),unmatched_clinical_records=nrow(sc)-nrow(d),available_scores=sum(vs),
    excluded_missing_or_invalid_status=sum(vs&!ve),excluded_missing_time=sum(vs&ve&!vt),
    excluded_nonpositive_time=sum(vs&ve&vt&!is.na(vp)&!vp),analysed=sum(keep),events=sum(ev[keep]),
    median_observed_months=median(tt[keep]),q1_months=quantile(tt[keep],.25),q3_months=quantile(tt[keep],.75),
    status_field=endpoints[[ep]][1],time_field=endpoints[[ep]][2],event_labels=paste(unique(status[keep]),collapse=";"))
}))
fwrite(flow,file.path(RESULT_DIR,"survival_endpoint_flow.csv"))
