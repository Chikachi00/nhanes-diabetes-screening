# Reconstruct the submitted F_ext analysis from local public NHANES files.
# No downloads, no original-file mutation. Arguments: workspace, raw, output, source, original-tables.
args <- commandArgs(trailingOnly=TRUE)
ws <- normalizePath(if(length(args)) args[1] else getwd(), winslash="/")
.libPaths(c(file.path(ws,"outputs/nhanes_nlr_depression_start_2026-10-20/runtime/library"),.libPaths()))
out <- if(length(args)>=3) args[3] else file.path(ws,"outputs/diabetes_manuscript_revision_20261104_v2/fext")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
src <- if(length(args)>=4) args[4] else file.path(ws,"outputs/diabetes_manuscript_revision_20261104/supplement/code/denoise_features_v5.R")
PROJECT_DIR <- out
DATA_DIR <- if(length(args)>=2) args[2] else file.path(ws,"outputs/diabetes_manuscript_revision_20261104_v2/analysis/raw")
OUT_DIR <- out
exprs <- parse(src,encoding="UTF-8")
# Evaluate supplied preprocessing and original feature-selection/fitting only.
# Stop before survey inference/plotting; all selected constants are preserved.
for(e in exprs){
  lhs <- if(is.call(e)&&identical(e[[1]],as.name("<-"))&&is.symbol(e[[2]])) as.character(e[[2]]) else ""
  if(lhs=="des") break
  if(lhs %in% c("PROJECT_DIR","DATA_DIR","OUT_DIR")) next
  if(is.call(e)&&identical(e[[1]],as.name("dir.create"))) next
  eval(e)
}
tabdir <- if(length(args)>=5) args[5] else file.path(ws,"outputs/diabetes_manuscript_revision_20261104/supplement/tables")
orig <- read.csv(file.path(tabdir,"S12_information_coefficients_by_model.csv"))
ck <- merge(coef_tab,orig,by=c("model","term"),suffixes=c("_reconstructed","_submitted"))
ck$difference <- ck$coef_reconstructed-ck$coef_submitted
write.csv(ck,file.path(out,"original_coefficient_reproduction.csv"),row.names=FALSE)
stopifnot(nrow(ck)==nrow(orig), max(abs(ck$difference))<1e-8)
write.csv(thr_tab,file.path(out,"original_thresholds_reproduced.csv"),row.names=FALSE)
full$period <- ifelse(full$dev,"development",ifelse(full$eval,"evaluation",NA))
sample <- full[full$in_sample,]
# All analytic participants are adults with DIQ010=2. Do not classify missing
# values as valid skips without documented questionnaire routing evidence.
status <- function(x) ifelse(is.na(x),"missing_unresolved",ifelse(x==1,"yes",ifelse(x==2,"no",ifelse(x==7,"refused",ifelse(x==9,"dont_know","other_code")))))
responses <- list()
for(pe in unique(sample$period)) for(ac in c("overall","high","low")) for(v in c("DIQ160","DIQ170")) {
 d <- sample[sample$period==pe & (ac=="overall" | sample$access==ac),]
 ss <- status(d[[v]])
 for(st in c("yes","no","refused","dont_know","missing_unresolved","other_code")) responses[[length(responses)+1]] <- data.frame(period=pe,access=ac,variable=v,status=st,N=nrow(d),n=sum(ss==st),unweighted_pct=100*mean(ss==st),weighted_pct=100*wmean(as.numeric(ss==st),d$wt))
}
write.csv(bind_rows(responses),file.path(out,"response_status_by_period_access.csv"),row.names=FALSE)
cytab <- bind_rows(lapply(split(sample,sample$cycle),function(d) bind_rows(lapply(c("DIQ160","DIQ170"),function(v){ss<-status(d[[v]]);data.frame(cycle=d$cycle[1],variable=v,status=names(table(ss)),n=as.integer(table(ss)),N=nrow(d))}))))
write.csv(cytab,file.path(out,"response_status_by_cycle.csv"),row.names=FALSE)
vars <- MODELS$F_ext
stopifnot("prediab" %in% vars, !"atrisk" %in% vars)
# Follow-up sensitivity: keep the submitted F_ext predictor structure fixed.
# (i) fit only known DIQ160 responses; (ii) code unknown DIQ160 as yes.
# A missingness indicator is not estimable when there are no unknown responses.
metric_rows<-list(); coef_rows<-list(); thr_rows<-list()
evaluate <- function(label,dtr,dte,alternate_yes=FALSE){
 if(alternate_yes){dtr$prediab<-as.numeric(dtr$DIQ160!=2 | is.na(dtr$DIQ160));dte$prediab<-as.numeric(dte$DIQ160!=2 | is.na(dte$DIQ160))}
 pdev<-loco_p(dtr,vars); old<-dtr$age>=SPLIT_AGE
 tar<-wmean(dtr$ref_bang,dtr$wt); ts<-wmean(dtr$ref_bang[old&dtr$case==1],dtr$wt[old&dtr$case==1])
 t0<-v5_threshold(pdev,dtr$wt,tar);t1<-min(cut_for_sens(pdev[old],dtr$case[old],dtr$wt[old],ts),t0)
 beta<-fit_l2(as.matrix(dtr[,vars]),dtr$case,dtr$wt)
 coef_rows[[length(coef_rows)+1]]<<-data.frame(handling=label,term=names(beta),coef=unname(beta))
 thr_rows[[length(thr_rows)+1]]<<-data.frame(handling=label,n_development=nrow(dtr),single=t0,older=t1,target_referral=tar)
 for(pop in c("all","known_DIQ160")){
  dx<-if(pop=="all") dte else dte[dte$DIQ160 %in% c(1,2),]
  pp<-plogis(lp_of(beta,dx))
  for(g in c("overall","under60","60plus")){
   k<-if(g=="overall") rep(TRUE,nrow(dx)) else if(g=="under60") dx$age<60 else dx$age>=60
   y<-dx$case[k];w<-dx$wt[k];p<-pp[k];older<-dx$age[k]>=60
   for(rule in c("S0","S1")){
    ref<-as.numeric(p>=if(rule=="S0") t0 else ifelse(older,t1,t0))
    metric_rows[[length(metric_rows)+1]]<<-data.frame(handling=label,evaluation_population=pop,group=g,rule=rule,N=sum(k),cases=sum(y),AUC=w_auc(p,y,w),referral=100*wmean(ref,w),sensitivity=100*sum(w*y*ref)/sum(w*y),specificity=100*sum(w*(1-y)*(1-ref))/sum(w*(1-y)),PPV=100*sum(w*y*ref)/sum(w*ref))
   }
  }
 }
 invisible(NULL)
}
ddev<-full[full$dev,];deval<-full[full$eval,]
evaluate("unknown_as_no_original",ddev,deval)
evaluate("known_only_fit",ddev[ddev$DIQ160 %in% c(1,2),],deval)
evaluate("unknown_as_yes",ddev,deval,TRUE)
mr<-bind_rows(metric_rows)
write.csv(mr,file.path(out,"fixed_structure_sensitivity_metrics.csv"),row.names=FALSE)
write.csv(bind_rows(coef_rows),file.path(out,"fixed_structure_sensitivity_coefficients.csv"),row.names=FALSE)
write.csv(bind_rows(thr_rows),file.path(out,"fixed_structure_sensitivity_thresholds.csv"),row.names=FALSE)
oa<-read.csv(file.path(tabdir,"S12_information_auc_by_model.csv"));pa<-mr[mr$handling=="unknown_as_no_original"&mr$evaluation_population=="all"&mr$group=="overall"&mr$rule=="S0",]
stopifnot(abs(pa$AUC-oa$AUC[oa$period=="eval"&oa$model=="F_ext"])<1e-10)
writeLines(c(capture.output(sessionInfo()),paste("source_md5",unname(tools::md5sum(src))),paste("max_coefficient_difference",format(max(abs(ck$difference)),digits=17))),file.path(out,"execution_receipt.txt"))
cat("COMPLETE; max coefficient difference",max(abs(ck$difference)),"\n")
print(bind_rows(responses)[bind_rows(responses)$access=="overall",])
print(mr[mr$group=="overall"&mr$rule=="S0",])
