# Reproduce the supplied implementation; only data, source and output paths change.
# Rscript rerun_s17.R DATA_DIR OUT_DIR SOURCE_R LIBRARY_DIR ORIGINAL_TABLE_DIR
args <- commandArgs(trailingOnly=TRUE)
here <- "outputs/diabetes_manuscript_revision_20261104_v2/analysis"
DATA_OVERRIDE <- if(length(args)>=1) args[1] else file.path(here,"raw")
OUT_OVERRIDE <- if(length(args)>=2) args[2] else file.path(here,"s17_rerun")
SOURCE <- if(length(args)>=3) args[3] else "outputs/diabetes_manuscript_revision_20261104/supplement/code/external_validation_v5_WTPH2YR.R"
LIB <- if(length(args)>=4) args[4] else "outputs/nhanes_nlr_depression_start_2026-10-20/runtime/library"
ORIG <- if(length(args)>=5) args[5] else "outputs/diabetes_manuscript_revision_20261104/supplement/tables"
.libPaths(c(LIB,.libPaths()));dir.create(OUT_OVERRIDE,recursive=TRUE,showWarnings=FALSE)
options(digits=17)
mult_records <- list(); boot_index <- 0L
expressions <- parse(SOURCE,encoding="UTF-8")
for(expr in expressions) {
  # Stop before original graphics/export; tables are written explicitly below.
  if(is.call(expr) && identical(expr[[1]],as.name("png"))) break
  lhs <- if(is.call(expr) && identical(expr[[1]],as.name("<-")) && is.symbol(expr[[2]])) as.character(expr[[2]]) else ""
  if(lhs=="PROJECT_DIR") {PROJECT_DIR <- here;next}
  if(lhs=="DATA_DIR") {DATA_DIR <- DATA_OVERRIDE;next}
  if(lhs=="OUT_DIR") {OUT_DIR <- OUT_OVERRIDE;next}
  if(lhs=="AUTO_DOWNLOAD") {AUTO_DOWNLOAD <- FALSE;next}
  eval(expr)
  if(lhs=="boot_weights") {
    boot_weights_original <- boot_weights
    boot_weights <- function(d) {
      wb <- boot_weights_original(d);boot_index <<- boot_index+1L
      ix <- !duplicated(paste(d$strata,d$psu))
      mult_records[[boot_index]] <<- data.frame(replicate=boot_index,strata=d$strata[ix],psu=d$psu[ix],multiplier=wb[ix]/d$wt[ix])
      wb
    }
  }
}
# Export newly generated source-implementation model bundles, not historical objects.
fit_export <- fit_l2
fb <- as.list(body(fit_export));fb[[length(fb)]] <- quote(list(beta=beta,standardized_parameters=op$par,weighted_mean=mu,weighted_sd=sg,optimization=op,C=C,n_norm=n_norm))
body(fit_export) <- as.call(fb)
bundles <- list();coef_rows <- list();pre_rows <- list()
for(m in c("V5nf","F3")) {
  v <- MODELS[[m]];obj <- fit_export(as.matrix(dd[,v]),dd$case,dd$wt)
  stopifnot(identical(obj$beta,FIT[[m]]))
  obj$model <- m;obj$variables <- v;obj$thresholds <- THR[[m]]
  obj$development_cycles <- DEV_CYCLES;obj$sample_n <- nrow(dd)
  obj$provenance <- "New rerun of the supplied R fitting implementation, not the historical fitted object"
  obj$source_md5 <- unname(tools::md5sum(SOURCE));obj$R_version <- R.version.string
  bundles[[m]] <- obj
  coef_rows[[m]] <- data.frame(model=m,term=names(obj$beta),coefficient=unname(obj$beta))
  pre_rows[[m]] <- data.frame(model=m,term=v,weighted_mean=unname(obj$weighted_mean),weighted_sd=unname(obj$weighted_sd))
  saveRDS(obj,file.path(OUT_DIR,paste0("S17_",m,"_source_rerun_model.rds")))
}
saveRDS(list(models=bundles,fit_l2=fit_l2,thresholds=THR[c("V5nf","F3")],Bang_nf_cut=BANG_NF_CUT,development_target_referral=target_ref,development_target_old_sensitivity=target_sens),file.path(OUT_DIR,"S17_derived_model_bundle.rds"))
co_export <- bind_rows(coef_rows);co_export$coefficient <- sprintf("%.17g",co_export$coefficient)
pre_export <- bind_rows(pre_rows);pre_export$weighted_mean <- sprintf("%.17g",pre_export$weighted_mean);pre_export$weighted_sd <- sprintf("%.17g",pre_export$weighted_sd)
write.csv(co_export,file.path(OUT_DIR,"S17_source_R_full_precision_coefficients.csv"),row.names=FALSE)
write.csv(pre_export,file.path(OUT_DIR,"S17_source_R_preprocessing.csv"),row.names=FALSE)
write.csv(data.frame(model=names(THR),cut_single=sapply(THR,`[`,"single"),cut_60plus=sapply(THR,`[`,"old")),file.path(OUT_DIR,"S17_thresholds_rerun.csv"),row.names=FALSE)
write.csv(data.frame(Bang_nf_cut=BANG_NF_CUT,target_referral=target_ref,target_old_sensitivity=target_sens),file.path(OUT_DIR,"S17_threshold_development_targets.csv"),row.names=FALSE)
# Bounded single-proportion intervals. Constant domains cannot identify design-based
# uncertainty; report point estimates with NA endpoints rather than a spurious 0-width CI.
bound <- function(var,idx,complement=FALSE) {
  idx <- idx %in% TRUE & !is.na(full[[var]])
  ddsg <- des[idx,];x <- full[[var]][idx];n <- sum(idx)
  p <- as.numeric(coef(svymean(as.formula(paste0("~",var)),ddsg,na.rm=TRUE)))
  df <- degf(ddsg)
  if(complement) p <- 1-p
  status <- if(n==0) "empty_domain" else if(length(unique(x))==1) "constant_domain_CI_not_estimable" else if(df<=0) "no_design_df" else "ok"
  ci <- c(NA_real_,NA_real_)
  if(status=="ok") {
    ci <- tryCatch(as.numeric(confint(svyciprop(as.formula(paste0("~",var)),ddsg,method="logit",df=df,na.rm=TRUE))),error=function(e) {status <<- paste0("error: ",conditionMessage(e));c(NA_real_,NA_real_)})
    if(complement) ci <- 1-rev(ci)
    if(any(!is.finite(ci))) {status <- "nonfinite_CI";ci <- c(NA_real_,NA_real_)}
  }
  data.frame(estimate=100*p,lower95=100*ci[1],upper95=100*ci[2],domain_n=n,domain_df=df,CI_method="survey_svyciprop_logit_t",CI_status=status)
}
bounded <- list()
for(pe in periods) for(gn in names(GROUPS)) {
  g <- if(isTRUE(GROUPS[[gn]])) rep(TRUE,nrow(full)) else eval(GROUPS[[gn]],full)%in%TRUE
  kg <- full[[pe]] & g
  for(r in RULES) {
    if(!avail_in(r,pe))next
    v <- refcol(r)
    for(me in c("referral","sensitivity","specificity","PPV")) {
      idx <- switch(me,referral=kg,sensitivity=kg & full$case==1,specificity=kg & full$case==0,PPV=kg & full[[v]]%in%1)
      z <- bound(if(me=="PPV") "case" else v,idx,me=="specificity")
      bounded[[length(bounded)+1L]] <- cbind(data.frame(period=pe,group=gn,strategy=r,measure=me),z)
    }
  }
}
bounded <- bind_rows(bounded)
write.csv(bounded,file.path(OUT_DIR,"S17_absolute_proportions_logit_CI.csv"),row.names=FALSE)
write.csv(results,file.path(OUT_DIR,"S17_Taylor_results_source_rerun.csv"),row.names=FALSE)
write.csv(perf_tab,file.path(OUT_DIR,"S17_auc_calibration_rerun.csv"),row.names=FALSE)
write.csv(gap_tab,file.path(OUT_DIR,"S17_access_gaps_rerun.csv"),row.names=FALSE)
write.csv(curve_tab,file.path(OUT_DIR,"S17_family_curves_rerun.csv"),row.names=FALSE)
write.csv(diff_tab,file.path(OUT_DIR,"S17_family_differences_rerun.csv"),row.names=FALSE)
write.csv(key_tab,file.path(OUT_DIR,"S17_key_budgets_rerun.csv"),row.names=FALSE)
write.csv(bind_rows(mult_records),file.path(OUT_DIR,"S17_bootstrap_PSU_multipliers.csv"),row.names=FALSE)
saveRDS(bt,file.path(OUT_DIR,"S17_bootstrap_curve_array.rds"))
long <- bind_rows(lapply(seq_len(B_BOOT),function(b)bind_rows(lapply(seq_len(ncol(bt)),function(j)data.frame(replicate=b,curve=dimnames(bt)[[2]][j],budget=BUDGETS,sensitivity=bt[,j,b])))))
long$status <- ifelse(is.na(long$sensitivity),"budget_outside_attainable_range","ok")
write.csv(long,file.path(OUT_DIR,"S17_bootstrap_curve_replicates.csv"),row.names=FALSE)
write.csv(long %>% group_by(curve,budget) %>% summarise(n_replicates=n(),n_valid=sum(is.finite(sensitivity)),.groups="drop"),file.path(OUT_DIR,"S17_bootstrap_effective_counts.csv"),row.names=FALSE)
write.csv(full[full$ext,c("SEQN","cycle","wt","strata","psu","age","age_grp","case","access","p_V5nf","p_F3","ref_bang_nf","ref_V5nf_R0","ref_V5nf_R1","ref_F3_R0","ref_F3_R1")],file.path(OUT_DIR,"S17_external_predictions_rerun.csv"),row.names=FALSE)
compare <- function(new,file,keys,cols) {
  orig <- read.csv(file,check.names=FALSE)
  mg <- merge(new,orig,by=keys,suffixes=c("_rerun","_supplied"))
  for(co in cols) mg[[paste0(co,"_difference")]] <- mg[[paste0(co,"_rerun")]]-mg[[paste0(co,"_supplied")]]
  mg
}
comparisons <- list(
  results=compare(results,file.path(ORIG,"S17_ext_results.csv"),c("period","group","strategy","measure"),c("estimate","lower95","upper95")),
  curves=compare(curve_tab,file.path(ORIG,"S17_ext_family_curves.csv"),c("curve","budget"),c("sens","lower95","upper95")),
  differences=compare(diff_tab,file.path(ORIG,"S17_ext_family_differences.csv"),c("comparison","budget"),c("diff","lower95","upper95")),
  key_budgets=compare(key_tab,file.path(ORIG,"S17_ext_differences_at_key_budgets.csv"),c("comparison","at"),c("budget","diff","lower95","upper95")))
summary_rows <- list()
for(nm in names(comparisons)) {
  z <- comparisons[[nm]];write.csv(z,file.path(OUT_DIR,paste0("S17_comparison_",nm,".csv")),row.names=FALSE)
  dc <- names(z)[grepl("_difference$",names(z))]
  summary_rows[[nm]] <- bind_rows(lapply(dc,function(co)data.frame(table=nm,column=co,n=nrow(z),max_absolute_difference=max(abs(z[[co]]),na.rm=TRUE))))
}
write.csv(bind_rows(summary_rows),file.path(OUT_DIR,"S17_reproduction_summary.csv"),row.names=FALSE)
write.csv(merge(bounded,results,by=c("period","group","strategy","measure"),suffixes=c("_logit","_Taylor")),file.path(OUT_DIR,"S17_logit_vs_original_Taylor.csv"),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(OUT_DIR,"sessionInfo.txt"))
write.csv(data.frame(file=list.files(DATA_DIR,full.names=TRUE),md5=unname(tools::md5sum(list.files(DATA_DIR,full.names=TRUE)))),file.path(OUT_DIR,"input_md5.csv"),row.names=FALSE)
cat("\nSOURCE RERUN FINISHED\n");print(bind_rows(summary_rows));print(table(bounded$CI_status));print(bounded[bounded$period=="ext" & bounded$measure=="sensitivity" & bounded$strategy%in%c("Bang_nf","V5nf_R1"),])
