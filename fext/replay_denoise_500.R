# Replay supplied denoise source through inference, without graphics.
# Arguments: workspace, raw-directory, output-directory, source-script, original-tables.
args<-commandArgs(trailingOnly=TRUE)
ws<-normalizePath(if(length(args))args[1] else getwd(),winslash="/")
.libPaths(c(file.path(ws,"outputs/nhanes_nlr_depression_start_2026-10-20/runtime/library"),.libPaths()))
replay_out<-if(length(args)>=3)args[3] else file.path(ws,"outputs/diabetes_manuscript_revision_20261104_v2/fext/denoise_500_replay")
dir.create(replay_out,recursive=TRUE,showWarnings=FALSE)
source_file<-if(length(args)>=4)args[4] else file.path(ws,"outputs/diabetes_manuscript_revision_20261104/supplement/code/denoise_features_v5.R")
PROJECT_DIR<-replay_out;OUT_DIR<-replay_out
DATA_DIR<-if(length(args)>=2)args[2] else file.path(ws,"outputs/diabetes_manuscript_revision_20261104_v2/analysis/raw")
start_time<-Sys.time();multiplier_records<-list()
for(e in parse(source_file,encoding="UTF-8")){
 lhs<-if(is.call(e)&&identical(e[[1]],as.name("<-"))&&is.symbol(e[[2]]))as.character(e[[2]]) else ""
 if(lhs=="band")break
 if(lhs %in% c("PROJECT_DIR","DATA_DIR","OUT_DIR"))next
 if(is.call(e)&&identical(e[[1]],as.name("dir.create")))next
 eval(e)
 if(lhs=="boot_weights"){
  supplied_boot_weights<-boot_weights
  boot_weights<-function(d){
   bw<-supplied_boot_weights(d)
   k<-!duplicated(paste(d$strata,d$psu))
   multiplier_records[[length(multiplier_records)+1]]<<-data.frame(replicate=length(multiplier_records)+1L,strata=d$strata[k],psu=d$psu[k],multiplier=bw[k]/d$wt[k])
   bw
  }
 }
}
stopifnot(B_BOOT==500,length(multiplier_records)==500,dim(CV$eval$boot)[3]==500)
write.csv(results,file.path(replay_out,"results_by_model.csv"),row.names=FALSE)
write.csv(auc_tab,file.path(replay_out,"auc_by_model.csv"),row.names=FALSE)
write.csv(coef_tab,file.path(replay_out,"coefficients_by_model.csv"),row.names=FALSE)
write.csv(thr_tab,file.path(replay_out,"thresholds_by_model.csv"),row.names=FALSE)
write.csv(curve_tab,file.path(replay_out,"family_curves.csv"),row.names=FALSE)
write.csv(diff_tab,file.path(replay_out,"family_differences.csv"),row.names=FALSE)
write.csv(key_tab,file.path(replay_out,"family_differences_at_key_budgets.csv"),row.names=FALSE)
write.csv(do.call(rbind,multiplier_records),file.path(replay_out,"PSU_bootstrap_multipliers.csv"),row.names=FALSE)
saveRDS(list(CV=CV,budgets=BUDGETS,comparisons=CMP,key_budgets=key_budgets),file.path(replay_out,"replicate_curve_arrays.rds"))
arr<-CV$eval$boot
long<-expand.grid(budget=BUDGETS,curve=dimnames(arr)[[2]],replicate=seq_len(B_BOOT),KEEP.OUT.ATTRS=FALSE,stringsAsFactors=FALSE)
long$sensitivity<-as.vector(arr)
write.csv(long,file.path(replay_out,"replicate_curve_metrics.csv"),row.names=FALSE)
status<-do.call(rbind,lapply(seq_len(B_BOOT),function(b)data.frame(replicate=b,curve=dimnames(arr)[[2]],finite_budget_points=colSums(is.finite(arr[,,b])),undefined_budget_points=colSums(is.na(arr[,,b])),total_budget_points=length(BUDGETS))))
write.csv(status,file.path(replay_out,"replicate_status.csv"),row.names=FALSE)
ref<-if(length(args)>=5)args[5] else file.path(ws,"outputs/diabetes_manuscript_revision_20261104/supplement/tables")
tables<-list(results_by_model=results,auc_by_model=auc_tab,coefficients_by_model=coef_tab,thresholds_by_model=thr_tab,family_curves=curve_tab,family_differences=diff_tab,family_differences_at_key_budgets=key_tab)
checks<-lapply(names(tables),function(nm){
 a<-tables[[nm]];b<-read.csv(file.path(ref,paste0("S12_information_",nm,".csv")),check.names=FALSE)
 stopifnot(identical(names(a),names(b)),nrow(a)==nrow(b))
 nums<-names(a)[vapply(a,is.numeric,logical(1))]
 nonnum<-setdiff(names(a),nums)
 stopifnot(all(vapply(nonnum,function(v)identical(as.character(a[[v]]),as.character(b[[v]])),logical(1))))
 do.call(rbind,lapply(nums,function(v){
  ok<-is.finite(a[[v]])&is.finite(b[[v]])
  data.frame(table=nm,field=v,n_rows=nrow(a),n_compared=sum(ok),NA_pattern_matches=identical(is.na(a[[v]]),is.na(b[[v]])),max_abs_difference=if(any(ok))max(abs(a[[v]][ok]-b[[v]][ok])) else NA_real_)
 }))
})
checks<-do.call(rbind,checks)
write.csv(checks,file.path(replay_out,"submitted_table_comparisons.csv"),row.names=FALSE)
writeLines(c(paste("started",start_time),paste("finished",Sys.time()),paste("replicates",B_BOOT),paste("seed",20261103),"Source computation and seed retained; only paths adapted, graphics omitted, and bootstrap returns observed for multiplier export.","Uncertainty is conditional on submitted fixed equations, selected features, and thresholds; no bootstrap refitting or feature reselection.","Undefined curve regions below the family's minimum feasible referral share remain NA. Original >5% undefined-iteration reporting rule retained.",paste("source_md5",unname(tools::md5sum(source_file))),capture.output(sessionInfo())),file.path(replay_out,"execution_receipt.txt"))
print(checks);cat("DENOISE_500_REPLAY_COMPLETE\n")
