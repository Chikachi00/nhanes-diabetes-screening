# Validate saved artifacts without rerunning fitting or bootstrap.
args <- commandArgs(trailingOnly=TRUE)
out <- if(length(args))args[1] else "outputs/diabetes_manuscript_revision_20261104_v2/analysis/s17_rerun"
b <- read.csv(file.path(out,"S17_absolute_proportions_logit_CI.csv"))
# A numerically failed interval has no reportable endpoint.
k <- b$CI_status=="nonfinite_CI";b[k,c("lower95","upper95")] <- NA_real_
write.csv(b,file.path(out,"S17_absolute_proportions_logit_CI.csv"),row.names=FALSE)
tg <- read.csv(file.path(out,"S17_Taylor_results_source_rerun.csv"))
write.csv(merge(b,tg,by=c("period","group","strategy","measure"),suffixes=c("_logit","_Taylor")),file.path(out,"S17_logit_vs_original_Taylor.csv"),row.names=FALSE)
bundle <- readRDS(file.path(out,"S17_derived_model_bundle.rds"))
thr17 <- do.call(rbind,lapply(names(bundle$models),function(m)data.frame(model=m,cut_single=sprintf("%.17g",bundle$models[[m]]$thresholds["single"]),cut_60plus=sprintf("%.17g",bundle$models[[m]]$thresholds["old"]))))
write.csv(thr17,file.path(out,"S17_derived_thresholds_17digits.csv"),row.names=FALSE)
co <- read.csv(file.path(out,"S17_source_R_full_precision_coefficients.csv"))
pre <- read.csv(file.path(out,"S17_source_R_preprocessing.csv"))
checks <- list()
for(m in c("V5nf","F3")) {
  z <- bundle$models[[m]]
  xx <- co[co$model==m,];pp <- pre[pre$model==m,]
  stopifnot(identical(xx$term,names(z$beta)))
  checks[[m]] <- data.frame(check=paste0(m,"_CSV_vs_RDS_coefficient"),value=max(abs(xx$coefficient-z$beta)))
  checks[[paste0(m,"_pre")]] <- data.frame(check=paste0(m,"_CSV_vs_RDS_preprocessing"),value=max(abs(pp$weighted_mean-z$weighted_mean),abs(pp$weighted_sd-z$weighted_sd)))
  checks[[paste0(m,"_opt")]] <- data.frame(check=paste0(m,"_optimizer_convergence_code"),value=z$optimization$convergence)
}
mp <- read.csv(file.path(out,"S17_bootstrap_PSU_multipliers.csv"))
ms <- aggregate(multiplier~replicate+strata,mp,sum)
ct <- aggregate(psu~replicate+strata,mp,length)
stopifnot(nrow(mp)==15000L,length(unique(mp$replicate))==500L)
checks[["multiplier"]] <- data.frame(check="PSU_multiplier_stratum_sum_minus_nPSU",value=max(abs(ms$multiplier-ct$psu)))
checks[["bounded"]] <- data.frame(check="bounded_CI_rows_outside_0_100",value=sum(b$lower95<0 | b$upper95>100,na.rm=TRUE))
checks[["point"]] <- data.frame(check="bounded_vs_Taylor_max_point_difference_pp",value=max(abs(merge(b,tg,by=c("period","group","strategy","measure"))$estimate.x-merge(b,tg,by=c("period","group","strategy","measure"))$estimate.y)))
checks[["effective"]] <- data.frame(check="curve_replicates_total",value=nrow(read.csv(file.path(out,"S17_bootstrap_curve_replicates.csv"))))
write.csv(do.call(rbind,checks),file.path(out,"S17_export_checks.csv"),row.names=FALSE)
print(do.call(rbind,checks));print(table(b$CI_status))
