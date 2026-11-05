# Usage: Rscript reproduce_main.R [output_directory] [--refit-all]
# Base R only. Checks the unchanged main equation and all archived intervals.
args <- commandArgs(trailingOnly=TRUE)
file_arg <- grep('^--file=',commandArgs(),value=TRUE)
here <- dirname(normalizePath(sub('^--file=','',file_arg[1]),winslash='/'))
output_arg <- args[!grepl('^--',args)]
base <- if(length(output_arg)) output_arg[1] else file.path(here,'reproduced_main')
dir.create(base,recursive=TRUE,showWarnings=FALSE)
old <- file.path(here,'main_inputs')
archive <- file.path(here,'main_archive')
source(file.path(here,'main_point_functions.R'),encoding='UTF-8')
locked <- read.csv(file.path(archive,'locked_strategy_parameters.csv'))
locked <- setNames(locked$value,locked$term)
stopifnot(max(abs(point$beta-locked[names(point$beta)]))<1e-10,
          abs(point$t0-locked['THRESHOLD'])<1e-12,
          abs(point$to-locked['THRESHOLD_OLD'])<1e-12)
expected_point <- read.csv(file.path(archive,'point_metrics.csv'))
stopifnot(setequal(expected_point$key,names(pointm)))
point_error <- max(abs(pointm[expected_point$key]-expected_point$estimate))
stopifnot(point_error<1e-10)
saved_predictions <- read.csv(file.path(archive,'evaluation_predictions.csv'))
prediction_error <- max(abs(plogis(point$lp)-saved_predictions$probability))
stopifnot(identical(as.integer(e$SEQN),as.integer(saved_predictions$SEQN)),prediction_error<1e-12)
stopifnot(all(pf$S1>=pf$S0),all(pf$S4>=pf$S1),all(pf$S0[e$age<60]==pf$S1[e$age<60]))
saveRDS(point,file.path(base,'main_point_fit_reproduced.rds'))
write.csv(data.frame(key=names(pointm),estimate=unname(pointm)),file.path(base,'point_metrics.csv'),row.names=FALSE)
expected <- read.csv(file.path(archive,'interval_summary.csv'))
intervals <- list()
for(label in c('full_pipeline','fixed_model_bootstrap')) {
    fname <- if(label=='full_pipeline') 'bootstrap_full_metrics.csv' else 'bootstrap_fixed_metrics.csv'
    z <- read.csv(file.path(archive,fname),check.names=FALSE)
    stopifnot(nrow(z)==1000,!anyDuplicated(z$replicate),setequal(z$replicate,1:1000))
    z <- z[order(z$replicate),];mat <- as.matrix(z[,-1])
    stopifnot(all(is.finite(mat)),setequal(colnames(mat),names(pointm)))
    intervals[[label]] <- do.call(rbind,lapply(seq_len(ncol(mat)),function(j){
        key<-strsplit(colnames(mat)[j],'::',fixed=TRUE)[[1]]
        value<-mat[,j];p<-unname(pointm[colnames(mat)[j]]);s<-sd(value)
        q<-quantile(value,c(.025,.975),names=FALSE)
        data.frame(analysis=label,strategy=key[1],group=key[2],measure=key[3],
                   estimate=p,bootstrap_se=s,lower95=p-1.96*s,upper95=p+1.96*s,
                   percentile_lower95=q[1],percentile_upper95=q[2],
                   bootstrap_mean=mean(value),valid_B=length(value))
    }))
}
got<-do.call(rbind,intervals)
key<-function(z)paste(z$analysis,z$strategy,z$group,z$measure,sep='::')
got<-got[match(key(expected),key(got)),]
cols<-c('estimate','bootstrap_se','lower95','upper95','percentile_lower95','percentile_upper95','bootstrap_mean','valid_B')
interval_error<-max(abs(as.matrix(got[,cols])-as.matrix(expected[,cols])))
stopifnot(nrow(got)==nrow(expected),interval_error<1e-10)
write.csv(got,file.path(base,'intervals_reproduced.csv'),row.names=FALSE)
refit_note<-'Not requested: archived replicate metrics were replayed; the development fit and thresholds were rerun once.'
if('--refit-all'%in%args){
    saved<-read.csv(file.path(archive,'bootstrap_full_metrics.csv'),check.names=FALSE)
    errors<-numeric(1000)
    for(b in 1:1000){
        mul<-obj$multipliers[,b];wd<-d$wt*mul[obj$idxdev];we<-e$wt*mul[obj$idxeval]
        f<-develop(wd);m<-flatten(flags(f),we)
        errors[b]<-max(abs(m-as.numeric(saved[b,names(m)])))
        if(b%%10==0){write.csv(data.frame(replicate=1:b,max_error=errors[1:b]),file.path(base,'refit_all_status.csv'),row.names=FALSE);cat('Refit',b,'\n');flush.console()}
    }
    stopifnot(max(errors)<1e-10)
    refit_note<-sprintf('All 1000 archived full-pipeline replicates refitted; maximum error %.17g.',max(errors))
}
writeLines(c(sprintf('Development n=%d; evaluation n=%d',nrow(d),nrow(e)),
             sprintf('Maximum point-metric error: %.17g',point_error),
             sprintf('Maximum probability error: %.17g',prediction_error),
             sprintf('Maximum interval replay error: %.17g',interval_error),
             'All 1000 unique archived replicate IDs, finite metrics, and nested referral sets checked.',
             refit_note),file.path(base,'verification.txt'))
writeLines(capture.output(sessionInfo()),file.path(base,'sessionInfo.txt'))
cat('PASS main point refit and complete archived-interval replay.\n')
