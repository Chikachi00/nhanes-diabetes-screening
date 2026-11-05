args<-commandArgs(trailingOnly=TRUE)
base<-normalizePath(if(length(args)>=1) args[1] else 'outputs/diabetes_manuscript_revision_20261105_v4',winslash='/')
out<-if(length(args)>=2) args[2] else 'outputs/diabetes_manuscript_revision_20261106_v5/analysis/local_calibration'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
.libPaths(c('outputs/nhanes_nlr_depression_start_2026-10-20/runtime/library',.libPaths()))
suppressPackageStartupMessages(library(survey))
options(survey.lonely.psu='adjust')
o<-readRDS(file.path(base,'reproducibility/main_inputs/prepared.rds'))
e<-o$eval;fr<-read.csv(file.path(base,'reproducibility/main_inputs/full_base.csv'))
idx<-match(e$SEQN,fr$SEQN)
stopifnot(!anyNA(idx),!anyDuplicated(idx))
stopifnot(all(fr$period[idx]=='eval'))
weight_scale<-median(e$wt/fr$wt[idx]);stopifnot(max(abs(e$wt/fr$wt[idx]/weight_scale-1))<1e-10)
print(names(fr));print(names(e));print(dim(fr));print(dim(e))
locked<-read.csv(file.path(base,'reproducibility/main_archive/locked_strategy_parameters.csv'))
cf<-setNames(locked$value,locked$term)
vars<-c('a','w','male','famhx','htn','inactive')
e$p<-plogis(cf['(Intercept)']+as.vector(as.matrix(e[,vars])%*%cf[vars]))
saved<-read.csv(file.path(base,'reproducibility/main_archive/evaluation_predictions.csv'))
stopifnot(identical(as.integer(e$SEQN),as.integer(saved$SEQN)),max(abs(e$p-saved$probability))<1e-12)
cut<-unname(cf['THRESHOLD'])
# All PSU/stratum observations remain in the design; non-domain values contribute zero to ratios.
fr$y<-0;fr$p<-0;fr$analytic<-FALSE;fr$access<-'outside'
fr$y[idx]<-e$case;fr$p[idx]<-e$p;fr$analytic[idx]<-TRUE;fr$access[idx]<-as.character(e$access)
design<-svydesign(ids=~psu,strata=~strata,weights=~wt,data=fr,nest=TRUE)
rows<-list();ii<-0
for(band in c('all','referred_R0','above_cutoff_band')){
    mask<-switch(band,all=rep(TRUE,nrow(e)),referred_R0=e$p>=cut,above_cutoff_band=e$p>=cut&e$p<.035)
    for(group in c('high','low','overall')){
        take<-mask & if(group=='overall') rep(TRUE,nrow(e)) else e$access==group
        k<-rep(FALSE,nrow(fr));k[idx]<-take
        subd<-subset(design,k)
        w<-e$wt[take];y<-e$case[take];pr<-e$p[take]
        obs<-sum(w*y)/sum(w);pred<-sum(w*pr)/sum(w)
        ci<-tryCatch(confint(svyciprop(~y,subd,method='logit',df=degf(subd))),error=function(z)c(NA,NA))
        ii<-ii+1;rows[[ii]]<-data.frame(domain=band,access=group,n=sum(take),positive_n=sum(y),weighted_predicted=pred,weighted_observed=obs,O_E=obs/pred,observed_lower95=ci[1],observed_upper95=ci[2],domain_psu=length(unique(paste(fr$strata[k],fr$psu[k]))),domain_strata=length(unique(fr$strata[k])),design_df=degf(subd),unweighted_observed=mean(y),cutoff_lower=cut,band_upper=if(band=='above_cutoff_band').035 else NA)
    }
}
res<-do.call(rbind,rows)
write.csv(res,file.path(out,'local_calibration_by_access.csv'),row.names=FALSE)
# Survey-linearised contrast of two domain ratios, conditional on the frozen equation and the inspected band.
kh<-fr$analytic&fr$access=='high'&fr$p>=cut&fr$p<.035
kl<-fr$analytic&fr$access=='low'&fr$p>=cut&fr$p<.035
design<-update(design,num_h=as.numeric(kh)*fr$y,den_h=as.numeric(kh),num_l=as.numeric(kl)*fr$y,den_l=as.numeric(kl))
tot<-svytotal(~num_h+den_h+num_l+den_l,design)
contrast<-svycontrast(tot,quote(num_l/den_l-num_h/den_h))
df<-degf(subset(design,kh|kl));estimate<-coef(contrast)[1];se<-SE(contrast)[1]
write.csv(data.frame(contrast='low minus high observed positivity in above-cutoff band',estimate=estimate,SE=se,lower95=estimate-qt(.975,df)*se,upper95=estimate+qt(.975,df)*se,df=df,scope='Conditional on frozen model and reviewer-selected single-sided band; no refit, no bootstrap, no selection adjustment'),file.path(out,'local_difference_Taylor.csv'),row.names=FALSE)
write.csv(e[,c('SEQN','cycle','access','case','wt','p')],file.path(out,'scoring_receipt.csv'),row.names=FALSE)
writeLines(c(sprintf('Frozen cutoff: %.17g',cut),sprintf('Prediction max error: %.17g',max(abs(e$p-saved$probability))),sprintf('Referral n: %d; positive n: %d',sum(e$p>=cut),sum(e$case[e$p>=cut])), 'The inspected interval is on the referred side only. It is not a point estimate at the threshold and does not include non-referred neighbours.', 'Band supplied in review; analysis is post hoc, conditional, and not evidence of access-specific threshold optimality.',capture.output(sessionInfo())),file.path(out,'README.txt'))
print(res);print(c(estimate=estimate,SE=se,lower95=estimate-qt(.975,df)*se,upper95=estimate+qt(.975,df)*se))
