# S15 check: age-standardised access gap with re-estimated versus fixed age standard across the 1,000 replicates
args<-commandArgs(TRUE);V<-args[1]
obj<-readRDS(file.path(V,'reproducibility/main_inputs/prepared.rds'));e<-obj$eval
pr<-read.csv(file.path(V,'reproducibility/main_archive/evaluation_predictions.csv'))
y<-e$case==1;ag<-e$age_grp;ac<-e$access
gap<-function(r,w,agew=NULL){if(is.null(agew)){agew<-tapply(w[y],ag[y],sum);agew<-agew/sum(agew)}
 s<-function(a)sapply(names(agew),function(g){k<-y&ag==g&ac==a;sum(w[k]*r[k])/sum(w[k])})
 sum(agew*(s('high')-s('low')))}
w0<-e$wt;a0<-tapply(w0[y],ag[y],sum);a0<-a0/sum(a0)
pt<-gap(e$ref_bang,w0)
re<-fx<-numeric(1000)
for(b in 1:1000){we<-e$wt*obj$multipliers[obj$idxeval,b];re[b]<-gap(e$ref_bang,we);fx[b]<-gap(e$ref_bang,we,a0)}
s6<-read.csv(file.path(V,'supplement/tables/S6_age_standardised_access.csv'))
print(s6[s6$strategy=='Bang',c('group','estimate','bootstrap_se','bootstrap_mean','percentile_lower95','percentile_upper95')])
cat(sprintf('point %.6f | re-estimated weights: se %.6f mean %.6f pct %.4f %.4f | fixed weights: se %.6f mean %.6f\n',pt,sd(re),mean(re),quantile(re,.025),quantile(re,.975),sd(fx),mean(fx)))
