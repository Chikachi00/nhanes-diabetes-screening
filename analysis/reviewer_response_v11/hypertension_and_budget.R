# v11 post hoc hypertension-definition sensitivity and public-health figure source.
# Usage: Rscript hypertension_and_budget.R INPUT_ROOT BPX_DIRECTORY OUTPUT_DIRECTORY [B=1000]
# INPUT_ROOT restores Data Sheet 2 + Data Sheet 3 and contains raw/ BPQ files.
# Existing archived main-equation fits and Rao-Wu multipliers are reused, not refitted.
args<-commandArgs(TRUE);V<-args[1];BPX<-args[2];out<-args[3]
B<-if(length(args)>3)as.integer(args[4]) else 1000L
suppressMessages(library(haven))
dir.create(out,recursive=TRUE,showWarnings=FALSE)
obj<-readRDS(file.path(V,'reproducibility/main_inputs/prepared.rds'))
e<-obj$eval;d<-obj$dev
pr<-read.csv(file.path(V,'reproducibility/main_archive/evaluation_predictions.csv'))
stopifnot(identical(e$SEQN,pr$SEQN),B<=ncol(obj$multipliers))
bp<-read.csv(file.path(V,'reproducibility/main_archive/bootstrap_parameters.csv'),check.names=FALSE)
lk<-read.csv(file.path(V,'reproducibility/main_archive/locked_strategy_parameters.csv'))
vars<-c('a','w','male','famhx','htn','inactive');terms<-c('(Intercept)',vars)
beta<-setNames(lk$value,lk$term)[terms]
X<-cbind('(Intercept)'=1,as.matrix(e[,vars]));p0<-plogis(as.vector(X%*%beta))
stopifnot(max(abs(p0-pr$probability))<1e-12)
probboot<-plogis(X%*%t(as.matrix(bp[seq_len(B),terms])))
rd<-function(cy){
 q<-as.data.frame(read_xpt(file.path(V,'raw',paste0('BPQ_',cy,'.xpt'))))
 z<-as.data.frame(read_xpt(file.path(BPX,paste0('BPX_',cy,'.xpt'))))
 stopifnot(!anyDuplicated(q$SEQN),!anyDuplicated(z$SEQN))
 q<-q[,c('SEQN','BPQ020','BPQ040A','BPQ050A')]
 z<-z[,c('SEQN',paste0('BPXSY',1:4),paste0('BPXDI',1:4))]
 m<-merge(q,z,by='SEQN',all=TRUE);m$cycle<-cy;m}
raw<-do.call(rbind,lapply(c('E','F','G','H','I','J'),rd))
stopifnot(!anyDuplicated(raw$SEQN))
# Mean of all valid available readings, 1-4; DBP zero is valid in the CDC codebook.
S<-as.matrix(raw[,paste0('BPXSY',1:4)]);D<-as.matrix(raw[,paste0('BPXDI',1:4)])
stopifnot(all(is.na(S)|S>0),all(is.na(D)|D>=0))
raw$n_sbp<-rowSums(!is.na(S));raw$n_dbp<-rowSums(!is.na(D))
raw$mean_sbp<-rowMeans(S,na.rm=TRUE);raw$mean_dbp<-rowMeans(D,na.rm=TRUE)
raw$mean_sbp[!is.finite(raw$mean_sbp)]<-NA;raw$mean_dbp[!is.finite(raw$mean_dbp)]<-NA
yn<-function(x)ifelse(x%in%c(1,2),as.numeric(x==1),NA_real_)
raw$self<-yn(raw$BPQ020);raw$med<-yn(raw$BPQ050A)
# Structurally skipped medication items are negative only after an explicit upstream no.
k<-is.na(raw$BPQ050A)&(raw$BPQ020==2|raw$BPQ040A==2);k[is.na(k)]<-FALSE
raw$med[k]<-0
or3<-function(...){A<-cbind(...);ifelse(rowSums(A==1,na.rm=TRUE)>0,1,
 ifelse(rowSums(is.na(A))==0,0,NA_real_))}
raw$elevated<-or3(ifelse(is.na(raw$mean_sbp),NA,as.numeric(raw$mean_sbp>=140)),
 ifelse(is.na(raw$mean_dbp),NA,as.numeric(raw$mean_dbp>=90)))
raw$history_med<-or3(raw$self,raw$med)
raw$composite<-or3(raw$self,raw$med,raw$elevated)
raw$bp_complete<-is.finite(raw$mean_sbp)&is.finite(raw$mean_dbp)
raw$single_reading<-as.numeric(raw$n_sbp==1|raw$n_dbp==1)
er<-raw[match(e$SEQN,raw$SEQN),];dr<-raw[match(d$SEQN,raw$SEQN),]
stopifnot(identical(er$SEQN,e$SEQN),identical(dr$SEQN,d$SEQN),all(er$self==e$htn),all(dr$self==d$htn))
keep<-!is.na(er$composite)&!is.na(er$history_med)
strict<-keep&er$bp_complete
stopifnot(all(er$composite[keep]>=er$self[keep]),all(er$history_med[keep]>=er$self[keep]))
groups<-list(overall=rep(TRUE,nrow(e)),young=e$age_grp=='20-39',mid=e$age_grp=='40-59',
 old=e$age_grp=='>=60',high=e$access=='high',low=e$access=='low')
# Record all analytic rules before any result is calculated.
writeLines(c('v11 analysis specified before viewing numerical output',
 'Post hoc evaluation-period sensitivity; no new external validation or strategy selection.',
 'Primary BP mean: all valid readings 1-4; a single available reading is retained.',
 'Composite: self-reported history OR current HBP medication OR mean SBP>=140 OR DBP>=90.',
 'Refused/dont know are missing. Medication skips set negative only after explicit upstream no.',
 'Common sample: all three hypertension definitions ascertainable. Strict BP sample also tabulated.',
 'Main equation and original development-derived R0/R1 cut-offs unchanged.',
 'Original 1000 Rao-Wu multiplier columns and stored refitted parameters reused in original order.',
 'Equal-budget equation threshold reselected in evaluation within every replicate; ties retained.',
 'Threshold minimises absolute distance to target referral among attainable complete tie groups.',
 'Performance intervals: estimate +/-1.96 bootstrap SD; percentile intervals also exported.',
 'Nested R1-R0 incremental metrics: percentile intervals shown, with zero-event replicate counts.',
 'Access gap: high-minus-low sensitivity; positive gap-change means widening.',
 'Age standard: combined evaluation HbA1c-positive age distribution re-estimated each replicate.'),file.path(out,'analysis_protocol.txt'))
covrow<-function(z,r,lab,w){data.frame(period=lab,n=nrow(z),cases=sum(z$case),
 common_n=sum(!is.na(r$composite)&!is.na(r$history_med)),
 common_cases=sum(z$case[!is.na(r$composite)&!is.na(r$history_med)]),
 common_weighted_pct=100*sum(w[!is.na(r$composite)&!is.na(r$history_med)])/sum(w),
 no_usable_BP_n=sum(!r$bp_complete),single_BP_n=sum(r$single_reading==1),
 history_med_added_n=sum(r$history_med==1&r$self==0,na.rm=TRUE),
 composite_added_n=sum(r$composite==1&r$self==0,na.rm=TRUE))}
coverage<-rbind(covrow(d,dr,'development',d$wt),covrow(e,er,'evaluation',e$wt))
write.csv(coverage,file.path(out,'S24_sample_coverage.csv'),row.names=FALSE)
getmet<-function(ref,w,subset){
 o<-numeric()
 for(g in names(groups)){k<-groups[[g]]&subset;ww<-w[k];y<-e$case[k];r<-ref[k]
  pos<-sum(ww*y);rf<-sum(ww*r);found<-sum(ww*y*r);tot<-sum(ww)
  v<-c(sensitivity=100*found/pos,referral=100*rf/tot,ppv=100*found/rf,
   found_per1000=1000*found/tot,referrals_per1000=1000*rf/tot)
  names(v)<-paste(g,names(v),sep='::');o<-c(o,v)}
 c(o,'access_gap::sensitivity'=unname(o['high::sensitivity']-o['low::sensitivity']))}
nearest<-function(p,w,target){k<-w>0;o<-order(p[k],decreasing=TRUE);pp<-p[k][o];ww<-w[k][o]
 en<-c(which(diff(pp)!=0),length(pp));a<-cumsum(ww)[en]/sum(ww)
 pp[en[which.min(abs(a-target))]]}
metall<-function(p,w,t0,to,subset){
 scores<-list(Bang_self=e$bang,Bang_history_med=e$bang-e$htn+er$history_med,
 Bang_composite=e$bang-e$htn+er$composite)
 R<-lapply(scores,function(s)as.numeric(s>=5));age<-e$age_grp=='>=60'
 R$R0<-as.numeric(p>=t0);R$R1<-as.numeric(p>=ifelse(age,to,t0))
 for(n in names(scores)){
  target<-sum(w[subset]*R[[n]][subset])/sum(w[subset])
  t<-nearest(p[subset],w[subset],target)
  R[[paste0('Eq_to_',n)]]<-as.numeric(p>=t)}
 st<-lapply(R,getmet,w=w,subset=subset);o<-numeric()
 for(n in names(st)){v<-st[[n]];names(v)<-paste(n,names(v),sep='::');o<-c(o,v)}
 cmp<-list(composite_minus_self=c('Bang_composite','Bang_self'),
  med_minus_self=c('Bang_history_med','Bang_self'),R0_minus_Bang_self=c('R0','Bang_self'),
  R0_minus_Bang_composite=c('R0','Bang_composite'),R1_minus_Bang_composite=c('R1','Bang_composite'),
  Eq_minus_Bang_self=c('Eq_to_Bang_self','Bang_self'),
  Eq_minus_Bang_composite=c('Eq_to_Bang_composite','Bang_composite'))
 for(n in names(cmp)){v<-st[[cmp[[n]][1]]]-st[[cmp[[n]][2]]];names(v)<-paste(n,names(v),sep='::');o<-c(o,v)}
 # Positive-case age standardisation, same age weights for both access groups and strategies.
 ag<-c('20-39','40-59','>=60');y<-e$case==1
 aw<-sapply(ag,function(a)sum(w[subset&y&e$age_grp==a]));aw<-aw/sum(aw)
 for(n in c('Bang_self','Bang_composite','R0','R1')){
  gap<-sapply(ag,function(a){f<-function(ac){k<-subset&y&e$age_grp==a&e$access==ac;sum(w[k]*R[[n]][k])/sum(w[k])};f('high')-f('low')})
  o<-c(o,setNames(100*sum(aw*gap),paste0(n,'::age_standardised_gap::sensitivity')))}
 for(n in c('Bang_composite','R0','R1'))o<-c(o,setNames(o[paste0(n,'::age_standardised_gap::sensitivity')]-o['Bang_self::age_standardised_gap::sensitivity'],paste0(n,'_minus_self::age_standardised_gap::sensitivity')))
 o}
desc<-function(w){o<-numeric()
 for(g in names(groups)){
  k<-groups[[g]]&keep
  vs<-list(self=er$self,history_med=er$history_med,composite=er$composite,
    unreported_elevated=as.numeric(er$self==0&er$elevated==1),
    newly_hypertensive=as.numeric(er$self==0&er$composite==1))
  for(n in names(vs)){v<-vs[[n]];kk<-k&!is.na(v);o<-c(o,setNames(100*sum(w[kk]*v[kk])/sum(w[kk]),paste(g,n,sep='::')))}}
 # Age-standardise descriptive prevalence to all common-sample participants.
 ag<-c('20-39','40-59','>=60');aw<-sapply(ag,function(a)sum(w[keep&e$age_grp==a]));aw<-aw/sum(aw)
 for(n in c('self','composite','unreported_elevated')){
  v<-switch(n,self=er$self,composite=er$composite,unreported_elevated=as.numeric(er$self==0&er$elevated==1))
  av<-sapply(c('high','low'),function(ac)sum(aw*sapply(ag,function(a){k<-keep&e$access==ac&e$age_grp==a&!is.na(v);sum(w[k]*v[k])/sum(w[k])})))
  o<-c(o,setNames(av*100,paste0(c('high','low'),'_age_standardised::',n)),setNames(100*(av['low']-av['high']),paste0('low_minus_high_age_standardised::',n)))}
 o}
t0<-lk$value[lk$term=='THRESHOLD'];to<-lk$value[lk$term=='THRESHOLD_OLD']
pt<-metall(p0,e$wt,t0,to,keep);ptstrict<-metall(p0,e$wt,t0,to,strict);pd<-desc(e$wt)
grid<-seq(.10,.60,by=.01)
curve<-function(p,w){k<-w>0;o<-order(p[k],decreasing=TRUE);pp<-p[k][o];ww<-w[k][o];y<-e$case[k][o]
 en<-c(which(diff(pp)!=0),length(pp));cw<-cumsum(ww)[en]/sum(ww);cp<-cumsum(ww*y)[en]/sum(ww*y)
 ii<-vapply(grid,function(t)which.min(abs(cw-t)),1L)
 c(cw[ii]*100,cp[ii]*100)}
pc<-curve(p0,e$wt)
M<-matrix(NA_real_,B,length(pt),dimnames=list(NULL,names(pt)))
MS<-matrix(NA_real_,B,length(ptstrict),dimnames=list(NULL,names(ptstrict)))
MD<-matrix(NA_real_,B,length(pd),dimnames=list(NULL,names(pd)))
MC<-matrix(NA_real_,B,length(pc))
for(b in seq_len(B)){
 w<-e$wt*obj$multipliers[obj$idxeval,b]
 M[b,]<-metall(probboot[,b],w,bp$THRESHOLD[b],bp$THRESHOLD_OLD[b],keep)
 MS[b,]<-metall(probboot[,b],w,bp$THRESHOLD[b],bp$THRESHOLD_OLD[b],strict)
 MD[b,]<-desc(w);MC[b,]<-curve(probboot[,b],w)
 if(b%%100==0){cat('v11 replicate',b,'of',B,'\n');flush.console()}}
summ<-function(pt,M){do.call(rbind,lapply(seq_along(pt),function(i){v<-M[,i];v<-v[is.finite(v)];se<-sd(v);q<-quantile(v,c(.025,.975),names=FALSE)
 data.frame(key=names(pt)[i],estimate=unname(pt[i]),se=se,lower95=unname(pt[i])-1.96*se,upper95=unname(pt[i])+1.96*se,
 percentile_lower95=q[1],percentile_upper95=q[2],valid_B=length(v))}))}
write.csv(summ(pt,M),file.path(out,'S25_hypertension_strategy_comparisons.csv'),row.names=FALSE)
write.csv(summ(ptstrict,MS),file.path(out,'S25_strict_BP_sample.csv'),row.names=FALSE)
write.csv(summ(pd,MD),file.path(out,'S24_hypertension_by_group.csv'),row.names=FALSE)
ng<-length(grid);curvept<-data.frame(target_referral_pct=grid*100,actual_referral_pct=pc[seq_len(ng)],sensitivity_pct=pc[ng+seq_len(ng)],
 lower95=apply(MC[,ng+seq_len(ng),drop=FALSE],2,quantile,.025),upper95=apply(MC[,ng+seq_len(ng),drop=FALSE],2,quantile,.975))
write.csv(curvept,file.path(out,'Figure5_equation_budget_curve.csv'),row.names=FALSE)
bangcurve<-do.call(rbind,lapply(sort(unique(e$bang),decreasing=TRUE),function(t){r<-e$bang>=t;v<-getmet(r,e$wt,rep(TRUE,nrow(e)));data.frame(threshold=t,referral_pct=v['overall::referral'],sensitivity_pct=v['overall::sensitivity'])}))
write.csv(bangcurve,file.path(out,'Figure5_Bang_attainable_points.csv'),row.names=FALSE)
write.csv(data.frame(replicate=seq_len(B),M,check.names=FALSE),file.path(out,'hypertension_bootstrap_replicates.csv'),row.names=FALSE)
write.csv(data.frame(replicate=seq_len(B),MD,check.names=FALSE),file.path(out,'hypertension_description_replicates.csv'),row.names=FALSE)
write.csv(data.frame(replicate=seq_len(B),MC,check.names=FALSE),file.path(out,'budget_curve_replicates.csv'),row.names=FALSE)
# These identifiers and BP values are companion-data inputs, not public aggregate code outputs.
write.csv(data.frame(SEQN=e$SEQN,cycle=e$cycle,common=keep,strict_BP=strict,er[,setdiff(names(er),c('SEQN','cycle'))]),file.path(out,'participant_mapping_PRIVATE.csv'),row.names=FALSE)
stopifnot(all(is.finite(M)),all(is.finite(MS)),all(is.finite(MD)),all(is.finite(MC)))
cat('COMPLETE',B,'replicates; common',sum(keep),'strict',sum(strict),'\n')
print(coverage)
print(summ(pt,M)[grepl('^(R0_minus_Bang_composite|composite_minus_self|Eq_minus_Bang_composite)::overall::(sensitivity|referral)$',names(pt)),])
