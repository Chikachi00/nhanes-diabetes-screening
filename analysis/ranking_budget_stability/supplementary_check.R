# Frozen-equation supplementary analysis. No model fitting or strategy selection.
# Rscript supplementary_check.R INPUT_ROOT V11_RESULTS OUTPUT_DIRECTORY [B=1000]
args <- commandArgs(TRUE)
stopifnot(length(args)>=3)
V <- args[1]; A <- args[2]; O <- args[3]
B <- if(length(args)>3) as.integer(args[4]) else 1000L
stopifnot(B>=2, B<=1000)
dir.create(O, recursive=TRUE, showWarnings=FALSE)
options(digits=17)
writeLines(c(
 'Material Passport: academic-research-suite / experiment-agent / user-authorised implementation',
 'Version: supplementary_check_20261111; post hoc; earlier v11 results already known.',
 'Rules written before this run produces new estimates; not study preregistration.',
 'Primary metric this round: paired weighted AUC difference, equation minus composite-hypertension Bang, common sample.',
 'Secondary: all integer Bang thresholds, nearest attainable equation budgets, fixed R0 and cycle stability.',
 'Point predictions, coefficients, self-report equation hypertension and development R0/R1 thresholds unchanged.',
 'Bootstrap reuses matched archived development parameter rows and Rao-Wu evaluation weights; no refitting.',
 'Main sample: v11 common; sensitivity: v11 strict BP sample. Cycle analyses use common sample.',
 'AUC uses weighted positive-negative pairs, awarding one half to ties; zero-weight rows excluded.',
 'All integer thresholds from minimum score through maximum+1 included, with max+1 yielding no referrals.',
 'Budget matching includes no-referral boundary and full tied probability groups.',
 'Closest referral fraction selected; equal distance within 1e-14 favours the lower referral fraction.',
 'Budget is reselected in each replicate and in each cycle/domain; matching does not use HbA1c labels.',
 'Equation grid: 10%-60% inclusive by 1 percentage point; pointwise intervals only.',
 'Normal intervals = point +/- 1.96 bootstrap SD; percentile intervals use R quantile type 7.',
 'Simultaneous family: both Bang definitions, both samples, all integer thresholds, paired sensitivity and positives per 1000.',
 'Simultaneous critical value = 95th percentile of maximum absolute centred replicate deviation / bootstrap SD.',
 'Only complete finite family columns with SD>1e-12 enter maximum; zero-SD contrasts have point-width intervals.',
 'Zero-referral PPV and any undefined denominator remain NA with explicit valid_B; never silently discarded.',
 'Core common-sample AUC and cutoff-5 comparisons must be finite in all requested replicates.',
 'No outcome improvement, causal effect, selection-adjusted confirmation or independent external validation claimed.',
 'Composite Bang uses more measurement information; activity mapping is still not identical to original Bang.'
), file.path(O,'analysis_protocol.txt'))

obj <- readRDS(file.path(V,'reproducibility/main_inputs/prepared.rds'))
e <- obj$eval
pr <- read.csv(file.path(V,'reproducibility/main_archive/evaluation_predictions.csv'))
mp <- read.csv(file.path(A,'participant_mapping_PRIVATE.csv'))
bp <- read.csv(file.path(V,'reproducibility/main_archive/bootstrap_parameters.csv'),check.names=FALSE)
lk <- read.csv(file.path(V,'reproducibility/main_archive/locked_strategy_parameters.csv'))
stopifnot(!anyDuplicated(e$SEQN), !anyDuplicated(mp$SEQN),
 identical(e$SEQN,pr$SEQN), identical(e$SEQN,mp$SEQN),
 identical(as.integer(bp$replicate),1:1000), B<=ncol(obj$multipliers),
 all(e$cycle==mp$cycle), all(e$htn==mp$self), all(e$wt>0))
terms <- c('(Intercept)','a','w','male','famhx','htn','inactive')
beta <- setNames(lk$value,lk$term)[terms]
X <- cbind('(Intercept)'=1,as.matrix(e[,terms[-1]]))
p <- plogis(as.vector(X%*%beta))
stopifnot(max(abs(p-pr$probability))<1e-12, all(is.finite(p)))
pb <- plogis(X%*%t(as.matrix(bp[seq_len(B),terms])))
t0 <- lk$value[lk$term=='THRESHOLD']
scores <- list(Bang_self=e$bang, Bang_composite=e$bang-e$htn+mp$composite)
keep <- mp$common; strict <- mp$strict_BP
stopifnot(sum(keep)==11624, sum(e$case[keep])==348,
 sum(strict)==11548, sum(e$case[strict])==344,
 all(strict<=keep), all(scores$Bang_composite[keep]>=scores$Bang_self[keep]))
masks <- list(common=keep, strict_BP=strict)
for(cy in c('H','I','J')) masks[[paste0('cycle_',cy)]] <- keep&e$cycle==cy
for(cy in c('H','I','J')) masks[[paste0('without_',cy)]] <- keep&e$cycle!=cy
cutoffs <- seq(min(unlist(lapply(scores,function(s)s[keep]))),
 max(unlist(lapply(scores,function(s)s[keep])))+1)
grid <- seq(.10,.60,.01)
ratio <- function(a,b) if(b>0) a/b else NA_real_
met <- function(r,y,w) {
 tot<-sum(w);pos<-sum(w*y);found<-sum(w*y*r);ref<-sum(w*r)
 c(sensitivity_pct=100*ratio(found,pos),referral_pct=100*ratio(ref,tot),
   ppv_pct=100*ratio(found,ref),found_per1000=1000*ratio(found,tot),
   referrals_per1000=1000*ratio(ref,tot))
}
auc <- function(s,y,w) {
 k<-w>0; s<-s[k];y<-y[k];w<-w[k]
 ix<-order(s);s<-s[ix];y<-y[ix];w<-w[ix]
 end<-c(which(diff(s)!=0),length(s));cp<-cumsum(w*y)[end];cn<-cumsum(w*(1-y))[end]
 wp<-diff(c(0,cp));wn<-diff(c(0,cn));below<-c(0,head(cn,-1))
 ratio(sum(wp*(below+.5*wn)),sum(wp)*sum(wn))
}
ranktable <- function(p,y,w) {
 k<-w>0;ix<-order(p[k],decreasing=TRUE);pp<-p[k][ix];ww<-w[k][ix];yy<-y[k][ix]
 end<-c(which(diff(pp)!=0),length(pp))
 list(threshold=c(Inf,pp[end]),ref=c(0,cumsum(ww)[end]/sum(ww)),
      found=c(0,cumsum(ww*yy)[end]),total=sum(ww),positive=sum(ww*yy))
}
choose <- function(tab,target) {
 distance<-abs(tab$ref-target); ii<-which(distance<=min(distance)+1e-14)
 ii[which.min(tab$ref[ii])]
}
rankmet <- function(tab,i) {
 ref<-tab$ref[i]*tab$total;found<-tab$found[i]
 c(sensitivity_pct=100*ratio(found,tab$positive),referral_pct=100*tab$ref[i],
 ppv_pct=100*ratio(found,ref),found_per1000=1000*ratio(found,tab$total),
 referrals_per1000=1000*tab$ref[i])
}
# Meaningful boundary checks, including no referrals and tied score groups.
stopifnot(abs(auc(c(1,1,2),c(1,0,1),c(1,1,1))-.75)<1e-14)
tabtest<-ranktable(c(.5,.5,.1),c(1,0,1),c(1,1,2))
stopifnot(choose(tabtest,0)==1, choose(tabtest,.25)==1,
 choose(tabtest,1)==3, length(tabtest$threshold)==3)

one <- function(p,w,t0) {
 ans<-numeric()
 add<-function(pre,v) {names(v)<-paste(pre,names(v),sep='::');ans<<-c(ans,v)}
 for(sn in names(masks)) {
  k<-masks[[sn]];yy<-e$case[k];ww<-w[k];pp<-p[k]
  ss<-lapply(scores,function(s)s[k]);rt<-ranktable(pp,yy,ww)
  aa<-c(equation=auc(pp,yy,ww),Bang_self=auc(ss$Bang_self,yy,ww),Bang_composite=auc(ss$Bang_composite,yy,ww))
  aa<-c(aa,equation_minus_self=unname(aa['equation']-aa['Bang_self']),
    equation_minus_composite=unname(aa['equation']-aa['Bang_composite']),
    composite_minus_self=unname(aa['Bang_composite']-aa['Bang_self']))
  add(paste(sn,'auc',sep='::'),aa)
  st<-list(R0=met(as.numeric(pp>=t0),yy,ww))
  for(n in names(ss)) {
   v<-met(as.numeric(ss[[n]]>=5),yy,ww);st[[n]]<-v
   eq<-rankmet(rt,choose(rt,v['referral_pct']/100));st[[paste0('Eq_to_',n)]]<-eq
   add(paste(sn,'matched5',paste0('difference_',n),sep='::'),eq-v)
   add(paste(sn,'fixed',paste0('R0_minus_',n),sep='::'),st$R0-v)
  }
  for(n in names(st))add(paste(sn,'fixed',n,sep='::'),st[[n]])
  if(sn %in% c('common','strict_BP')) {
   for(n in names(ss)) for(t in cutoffs) {
    v<-met(as.numeric(ss[[n]]>=t),yy,ww)
    eq<-rankmet(rt,choose(rt,v['referral_pct']/100))
    pre<-paste(sn,'budget',n,paste0('k',t),sep='::')
    add(paste(pre,'Bang',sep='::'),v);add(paste(pre,'Equation',sep='::'),eq)
    add(paste(pre,'difference',sep='::'),eq-v)
   }
  }
  if(sn=='common') for(g in grid) add(paste(sn,'grid',sprintf('g%02d',round(100*g)),sep='::'),rankmet(rt,choose(rt,g)))
 }
 ans
}
coverage<-do.call(rbind,lapply(names(masks),function(n){k<-masks[[n]]
 data.frame(sample=n,n=sum(k),cases=sum(e$case[k]),weighted_positivity_pct=100*sum(e$wt[k]*e$case[k])/sum(e$wt[k]))}))
write.csv(coverage,file.path(O,'sample_coverage.csv'),row.names=FALSE)
pt<-one(p,e$wt,t0)
old<-read.csv(file.path(A,'S25_hypertension_strategy_comparisons.csv'))
old<-setNames(old$estimate,old$key)
check<-numeric()
for(n in c('Bang_self','Bang_composite','R0')) for(v in c('sensitivity','referral','ppv','found_per1000','referrals_per1000')) {
 newv<-if(v %in% c('sensitivity','referral','ppv'))paste0(v,'_pct') else v
 check<-c(check,abs(pt[paste('common','fixed',n,newv,sep='::')]-old[paste(n,'overall',v,sep='::')]))
}
stopifnot(all(is.finite(check)),max(check)<1e-9)
write.csv(data.frame(key=names(pt),estimate=unname(pt)),file.path(O,'point_estimates.csv'),row.names=FALSE)
write.csv(data.frame(check='v11_common_fixed_metrics',comparisons=length(check),maximum_error=max(check)),file.path(O,'v11_reproduction.csv'),row.names=FALSE)
M<-matrix(NA_real_,B,length(pt),dimnames=list(NULL,names(pt)))
start<-proc.time()[3]
for(b in seq_len(B)) {
 w<-e$wt*obj$multipliers[obj$idxeval,b]
 M[b,]<-one(pb[,b],w,bp$THRESHOLD[b])
 if(b%%50==0) {cat('replicate',b,'/',B,'elapsed_seconds',round(proc.time()[3]-start,1),'\n');flush.console()}
}
write.csv(data.frame(replicate=seq_len(B),M,check.names=FALSE),file.path(O,'bootstrap_replicates.csv'),row.names=FALSE)
S<-do.call(rbind,lapply(seq_along(pt),function(i){v<-M[,i];ok<-is.finite(v);n<-sum(ok)
 se<-if(n>=2)sd(v[ok]) else NA_real_;q<-if(n)quantile(v[ok],c(.025,.975),names=FALSE,type=7) else c(NA,NA)
 data.frame(key=names(pt)[i],estimate=unname(pt[i]),se=se,normal_lower=unname(pt[i])-1.96*se,
 normal_upper=unname(pt[i])+1.96*se,percentile_lower=q[1],percentile_upper=q[2],valid_B=n,requested_B=B,
 status=if(n==B)'complete' else if(n==0&&!is.finite(pt[i]))'undefined_denominator' else 'incomplete')
}))
write.csv(S,file.path(O,'summary_all.csv'),row.names=FALSE)
write.csv(S[grepl('::auc::',S$key),],file.path(O,'auc_summary.csv'),row.names=FALSE)
write.csv(S[grepl('::(fixed|matched5)::',S$key),],file.path(O,'strategy_summary.csv'),row.names=FALSE)
write.csv(S[grepl('::budget::',S$key),],file.path(O,'budget_summary.csv'),row.names=FALSE)
write.csv(S[grepl('::grid::',S$key),],file.path(O,'equation_grid_summary.csv'),row.names=FALSE)
family<-grepl('::budget::.*::difference::(sensitivity_pct|found_per1000)$',S$key)
stopifnot(all(S$valid_B[family]==B))
idx<-which(family&S$se>1e-12);zero<-which(family&S$se<=1e-12)
Z<-sweep(sweep(M[,idx,drop=FALSE],2,pt[idx],'-'),2,S$se[idx],'/')
mx<-apply(abs(Z),1,max);critical<-unname(quantile(mx,.95,type=7))
sim<-S[family,];sim$simultaneous_critical<-critical
sim$simultaneous_lower<-sim$estimate-critical*sim$se;sim$simultaneous_upper<-sim$estimate+critical*sim$se
sim$family<-'both_samples_both_Bang_definitions_all_integer_cutoffs_sensitivity_and_found'
write.csv(sim,file.path(O,'budget_simultaneous_intervals.csv'),row.names=FALSE)
write.csv(data.frame(replicate=seq_len(B),maximum_absolute_standardised_deviation=mx),file.path(O,'simultaneous_maxima.csv'),row.names=FALSE)
core<-grepl('^common::auc::|^common::matched5::|^common::fixed::',S$key)
coreok<-all(S$valid_B[core]==B)
write.csv(data.frame(requested_B=B,core_complete=coreok,metric_columns=ncol(M),
 complete_columns=sum(S$valid_B==B),undefined_columns=sum(S$valid_B==0),incomplete_columns=sum(S$valid_B>0&S$valid_B<B),
 simultaneous_family_columns=sum(family),positive_SD_columns=length(idx),zero_SD_columns=length(zero),
 simultaneous_critical=critical,seconds=proc.time()[3]-start),file.path(O,'run_status.csv'),row.names=FALSE)
if(!coreok)stop('Core comparisons incomplete: no inferential headline permitted. See summary_all.csv.')
cat('COMPLETE',B,'replicates;',length(pt),'metrics\n')
print(S[S$key=='common::auc::equation_minus_composite',])
print(S[S$key=='common::matched5::difference_Bang_composite::sensitivity_pct',])
