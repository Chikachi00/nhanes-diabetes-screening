# Reviewer-response analyses C (subgroup and near-cut-off calibration) and D (self-report completeness by access)
args<-commandArgs(TRUE);V<-args[1];out<-args[2];dir.create(out,showWarnings=FALSE,recursive=TRUE)
suppressMessages({library(survey);library(haven)})
options(survey.lonely.psu='adjust')
obj<-readRDS(file.path(V,'reproducibility/main_inputs/prepared.rds'));e<-obj$eval
pr<-read.csv(file.path(V,'reproducibility/main_archive/evaluation_predictions.csv'))
stopifnot(all(pr$SEQN==e$SEQN));e$p<-pr$probability;e$age<-round((e$a+5)*10)
lk<-read.csv(file.path(V,'reproducibility/main_archive/locked_strategy_parameters.csv'));t0<-lk$value[lk$term=='THRESHOLD'];to<-lk$value[lk$term=='THRESHOLD_OLD']
# full examined PSU frame: keep all evaluation-period PSUs via a zero-weight-free design on the analytic sample;
# domain estimation by subset() keeps the design's PSUs.
des<-svydesign(ids=~psu,strata=~strata,weights=~wt,nest=TRUE,data=e)
grp<-list(overall=rep(TRUE,nrow(e)),age20_39=e$age_grp=='20-39',age40_59=e$age_grp=='40-59',age60=e$age_grp=='>=60',high=e$access=='high',low=e$access=='low')
oe<-function(k){ds<-subset(des,k);r<-svyratio(~case,~p,ds);ci<-confint(r,df=degf(des));
 c(n=sum(k),pos=sum(e$case[k]),mean_pred=100*sum(e$wt[k]*e$p[k])/sum(e$wt[k]),obs=100*sum(e$wt[k]*e$case[k])/sum(e$wt[k]),OE=coef(r)[[1]],lo=ci[1],hi=ci[2])}
O<-do.call(rbind,lapply(names(grp),function(g)data.frame(group=g,t(oe(grp[[g]])))))
write.csv(O,file.path(out,'C1_subgroup_OE_Taylor.csv'),row.names=FALSE);print(O,digits=3)
# near-cut-off bands: predictions within a factor 1.5 of the cut-off, below and above
band<-function(lo,hi,k,label){kk<-k&e$p>=lo&e$p<hi;if(sum(kk)<5)return(NULL);ds<-subset(des,kk)
 m<-svyciprop(~case,ds,method='logit',df=degf(des));ci<-attr(m,'ci')
 data.frame(band=label,n=sum(kk),pos=sum(e$case[kk]),mean_pred=100*sum(e$wt[kk]*e$p[kk])/sum(e$wt[kk]),obs=100*coef(m)[[1]],lo=100*ci[[1]],hi=100*ci[[2]])}
rows<-list()
for(g in names(grp)){k<-grp[[g]]
 rows[[length(rows)+1]]<-cbind(cutoff='R0 0.0226',group=g,band(t0/1.5,t0,k,'below (t/1.5 to t)'))
 rows[[length(rows)+1]]<-cbind(cutoff='R0 0.0226',group=g,band(t0,t0*1.5,k,'above (t to 1.5t)'))
 rows[[length(rows)+1]]<-cbind(cutoff='R0 0.0226',group=g,band(t0/1.5,t0*1.5,k,'both (t/1.5 to 1.5t)'))}
for(g in c('age60')){k<-grp[[g]]
 for(a in c('all','high','low')){kk<-k&(a=='all'|e$access==a)
 rows[[length(rows)+1]]<-cbind(cutoff='R1 older 0.0147',group=paste(g,a,sep='_'),band(to/1.5,to,kk,'below (t/1.5 to t)'))
 rows[[length(rows)+1]]<-cbind(cutoff='R1 older 0.0147',group=paste(g,a,sep='_'),band(to,to*1.5,kk,'above (t to 1.5t)'))
 rows[[length(rows)+1]]<-cbind(cutoff='R1 older 0.0147',group=paste(g,a,sep='_'),band(to/1.5,to*1.5,kk,'both (t/1.5 to 1.5t)'))}}
Bd<-do.call(rbind,rows);write.csv(Bd,file.path(out,'C2_near_cutoff_calibration.csv'),row.names=FALSE);print(Bd,digits=3)
# D: completeness of self-reported items by access, evaluation target population with HbA1c and known access
fb<-read.csv(file.path(V,'reproducibility/main_inputs/full_base.csv'));fb<-fb[fb$period=='eval'&fb$scope==1&!is.na(fb$hba1c)&!is.na(fb$ins)&!is.na(fb$usual),]
rd<-function(p,cols)do.call(rbind,lapply(c('H','I','J'),function(cy){x<-as.data.frame(read_xpt(file.path(V,'raw',paste0(p,'_',cy,'.xpt'))));x[,cols]}))
mcq<-rd('MCQ',c('SEQN','MCQ300C'));bpq<-rd('BPQ',c('SEQN','BPQ020'))
fb<-merge(fb,mcq,by='SEQN',all.x=TRUE);fb<-merge(fb,bpq,by='SEQN',all.x=TRUE)
fb$access<-ifelse(fb$ins==1&fb$usual==1,'high','low')
fb$famhx_unknown<-as.numeric(!(fb$MCQ300C%in%c(1,2)))
fb$htn_unknown<-as.numeric(!(fb$BPQ020%in%c(1,2)))
fb$in_common<-as.numeric(fb$SEQN%in%e$SEQN)
fb$famhx_yes<-ifelse(fb$MCQ300C%in%c(1,2),as.numeric(fb$MCQ300C==1),NA)
fb$htn_yes<-ifelse(fb$BPQ020%in%c(1,2),as.numeric(fb$BPQ020==1),NA)
dfb<-svydesign(ids=~psu,strata=~strata,weights=~wt,nest=TRUE,data=fb)
res<-list()
for(v in c('famhx_unknown','htn_unknown','in_common','famhx_yes','htn_yes')){
 f<-as.formula(paste0('~',v));b<-svyby(f,~access,dfb,svymean,na.rm=TRUE,covmat=TRUE);ct<-svycontrast(b,c(high=1,low=-1))
 res[[v]]<-data.frame(item=v,high=100*coef(b)[['high']],low=100*coef(b)[['low']],diff_high_minus_low=100*coef(ct)[[1]],lo=100*(coef(ct)[[1]]-qt(.975,degf(dfb))*SE(ct)[[1]]),hi=100*(coef(ct)[[1]]+qt(.975,degf(dfb))*SE(ct)[[1]]),
   n_high=sum(fb$access=='high'&!is.na(fb[[v]])),n_low=sum(fb$access=='low'&!is.na(fb[[v]])),unknown_n_high=if(grepl('unknown',v))sum(fb[[v]][fb$access=='high']) else NA,unknown_n_low=if(grepl('unknown',v))sum(fb[[v]][fb$access=='low']) else NA)}
D<-do.call(rbind,res);write.csv(D,file.path(out,'D_selfreport_by_access.csv'),row.names=FALSE);print(D,digits=3)
# family history yes by access among positives in common sample (awareness among cases)
dc<-subset(des,e$case==1);print(svyby(~famhx+htn,~access,dc,svymean))
write.csv(as.data.frame(svyby(~famhx+htn,~access,dc,svymean)),file.path(out,'D2_items_among_positives.csv'),row.names=FALSE)
cat('n target',nrow(fb),'\n')
