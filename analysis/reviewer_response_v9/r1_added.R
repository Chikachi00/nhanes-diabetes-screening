# Reviewer-response analysis F: positivity among older adults referred by R1 but not by R0 (evaluation 2013-2018)
args<-commandArgs(TRUE);V<-args[1];out<-args[2];dir.create(out,showWarnings=FALSE,recursive=TRUE)
suppressMessages(library(survey));options(survey.lonely.psu='adjust')
e<-readRDS(file.path(V,'reproducibility/main_inputs/prepared.rds'))$eval
pr<-read.csv(file.path(V,'reproducibility/main_archive/evaluation_predictions.csv'))
stopifnot(all(pr$SEQN==e$SEQN));e$p<-pr$probability
lk<-read.csv(file.path(V,'reproducibility/main_archive/locked_strategy_parameters.csv'))
t0<-lk$value[lk$term=='THRESHOLD'];to<-lk$value[lk$term=='THRESHOLD_OLD']
e$added<-as.numeric(pr$S1==1&pr$S0==0)
stopifnot(all(e$added==as.numeric(e$age_grp=='>=60'&e$p>=to&e$p<t0)))
des<-svydesign(ids=~psu,strata=~strata,weights=~wt,nest=TRUE,data=e)
f<-function(k,lab){ds<-subset(des,k);m<-svyciprop(~case,ds,method='logit',df=degf(des));ci<-attr(m,'ci')
 data.frame(group=lab,n=sum(k),positives=sum(e$case[k]),pop_share_pct=100*sum(e$wt[k])/sum(e$wt),
  positivity_pct=100*coef(m)[[1]],lower95=100*ci[[1]],upper95=100*ci[[2]])}
R<-rbind(f(e$added==1,'R1-added (age>=60, R1 yes, R0 no)'),
 f(e$age_grp=='>=60'&pr$S0==1,'age>=60 referred by R0'),
 f(e$age_grp=='>=60'&pr$S1==0,'age>=60 not referred by R1'))
write.csv(R,file.path(out,'F1_R1_added_positivity.csv'),row.names=FALSE);print(R,digits=3)
