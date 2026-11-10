# D3: self-reported hypertension and family history by age group and access (evaluation 2013-2018)
args<-commandArgs(TRUE);V<-args[1];out<-args[2];dir.create(out,showWarnings=FALSE,recursive=TRUE)
suppressMessages(library(survey));options(survey.lonely.psu='adjust')
obj<-readRDS(file.path(V,'reproducibility/main_inputs/prepared.rds'));e<-obj$eval
des<-svydesign(ids=~psu,strata=~strata,weights=~wt,nest=TRUE,data=e)
r<-svyby(~htn+famhx,~age_grp+access,des,svymean);print(r,digits=3)
write.csv(as.data.frame(r),file.path(out,'D3_htn_famhx_by_age_access.csv'),row.names=FALSE)
