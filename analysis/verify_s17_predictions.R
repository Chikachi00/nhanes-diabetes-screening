args <- commandArgs(trailingOnly=TRUE)
base <- "outputs/diabetes_manuscript_revision_20261104_v2/analysis"
raw <- if(length(args)>=1)args[1] else file.path(base,"raw")
out <- if(length(args)>=2)args[2] else file.path(base,"s17_rerun")
source(if(length(args)>=3)args[3] else file.path(base,"predict_s17_derived.R"))
rd <- function(m)foreign::read.xport(file.path(raw,paste0(m,"_L.xpt")))
d <- read.csv(file.path(out,"S17_external_predictions_rerun.csv"))
de <- rd("DEMO");bm <- rd("BMX");bp <- rd("BPQ");pa <- rd("PAQ")
de <- de[match(d$SEQN,de$SEQN),];bm <- bm[match(d$SEQN,bm$SEQN),];bp <- bp[match(d$SEQN,bp$SEQN),];pa <- pa[match(d$SEQN,pa$SEQN),]
pw <- function(q,u) {
  u <- toupper(trimws(as.character(u)));q <- as.numeric(q)
  ans <- rep(NA_real_,length(q))
  for(un in c("D","W","M","Y")) {k <- (u==un & q<7777)%in%TRUE;ans[k] <- q[k]*switch(un,D=7,W=1,M=1/4.345,Y=1/52.18)}
  ans[q%in%0] <- 0;ans
}
pred <- data.frame(SEQN=d$SEQN,age=de$RIDAGEYR,whtr=bm$BMXWAIST/bm$BMXHT,male=as.numeric(de$RIAGENDR==1),htn=ifelse(bp$BPQ020==1,1,ifelse(bp$BPQ020==2,0,NA)),inactive=as.numeric(pw(pa$PAD790Q,pa$PAD790U)+pw(pa$PAD810Q,pa$PAD810U)<5))
stopifnot(complete.cases(pred))
checks <- list()
for(m in c("V5nf","F3")) {
  z <- predict_s17_derived(file.path(out,paste0("S17_",m,"_source_rerun_model.rds")),pred)
  checks[[m]] <- data.frame(model=m,n=nrow(d),max_probability_difference=max(abs(z$probability-d[[paste0("p_",m)]])),R0_mismatches=sum(z$R0!=d[[paste0("ref_",m,"_R0")]]),R1_mismatches=sum(z$R1!=d[[paste0("ref_",m,"_R1")]]))
}
write.csv(pred,file.path(out,"S17_external_scoring_inputs.csv"),row.names=FALSE)
write.csv(do.call(rbind,checks),file.path(out,"S17_saved_model_prediction_checks.csv"),row.names=FALSE)
print(do.call(rbind,checks))
