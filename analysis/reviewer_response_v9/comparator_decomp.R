# Reviewer-response analyses A (Bang mapping stability), B (same-information model), E (equal budget)
# Uses frozen prepared.rds, original fit_l2 and threshold rules, and the original Rao-Wu multipliers.
args<-commandArgs(TRUE); V<-args[1]; out<-args[2]; B<-as.integer(ifelse(length(args)>2,args[3],1000))
suppressMessages({library(parallel);library(haven)})
snap<-file.path(V,'reproducibility/main_inputs/source_snapshot/refine_v5_equation.R')
C_L2<-1;DEV_CYCLES<-c('E','F','G')
for(ex in parse(snap,encoding='UTF-8')) if(is.call(ex)&&identical(ex[[1]],as.name('<-'))&&is.symbol(ex[[2]])&&as.character(ex[[2]])=='fit_l2')eval(ex)
fit_orig<-fit_l2; fit_l2<-function(...)withCallingHandlers(fit_orig(...),warning=function(w)stop(conditionMessage(w)))
obj<-readRDS(file.path(V,'reproducibility/main_inputs/prepared.rds'));d<-obj$dev;e<-obj$eval
d$age<-round((d$a+5)*10);e$age<-round((e$a+5)*10)
thr_ref<-function(p,w,t){o<-order(p,decreasing=TRUE);p<-p[o];w<-w[o];en<-c(which(diff(p)!=0),length(p));p[en[which.min(abs(cumsum(w)[en]/sum(w)-t))]]}
thr_sens<-function(p,y,w,t){k<-y==1&w>0;o<-order(p[k],decreasing=TRUE);pp<-p[k][o];ww<-w[k][o];en<-c(which(diff(pp)!=0),length(pp));pp[en[which(cumsum(ww)[en]/sum(ww)>=t-1e-12)[1]]]}
# ---- rebuild Bang items from source inputs ----
fb<-read.csv(file.path(V,'reproducibility/main_inputs/full_base.csv'))
rd<-function(p,cols)do.call(rbind,lapply(c('E','F','G','H','I','J'),function(cy){x<-as.data.frame(read_xpt(file.path(V,'raw',paste0(p,'_',cy,'.xpt'))));x[,cols]}))
paq<-rd('PAQ',c('SEQN','PAQ605','PAQ610','PAQ620','PAQ625','PAQ635','PAQ640','PAQ650','PAQ655','PAQ665','PAQ670'))
dd<-function(yn,days)ifelse(yn==2,0,ifelse(yn==1&days%in%1:7,days,NA))
paq$d1<-dd(paq$PAQ605,paq$PAQ610);paq$d2<-dd(paq$PAQ620,paq$PAQ625);paq$d3<-dd(paq$PAQ635,paq$PAQ640);paq$d4<-dd(paq$PAQ650,paq$PAQ655);paq$d5<-dd(paq$PAQ665,paq$PAQ670)
items<-function(z){
 m<-merge(z[,c('SEQN','cycle','bang')],fb[,c('SEQN','age','male','bmi','waist','famhx','htn','inactive')],by='SEQN',sort=FALSE,suffixes=c('','.fb'))
 m<-merge(m,paq[,c('SEQN','d1','d2','d3','d4','d5')],by='SEQN',sort=FALSE)
 m<-m[match(z$SEQN,m$SEQN),]
 m$age_pts<-ifelse(m$age<40,0,ifelse(m$age<50,1,ifelse(m$age<60,2,3)))
 m$bmi_pts<-ifelse(m$bmi<25,0,ifelse(m$bmi<30,1,ifelse(m$bmi<40,2,3)))
 wi<-m$waist/2.54
 m$waist_pts<-ifelse(m$male==1,ifelse(wi<37,0,ifelse(wi<40,1,ifelse(wi<50,2,3))),ifelse(wi<31.5,0,ifelse(wi<35,1,ifelse(wi<49,2,3))))
 m$obese_pts<-pmax(m$bmi_pts,m$waist_pts); m$active<-1-m$inactive
 D<-as.matrix(m[,c('d1','d2','d3','d4','d5')])
 m$active_any<-ifelse(rowSums(D>=1,na.rm=TRUE)>0,1,ifelse(rowSums(is.na(D))==0,0,NA))
 base<-m$age_pts+m$male+m$famhx+m$htn
 m$bang_main<-base+m$obese_pts-m$active
 m$bang_any<-base+m$obese_pts-m$active_any
 m$bang_bmi<-base+m$bmi_pts-m$active
 m$bang_waist<-base+m$waist_pts-m$active
 stopifnot(all(m$SEQN==z$SEQN),max(abs(m$bang_main-z$bang))==0,!anyNA(m$bang_any),!anyNA(m$bang_bmi),!anyNA(m$bang_waist))
 m}
Id<-items(d);Ie<-items(e)
# same-information design: Bang items in their categorical form
sx<-function(m)cbind(age1=as.numeric(m$age_pts==1),age2=as.numeric(m$age_pts==2),age3=as.numeric(m$age_pts==3),
  ob1=as.numeric(m$obese_pts==1),ob2=as.numeric(m$obese_pts==2),ob3=as.numeric(m$obese_pts==3),
  male=m$male,famhx=m$famhx,htn=m$htn,active=m$active)
vars<-c('a','w','male','famhx','htn','inactive')
XS<-list(main=list(dev=as.matrix(d[,vars]),eval=as.matrix(e[,vars])),same=list(dev=sx(Id),eval=sx(Ie)))
develop<-function(X,XE,w){
 lp<-rep(NA_real_,nrow(d))
 for(cy in DEV_CYCLES){tr<-d$cycle!=cy;cf<-fit_l2(X[tr,,drop=FALSE],d$case[tr],w[tr]);lp[!tr]<-as.vector(cf[1]+X[!tr,,drop=FALSE]%*%cf[-1])}
 k<-w>0;target<-sum(w*d$ref_bang)/sum(w);ko<-k&d$age>=60&d$case==1
 ts<-sum(w[ko]*d$ref_bang[ko])/sum(w[ko]);t0<-thr_ref(plogis(lp[k]),w[k],target)
 to<-min(t0,thr_sens(plogis(lp[k&d$age>=60]),d$case[k&d$age>=60],w[k&d$age>=60],ts))
 cf<-fit_l2(X,d$case,w)
 list(beta=cf,t0=t0,to=to,lp=as.vector(cf[1]+XE%*%cf[-1]))}
masks<-obj$masks
auc_w<-function(s,y,w){k<-w>0;s<-s[k];y<-y[k];w<-w[k];u<-sort(unique(s));g<-match(s,u)
 wp<-tapply(w*y,factor(g,levels=seq_along(u)),sum);wn<-tapply(w*(1-y),factor(g,levels=seq_along(u)),sum);wp[is.na(wp)]<-0;wn[is.na(wn)]<-0
 below<-cumsum(wn)-wn; sum(wp*(below+0.5*wn))/(sum(wp)*sum(wn))}
met<-function(r,w){y<-e$case;o<-numeric()
 for(g in names(masks)){k<-masks[[g]];ww<-w[k];yy<-y[k];rr<-r[k];tot<-sum(ww);pos<-sum(ww*yy);f<-sum(ww*yy*rr);rf<-sum(ww*rr)
  z<-c(sensitivity=100*f/pos,referral=100*rf/tot,ppv=100*f/rf,found_per1000=1000*f/tot,referrals_per1000=1000*rf/tot);names(z)<-paste(g,names(z),sep='::');o<-c(o,z)}
 c(o,'access_gap::sensitivity'=unname(o['high::sensitivity']-o['low::sensitivity'])) }
one<-function(wd,we){
 fm<-develop(XS$main$dev,XS$main$eval,wd); fs<-develop(XS$same$dev,XS$same$eval,wd)
 pm<-plogis(fm$lp); ps<-plogis(fs$lp)
 R<-list(Bang=Ie$bang_main>=5,Bang_any=Ie$bang_any>=5,Bang_bmi=Ie$bang_bmi>=5,Bang_waist=Ie$bang_waist>=5,
  R0=pm>=fm$t0,R1=pm>=ifelse(e$age>=60,fm$to,fm$t0),Same0=ps>=fs$t0,Same1=ps>=ifelse(e$age>=60,fs$to,fs$t0))
 bref<-sum(we*R$Bang)/sum(we)
 R$EqMain<-pm>=thr_ref(pm[we>0],we[we>0],bref); R$EqSame<-ps>=thr_ref(ps[we>0],we[we>0],bref)
 R<-lapply(R,as.numeric)
 st<-lapply(R,met,w=we); o<-numeric()
 for(n in names(st)){z<-st[[n]];names(z)<-paste(n,names(z),sep='::');o<-c(o,z)}
 cmp<-list(R0_minus_Bang=c('R0','Bang'),R0_minus_Bang_any=c('R0','Bang_any'),R0_minus_Bang_bmi=c('R0','Bang_bmi'),R0_minus_Bang_waist=c('R0','Bang_waist'),
  R1_minus_Bang=c('R1','Bang'),R1_minus_Bang_any=c('R1','Bang_any'),R1_minus_Bang_bmi=c('R1','Bang_bmi'),R1_minus_Bang_waist=c('R1','Bang_waist'),
  Same0_minus_Bang=c('Same0','Bang'),R0_minus_Same0=c('R0','Same0'),Same1_minus_Bang=c('Same1','Bang'),R1_minus_Same1=c('R1','Same1'),
  EqMain_minus_Bang=c('EqMain','Bang'),EqSame_minus_Bang=c('EqSame','Bang'),EqMain_minus_EqSame=c('EqMain','EqSame'),R0_minus_EqMain=c('R0','EqMain'),R1_minus_R0=c('R1','R0'))
 for(n in names(cmp)){z<-st[[cmp[[n]][1]]]-st[[cmp[[n]][2]]];names(z)<-paste(n,names(z),sep='::');o<-c(o,z)}
 a<-c(main=auc_w(pm,e$case,we),same=auc_w(ps,e$case,we),bang=auc_w(Ie$bang_main,e$case,we))
 o<-c(o,'AUC::main'=a[['main']],'AUC::same'=a[['same']],'AUC::bang'=a[['bang']],'AUC::main_minus_same'=a[['main']]-a[['same']],'AUC::same_minus_bang'=a[['same']]-a[['bang']],
  'thr::main_t0'=fm$t0,'thr::main_to'=fm$to,'thr::same_t0'=fs$t0,'thr::same_to'=fs$to)
 list(o=o,fm=fm,fs=fs)}
pt<-one(d$wt,e$wt)
locked<-read.csv(file.path(V,'reproducibility/main_archive/locked_strategy_parameters.csv'))
stopifnot(abs(pt$fm$t0-locked$value[locked$term=='THRESHOLD'])<1e-12,abs(pt$fm$to-locked$value[locked$term=='THRESHOLD_OLD'])<1e-12)
dir.create(out,showWarnings=FALSE,recursive=TRUE)
write.csv(data.frame(term=names(pt$fs$beta),value=pt$fs$beta),file.path(out,'same_info_coefficients.csv'),row.names=FALSE)
write.csv(data.frame(key=names(pt$o),estimate=pt$o),file.path(out,'point.csv'),row.names=FALSE)
cat('point done\n');flush.console()
arch<-read.csv(file.path(V,'reproducibility/main_archive/bootstrap_full_metrics.csv'),check.names=FALSE)
RNGkind("L'Ecuyer-CMRG");set.seed(1)
res<-mclapply(seq_len(B),function(b){
  if(b %% 50L == 0L){cat("Replicate",b,"of",B,"\n");flush.console()}
  mul<-obj$multipliers[,b];wd<-d$wt*mul[obj$idxdev];we<-e$wt*mul[obj$idxeval]
  r<-tryCatch(one(wd,we),error=function(x)NULL);if(is.null(r))return(NULL)
  chk<-max(abs(c(r$o['R0::overall::sensitivity']/100-arch[b,'S0::overall::sensitivity'],r$o['R1::overall::sensitivity']/100-arch[b,'S1::overall::sensitivity'],r$o['Bang::overall::sensitivity']/100-arch[b,'Bang::overall::sensitivity'])))
  c(replicate=b,r$o,archive_check=chk)},mc.cores=if(.Platform$OS.type == "windows") 1L else 4L)
ok<-!vapply(res,is.null,TRUE);M<-do.call(rbind,res[ok])
write.csv(M,file.path(out,'replicates.csv'),row.names=FALSE)
cat('valid',sum(ok),'max archive mismatch',max(M[,'archive_check']),'\n')
keys<-names(pt$o)
S<-do.call(rbind,lapply(keys,function(k){v<-M[,k];v<-v[is.finite(v)];s<-sd(v);q<-quantile(v,c(.025,.975),names=FALSE)
 data.frame(key=k,estimate=pt$o[[k]],se=s,normal_lo=pt$o[[k]]-1.96*s,normal_hi=pt$o[[k]]+1.96*s,pct_lo=q[1],pct_hi=q[2],share_le0=mean(v<=1e-12),valid_B=length(v))}))
write.csv(S,file.path(out,'summary.csv'),row.names=FALSE)
cat('done\n')
