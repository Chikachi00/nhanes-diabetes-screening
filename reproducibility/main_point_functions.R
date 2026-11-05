# Adapted from the frozen strategies.R: paths only; fitting and rule functions unchanged.
C_L2<-1;DEV_CYCLES<-c('E','F','G')
for(ex in parse(file.path(old,'source_snapshot/refine_v5_equation.R'),encoding='UTF-8'))
 if(is.call(ex)&&identical(ex[[1]],as.name('<-'))&&is.symbol(ex[[2]])&&as.character(ex[[2]])=='fit_l2')eval(ex)
fit_orig<-fit_l2
fit_l2<-function(...)withCallingHandlers(fit_orig(...),warning=function(w)stop(conditionMessage(w)))
obj<-readRDS(file.path(old,'prepared.rds'));d<-obj$dev;e<-obj$eval
d$age<-round((d$a+5)*10);e$age<-round((e$a+5)*10)
vars<-c('a','w','male','famhx','htn','inactive');X<-as.matrix(d[,vars]);XE<-as.matrix(e[,vars])
thr_ref<-function(p,w,t){o<-order(p,decreasing=TRUE);p<-p[o];w<-w[o];en<-c(which(diff(p)!=0),length(p));p[en[which.min(abs(cumsum(w)[en]/sum(w)-t))]]}
thr_sens<-function(p,y,w,t){k<-y==1&w>0;o<-order(p[k],decreasing=TRUE);pp<-p[k][o];ww<-w[k][o];en<-c(which(diff(pp)!=0),length(pp));pp[en[which(cumsum(ww)[en]/sum(ww)>=t-1e-12)[1]]]}
develop<-function(w){
 lp<-rep(NA_real_,nrow(d));folds<-list()
 for(cy in DEV_CYCLES){tr<-d$cycle!=cy;cf<-fit_l2(X[tr,,drop=FALSE],d$case[tr],w[tr]);folds[[cy]]<-cf;lp[!tr]<-as.vector(cf[1]+X[!tr,,drop=FALSE]%*%cf[-1])}
 k<-w>0;target<-sum(w*d$ref_bang)/sum(w);ko<-k&d$age>=60&d$case==1
 ts<-sum(w[ko]*d$ref_bang[ko])/sum(w[ko]);t0<-thr_ref(plogis(lp[k]),w[k],target)
 to<-min(t0,thr_sens(plogis(lp[k&d$age>=60]),d$case[k&d$age>=60],w[k&d$age>=60],ts))
 cf<-fit_l2(X,d$case,w)
 list(beta=cf,t0=t0,to=to,target=target,target_old_sens=ts,folds=folds,lp_oof=lp,lp=as.vector(cf[1]+XE%*%cf[-1]))
}
flags<-function(f){p<-plogis(f$lp);r0<-as.numeric(p>=f$t0);r1<-as.numeric(p>=ifelse(e$age>=60,f$to,f$t0));r4<-ifelse(e$age>=60,1,r0);stopifnot(all(r1>=r0),all(r1[e$age<60]==r0[e$age<60]));list(Bang=e$ref_bang,S0=r0,S1=r1,S4=r4)}
masks<-obj$masks
metrics<-function(ff,w){out<-numeric();y<-e$case
 for(g in names(masks)){k<-masks[[g]];ww<-w[k];yy<-y[k];rr<-ff[k];tot<-sum(ww);pos<-sum(ww*yy);neg<-tot-pos;found<-sum(ww*yy*rr);ref<-sum(ww*rr)
  z<-c(sensitivity=found/pos,specificity=(neg-ref+found)/neg,ppv=found/ref,referral=ref/tot,found_per1000=1000*found/tot,missed_per1000=1000*(pos-found)/tot,cases_per1000=1000*pos/tot,referrals_per1000=1000*ref/tot)
  names(z)<-paste(g,names(z),sep='::');out<-c(out,z)}
 c(out,'access_gap::sensitivity'=unname(out['high::sensitivity']-out['low::sensitivity']))}
comparisons<-list(S0_minus_Bang=c('S0','Bang'),S1_minus_Bang=c('S1','Bang'),S1_minus_S0=c('S1','S0'),S4_minus_S1=c('S4','S1'),S4_minus_Bang=c('S4','Bang'))
flatten<-function(ff,w){st<-lapply(ff,metrics,w=w);out<-numeric();for(n in names(st)){z<-st[[n]];names(z)<-paste(n,names(z),sep='::');out<-c(out,z)};for(n in names(comparisons)){cc<-comparisons[[n]];z<-st[[cc[1]]]-st[[cc[2]]];names(z)<-paste(n,names(z),sep='::');out<-c(out,z)};out}
point<-develop(d$wt);pf<-flags(point);pointm<-flatten(pf,e$wt)
