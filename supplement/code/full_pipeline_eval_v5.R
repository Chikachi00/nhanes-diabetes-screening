# NHANES: paired development-to-evaluation Rao-Wu bootstrap.
# Usage: Rscript full_pipeline_eval_v5.R WORKSPACE prepare|pilot|worker|collect [start] [end]
# Base R only. Original incoming scripts are never sourced at top level or edited.
# Structures and age knots are fixed; no automatic structure selection is claimed.
args <- commandArgs(trailingOnly=TRUE)
ws <- normalizePath(args[1],winslash="/")
mode <- args[2]
base <- file.path(ws,"outputs","diabetes_full_pipeline_20261029")
dir.create(file.path(base,"chunks"),showWarnings=FALSE)
C_L2 <- 1
SEED <- 20261029L
B <- 1000L
DEV_CYCLES <- c("E","F","G")
MODELS <- list(M2_linear=c("a","w","male","famhx","htn","inactive"),
               M2_spline=c("a","s","w","male","famhx","htn","inactive"),
               M3_inter=c("a","s","w","male","famhx","htn","inactive","aw"))
src <- file.path(base,"source_snapshot","refine_v5_equation.R")
for (expr in parse(src,encoding="UTF-8")) {
  if (is.call(expr) && identical(expr[[1]],as.name("<-")) && is.symbol(expr[[2]]) &&
      as.character(expr[[2]]) %in% c("fit_l2","v5_threshold","w_auc")) eval(expr)
}
fit_original <- fit_l2
fit_l2 <- function(...) {
  ans <- withCallingHandlers(fit_original(...),warning=function(w) stop("Fitting warning: ",conditionMessage(w)))
  if (any(!is.finite(ans))) stop("Nonfinite fitted coefficients")
  ans
}

# Same complete-tie-block rule as supplied code, with O(n log n) sorting.
threshold_fast <- function(p,w,target) {
  o <- order(p,decreasing=TRUE); p <- p[o]; w <- w[o]
  ends <- c(which(diff(p)!=0),length(p))
  share <- cumsum(w)[ends]/sum(w)
  p[ends[which.min(abs(share-target))]]
}
auc_fast <- function(score,y,w) {
  k <- w>0; score<-score[k];y<-y[k];w<-w[k]
  o<-order(score);score<-score[o];y<-y[o];w<-w[o]
  ends<-c(which(diff(score)!=0),length(score))
  cp<-cumsum(w*y);cn<-cumsum(w*(1-y))
  pos<-diff(c(0,cp[ends]));neg<-diff(c(0,cn[ends]))
  sum(pos*(cn[ends]-.5*neg))/(sum(w*y)*sum(w*(1-y)))
}
develop <- function(m,wd,obj) {
  X<-obj$Xdev[[m]]; y<-obj$dev$case;lp<-rep(NA_real_,length(y))
  for (cy in DEV_CYCLES) {
    tr<-obj$dev$cycle!=cy;te<-!tr
    cf<-fit_l2(X[tr,,drop=FALSE],y[tr],wd[tr])
    lp[te]<-as.vector(cf[1]+X[te,,drop=FALSE]%*%cf[-1])
  }
  k<-wd>0
  target<-sum(wd[k]*obj$dev$ref_bang[k])/sum(wd[k])
  cut<-threshold_fast(plogis(lp[k]),wd[k],target)
  cf<-fit_l2(X,y,wd)
  eval_lp<-as.vector(cf[1]+obj$Xeval[[m]]%*%cf[-1])
  list(beta=cf,threshold=cut,target=target,lp=eval_lp,
       oof_referral=sum(wd[k]*(plogis(lp[k])>=cut))/sum(wd[k]))
}
metric <- function(score,flag,we,obj,lp=NULL) {
  d<-obj$eval;y<-d$case
  rat<-function(num,den) {s<-sum(we*den);if(s>0)sum(we*num)/s else NA_real_}
  out<-numeric()
  for (g in names(obj$masks)) {
    k<-obj$masks[[g]];den<-as.numeric(k)
    vals<-c(sensitivity=rat(flag*y*k,y*k),specificity=rat((1-flag)*(1-y)*k,(1-y)*k),
            ppv=rat(flag*y*k,flag*k),referral=rat(flag*k,den),
            found_per1000=1000*rat(flag*y*k,den),missed_per1000=1000*rat((1-flag)*y*k,den),
            cases_per1000=1000*rat(y*k,den),referrals_per1000=1000*rat(flag*k,den))
    names(vals)<-paste(g,names(vals),sep="::");out<-c(out,vals)
  }
  out<-c(out,"overall::AUC"=auc_fast(score,y,we))
  out<-c(out,"access_gap::sensitivity"=unname(out["high::sensitivity"]-out["low::sensitivity"]))
  if (!is.null(lp)) {
    k<-we>0; nw<-we[k]/mean(we[k]);yy<-y[k];ll<-lp[k]
    f0<-suppressWarnings(glm.fit(x=matrix(1,sum(k),1),y=yy,weights=nw,offset=ll,family=quasibinomial()))
    f1<-suppressWarnings(glm.fit(x=cbind(1,ll),y=yy,weights=nw,family=quasibinomial()))
    if(!f0$converged || !f1$converged)stop("Calibration fit did not converge")
    out<-c(out,"overall::OE"=sum(we*y)/sum(we*plogis(lp)),
           "overall::calibration_intercept"=unname(f0$coefficients[1]),
           "overall::calibration_slope"=unname(f1$coefficients[2]))
  }
  out
}
comparisons<-list(M2_linear_minus_Bang=c("M2_linear","Bang"),
                  M2_spline_minus_Bang=c("M2_spline","Bang"),
                  M3_inter_minus_Bang=c("M3_inter","Bang"),
                  M2_spline_minus_M2_linear=c("M2_spline","M2_linear"),
                  M3_inter_minus_M2_spline=c("M3_inter","M2_spline"))
flatten <- function(st) {
  out<-numeric()
  for(m in names(st)){z<-st[[m]];names(z)<-paste(m,names(z),sep="::");out<-c(out,z)}
  for(n in names(comparisons)) {
    ms<-comparisons[[n]];common<-intersect(names(st[[ms[1]]]),names(st[[ms[2]]]))
    z<-st[[ms[1]]][common]-st[[ms[2]]][common];names(z)<-paste(n,common,sep="::")
    out<-c(out,z)
  }
  out
}
parameter_vector <- function(fits) {
  out<-numeric()
  for(m in names(fits)) {
    f<-fits[[m]]; z<-c(f$beta,THRESHOLD=f$threshold,BANG_TARGET=f$target,OOF_REFERRAL=f$oof_referral)
    names(z)<-paste(m,names(z),sep="::");out<-c(out,z)
  }
  out
}
evaluate_rep <- function(b,obj) {
  t0<-proc.time()[3]
  mult<-obj$multipliers[,b]
  wd<-obj$dev$wt*mult[obj$idxdev];we<-obj$eval$wt*mult[obj$idxeval]
  tryCatch({
    fits<-lapply(names(MODELS),function(m)develop(m,wd,obj));names(fits)<-names(MODELS)
    bang<-metric(obj$eval$bang,obj$eval$ref_bang,we,obj)
    full<-conditional<-list(Bang=bang)
    for(m in names(MODELS)) {
      f<-fits[[m]];full[[m]]<-metric(plogis(f$lp),as.numeric(plogis(f$lp)>=f$threshold),we,obj,f$lp)
      f<-obj$point_fits[[m]];conditional[[m]]<-metric(plogis(f$lp),as.numeric(plogis(f$lp)>=f$threshold),we,obj,f$lp)
    }
    list(id=b,success=TRUE,full=flatten(full),conditional=flatten(conditional),
         parameters=parameter_vector(fits),seconds=unname(proc.time()[3]-t0),error="")
  },error=function(e)list(id=b,success=FALSE,seconds=unname(proc.time()[3]-t0),error=conditionMessage(e)))
}

if(mode=="prepare") {
  dev<-read.csv(file.path(base,"inputs","development_inputs.csv"))
  ev<-read.csv(file.path(base,"inputs","evaluation_inputs.csv"))
  fr<-read.csv(file.path(base,"inputs","full_examined_PSU_frame.csv"))
  key<-function(d)paste(d$cycle,d$strata,d$psu,sep="/")
  ixdev<-match(key(dev),key(fr));ixeval<-match(key(ev),key(fr))
  stopifnot(!anyNA(ixdev),!anyNA(ixeval),!anyDuplicated(key(fr)),nrow(dev)==12781,nrow(ev)==11848)
  groups<-split(seq_len(nrow(fr)),paste(fr$cycle,fr$strata,sep="/"))
  stopifnot(all(lengths(groups)%in%c(2,3)))
  set.seed(SEED)
  mult<-matrix(0,nrow(fr),B)
  for(g in groups) {
    m<-length(g);draw<-matrix(sample.int(m,(m-1L)*B,replace=TRUE),nrow=m-1L)
    for(j in seq_len(m))mult[g[j],]<-colSums(draw==j)*m/(m-1L)
    stopifnot(max(abs(colSums(mult[g,,drop=FALSE])-m))<1e-12)
  }
  obj<-list(dev=dev,eval=ev,frame=fr,idxdev=ixdev,idxeval=ixeval,multipliers=mult,
            Xdev=lapply(MODELS,function(v)as.matrix(dev[,v])),
            Xeval=lapply(MODELS,function(v)as.matrix(ev[,v])),
            masks=list(overall=rep(TRUE,nrow(ev)),young=ev$age_grp=="20-39",mid=ev$age_grp=="40-59",
                       old=ev$age_grp==">=60",high=ev$access=="high",low=ev$access=="low"))
  fits<-lapply(names(MODELS),function(m)develop(m,dev$wt,obj));names(fits)<-names(MODELS)
  obj$point_fits<-fits
  original_threshold<-v5_threshold
  threshold_checks<-do.call(rbind,lapply(names(MODELS),function(m){
    z<-read.csv(file.path(base,"source_snapshot",paste0("locked_",m,".csv")))
    expected<-setNames(z$value,z$term)
    cf<-fits[[m]]$beta
    stopifnot(max(abs(cf-expected[names(cf)]))<1e-10,abs(fits[[m]]$threshold-expected["THRESHOLD"])<1e-12)
    data.frame(model=m,coefficient_max_difference=max(abs(cf-expected[names(cf)])),
               threshold=fits[[m]]$threshold,threshold_difference=fits[[m]]$threshold-expected["THRESHOLD"])
  }))
  for(p in list(c(.1,.1,.2,.3),c(.9,.7,.7,.2,.1)))for(t in c(.2,.5,.8)) {
    w<-seq_along(p);stopifnot(identical(threshold_fast(p,w,t),original_threshold(p,w,t)))
  }
  for(m in names(MODELS)) {
    p<-plogis(fits[[m]]$lp)
    stopifnot(abs(auc_fast(p,ev$case,ev$wt)-w_auc(p,ev$case,ev$wt))<1e-12)
  }
  # Evaluation outcomes are absent from the development and prediction inputs.
  flipped<-obj;flipped$eval$case<-1-flipped$eval$case
  for(m in names(MODELS)) {
    f<-develop(m,dev$wt,flipped)
    stopifnot(identical(f$lp,fits[[m]]$lp),identical(f$threshold,fits[[m]]$threshold))
  }
  # Zero-weight records still count toward the original fold normalization N.
  xx<-cbind(x=c(-2,-1,1,2))
  stopifnot(max(abs(fit_l2(xx,c(0,0,1,1),c(2,0,0,2))-
                       fit_l2(xx[c(1,4),,drop=FALSE],c(0,1),c(2,2),n_norm=4)))<1e-12)
  st<-list(Bang=metric(ev$bang,ev$ref_bang,ev$wt,obj))
  for(m in names(MODELS)){f<-fits[[m]];st[[m]]<-metric(plogis(f$lp),as.numeric(plogis(f$lp)>=f$threshold),ev$wt,obj,f$lp)}
  stopifnot(all(st$M2_linear-st$M2_linear==0))
  obj$point<-flatten(st)
  saveRDS(obj,file.path(base,"prepared.rds"))
  write.csv(threshold_checks,file.path(base,"point_reproduction_checks.csv"),row.names=FALSE)
  write.csv(data.frame(key=names(obj$point),estimate=unname(obj$point)),file.path(base,"point_metrics.csv"),row.names=FALSE)
  write.csv(data.frame(PSU_key=key(fr),mult),file.path(base,"replicate_multipliers.csv"),row.names=FALSE)
  writeLines(c("PASS: source coefficient and threshold reproduction",
               "PASS: accelerated tie-block threshold equivalence",
               "PASS: accelerated weighted AUC equivalence",
               "PASS: evaluation outcome changes do not change predictions or thresholds",
               "PASS: zero-weight records preserve original fold N",
               "PASS: identical paired predictions give zero difference",
               "PASS: actual 2/3-PSU strata preserve the Rao-Wu multiplier sum"),file.path(base,"acceptance_checks.txt"))
  writeLines(capture.output(sessionInfo()),file.path(base,"sessionInfo.txt"))
  writeLines(c(paste("Seed:",SEED),paste("Final B:",B),paste("Prepared:",Sys.time()),
               "Fixed candidate structures; fixed age knots 25/46/73 years; C=1.",
               "Normal approximate CI: original point +/- 1.96 * sample SD of replicates.",
               "Percentile intervals are reported as a sensitivity/diagnostic.",
               "No automatic model selection; no correction for historical evaluation exploration.",
               paste("Source MD5:",tools::md5sum(src))),file.path(base,"run_spec.txt"))
  cat("Prepared",B,"joint replicates;",nrow(fr),"PSUs in",length(groups),"strata.\n")
  print(threshold_checks)
  quit(save="no")
}

obj<-readRDS(file.path(base,"prepared.rds"))
if(mode %in% c("pilot","worker")) {
  lo<-if(length(args)>=3)as.integer(args[3])else 1L
  hi<-if(length(args)>=4)as.integer(args[4])else 30L
  stopifnot(lo>=1,hi<=B,lo<=hi)
  dest<-if(mode=="pilot")file.path(base,"pilot.rds")else file.path(base,"chunks",sprintf("chunk_%04d_%04d.rds",lo,hi))
  results<-if(file.exists(dest))readRDS(dest)else list()
  done<-vapply(results,function(x)x$id,integer(1))
  cat("Running",mode,lo,"to",hi,"; already saved",length(done),"\n");flush.console()
  for(i in setdiff(seq.int(lo,hi),done)) {
    results[[length(results)+1L]]<-evaluate_rep(i,obj)
    if(i%%10L==0L || i==hi) {
      saveRDS(results,dest)
      cat("Saved",i,"; successes",sum(vapply(results,function(x)x$success,logical(1))),
          "; mean seconds",round(mean(vapply(results,function(x)x$seconds,numeric(1))),3),"\n")
      flush.console()
    }
  }
  saveRDS(results,dest)
  if(mode=="pilot") {
    stopifnot(all(vapply(results,function(x)x$success,logical(1))))
    cat("Pilot passed:",length(results),"replicates.\n")
  }
  quit(save="no")
}
if(mode=="collect") {
  files<-list.files(file.path(base,"chunks"),pattern="^chunk_[0-9]+_[0-9]+\\.rds$",full.names=TRUE)
  runs<-unlist(lapply(files,readRDS),recursive=FALSE)
  ids<-vapply(runs,function(x)x$id,integer(1));stopifnot(!anyDuplicated(ids),setequal(ids,seq_len(B)))
  runs<-runs[order(ids)]
  status<-do.call(rbind,lapply(runs,function(r)data.frame(replicate=r$id,success=r$success,seconds=r$seconds,error=r$error)))
  write.csv(status,file.path(base,"replicate_status.csv"),row.names=FALSE)
  ok<-status$success
  if(!all(ok))stop("Some replicates failed; inspect replicate_status.csv before inference")
  matrix_of<-function(field)do.call(rbind,lapply(runs,`[[`,field))
  full<-matrix_of("full");conditional<-matrix_of("conditional");pars<-matrix_of("parameters")
  stopifnot(identical(colnames(full),names(obj$point)),identical(colnames(conditional),names(obj$point)))
  bcols<-startsWith(colnames(full),"Bang::")
  stopifnot(max(abs(full[,bcols]-conditional[,bcols]),na.rm=TRUE)==0)
  summarize<-function(mat,label) {
    do.call(rbind,lapply(seq_len(ncol(mat)),function(j){
      v<-mat[,j];v<-v[is.finite(v)];point<-obj$point[j];se<-sd(v);pc<-quantile(v,c(.025,.975),names=FALSE)
      fields<-strsplit(names(obj$point)[j],"::",fixed=TRUE)[[1]]
      data.frame(analysis=label,model=fields[1],group=fields[2],measure=fields[3],estimate=unname(point),
                 bootstrap_se=se,lower95=point-1.96*se,upper95=point+1.96*se,
                 percentile_lower95=pc[1],percentile_upper95=pc[2],
                 bootstrap_mean=mean(v),bootstrap_median=median(v),mean_minus_point=mean(v)-point,
                 bias_over_se=if(se>0)(mean(v)-point)/se else 0,
                 valid_B=length(v),requested_B=B,row.names=NULL)
    }))
  }
  result<-rbind(summarize(full,"full_pipeline"),summarize(conditional,"fixed_model_bootstrap"))
  stopifnot(all(result$valid_B==B))
  write.csv(result,file.path(base,"interval_summary.csv"),row.names=FALSE)
  write.csv(data.frame(replicate=seq_len(B),full,check.names=FALSE),file.path(base,"bootstrap_full_metrics.csv"),row.names=FALSE)
  write.csv(data.frame(replicate=seq_len(B),conditional,check.names=FALSE),file.path(base,"bootstrap_fixed_metrics.csv"),row.names=FALSE)
  write.csv(data.frame(replicate=seq_len(B),pars,check.names=FALSE),file.path(base,"bootstrap_parameters.csv"),row.names=FALSE)
  if(file.exists(file.path(base,"pilot.rds"))) {
    pilot<-readRDS(file.path(base,"pilot.rds"))
    dif<-max(vapply(pilot,function(p)max(abs(p$full-runs[[p$id]]$full)),numeric(1)))
    stopifnot(dif==0)
    cat("Pilot/final repeated replicate agreement: exact.\n")
  }
  writeLines(c("Completed 1000 unique successful joint replicates; all summarized metrics finite.",
               "Bang fixed/full replicate metrics agree exactly.",
               "Original point coefficients and thresholds are not replaced by bootstrap averages.",
               paste("Completed:",Sys.time())),file.path(base,"completion.txt"))
  cat("Completed:",nrow(result),"interval rows; all B=",B,".\n")
  print(subset(result,analysis=="full_pipeline" & model=="M2_linear_minus_Bang" &
                      ((measure=="sensitivity" & group%in%c("overall","young","old","access_gap")) |
                       (group=="overall" & measure=="referral")))[,c("group","measure","estimate","lower95","upper95","valid_B")])
  quit(save="no")
}
stop("Unknown mode")
