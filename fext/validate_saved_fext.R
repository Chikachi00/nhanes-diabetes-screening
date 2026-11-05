# Validate saved point estimates and write matched-population differences.
args<-commandArgs(trailingOnly=TRUE)
ws<-normalizePath(if(length(args))args[1] else getwd(),winslash="/")
out<-if(length(args)>1)args[2] else file.path(ws,"outputs/diabetes_manuscript_revision_20261104_v2/fext")
ref<-if(length(args)>2)args[3] else file.path(ws,"outputs/diabetes_manuscript_revision_20261104/supplement/tables")
a<-read.csv(file.path(out,"original_thresholds_reproduced.csv"));b<-read.csv(file.path(ref,"S12_information_thresholds_by_model.csv"))
z<-merge(a,b,by="model",suffixes=c("_reproduced","_submitted"))
z$single_difference<-z$cut_single_reproduced-z$cut_single_submitted
z$older_difference<-z$cut_60plus_S1_reproduced-z$cut_60plus_S1_submitted
write.csv(z,file.path(out,"original_threshold_reproduction_checks.csv"),row.names=FALSE)
stopifnot(max(abs(z$single_difference))<1e-10,max(abs(z$older_difference))<1e-10)
m<-read.csv(file.path(out,"fixed_structure_sensitivity_metrics.csv"))
r<-read.csv(file.path(ref,"S12_information_results_by_model.csv"))
rows<-list()
for(i in which(m$handling=="unknown_as_no_original"&m$evaluation_population=="all")) for(v in c("referral","sensitivity","specificity")){
 gr<-if(m$group[i]=="60plus")">=60" else m$group[i]
 if(gr=="under60")next
 j<-which(r$period=="eval"&r$group==gr&r$strategy==paste0("F_ext_",m$rule[i])&r$measure==v)
 stopifnot(length(j)==1)
 rows[[length(rows)+1]]<-data.frame(group=gr,rule=m$rule[i],measure=v,reproduced=m[[v]][i],submitted=r$estimate[j],difference=m[[v]][i]-r$estimate[j])
}
ck<-do.call(rbind,rows);stopifnot(max(abs(ck$difference))<1e-10)
write.csv(ck,file.path(out,"original_metric_reproduction_checks.csv"),row.names=FALSE)
deltas<-list()
for(h in c("known_only_fit","unknown_as_yes")){
 pop<-if(h=="known_only_fit")"known_DIQ160" else "all"
 baseline<-m[m$handling=="unknown_as_no_original"&m$evaluation_population==pop,]
 alt<-m[m$handling==h&m$evaluation_population==pop,]
 dz<-merge(baseline,alt,by=c("evaluation_population","group","rule","N","cases"),suffixes=c("_original","_alternative"))
 for(v in c("AUC","referral","sensitivity","specificity","PPV")){
  deltas[[length(deltas)+1]]<-data.frame(alternative=h,dz[,c("evaluation_population","group","rule","N","cases")],measure=v,original=dz[[paste0(v,"_original")]],alternative_estimate=dz[[paste0(v,"_alternative")]],difference=dz[[paste0(v,"_alternative")]]-dz[[paste0(v,"_original")]])
 }
}
write.csv(do.call(rbind,deltas),file.path(out,"matched_population_differences.csv"),row.names=FALSE)
raw<-if(length(args)>3)args[4] else file.path(ws,"outputs/diabetes_manuscript_revision_20261104_v2/analysis/raw")
files<-list.files(raw,pattern="_(E|F|G|H|I|J)\\.xpt$",full.names=TRUE)
write.csv(data.frame(file=basename(files),md5=unname(tools::md5sum(files))),file.path(out,"raw_input_md5.csv"),row.names=FALSE)
cat("Saved results validated. Maximum original point-metric difference:",max(abs(ck$difference)),"\n")
