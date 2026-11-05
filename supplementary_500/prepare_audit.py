from pathlib import Path
import hashlib,json

out=Path(__file__).resolve().parent
workspace=out.parents[2]
src=workspace/'outputs/diabetes_manuscript_revision_20261104/supplement/code'
receipts=[]
for kind in ['joint_cutoff','access_gap']:
    path=src/f'{kind}_v5.R'
    text=path.read_text(encoding='utf-8-sig')
    start=text.index('PROJECT_DIR <-')
    end=text.index('dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)',start)
    header=f'''# Portable execution adapter. Optional positional arguments: raw_directory output_directory.
.audit_file_arg <- grep("^--file=", commandArgs(), value = TRUE)[1]
.audit_script_dir <- dirname(normalizePath(sub("^--file=", "", .audit_file_arg), winslash = "/", mustWork = TRUE))
.audit_args <- commandArgs(trailingOnly = TRUE)
PROJECT_DIR <- .audit_script_dir # optional original locked_main.csv diagnostic only
DATA_DIR <- if (length(.audit_args) >= 1) normalizePath(.audit_args[1], winslash="/", mustWork=TRUE) else file.path(.audit_script_dir, "..", "analysis", "raw")
OUT_DIR <- if (length(.audit_args) >= 2) .audit_args[2] else file.path(.audit_script_dir, "{kind}_rerun")
'''
    text=text[:start]+header+text[end:]
    edits=['Replaced author-specific PROJECT_DIR, DATA_DIR and dated OUT_DIR with portable positional arguments and script-relative defaults.']
    if kind=='access_gap':
        target='  d$wt * f[key]\n}'
        assert text.count(target)==1
        text=text.replace(target,'  .audit_last_f <<- f\n  d$wt * f[key]\n}',1)
        target='for (b in seq_len(B_BOOT)) bt[, , b] <- curves_for(dC, boot_weights(dC))'
        replacement='''audit_status <- vector("list", B_BOOT)
audit_multipliers <- vector("list", B_BOOT)
audit_rng_states <- vector("list", B_BOOT)
for (b in seq_len(B_BOOT)) {
  audit_rng_states[[b]] <- .Random.seed
  .audit_start <- proc.time()[3]
  .audit_weights <- boot_weights(dC)
  bt[, , b] <- curves_for(dC, .audit_weights)
  audit_multipliers[[b]] <- data.frame(replicate=b, stratum_psu=names(.audit_last_f), multiplier=as.numeric(.audit_last_f))
  audit_status[[b]] <- data.frame(replicate=b, status="completed", finite_cells=sum(is.finite(bt[,,b])), missing_cells=sum(is.na(bt[,,b])), total_cells=length(bt[,,b]), elapsed_seconds=unname(proc.time()[3]-.audit_start))
  if (b %% 50 == 0) cat("Audit completed", b, "of", B_BOOT, "replicates\\n")
}
saveRDS(list(estimate=est, replicate_values=bt, budgets=BUDGETS, pre_draw_rng_states=audit_rng_states, seed=20261103, B=B_BOOT, conditional_on_fitted_models=TRUE), file.path(OUT_DIR,"replicate_archive.rds"))
write.csv(bind_rows(audit_multipliers), file.path(OUT_DIR,"replicate_psu_multipliers.csv"),row.names=FALSE)
write.csv(bind_rows(audit_status), file.path(OUT_DIR,"replicate_status.csv"),row.names=FALSE)
audit_valid_matrix <- apply(is.finite(bt),c(1,2),sum)
dimnames(audit_valid_matrix) <- list(as.character(seq_along(BUDGETS)),colnames(est))
audit_valid <- as.data.frame(as.table(audit_valid_matrix), stringsAsFactors=FALSE)
names(audit_valid) <- c("budget_index","curve","valid_replicates")
audit_valid$budget <- BUDGETS[as.integer(audit_valid$budget_index)]
write.csv(audit_valid,file.path(OUT_DIR,"replicate_valid_counts.csv"),row.names=FALSE)
audit_bt <- bt
dimnames(audit_bt) <- list(as.character(seq_along(BUDGETS)),colnames(est),as.character(seq_len(B_BOOT)))
audit_long <- as.data.frame(as.table(audit_bt),stringsAsFactors=FALSE)
names(audit_long) <- c("budget_index","curve","replicate","value")
audit_long$budget <- BUDGETS[as.integer(audit_long$budget_index)]
audit_long$replicate <- as.integer(audit_long$replicate)
write.csv(audit_long,file.path(OUT_DIR,"replicate_values.csv"),row.names=FALSE)
'''
        assert text.count(target)==1
        text=text.replace(target,replacement,1)
        edits.append('Instrumented original 500-replicate loop with pre-draw RNG state, PSU multipliers, replicate values, finite/missing counts and timing; no random draws added, no inferential calculation changed.')
    else:
        edits.append('Original script uses fixed-threshold Taylor intervals and has no bootstrap. No new bootstrap was added or represented as original analysis.')
    newpath=out/f'{kind}_v5_audit.R'
    newpath.write_text(text,encoding='utf-8')
    receipts.append({'source':str(path),'source_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'adapted':str(newpath),'adapted_sha256':hashlib.sha256(newpath.read_bytes()).hexdigest(),'changes':edits})
(out/'adaptation_manifest.json').write_text(json.dumps(receipts,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps({'scripts_created':len(receipts),'source_originals_modified':False}))
