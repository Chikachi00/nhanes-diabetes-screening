# =====================================================================
# 锁定方程的一次性评价：NHANES 2013–2018（V5 探索性时间分段评价期）
# ---------------------------------------------------------------------
# 只读取 refine_v5_equation.R 生成的 locked_main.csv，系数和阈值一律不改。
# 内容：区分度、转检负担、敏感度（总体、年龄组、可及性组）及与 Bang 的配对差，
#       校准（O:E、截距、斜率、校准图）、截距更新（属于"更新"，不是验证）、
#       每千人表、决策曲线。区间用调查设计 Taylor 线性化（权重、层、PSU）。
#
# 用法：先跑完 refine_v5_equation.R；在第 0 部分填 LOCKED_FILE 的路径；点 Source。
# 重要：这一步只跑一次。看到结果后如果修改方程，就属于"继续开发"。
# =====================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(survey)
})
options(survey.lonely.psu = "adjust")
set.seed(20261028)


# =====================================================================
# 第 0 部分：设置（只改这里）
# =====================================================================
PROJECT_DIR <- "D:/Nottingham/Diabetes_Screening_NHANES"
DATA_DIR    <- file.path(PROJECT_DIR, "data", "raw")
# 改成 refine_v5_equation.R 输出文件夹里的 locked_main.csv，例如：
LOCKED_FILE <- list.files(file.path(PROJECT_DIR, "results"), pattern = "^locked_main\\.csv$",
                          recursive = TRUE, full.names = TRUE)
LOCKED_FILE <- LOCKED_FILE[which.max(file.mtime(LOCKED_FILE))]
cat("使用的锁定文件：", LOCKED_FILE, "\n")
OUT_DIR     <- file.path(PROJECT_DIR, "results",
                         paste0("evaluate_locked_", format(Sys.Date(), "%Y%m%d")))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

DEV_CYCLES  <- character(0)       # 评价时不读开发期
EVAL_CYCLES <- c("H", "I", "J")   # 2013–2018
N_BOOT      <- 1000               # 只用于 AUC 的区间（方程固定，不重拟合）

# 读锁定文件
locked <- read.csv(LOCKED_FILE)
val <- setNames(locked$value, locked$term)
KNOTS <- unname(val[c("KNOT1", "KNOT2", "KNOT3")])
THRESHOLD <- unname(val["THRESHOLD"])
BETA <- val[!names(val) %in% c("THRESHOLD", "KNOT1", "KNOT2", "KNOT3", "C_L2")]
cat("已读取锁定方程：\n"); print(BETA, digits = 10)
cat("锁定阈值：", format(THRESHOLD, digits = 10), "\n")
cat("锁定文件 MD5：", tools::md5sum(LOCKED_FILE), "（写进研究记录）\n")

# 第 1 部分：读数据并整理变量
# =====================================================================
read_nhanes <- function(module, suffix) {
  f <- list.files(DATA_DIR, pattern = paste0("^", module, "_", suffix, "\\.xpt$"),
                  recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  if (length(f) == 0) stop("找不到文件：", module, "_", suffix, ".xpt（检查 DATA_DIR）")
  foreign::read.xport(f[1])
}

read_cycle <- function(suffix) {
  m <- function(x) read_nhanes(x, suffix)
  m("DEMO") %>%
    select(SEQN, RIDSTATR, RIDAGEYR, RIAGENDR, RIDEXPRG, WTMEC2YR, SDMVSTRA, SDMVPSU) %>%
    left_join(m("DIQ") %>% select(SEQN, DIQ010), by = "SEQN") %>%
    left_join(m("BMX") %>% select(SEQN, BMXWAIST, BMXHT, BMXBMI), by = "SEQN") %>%
    left_join(m("BPQ") %>% select(SEQN, BPQ020), by = "SEQN") %>%
    left_join(m("MCQ") %>% select(SEQN, MCQ300C), by = "SEQN") %>%
    left_join(m("PAQ") %>% select(SEQN, PAQ605, PAQ610, PAQ620, PAQ625, PAQ635,
                                  PAQ640, PAQ650, PAQ655, PAQ665, PAQ670), by = "SEQN") %>%
    left_join(m("GHB") %>% select(SEQN, LBXGH), by = "SEQN") %>%
    left_join(m("HIQ") %>% select(SEQN, HIQ011), by = "SEQN") %>%
    left_join(m("HUQ") %>% select(SEQN, HUQ030), by = "SEQN") %>%
    mutate(cycle = suffix)
}

domain_days <- function(yesno, days) {
  case_when(yesno == 2 ~ 0,
            yesno == 1 & days %in% 1:7 ~ as.numeric(days),
            TRUE ~ NA_real_)
}
yn <- function(x) case_when(x == 1 ~ 1, x == 2 ~ 0, TRUE ~ NA_real_)

# 限制性立方样条（与 V5 附录 B 一致）
rcs_term <- function(A, k) {
  pos <- function(u) pmax(u, 0)^3
  (pos(A - k[1]) -
      pos(A - k[2]) * (k[3] - k[1]) / (k[3] - k[2]) +
      pos(A - k[3]) * (k[2] - k[1]) / (k[3] - k[2])) / (k[3] - k[1])^2
}

all_cycles <- c(DEV_CYCLES, EVAL_CYCLES)
raw <- bind_rows(lapply(all_cycles, read_cycle))

full <- raw %>%
  filter(RIDSTATR == 2) %>%
  mutate(
    wt     = WTMEC2YR / 3,            # 每段 3 个周期，权重除以 3（与 V5 一致）
    strata = SDMVSTRA, psu = SDMVPSU,
    age    = RIDAGEYR,
    male   = as.numeric(RIAGENDR == 1),
    waist  = BMXWAIST, height = BMXHT, bmi = BMXBMI,
    whtr   = BMXWAIST / BMXHT,
    famhx  = yn(MCQ300C),
    htn    = yn(BPQ020),
    d1 = domain_days(PAQ605, PAQ610), d2 = domain_days(PAQ620, PAQ625),
    d3 = domain_days(PAQ635, PAQ640), d4 = domain_days(PAQ650, PAQ655),
    d5 = domain_days(PAQ665, PAQ670),
    any_active = (d1 >= 5) | (d2 >= 5) | (d3 >= 5) | (d4 >= 5) | (d5 >= 5),
    inactive = case_when(any_active %in% TRUE ~ 0,
                         rowSums(is.na(cbind(d1, d2, d3, d4, d5))) == 0 ~ 1,
                         TRUE ~ NA_real_),
    case   = ifelse(is.na(LBXGH), NA, as.numeric(LBXGH >= 6.5)),
    ins    = yn(HIQ011),
    usual  = case_when(HUQ030 %in% c(1, 3) ~ 1, HUQ030 == 2 ~ 0, TRUE ~ NA_real_),
    access = case_when(ins == 1 & usual == 1 ~ "high", ins == 0 | usual == 0 ~ "low"),
    access = ifelse(is.na(ins) | is.na(usual), NA, access),
    age_grp = cut(age, c(19, 39, 59, Inf), labels = c("20-39", "40-59", ">=60")),
    
    # ---- 方程用到的变换（与 V5 附录 B 一致）----
    A  = age / 10,
    a  = A - 5,                          # 以 50 岁为中心
    s  = rcs_term(A, KNOTS) - rcs_term(5, KNOTS),
    w  = 10 * (whtr - 0.55),             # 以 WHtR 0.55 为中心
    aw = a * w,                          # 年龄 × WHtR 交互
    
    # ---- Bang 积分 ----
    age_pts   = case_when(age < 40 ~ 0, age < 50 ~ 1, age < 60 ~ 2, TRUE ~ 3),
    bmi_pts   = case_when(bmi < 25 ~ 0, bmi < 30 ~ 1, bmi < 40 ~ 2, TRUE ~ 3),
    w_in      = waist / 2.54,
    waist_pts = ifelse(male == 1,
                       case_when(w_in < 37 ~ 0, w_in < 40 ~ 1, w_in < 50 ~ 2, TRUE ~ 3),
                       case_when(w_in < 31.5 ~ 0, w_in < 35 ~ 1, w_in < 49 ~ 2, TRUE ~ 3)),
    obese_pts = pmax(bmi_pts, waist_pts),
    active    = 1 - inactive,
    bang      = age_pts + male + famhx + htn + obese_pts - active,
    ref_bang  = as.numeric(bang >= 5)
  )

need <- c("age", "male", "waist", "height", "bmi", "whtr", "famhx", "htn",
          "inactive", "case", "access")
full <- full %>%
  mutate(in_sample = age >= 20 & DIQ010 %in% 2 & !(RIDEXPRG %in% 1) &
           complete.cases(across(all_of(need))))
full$in_sample[is.na(full$in_sample)] <- FALSE
full$dev  <- full$in_sample & full$cycle %in% DEV_CYCLES
full$eval <- full$in_sample & full$cycle %in% EVAL_CYCLES

cat("\n开发期样本：", sum(full$dev), "人，病例", sum(full$case[full$dev]), "人\n")
if (any(full$eval)) cat("评价期样本：", sum(full$eval), "人，病例", sum(full$case[full$eval]), "人\n")
print(table(cycle = full$cycle[full$dev], case = full$case[full$dev]))




# =====================================================================
# 第 2 部分：原样套用锁定方程
# =====================================================================
lp_of <- function(beta, d) as.vector(beta[1] + as.matrix(d[, names(beta)[-1]]) %*% beta[-1])
full$lp  <- lp_of(BETA, full)
full$p   <- plogis(full$lp)
full$ref_new <- as.numeric(full$p >= THRESHOLD)
full$both <- full$ref_new - full$ref_bang        # 配对差：新方程 − Bang

des   <- svydesign(ids = ~psu, strata = ~strata, weights = ~wt, nest = TRUE, data = full)
des_e <- subset(des, eval)
des_c <- subset(des_e, case == 1)                # 病例（敏感度）
des_n <- subset(des_e, case == 0)                # 非病例（特异度）
ev    <- full[full$eval, ]

ci_row <- function(est, label, group = "overall") {
  ci <- confint(est)
  data.frame(measure = label, group = group, estimate = as.numeric(coef(est)),
             lower95 = ci[, 1], upper95 = ci[, 2], row.names = NULL)
}


# =====================================================================
# 第 3 部分：转检负担、敏感度、特异度、PPV（Taylor 区间）
# =====================================================================
tab <- bind_rows(
  ci_row(svymean(~ref_new,  des_e), "refer_new"),
  ci_row(svymean(~ref_bang, des_e), "refer_bang"),
  ci_row(svymean(~both,     des_e), "refer_diff"),
  ci_row(svymean(~ref_new,  des_c), "sens_new"),
  ci_row(svymean(~ref_bang, des_c), "sens_bang"),
  ci_row(svymean(~both,     des_c), "sens_diff"),
  ci_row(svymean(~I(1 - ref_new),  des_n), "spec_new"),
  ci_row(svymean(~I(1 - ref_bang), des_n), "spec_bang"),
  ci_row(svymean(~case, subset(des_e, ref_new == 1)),  "ppv_new"),
  ci_row(svymean(~case, subset(des_e, ref_bang == 1)), "ppv_bang")
)

by_grp <- function(var) {
  bind_rows(lapply(c("ref_new", "ref_bang", "both"), function(r) {
    b <- svyby(as.formula(paste0("~", r)), as.formula(paste0("~", var)), des_c, svymean)
    ci <- confint(b)
    data.frame(measure = c(ref_new = "sens_new", ref_bang = "sens_bang", both = "sens_diff")[r],
               group = paste(var, b[[var]]), estimate = b[[r]],
               lower95 = ci[, 1], upper95 = ci[, 2], row.names = NULL)
  }))
}
tab <- bind_rows(tab, by_grp("age_grp"), by_grp("access"))

# 可及性差距（较高 − 低），Taylor 区间
gap <- function(r) {
  b <- svyby(as.formula(paste0("~", r)), ~access, des_c, svymean, covmat = TRUE)
  g <- svycontrast(b, c(high = 1, low = -1))
  data.frame(measure = paste0("gap_access_", sub("ref_", "", r)), group = "high - low",
             estimate = as.numeric(coef(g)), lower95 = confint(g)[1], upper95 = confint(g)[2])
}
tab <- bind_rows(tab, gap("ref_new"), gap("ref_bang"))


# =====================================================================
# 第 4 部分：区分度（加权 AUC，bootstrap 标准误）
# =====================================================================
w_auc <- function(score, y, w) {
  k <- w > 0; score <- score[k]; y <- y[k]; w <- w[k]
  o <- order(score); score <- score[o]; y <- y[o]; w <- w[o]
  w0 <- ifelse(y == 0, w, 0); w1 <- ifelse(y == 1, w, 0)
  grp <- match(score, unique(score))
  W0g <- tapply(w0, grp, sum); W1g <- tapply(w1, grp, sum)
  below <- cumsum(W0g) - W0g
  sum(W1g * (below + 0.5 * W0g)) / (sum(w1) * sum(w0))
}
rep_e <- subset(as.svrepdesign(des, type = "subbootstrap", replicates = N_BOOT), eval)
Wr <- weights(rep_e, type = "analysis")
auc_fun <- function(w) c(new = w_auc(ev$p, ev$case, w), bang = w_auc(ev$bang, ev$case, w))
auc_pt <- auc_fun(ev$wt)
auc_bt <- apply(Wr, 2, auc_fun)
auc_se <- c(apply(auc_bt, 1, sd), diff = sd(auc_bt[1, ] - auc_bt[2, ]))
auc_est <- c(auc_pt, diff = unname(auc_pt[1] - auc_pt[2]))
tab <- bind_rows(tab, data.frame(measure = c("AUC_new", "AUC_bang", "AUC_diff"), group = "overall",
                                 estimate = auc_est, lower95 = auc_est - 1.96 * auc_se,
                                 upper95 = auc_est + 1.96 * auc_se, row.names = NULL))


# =====================================================================
# 第 5 部分：校准
# ---------------------------------------------------------------------
# O:E < 1 表示预测偏高；截距 < 0 同义；斜率 < 1 表示预测太极端
# 截距更新：只改截距，让平均预测等于实际比例（这是"更新"，要单独报告）
# =====================================================================
oe  <- svyratio(~case, ~p, des_e)
fi  <- svyglm(case ~ offset(lp), design = des_e, family = quasibinomial())
fs  <- svyglm(case ~ lp,         design = des_e, family = quasibinomial())
calib <- data.frame(
  measure = c("O:E", "calibration intercept", "calibration slope"),
  estimate = c(coef(oe), coef(fi)[1], coef(fs)[2]),
  lower95 = c(confint(oe)[1], confint(fi)[1, 1], confint(fs)[2, 1]),
  upper95 = c(confint(oe)[2], confint(fi)[1, 2], confint(fs)[2, 2]))
updated_intercept <- unname(BETA[1] + coef(fi)[1])
cat("\n截距更新后的截距（仅作更新报告）：", format(updated_intercept, digits = 10), "\n")

ev$dec <- cut(ev$p, quantile(ev$p, 0:10 / 10), include.lowest = TRUE, labels = FALSE)
full$dec <- NA; full$dec[full$eval] <- ev$dec
dz <- subset(svydesign(ids = ~psu, strata = ~strata, weights = ~wt, nest = TRUE, data = full), eval)
ob <- svyby(~case, ~dec, dz, svymean); pr <- svyby(~p, ~dec, dz, svymean)
cal_tab <- data.frame(decile = ob$dec, predicted = pr$p, observed = ob$case, se = ob$se)
png(file.path(OUT_DIR, "calibration_eval.png"), width = 1600, height = 1600, res = 250)
lim <- c(0, max(ob$case + 1.96 * ob$se, pr$p) * 1.05)
plot(pr$p, ob$case, xlim = lim, ylim = lim, pch = 19,
     xlab = "Mean predicted risk", ylab = "Observed (HbA1c >= 6.5%)", main = "Calibration, 2013-2018")
arrows(pr$p, pmax(ob$case - 1.96 * ob$se, 0), pr$p, ob$case + 1.96 * ob$se, angle = 90, code = 3, length = 0.03)
abline(0, 1, lty = 2, col = "grey40")
invisible(dev.off())


# =====================================================================
# 第 6 部分：次要分析——在评价期严格等转检比例时的敏感度
# =====================================================================
wmean <- function(x, w) sum(w * x) / sum(w)
o <- order(ev$p, decreasing = TRUE)
share <- cumsum(ev$wt[o]) / sum(ev$wt)
thr_eq <- ev$p[o][which(share >= wmean(ev$ref_bang, ev$wt))[1]]
full$ref_eq <- as.numeric(full$p >= thr_eq); full$both_eq <- full$ref_eq - full$ref_bang
des_c2 <- subset(subset(svydesign(ids = ~psu, strata = ~strata, weights = ~wt, nest = TRUE, data = full), eval), case == 1)
tab <- bind_rows(tab,
                 ci_row(svymean(~ref_eq,  des_c2), "sens_new_equal_referral"),
                 ci_row(svymean(~both_eq, des_c2), "sens_diff_equal_referral"))


# =====================================================================
# 第 7 部分：每千人表（按年龄组）
# =====================================================================
per1000 <- ev %>% group_by(age_group = age_grp) %>%
  summarise(n = n(), cases = sum(case),
            cases_per1000 = 1000 * wmean(case, wt),
            found_bang  = 1000 * wmean(case * ref_bang, wt),
            found_new   = 1000 * wmean(case * ref_new, wt),
            extra_found = found_new - found_bang,
            refer_bang  = 1000 * wmean(ref_bang, wt),
            refer_new   = 1000 * wmean(ref_new, wt), .groups = "drop")


# =====================================================================
# 第 8 部分：决策曲线（加权净获益）
# =====================================================================
ts <- seq(0.01, 0.15, by = 0.005)
nb <- function(ref, t) wmean(ref * ev$case, ev$wt) - wmean(ref * (1 - ev$case), ev$wt) * t / (1 - t)
dca <- data.frame(threshold = ts,
                  new_eq = sapply(ts, function(t) nb(as.numeric(ev$p >= t), t)),
                  bang = sapply(ts, function(t) nb(ev$ref_bang, t)),
                  refer_all = sapply(ts, function(t) nb(rep(1, nrow(ev)), t)), refer_none = 0)
png(file.path(OUT_DIR, "decision_curve_eval.png"), width = 2000, height = 1500, res = 250)
matplot(dca$threshold, dca[, 2:5], type = "l", lty = c(1, 1, 2, 3), lwd = 2,
        col = c("#1f77b4", "#d62728", "grey50", "black"), ylim = c(-0.01, max(dca[, 2:4]) * 1.1),
        xlab = "Threshold probability", ylab = "Net benefit", main = "Decision curve, 2013-2018")
legend("topright", c("Locked equation", "Bang >=5", "Refer all", "Refer none"),
       col = c("#1f77b4", "#d62728", "grey50", "black"), lty = c(1, 1, 2, 3), lwd = 2, bty = "n")
invisible(dev.off())


# =====================================================================
# 保存
# =====================================================================
write.csv(tab,     file.path(OUT_DIR, "eval_main_results.csv"), row.names = FALSE)
write.csv(calib,   file.path(OUT_DIR, "eval_calibration.csv"), row.names = FALSE)
write.csv(cal_tab, file.path(OUT_DIR, "eval_calibration_deciles.csv"), row.names = FALSE)
write.csv(per1000, file.path(OUT_DIR, "eval_per1000.csv"), row.names = FALSE)
write.csv(dca,     file.path(OUT_DIR, "eval_decision_curve.csv"), row.names = FALSE)
writeLines(c(paste("锁定文件：", LOCKED_FILE), paste("锁定文件 MD5：", tools::md5sum(LOCKED_FILE)),
             paste("截距更新后的截距（更新，不是验证）：", format(updated_intercept, digits = 10)),
             paste("评价期等转检比例阈值（次要分析）：", format(thr_eq, digits = 10)),
             paste("运行时间：", Sys.time())),
           file.path(OUT_DIR, "eval_info.txt"))

cat("\n===== 评价期主要结果 =====\n"); print(tab, digits = 3)
cat("\n===== 校准 =====\n"); print(calib, digits = 3)
cat("\n===== 每千人表 =====\n"); print(as.data.frame(per1000), digits = 3)
cat("\n结果已保存到：", OUT_DIR, "\n")