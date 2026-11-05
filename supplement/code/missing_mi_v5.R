# =====================================================================
# 缺失数据描述 + 多重插补敏感性分析（V5 锁定方程，2013–2018）
# ---------------------------------------------------------------------
# 主分析仍是完整病例（与 V5 一致）。本脚本做两件事：
#   1. 缺失描述：各变量加权缺失比例；纳入者与被排除者的特征比较
#   2. 多重插补（mice，m = 20）：插补预测因素和可及性变量；HbA1c 缺失者只帮助插补、
#      不进入评价（von Hippel 的“插补后删除”）。每套数据原样套用锁定方程
#      （若找到 seg_thresholds.csv，也套用 S1 老年组单独切点），
#      Taylor 方差 + Rubin 规则合并，与完整病例结果并排。
# 方程和切点一律不改，不重新拟合。
#
# 第一次使用先安装：install.packages("mice")
# 用法：放到 scripts 文件夹，只改第 0 部分，点 Source。m = 20 约需 3–10 分钟。
# 输出：results\missing_mi_日期\
# =====================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(survey)
  library(mice)
})
options(survey.lonely.psu = "adjust")
set.seed(20261031)


# =====================================================================
# 第 0 部分：设置（只改这里）
# =====================================================================
PROJECT_DIR <- "D:/Nottingham/Diabetes_Screening_NHANES"
DATA_DIR    <- file.path(PROJECT_DIR, "data", "raw")
OUT_DIR     <- file.path(PROJECT_DIR, "results", paste0("missing_mi_", format(Sys.Date(), "%Y%m%d")))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

DEV_CYCLES  <- c("E", "F", "G")
EVAL_CYCLES <- c("H", "I", "J")
M_IMP  <- 20     # 插补套数；想先试跑可改成 5
MAXIT  <- 10     # 每套迭代次数

newest <- function(pattern) {
  f <- list.files(file.path(PROJECT_DIR, "results"), pattern = pattern, recursive = TRUE, full.names = TRUE)
  if (length(f) == 0) NULL else f[which.max(file.mtime(f))]
}
LOCKED_FILE <- newest("^locked_main\\.csv$")
if (is.null(LOCKED_FILE)) stop("results 里找不到 locked_main.csv，请先运行 refine_v5_equation.R")
locked <- read.csv(LOCKED_FILE); val <- setNames(locked$value, locked$term)
THRESHOLD <- unname(val["THRESHOLD"])
BETA <- val[!names(val) %in% c("THRESHOLD", "KNOT1", "KNOT2", "KNOT3", "C_L2")]
cat("锁定文件：", LOCKED_FILE, "\nMD5：", tools::md5sum(LOCKED_FILE), "\n")
if (any(!names(BETA)[-1] %in% c("a", "w", "male", "famhx", "htn", "inactive")))
  stop("本脚本只支持线性年龄主方程（M2_linear）")

SEG_FILE <- newest("^seg_thresholds\\.csv$")   # 可选：S1 老年组单独切点
THR_OLD <- NA
if (!is.null(SEG_FILE)) {
  st <- read.csv(SEG_FILE)
  THR_OLD <- st$cut_60plus[st$strategy == "S1_single_agecut"]
  cat("找到 S1 切点文件：", SEG_FILE, "；≥60 岁切点 =", signif(THR_OLD, 6), "\n")
} else cat("没有找到 seg_thresholds.csv，只评价 S0（单一切点）\n")


# =====================================================================
# 第 1 部分：读数据，整理“原始变量”（缺失保持为 NA）
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

raw <- bind_rows(lapply(c(DEV_CYCLES, EVAL_CYCLES), read_cycle))

base <- raw %>%
  filter(RIDSTATR == 2) %>%
  mutate(
    wt = WTMEC2YR / 3, strata = SDMVSTRA, psu = SDMVPSU,
    period = ifelse(cycle %in% DEV_CYCLES, "dev", "eval"),
    age = RIDAGEYR, male = as.numeric(RIAGENDR == 1),
    waist = BMXWAIST, height = BMXHT, bmi = BMXBMI,
    famhx = yn(MCQ300C), htn = yn(BPQ020),
    d1 = domain_days(PAQ605, PAQ610), d2 = domain_days(PAQ620, PAQ625),
    d3 = domain_days(PAQ635, PAQ640), d4 = domain_days(PAQ650, PAQ655),
    d5 = domain_days(PAQ665, PAQ670),
    any_active = (d1 >= 5) | (d2 >= 5) | (d3 >= 5) | (d4 >= 5) | (d5 >= 5),
    inactive = case_when(any_active %in% TRUE ~ 0,
                         rowSums(is.na(cbind(d1, d2, d3, d4, d5))) == 0 ~ 1,
                         TRUE ~ NA_real_),
    hba1c = LBXGH,
    ins = yn(HIQ011),
    usual = case_when(HUQ030 %in% c(1, 3) ~ 1, HUQ030 == 2 ~ 0, TRUE ~ NA_real_),
    # 目标人群（与主分析相同，但不要求变量完整）
    scope = age >= 20 & DIQ010 %in% 2 & !(RIDEXPRG %in% 1)
  ) %>%
  select(SEQN, cycle, period, wt, strata, psu, scope, age, male, waist, height, bmi,
         famhx, htn, inactive, hba1c, ins, usual)
base$scope[is.na(base$scope)] <- FALSE

# 由原始变量生成方程、Bang、USPSTF 等派生变量（完整病例和每套插补数据都用这一个函数）
derive <- function(d) {
  d %>% mutate(
    whtr = waist / height,
    case = as.numeric(hba1c >= 6.5),
    access = ifelse(is.na(ins) | is.na(usual), NA,
                    ifelse(ins == 1 & usual == 1, "high", "low")),
    age_grp = cut(age, c(19, 39, 59, Inf), labels = c("20-39", "40-59", ">=60")),
    a = age / 10 - 5, w = 10 * (whtr - 0.55),
    age_pts = case_when(age < 40 ~ 0, age < 50 ~ 1, age < 60 ~ 2, TRUE ~ 3),
    bmi_pts = case_when(bmi < 25 ~ 0, bmi < 30 ~ 1, bmi < 40 ~ 2, TRUE ~ 3),
    w_in = waist / 2.54,
    waist_pts = ifelse(male == 1,
                       case_when(w_in < 37 ~ 0, w_in < 40 ~ 1, w_in < 50 ~ 2, TRUE ~ 3),
                       case_when(w_in < 31.5 ~ 0, w_in < 35 ~ 1, w_in < 49 ~ 2, TRUE ~ 3)),
    bang = age_pts + male + famhx + htn + pmax(bmi_pts, waist_pts) - (1 - inactive),
    ref_bang = as.numeric(bang >= 5),
    ref_uspstf = as.numeric(age >= 35 & age <= 70 & bmi >= 25),
    p = plogis(as.vector(BETA[1] + as.matrix(cbind(a, w, male, famhx, htn, inactive)[, names(BETA)[-1]]) %*% BETA[-1])),
    ref_S0 = as.numeric(p >= THRESHOLD),
    ref_S1 = if (is.na(THR_OLD)) NA_real_ else as.numeric(p >= ifelse(age >= 60, THR_OLD, THRESHOLD)),
    dB_S0 = ref_S0 - ref_bang, dB_S1 = ref_S1 - ref_bang
  )
}

VARS_NEED <- c("waist", "height", "bmi", "famhx", "htn", "inactive", "hba1c", "ins", "usual")
cc <- derive(base)
cc$complete <- cc$scope & complete.cases(cc[, VARS_NEED])
cat("\n目标人群：", sum(base$scope), "人；完整病例：", sum(cc$complete), "人\n")
print(table(period = cc$period[cc$scope], complete = cc$complete[cc$scope]))


# =====================================================================
# 第 2 部分：缺失描述
# =====================================================================
des0 <- svydesign(ids = ~psu, strata = ~strata, weights = ~wt, nest = TRUE, data = cc)
miss_tab <- bind_rows(lapply(c("dev", "eval"), function(pe) {
  k <- cc$scope & cc$period == pe
  bind_rows(lapply(c(VARS_NEED, "any_missing"), function(v) {
    x <- if (v == "any_missing") !cc$complete[k] else is.na(cc[[v]][k])
    data.frame(period = pe, variable = v, n_missing = sum(x), n = sum(k),
               pct_unweighted = 100 * mean(x),
               pct_weighted = 100 * sum(cc$wt[k] * x) / sum(cc$wt[k]))
  }))
}))

# 纳入 vs 被排除：加权均值 / 比例，及调查设计检验 p 值
cmp_vars <- c("age", "male", "bmi", "waist", "famhx", "htn", "inactive", "hba1c")
cmp_tab <- bind_rows(lapply(c("dev", "eval"), function(pe) {
  bind_rows(lapply(cmp_vars, function(v) {
    d_s <- des0[cc$scope & cc$period == pe & !is.na(cc[[v]]), ]
    f <- as.formula(paste0("~", v))
    m <- svyby(f, ~complete, d_s, svymean, na.rm = TRUE)
    pv <- tryCatch(svyttest(as.formula(paste0(v, " ~ complete")), d_s)$p.value, error = function(e) NA)
    data.frame(period = pe, variable = v,
               included = m[m$complete == TRUE, v], excluded = m[m$complete == FALSE, v],
               p_value = unname(pv))
  }))
}))


# =====================================================================
# 第 3 部分：多重插补（评价期 2013–2018 目标人群）
# 插补模型：全部预测因素、HbA1c、可及性、年龄、性别、周期、log(权重)
#           （层数过多不放入；连续变量用 pmm，二分类用 logreg）
# =====================================================================
ev0 <- base %>% filter(scope, period == "eval") %>%
  mutate(across(c(famhx, htn, inactive, ins, usual), ~ factor(.x, levels = c(0, 1))),
         cycle = factor(cycle), logwt = log(wt))
imp_vars <- c("age", "male", "cycle", "logwt", "waist", "height", "bmi",
              "famhx", "htn", "inactive", "hba1c", "ins", "usual")
ini  <- mice(ev0[, imp_vars], maxit = 0, printFlag = FALSE)
meth <- ini$method
meth[c("waist", "height", "bmi", "hba1c")] <- ifelse(meth[c("waist", "height", "bmi", "hba1c")] == "", "", "pmm")
meth[c("famhx", "htn", "inactive", "ins", "usual")] <-
  ifelse(meth[c("famhx", "htn", "inactive", "ins", "usual")] == "", "", "logreg")
cat("\n开始插补：m =", M_IMP, "，每套", MAXIT, "次迭代……\n"); flush.console()
imp <- mice(ev0[, imp_vars], m = M_IMP, maxit = MAXIT, method = meth, printFlag = FALSE, seed = 20261031)


# =====================================================================
# 第 4 部分：每套数据的估计（Taylor）+ Rubin 合并；完整病例并排
# =====================================================================
num01 <- function(x) as.numeric(as.character(x))
est_one <- function(d, hba1c_obs) {
  # d：评价期全部受检者（含目标人群外的，用来保持 PSU 结构）
  dd <- derive(d)
  des <- svydesign(ids = ~psu, strata = ~strata, weights = ~wt, nest = TRUE, data = dd)
  k_all <- dd$period == "eval" & dd$scope & hba1c_obs   # HbA1c 实测者才进入评价
  out <- list()
  one <- function(label, var, idx) {
    if (all(is.na(dd[[var]][idx]))) return(NULL)
    m <- svymean(as.formula(paste0("~", var)), des[idx, ], na.rm = TRUE)
    data.frame(measure = label, est = 100 * coef(m)[[1]], se = 100 * SE(m)[[1]])
  }
  cases <- k_all & dd$case %in% 1; nonc <- k_all & dd$case %in% 0
  rows <- list(
    one("referral_Bang", "ref_bang", k_all), one("referral_S0", "ref_S0", k_all),
    one("referral_S1", "ref_S1", k_all), one("referral_USPSTF", "ref_uspstf", k_all),
    one("sens_Bang", "ref_bang", cases), one("sens_S0", "ref_S0", cases),
    one("sens_S1", "ref_S1", cases), one("sens_USPSTF", "ref_uspstf", cases),
    one("sens_diff_S0_vs_Bang", "dB_S0", cases), one("sens_diff_S1_vs_Bang", "dB_S1", cases),
    one("nonref_nonc_S0 (=specificity)", "ref_S0", nonc)
  )
  for (g in levels(dd$age_grp)) {
    kg <- cases & dd$age_grp %in% g
    rows <- c(rows, list(one(paste0("sens_Bang_", g), "ref_bang", kg), one(paste0("sens_S0_", g), "ref_S0", kg),
                         one(paste0("sens_S1_", g), "ref_S1", kg),
                         one(paste0("sens_diff_S0_vs_Bang_", g), "dB_S0", kg),
                         one(paste0("sens_diff_S1_vs_Bang_", g), "dB_S1", kg)))
  }
  r <- bind_rows(rows)
  r$est[r$measure == "nonref_nonc_S0 (=specificity)"] <- 100 - r$est[r$measure == "nonref_nonc_S0 (=specificity)"]
  kk <- k_all & !is.na(dd$p) & !is.na(dd$case)
  auc <- function(score, y, w) {
    o <- order(score); score <- score[o]; y <- y[o]; w <- w[o]
    w0 <- ifelse(y == 0, w, 0); w1 <- ifelse(y == 1, w, 0)
    grp <- match(score, unique(score)); W0g <- tapply(w0, grp, sum); W1g <- tapply(w1, grp, sum)
    sum(W1g * (cumsum(W0g) - W0g + 0.5 * W0g)) / (sum(w1) * sum(w0))
  }
  bind_rows(r, data.frame(measure = c("AUC_equation", "AUC_Bang"),
                          est = c(auc(dd$p[kk], dd$case[kk], dd$wt[kk]), auc(dd$bang[kk], dd$case[kk], dd$wt[kk])),
                          se = NA))
}

ev_all <- base %>% filter(period == "eval")                 # 评价期全部受检者
pos <- match(ev0$SEQN, ev_all$SEQN)
hba1c_obs_all <- !is.na(ev_all$hba1c)

# 完整病例（主分析口径）
cc_ev <- ev_all; cc_ev$scope <- cc_ev$scope & complete.cases(cc_ev[, VARS_NEED])
res_cc <- est_one(cc_ev, hba1c_obs_all) %>% mutate(lower95 = est - 1.96 * se, upper95 = est + 1.96 * se)

# 每套插补
res_m <- lapply(seq_len(M_IMP), function(i) {
  ci <- complete(imp, i)
  d <- ev_all
  for (v in c("waist", "height", "bmi", "hba1c")) d[[v]][pos] <- ci[[v]]
  for (v in c("famhx", "htn", "inactive", "ins", "usual")) d[[v]][pos] <- num01(ci[[v]])
  if (i %% 5 == 0) { cat("  已完成", i, "/", M_IMP, "套\n"); flush.console() }
  est_one(d, hba1c_obs_all) %>% mutate(imp = i)
})
res_m <- bind_rows(res_m)

# Rubin 规则
pooled <- res_m %>% group_by(measure) %>%
  summarise(m = n(), qbar = mean(est), W = mean(se^2), B = var(est), .groups = "drop") %>%
  mutate(T = W + (1 + 1 / m) * B,
         df = ifelse(B > 0, (m - 1) * (1 + W / ((1 + 1 / m) * B))^2, Inf),
         lower95 = qbar - qt(0.975, pmin(df, 1e6)) * sqrt(T),
         upper95 = qbar + qt(0.975, pmin(df, 1e6)) * sqrt(T),
         frac_missing_info = ifelse(T > 0, (1 + 1 / m) * B / T, NA))

side <- res_cc %>% select(measure, cc_est = est, cc_lower95 = lower95, cc_upper95 = upper95) %>%
  full_join(pooled %>% select(measure, mi_est = qbar, mi_lower95 = lower95, mi_upper95 = upper95, frac_missing_info),
            by = "measure") %>%
  mutate(difference_mi_minus_cc = mi_est - cc_est)


# =====================================================================
# 第 5 部分：保存和打印
# =====================================================================
write.csv(miss_tab, file.path(OUT_DIR, "missing_by_variable.csv"), row.names = FALSE)
write.csv(cmp_tab,  file.path(OUT_DIR, "included_vs_excluded.csv"), row.names = FALSE)
write.csv(res_m,    file.path(OUT_DIR, "mi_each_imputation.csv"), row.names = FALSE)
write.csv(side,     file.path(OUT_DIR, "mi_vs_complete_case.csv"), row.names = FALSE)
png(file.path(OUT_DIR, "mi_convergence.png"), width = 2000, height = 1400, res = 160)
print(plot(imp)); invisible(dev.off())
writeLines(c(paste("锁定文件：", LOCKED_FILE), paste("MD5：", tools::md5sum(LOCKED_FILE)),
             paste("S1 切点文件：", ifelse(is.null(SEG_FILE), "无", SEG_FILE)),
             paste("插补套数 m =", M_IMP, "；迭代", MAXIT, "；种子 20261031"),
             "主分析为完整病例；本结果为敏感性分析。HbA1c 缺失者只参与插补，不进入评价。",
             paste("运行时间：", Sys.time())),
           file.path(OUT_DIR, "mi_info.txt"))
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))

cat("\n===== 各变量缺失比例（%）=====\n"); print(miss_tab, digits = 3, row.names = FALSE)
cat("\n===== 纳入 vs 被排除（加权）=====\n"); print(cmp_tab, digits = 3, row.names = FALSE)
cat("\n===== 完整病例 vs 多重插补（2013–2018，%）=====\n"); print(as.data.frame(side), digits = 3, row.names = FALSE)
cat("\n结果已保存到：", OUT_DIR, "\n")
