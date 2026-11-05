# =====================================================================
# 锁定方程的重校准（"更新"，补充分析）：NHANES 2013–2018
# ---------------------------------------------------------------------
# logistic 重校准：logit(p_new) = a + b × lp，lp 为锁定方程的线性预测值
#   a 修正整体偏高/偏低，b 修正"太极端"（b<1 时把高的压低、低的抬高）
# 重校准不改变谁高谁低，所以 AUC、按名额转检的敏感度都不变；只改变概率数值。
# 步骤：
#   ① 三个周期轮流留一：两个周期估计 a、b，在第三个周期检查校准（诚实的评估）
#   ② 用全部 2013–2018 估计最终 a、b，给出更新后的方程（仍需外部数据验证）
# 用法：先跑过 refine_v5_equation.R；放到 scripts 文件夹，点 Source。
# =====================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(survey)
})
options(survey.lonely.psu = "adjust")


# =====================================================================
# 第 0 部分：设置（只改这里）
# =====================================================================
PROJECT_DIR <- "D:/Nottingham/Diabetes_Screening_NHANES"
DATA_DIR    <- file.path(PROJECT_DIR, "data", "raw")
OUT_DIR     <- file.path(PROJECT_DIR, "results",
                         paste0("recalibrate_", format(Sys.Date(), "%Y%m%d")))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

DEV_CYCLES  <- character(0)       # 不读开发期
EVAL_CYCLES <- c("H", "I", "J")   # 2013–2018

# 自动找到 results 里最新的 locked_main.csv（也可以直接写完整路径）
LOCKED_FILE <- list.files(file.path(PROJECT_DIR, "results"), pattern = "^locked_main\\.csv$",
                          recursive = TRUE, full.names = TRUE)
if (length(LOCKED_FILE) == 0) stop("results 里找不到 locked_main.csv，请先运行 refine_v5_equation.R")
LOCKED_FILE <- LOCKED_FILE[which.max(file.mtime(LOCKED_FILE))]
cat("使用的锁定文件：", LOCKED_FILE, "\n")

locked <- read.csv(LOCKED_FILE)
val <- setNames(locked$value, locked$term)
KNOTS <- unname(val[c("KNOT1", "KNOT2", "KNOT3")])
THRESHOLD <- unname(val["THRESHOLD"])
BETA <- val[!names(val) %in% c("THRESHOLD", "KNOT1", "KNOT2", "KNOT3", "C_L2")]


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
# 第 2 部分：锁定方程的线性预测值
# =====================================================================
lp_of <- function(beta, d) as.vector(beta[1] + as.matrix(d[, names(beta)[-1]]) %*% beta[-1])
full$lp <- lp_of(BETA, full)
ev <- full[full$eval, ]

wmean <- function(x, w) sum(w * x) / sum(w)

# 加权重校准拟合：返回 a、b
fit_recal <- function(d) {
  nw <- d$wt / mean(d$wt)
  cf <- coef(suppressWarnings(glm(case ~ lp, data = d, family = binomial, weights = nw)))
  c(a = unname(cf[1]), b = unname(cf[2]))
}
# 加权校准指标
calib_stats <- function(lp, d) {
  p <- plogis(lp); nw <- d$wt / mean(d$wt)
  c(OE = wmean(d$case, d$wt) / wmean(p, d$wt),
    cal_int = unname(coef(suppressWarnings(glm(d$case ~ offset(lp), family = binomial, weights = nw)))[1]),
    cal_slope = unname(coef(suppressWarnings(glm(d$case ~ lp, family = binomial, weights = nw)))[2]))
}


# =====================================================================
# 第 3 部分：三个周期轮流留一（在没参与估计 a、b 的周期上检查）
# =====================================================================
loco <- bind_rows(lapply(EVAL_CYCLES, function(cy) {
  tr <- ev[ev$cycle != cy, ]; te <- ev[ev$cycle == cy, ]
  ab <- fit_recal(tr)
  before <- calib_stats(te$lp, te)
  after  <- calib_stats(ab["a"] + ab["b"] * te$lp, te)
  data.frame(heldout_cycle = cy, a = ab["a"], b = ab["b"],
             OE_before = before["OE"], OE_after = after["OE"],
             int_before = before["cal_int"], int_after = after["cal_int"],
             slope_before = before["cal_slope"], slope_after = after["cal_slope"], row.names = NULL)
}))
cat("\n===== 轮流留一：重校准前后（留出周期上）=====\n"); print(loco, digits = 3)


# =====================================================================
# 第 4 部分：用全部 2013–2018 估计最终 a、b（调查设计给区间）
# =====================================================================
des   <- svydesign(ids = ~psu, strata = ~strata, weights = ~wt, nest = TRUE, data = full)
des_e <- subset(des, eval)
fit_ab <- svyglm(case ~ lp, design = des_e, family = quasibinomial())
ab <- coef(fit_ab); ab_ci <- confint(fit_ab)
recal_tab <- data.frame(parameter = c("a (截距修正)", "b (斜率修正)"), estimate = ab,
                        lower95 = ab_ci[, 1], upper95 = ab_ci[, 2], row.names = NULL)
cat("\n===== 最终重校准参数 =====\n"); print(recal_tab, digits = 4)

# 更新后的方程：β0' = a + b·β0，βj' = b·βj
beta_new <- c(ab[1] + ab[2] * BETA[1], ab[2] * BETA[-1]); names(beta_new) <- names(BETA)
coef_tab <- data.frame(term = names(BETA), locked = BETA, recalibrated = beta_new, row.names = NULL)
cat("\n===== 锁定方程 vs 重校准后方程 =====\n"); print(coef_tab, digits = 6)

# 同一个转检切点在新概率尺度上对应的数值（谁被转检完全不变）
thr_new <- plogis(ab[1] + ab[2] * qlogis(THRESHOLD))
cat("\n原切点", signif(THRESHOLD, 6), "在重校准后对应", signif(thr_new, 6), "（转检的人完全相同）\n")


# =====================================================================
# 第 5 部分：校准图（重校准前 vs 后）
# =====================================================================
ev$p_before <- plogis(ev$lp)
ev$p_after  <- plogis(ab[1] + ab[2] * ev$lp)
ev$dec <- cut(ev$p_before, quantile(ev$p_before, 0:10 / 10), include.lowest = TRUE, labels = FALSE)
dec_tab <- ev %>% group_by(dec) %>%
  summarise(pred_before = wmean(p_before, wt), pred_after = wmean(p_after, wt),
            observed = wmean(case, wt), n = n(), cases = sum(case), .groups = "drop")
png(file.path(OUT_DIR, "calibration_before_after.png"), width = 2400, height = 1200, res = 220)
par(mfrow = c(1, 2))
lim <- c(0, max(dec_tab$pred_before, dec_tab$observed) * 1.05)
plot(dec_tab$pred_before, dec_tab$observed, xlim = lim, ylim = lim, pch = 19,
     xlab = "Mean predicted risk", ylab = "Observed", main = "Locked equation, 2013-2018")
abline(0, 1, lty = 2, col = "grey40")
plot(dec_tab$pred_after, dec_tab$observed, xlim = lim, ylim = lim, pch = 19, col = "#1f77b4",
     xlab = "Mean predicted risk", ylab = "Observed", main = "After recalibration (update)")
abline(0, 1, lty = 2, col = "grey40")
invisible(dev.off())


# =====================================================================
# 保存
# =====================================================================
write.csv(loco,      file.path(OUT_DIR, "recal_leave_one_cycle_out.csv"), row.names = FALSE)
write.csv(recal_tab, file.path(OUT_DIR, "recal_parameters.csv"), row.names = FALSE)
write.csv(coef_tab,  file.path(OUT_DIR, "recal_coefficients.csv"), row.names = FALSE)
write.csv(dec_tab,   file.path(OUT_DIR, "recal_deciles.csv"), row.names = FALSE)
writeLines(c(paste("锁定文件：", LOCKED_FILE),
             paste("重校准 a =", format(ab[1], digits = 10), "；b =", format(ab[2], digits = 10)),
             paste("原切点", format(THRESHOLD, digits = 10), "→ 重校准后", format(thr_new, digits = 10)),
             "这是用 2013–2018 估计的更新，需要外部数据验证。",
             paste("运行时间：", Sys.time())),
           file.path(OUT_DIR, "recal_info.txt"))
cat("\n结果已保存到：", OUT_DIR, "\n")
