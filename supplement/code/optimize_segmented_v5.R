# =====================================================================
# 分段方程的系统优化（只在 2007–2012 内选择）+ 2013–2018 一次性比较
# ---------------------------------------------------------------------
# 候选方案（全部用 V5 的六项变量、调查加权 L2，截距不惩罚）：
#   A. 硬分段：分界年龄 40/45/50/55/60/65 × 惩罚 C 0.1/0.3/1/3/不惩罚
#              × 年龄写法（线性 / 三结点样条）
#   B. 平滑分段：一个方程 + 年龄与其余五项的交互项（系数随年龄连续变化，
#              没有分界处的跳变）× 同样的 C
#   参照 S0：已锁定的单一方程（线性年龄，C=1）
# 选择规则（事先规定）：开发期留一周期预测下，转检比例与 Bang 相同时的
#   总体敏感度最高（V5 主要终点）；并列时取 AUC 高、再取更简单者。
# 选定后锁定，在 2013–2018 上与 S0、Bang 比较一次（Taylor 区间）。
# 注意：从几十个方案里挑最好的，开发期成绩会偏乐观；2013–2018 是事后再评价。
#
# 用法：放到 scripts 文件夹，只改第 0 部分，点 Source。约 3–10 分钟。
# 输出：results\optimize_segmented_日期\
# =====================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(survey)
})
options(survey.lonely.psu = "adjust")
set.seed(20261031)


# =====================================================================
# 第 0 部分：设置（只改这里）
# =====================================================================
PROJECT_DIR <- "D:/Nottingham/Diabetes_Screening_NHANES"
DATA_DIR    <- file.path(PROJECT_DIR, "data", "raw")
OUT_DIR     <- file.path(PROJECT_DIR, "results",
                         paste0("optimize_segmented_", format(Sys.Date(), "%Y%m%d")))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

DEV_CYCLES  <- c("E", "F", "G")
EVAL_CYCLES <- c("H", "I", "J")
C_L2   <- 1
KNOTS  <- c(2.5, 4.6, 7.3)
SPLITS <- c(40, 45, 50, 55, 60, 65)
C_GRID <- c(0.1, 0.3, 1, 3, Inf)
BASE   <- c("w", "male", "famhx", "htn", "inactive")   # 年龄以外的五项

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
    ref_bang  = as.numeric(bang >= 5),
    # USPSTF 2021：35–70 岁且超重/肥胖（BMI≥25；NHANES 无法单独处理亚裔 BMI≥23 的标准）
    ref_uspstf = as.numeric(age >= 35 & age <= 70 & bmi >= 25)
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
# 第 3 部分：L2 调查加权 Logistic（复刻 V5 附录 B 的训练实现）
# ---------------------------------------------------------------------
# 目标函数（与 scikit-learn 相同）：
#   0.5 × Σγ² + C × Σ w_i × [log(1+e^η) − y·η]，η = α + Zγ，截距 α 不惩罚
# Z 是训练折内按加权均值 μ、加权标准差 σ 标准化后的变量
# 拟合后还原：β_j = γ_j / σ_j，β0 = α − Σ γ_j μ_j / σ_j
# 权重先归一化到"该折原始行数"（含零权重行），再删去零权重行拟合
# =====================================================================
fit_l2 <- function(X, y, w, C = C_L2, n_norm = length(w)) {
  w <- w * n_norm / sum(w)
  keep <- w > 0; X <- X[keep, , drop = FALSE]; y <- y[keep]; w <- w[keep]
  mu <- colSums(X * w) / sum(w)
  sg <- sqrt(colSums(w * sweep(X, 2, mu)^2) / sum(w))
  Z  <- sweep(sweep(X, 2, mu), 2, sg, "/")
  lam <- if (is.infinite(C)) 0 else 1 / C        # 把 C 换成惩罚强度
  obj <- function(par) {
    eta <- par[1] + Z %*% par[-1]
    0.5 * lam * sum(par[-1]^2) + sum(w * (log1p(exp(-abs(eta))) + pmax(eta, 0) - y * eta))
  }
  grad <- function(par) {
    p <- plogis(par[1] + Z %*% par[-1]); r <- w * (p - y)
    c(sum(r), lam * par[-1] + colSums(Z * as.vector(r)))
  }
  op <- optim(rep(0, ncol(Z) + 1), obj, grad, method = "BFGS",
              control = list(maxit = 2000, reltol = 1e-12))
  if (op$convergence != 0) warning("优化未完全收敛")
  g <- op$par[-1]
  beta <- c("(Intercept)" = op$par[1] - sum(g * mu / sg), g / sg)
  names(beta)[-1] <- colnames(X)
  beta
}
lp_of <- function(beta, d) as.vector(beta[1] + as.matrix(d[, names(beta)[-1]]) %*% beta[-1])


# =====================================================================
# 第 4 部分：工具函数
# =====================================================================
wmean <- function(x, w) sum(w * x) / sum(w)

w_auc <- function(score, y, w) {
  o <- order(score); score <- score[o]; y <- y[o]; w <- w[o]
  w0 <- ifelse(y == 0, w, 0); w1 <- ifelse(y == 1, w, 0)
  grp <- match(score, unique(score))
  W0g <- tapply(w0, grp, sum); W1g <- tapply(w1, grp, sum)
  below <- cumsum(W0g) - W0g
  sum(W1g * (below + 0.5 * W0g)) / (sum(w1) * sum(w0))
}

# V5 阈值规则：在留一周期预测中，选转检比例最接近目标的完整并列块阈值；偏差相同取较高阈值
v5_threshold <- function(p, w, target) {
  u <- sort(unique(p), decreasing = TRUE)
  share <- sapply(u, function(t) sum(w[p >= t]) / sum(w))
  dev <- abs(share - target)
  max(u[dev == min(dev)])
}

calib_pts <- function(lp, y, w) {   # 加权校准截距、斜率（点估计）
  nw <- w / mean(w)
  c(cal_int   = unname(coef(suppressWarnings(glm(y ~ offset(lp), family = binomial, weights = nw)))[1]),
    cal_slope = unname(coef(suppressWarnings(glm(y ~ lp, family = binomial, weights = nw)))[2]))
}

# =====================================================================
# 第 5 部分：候选方案
# =====================================================================
for (v in BASE) full[[paste0("a_", v)]] <- full$a * full[[v]]   # 平滑分段用的交互项
INTER <- paste0("a_", BASE)

specs <- list(S0_reference = list(type = "single", vars = c("a", BASE), C = 1, split = NA, age = "linear"))
for (sp in SPLITS) for (C in C_GRID) for (af in c("linear", "spline")) {
  nm <- sprintf("split%d_C%s_%s", sp, ifelse(is.infinite(C), "Inf", C), af)
  specs[[nm]] <- list(type = "split", split = sp, C = C, age = af,
                      vars = c(if (af == "linear") "a" else c("a", "s"), BASE))
}
for (C in C_GRID) {
  nm <- sprintf("smooth_C%s", ifelse(is.infinite(C), "Inf", C))
  specs[[nm]] <- list(type = "single", vars = c("a", BASE, INTER), C = C, split = NA, age = "linear+inter")
}
cat("\n候选方案共", length(specs), "个\n")

fit_spec <- function(sp, d, w) {
  if (sp$type == "single") return(list(all = fit_l2(as.matrix(d[, sp$vars]), d$case, w, C = sp$C)))
  g <- ifelse(d$age >= sp$split, "old", "young")
  lapply(c(young = "young", old = "old"), function(s) {
    k <- g == s
    fit_l2(as.matrix(d[k, sp$vars]), d$case[k], w[k], C = sp$C)
  })
}
pred_spec <- function(sp, fit, d) {
  if (sp$type == "single") return(lp_of(fit$all, d))
  g <- ifelse(d$age >= sp$split, "old", "young"); lp <- rep(NA_real_, nrow(d))
  for (s in c("young", "old")) if (any(g == s)) lp[g == s] <- lp_of(fit[[s]], d[g == s, ])
  lp
}
loco_spec <- function(sp, d) {
  lp <- rep(NA_real_, nrow(d))
  for (cy in DEV_CYCLES) {
    tr <- d$cycle != cy; te <- d$cycle == cy
    lp[te] <- pred_spec(sp, suppressWarnings(fit_spec(sp, d[tr, ], d$wt[tr])), d[te, ])
  }
  lp
}
wsens <- function(ref, d, keep = TRUE) { k <- d$case == 1 & keep; sum(d$wt[k] * ref[k]) / sum(d$wt[k]) }


# =====================================================================
# 第 6 部分：开发期留一周期比较，按事先规定的规则选出最佳方案
# =====================================================================
dd <- full[full$dev, ]
target_ref <- wmean(dd$ref_bang, dd$wt)
t0 <- Sys.time()
grid <- bind_rows(lapply(seq_along(specs), function(i) {
  sp <- specs[[i]]
  if (i %% 10 == 0) {
    cat(sprintf("  已完成 %d / %d 个方案，用时 %.1f 分钟\n", i, length(specs),
                as.numeric(difftime(Sys.time(), t0, units = "mins")))); flush.console()
  }
  p <- plogis(loco_spec(sp, dd))
  thr <- v5_threshold(p, dd$wt, target_ref)
  ref <- as.numeric(p >= thr)
  data.frame(spec = names(specs)[i], type = sp$type, split = sp$split, C = sp$C, age = sp$age,
             n_params = length(sp$vars) * ifelse(sp$type == "split", 2, 1) + ifelse(sp$type == "split", 2, 1),
             AUC = w_auc(p, dd$case, dd$wt), referral = 100 * wmean(ref, dd$wt),
             sens = 100 * wsens(ref, dd),
             sens_20_39 = 100 * wsens(ref, dd, dd$age_grp == "20-39"),
             sens_40_59 = 100 * wsens(ref, dd, dd$age_grp == "40-59"),
             sens_60plus = 100 * wsens(ref, dd, dd$age_grp == ">=60"))
}))
s0 <- grid[grid$spec == "S0_reference", ]
grid <- grid %>% mutate(sens_vs_S0 = sens - s0$sens, AUC_vs_S0 = AUC - s0$AUC,
                        sens60_vs_S0 = sens_60plus - s0$sens_60plus) %>%
  arrange(desc(round(sens, 6)), desc(AUC), n_params)
best_name <- grid$spec[grid$spec != "S0_reference"][1]
best <- specs[[best_name]]
cat("\n开发期选出的最佳方案：", best_name, "\n")


# =====================================================================
# 第 7 部分：锁定最佳方案，2013–2018 一次性比较
# =====================================================================
lock_one <- function(sp) {
  p <- plogis(loco_spec(sp, dd))
  list(fit = suppressWarnings(fit_spec(sp, dd, dd$wt)), thr = v5_threshold(p, dd$wt, target_ref))
}
L0 <- lock_one(specs$S0_reference); LB <- lock_one(best)

coef_tab <- bind_rows(lapply(names(LB$fit), function(s)
  data.frame(spec = best_name, part = s, term = names(LB$fit[[s]]), coef = unname(LB$fit[[s]]))))
coef_tab$OR <- ifelse(coef_tab$term == "(Intercept)", NA, exp(coef_tab$coef))

eval_tab <- NULL
if (any(full$eval)) {
  ee <- full$eval
  full$ref_S0 <- NA_real_; full$ref_best <- NA_real_
  full$ref_S0[ee]   <- as.numeric(plogis(pred_spec(specs$S0_reference, L0$fit, full[ee, ])) >= L0$thr)
  full$ref_best[ee] <- as.numeric(plogis(pred_spec(best, LB$fit, full[ee, ])) >= LB$thr)
  full$p_best <- NA_real_; full$p_best[ee] <- plogis(pred_spec(best, LB$fit, full[ee, ]))
  full$p_S0 <- NA_real_;   full$p_S0[ee]   <- plogis(pred_spec(specs$S0_reference, L0$fit, full[ee, ]))
  full$d_best_S0   <- full$ref_best - full$ref_S0
  full$d_best_Bang <- full$ref_best - full$ref_bang
  des <- svydesign(ids = ~psu, strata = ~strata, weights = ~wt, nest = TRUE, data = full)
  tay <- function(var, idx) {
    m <- svymean(as.formula(paste0("~", var)), des[idx, ], na.rm = TRUE); ci <- confint(m)
    100 * c(estimate = unname(coef(m)), lower95 = ci[1], upper95 = ci[2])
  }
  rows <- list()
  for (g in c("overall", "20-39", "40-59", ">=60")) {
    gk <- if (g == "overall") rep(TRUE, nrow(full)) else full$age_grp %in% g
    a <- ee & gk; c1 <- a & full$case == 1
    for (x in list(c("referral", "Bang", "ref_bang", "a"), c("referral", "S0", "ref_S0", "a"),
                   c("referral", "best", "ref_best", "a"),
                   c("sensitivity", "Bang", "ref_bang", "c"), c("sensitivity", "S0", "ref_S0", "c"),
                   c("sensitivity", "best", "ref_best", "c"),
                   c("sens_diff_best_vs_S0", "best", "d_best_S0", "c"),
                   c("sens_diff_best_vs_Bang", "best", "d_best_Bang", "c"),
                   c("refer_diff_best_vs_S0", "best", "d_best_S0", "a"))) {
      v <- tay(x[3], if (x[4] == "a") a else c1)
      rows[[length(rows) + 1]] <- data.frame(group = g, measure = x[1], strategy = x[2],
                                             estimate = v[1], lower95 = v[2], upper95 = v[3], row.names = NULL)
    }
  }
  eval_tab <- bind_rows(rows)
  auc_eval <- data.frame(AUC_S0 = w_auc(full$p_S0[ee], full$case[ee], full$wt[ee]),
                         AUC_best = w_auc(full$p_best[ee], full$case[ee], full$wt[ee]),
                         AUC_Bang = w_auc(full$bang[ee], full$case[ee], full$wt[ee]))
}


# =====================================================================
# 第 8 部分：图、保存和打印
# =====================================================================
png(file.path(OUT_DIR, "dev_grid_split.png"), width = 2400, height = 900, res = 200)
par(mfrow = c(1, 2), mar = c(4, 4.5, 3, 1))
for (af in c("linear", "spline")) {
  gg <- grid[grid$type == "split" & grid$age == af, ]
  plot(NA, xlim = range(SPLITS), ylim = range(c(grid$sens, s0$sens)), xlab = "Split age",
       ylab = "Sensitivity at Bang-matched referral (%)", main = paste0("Hard split, age ", af, " (2007-2012 LOCO)"))
  abline(h = s0$sens, lty = 2, col = "grey40")
  cols <- c("#d73027", "#fc8d59", "#4575b4", "#91bfdb", "black")
  for (j in seq_along(C_GRID)) {
    r <- gg[gg$C == C_GRID[j], ]; r <- r[order(r$split), ]
    lines(r$split, r$sens, type = "b", pch = 19, col = cols[j])
  }
  if (af == "linear") legend("bottomright", c(paste("C =", C_GRID), "S0"), col = c(cols, "grey40"),
                             lty = c(rep(1, 5), 2), pch = c(rep(19, 5), NA), bty = "n", cex = 0.75)
}
invisible(dev.off())

write.csv(grid, file.path(OUT_DIR, "dev_grid_all_specs.csv"), row.names = FALSE)
write.csv(coef_tab, file.path(OUT_DIR, "best_spec_coefficients.csv"), row.names = FALSE)
write.csv(data.frame(spec = c("S0_reference", best_name), threshold = c(L0$thr, LB$thr)),
          file.path(OUT_DIR, "best_spec_thresholds.csv"), row.names = FALSE)
if (!is.null(eval_tab)) {
  write.csv(eval_tab, file.path(OUT_DIR, "eval_best_vs_S0.csv"), row.names = FALSE)
  write.csv(auc_eval, file.path(OUT_DIR, "eval_auc.csv"), row.names = FALSE)
}
writeLines(c(paste("候选方案数：", length(specs)), paste("开发期最佳：", best_name),
             "选择规则：开发期留一周期、与 Bang 等转检比例下的总体敏感度最高（并列取 AUC、再取简单）。",
             "2013–2018 为事后再评价，不是独立验证。", paste("运行时间：", Sys.time())),
           file.path(OUT_DIR, "optimize_info.txt"))
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))

cat("\n===== 开发期前 10 名（及 S0）=====\n")
print(rbind(head(grid[grid$spec != "S0_reference", ], 10), grid[grid$spec == "S0_reference", ])[, c("spec", "AUC", "referral", "sens",
                                                                                                    "sens_20_39", "sens_40_59", "sens_60plus", "sens_vs_S0")], digits = 3, row.names = FALSE)
cat("\n===== 各类最佳（开发期）=====\n")
print(grid %>% group_by(type, age) %>% slice(1) %>% ungroup() %>%
        select(spec, AUC, sens, sens_60plus, sens_vs_S0) %>% as.data.frame(), digits = 3, row.names = FALSE)
if (!is.null(eval_tab)) {
  cat("\n===== 2013–2018：最佳方案 vs S0 vs Bang =====\n")
  et <- eval_tab; et$value <- sprintf("%.1f (%.1f to %.1f)", et$estimate, et$lower95, et$upper95)
  print(reshape(et[, c("group", "measure", "strategy", "value")] %>% mutate(key = paste(measure, strategy)) %>%
                  select(group, key, value), idvar = "key", timevar = "group", direction = "wide"), row.names = FALSE)
  cat("\nAUC：\n"); print(auc_eval, digits = 3)
}
cat("\n结果已保存到：", OUT_DIR, "\n")
