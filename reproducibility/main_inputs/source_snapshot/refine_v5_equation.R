# =====================================================================
# V5 附录 B 方程完善脚本：L2 调查加权拟合 + 调查设计分析 + 校准分析
# ---------------------------------------------------------------------
# 严格按 V5 附录 B：
#   η = β0 + f_age(A) + βW·w + βM·M + βF·F + βH·H + βI·I ，p = 1/(1+exp(−η))
#   A = 年龄/10，a = A − 5，w = 10×(WHtR − 0.55)
#   f_age 两种：线性 βA·a；或三结点样条 βA·a + βS·s，s = g(A) − g(5)，结点 25/46/73 岁
#   M3 扩展：样条 M2 + βAW·a·w
#   六项全部保留，不逐步删项
#   拟合：调查加权 Logistic，L2 C=1，截距不惩罚；训练折内标准化，
#         权重按该折原始行数归一化；系数还原回原单位
#
# 用法：放到 D:\Nottingham\Diabetes_Screening_NHANES\scripts\，只改第 0 部分，点 Source。
#       中文乱码时：File → Reopen with Encoding → UTF-8
# 输出：results\refine_v5_日期\
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
OUT_DIR     <- file.path(PROJECT_DIR, "results",
                         paste0("refine_v5_", format(Sys.Date(), "%Y%m%d")))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

DEV_CYCLES  <- c("E", "F", "G")   # V5：2007–2012 开发
EVAL_CYCLES <- c("H", "I", "J")   # V5：2013–2018 评价
RUN_EVAL    <- FALSE              # 锁定之前保持 FALSE

C_L2   <- 1      # V5 固定 C=1；改成 Inf 即不惩罚（普通调查加权 Logistic）
N_BOOT <- 1000   # 正式结果 1000；想先快速试跑可改成 200
KNOTS  <- c(2.5, 4.6, 7.3)
MAIN_MODEL <- "M2_linear"   # 10-28 根据开发期结果选定：线性年龄版为主方程

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
# 第 2 部分：V5 附录 B 的三个方程
# =====================================================================
MODELS <- list(
  M2_linear = c("a", "w", "male", "famhx", "htn", "inactive"),          # 6 个斜率
  M2_spline = c("a", "s", "w", "male", "famhx", "htn", "inactive"),     # 7 个斜率
  M3_inter  = c("a", "s", "w", "male", "famhx", "htn", "inactive", "aw") # 8 个斜率
)


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

# 一次完整开发流程：留一周期预测 → V5 阈值 → 全开发期拟合
# w = 开发期每人的权重（主分析用调查权重；bootstrap 时用重抽权重）
run_pipeline <- function(vars, d, w) {
  X <- as.matrix(d[, vars]); lp_loco <- rep(NA_real_, nrow(d))
  for (cy in DEV_CYCLES) {
    tr <- d$cycle != cy; te <- d$cycle == cy
    b  <- fit_l2(X[tr, , drop = FALSE], d$case[tr], w[tr])
    lp_loco[te] <- lp_of(b, d[te, ])
  }
  k <- w > 0
  target <- wmean(d$ref_bang[k], w[k])
  thr <- v5_threshold(plogis(lp_loco[k]), w[k], target)
  beta <- fit_l2(X, d$case, w)
  list(beta = beta, lp_loco = lp_loco, thr = thr, target = target)
}

summ <- function(res, d, w) {   # 由留一周期预测算主要指标
  k <- w > 0; d <- d[k, ]; w <- w[k]; lp <- res$lp_loco[k]; p <- plogis(lp)
  ref <- as.numeric(p >= res$thr)
  s <- function(r, keep = TRUE) { kk <- d$case == 1 & keep; sum(w[kk] * r[kk]) / sum(w[kk]) }
  c(AUC = w_auc(p, d$case, w), AUC_bang = w_auc(d$bang, d$case, w),
    OE = wmean(d$case, w) / wmean(p, w), calib_pts(lp, d$case, w),
    refer = wmean(ref, w), refer_bang = res$target,
    sens = s(ref), sens_bang = s(d$ref_bang), sens_diff = s(ref) - s(d$ref_bang),
    young_diff = s(ref, d$age_grp == "20-39") - s(d$ref_bang, d$age_grp == "20-39"),
    mid_diff   = s(ref, d$age_grp == "40-59") - s(d$ref_bang, d$age_grp == "40-59"),
    old_diff   = s(ref, d$age_grp == ">=60")  - s(d$ref_bang, d$age_grp == ">=60"),
    gap_access = s(ref, d$access == "high") - s(ref, d$access == "low"),
    gap_bang   = s(d$ref_bang, d$access == "high") - s(d$ref_bang, d$access == "low"))
}


# =====================================================================
# 第 5 部分：主分析（开发期，调查权重）
# =====================================================================
dd <- full[full$dev, ]
main <- lapply(MODELS, function(v) run_pipeline(v, dd, dd$wt))
main_tab <- sapply(names(MODELS), function(m) summ(main[[m]], dd, dd$wt))


# =====================================================================
# 第 6 部分：PSU 分层重抽 bootstrap（V5 附录 C：层内抽 m−1 个 PSU 并重缩放）
# 每次重抽对三个方程重做整个流程，得到系数和配对差值的 95% 区间
# =====================================================================
des <- svydesign(ids = ~psu, strata = ~strata, weights = ~wt, nest = TRUE, data = full)
rep_des <- as.svrepdesign(des, type = "subbootstrap", replicates = N_BOOT)
W <- weights(rep_des, type = "analysis")[full$dev, ]
cat("\n正在重抽", N_BOOT, "次，每次重做 3 个方程的完整流程，请稍候……\n")
t0 <- Sys.time()
boot <- lapply(seq_len(N_BOOT), function(b) {
  if (b %% 20 == 0) {   # 每 20 次报告一次进度和预计剩余时间
    el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    cat(sprintf("  已完成 %d / %d 次，已用 %.1f 分钟，预计还需 %.1f 分钟\n", b, N_BOOT, el, el / b * (N_BOOT - b)))
    flush.console()
  }
  lapply(MODELS, function(v) {
    r <- run_pipeline(v, dd, W[, b])
    list(beta = r$beta, stats = summ(r, dd, W[, b]))
  })
})
# 区间算法（10-28 修正）：点估计 ± 1.96 × bootstrap 标准误
# 原因：NHANES 每层只有 2 个 PSU，层内抽 m−1=1 个，重抽分布会偏到一侧，
#       百分位区间可能不包含点估计；用标准误法更稳。百分位区间仍另存供对照。
se_of <- function(m) apply(m, 1, sd, na.rm = TRUE)
pct_of <- function(m) t(apply(m, 1, quantile, c(0.025, 0.975), na.rm = TRUE))

perf <- bind_rows(lapply(names(MODELS), function(m) {
  bs <- sapply(boot, function(b) b[[m]]$stats); se <- se_of(bs); pc <- pct_of(bs)
  data.frame(model = m, measure = rownames(main_tab), estimate = main_tab[, m], se = se,
             lower95 = main_tab[, m] - 1.96 * se, upper95 = main_tab[, m] + 1.96 * se,
             pct_lo = pc[, 1], pct_hi = pc[, 2], row.names = NULL)
}))
coefs <- bind_rows(lapply(names(MODELS), function(m) {
  bb <- sapply(boot, function(b) b[[m]]$beta); se <- se_of(bb); b0 <- main[[m]]$beta
  data.frame(model = m, term = names(b0), coef_L2 = b0, se = se,
             lower95 = b0 - 1.96 * se, upper95 = b0 + 1.96 * se,
             OR = exp(b0), OR_lo = exp(b0 - 1.96 * se), OR_hi = exp(b0 + 1.96 * se), row.names = NULL)
}))
coefs[coefs$term == "(Intercept)", c("OR", "OR_lo", "OR_hi")] <- NA

# 配对比较：样条 vs 线性、交互 vs 样条（同一次重抽内相减）
pair <- function(m1, m0, st) {
  dif <- sapply(boot, function(b) b[[m1]]$stats[st] - b[[m0]]$stats[st])
  est <- main_tab[st, m1] - main_tab[st, m0]; se <- sd(dif, na.rm = TRUE)
  data.frame(compare = paste(m1, "-", m0), measure = st, estimate = est, se = se,
             lower95 = est - 1.96 * se, upper95 = est + 1.96 * se, row.names = NULL)
}
pairs <- bind_rows(lapply(c("sens", "young_diff", "old_diff", "AUC"), function(st)
  rbind(pair("M2_spline", "M2_linear", st), pair("M3_inter", "M2_spline", st))))


# =====================================================================
# 第 7 部分：调查设计分析（Taylor 线性化，含权重、层、PSU）
# ---------------------------------------------------------------------
# ① 不惩罚的调查加权 Logistic（svyglm），给出设计校正的标准误和区间，
#    用来和 L2 系数对照；两者接近说明 C=1 的惩罚影响很小
# ② 留一周期预测下的敏感度，用调查设计算区间
# =====================================================================
des_dev <- subset(des, dev)
svy_coefs <- bind_rows(lapply(names(MODELS), function(m) {
  f <- reformulate(MODELS[[m]], response = "case")
  fit <- svyglm(f, design = des_dev, family = quasibinomial())
  data.frame(model = m, term = names(coef(fit)), coef_svyglm = coef(fit),
             se = SE(fit), lower95 = confint(fit)[, 1], upper95 = confint(fit)[, 2], row.names = NULL)
}))

for (m in names(MODELS)) {
  full[[paste0("lp_", m)]] <- NA_real_
  full[[paste0("lp_", m)]][full$dev] <- main[[m]]$lp_loco
  full[[paste0("ref_", m)]] <- as.numeric(plogis(full[[paste0("lp_", m)]]) >= main[[m]]$thr)
}
des <- svydesign(ids = ~psu, strata = ~strata, weights = ~wt, nest = TRUE, data = full)
des_dev <- subset(des, dev)
des_cases <- subset(des_dev, case == 1)

taylor_sens <- bind_rows(lapply(c("ref_bang", paste0("ref_", names(MODELS))), function(r) {
  f <- as.formula(paste0("~", r))
  by_age <- svyby(f, ~age_grp, des_cases, svymean)
  by_acc <- svyby(f, ~access, des_cases, svymean)
  tot <- svymean(f, des_cases)
  data.frame(rule = sub("ref_", "", r),
             group = c("overall", paste0("age ", by_age$age_grp), paste0("access ", by_acc$access)),
             sens = c(coef(tot), by_age[[r]], by_acc[[r]]),
             se = c(SE(tot), SE(by_age), SE(by_acc)), row.names = NULL) %>%
    mutate(lower95 = sens - 1.96 * se, upper95 = sens + 1.96 * se)
}))


# =====================================================================
# 第 8 部分：校准分析（留一周期预测，调查设计）
# ---------------------------------------------------------------------
# O:E = 实际病例比例 / 平均预测风险（理想 1；<1 表示预测偏高）
# 校准截距：固定斜率为 1 时的偏移（理想 0）
# 校准斜率：理想 1；<1 表示预测太极端（高的太高、低的太低）
# 重校准：用 logit(p_new) = a + b × lp 把概率拉回，a、b 由留一周期预测估计
# =====================================================================
calib <- bind_rows(lapply(names(MODELS), function(m) {
  f_int <- as.formula(paste0("case ~ offset(lp_", m, ")"))
  f_slp <- as.formula(paste0("case ~ lp_", m))
  ci_fit <- svyglm(f_int, design = des_dev, family = quasibinomial())
  cs_fit <- svyglm(f_slp, design = des_dev, family = quasibinomial())
  oe <- svyratio(~case, as.formula(paste0("~I(plogis(lp_", m, "))")), des_dev)
  ab <- coef(cs_fit)
  b  <- main[[m]]$beta
  recal <- c(ab[1] + ab[2] * b[1], ab[2] * b[-1])    # 重校准后的系数
  data.frame(model = m,
             OE = coef(oe), OE_lo = confint(oe)[1], OE_hi = confint(oe)[2],
             cal_int = coef(ci_fit)[1], int_lo = confint(ci_fit)[1, 1], int_hi = confint(ci_fit)[1, 2],
             cal_slope = ab[2], slope_lo = confint(cs_fit)[2, 1], slope_hi = confint(cs_fit)[2, 2],
             recal_a = ab[1], recal_b = ab[2],
             recal_coefs = paste(sprintf("%s=%.6f", names(b), recal), collapse = "; "),
             row.names = NULL)
}))

# 校准图：十分位 + 调查设计 95% 区间
png(file.path(OUT_DIR, "calibration_plot.png"), width = 2400, height = 900, res = 220)
par(mfrow = c(1, length(MODELS)), mar = c(4.5, 4.5, 3, 1))
for (m in names(MODELS)) {
  full$p_tmp <- plogis(full[[paste0("lp_", m)]])
  full$dec <- NA
  full$dec[full$dev] <- cut(full$p_tmp[full$dev], quantile(full$p_tmp[full$dev], 0:10 / 10),
                            include.lowest = TRUE, labels = FALSE)
  dz <- subset(svydesign(ids = ~psu, strata = ~strata, weights = ~wt, nest = TRUE, data = full), dev)
  ob <- svyby(~case, ~dec, dz, svymean); pr <- svyby(~p_tmp, ~dec, dz, svymean)
  lim <- c(0, max(ob$case + 1.96 * ob$se, pr$p_tmp) * 1.05)
  plot(pr$p_tmp, ob$case, xlim = lim, ylim = lim, pch = 19,
       xlab = "Mean predicted risk", ylab = "Observed (HbA1c >= 6.5%)", main = m)
  arrows(pr$p_tmp, pmax(ob$case - 1.96 * ob$se, 0), pr$p_tmp, ob$case + 1.96 * ob$se,
         angle = 90, code = 3, length = 0.03)
  abline(0, 1, lty = 2, col = "grey40")
}
invisible(dev.off())


# =====================================================================
# 第 9 部分（锁定后才运行）：评价期一次性检查 + 截距更新
# =====================================================================
if (RUN_EVAL && any(full$eval)) {
  ev <- full[full$eval, ]
  eval_tab <- bind_rows(lapply(names(MODELS), function(m) {
    lp <- lp_of(main[[m]]$beta, ev); p <- plogis(lp); ref <- as.numeric(p >= main[[m]]$thr)
    s <- function(r, keep = TRUE) { kk <- ev$case == 1 & keep; sum(ev$wt[kk] * r[kk]) / sum(ev$wt[kk]) }
    nw <- ev$wt / mean(ev$wt)
    upd <- coef(suppressWarnings(glm(ev$case ~ offset(lp), family = binomial, weights = nw)))[1]
    data.frame(model = m, AUC = w_auc(p, ev$case, ev$wt), OE = wmean(ev$case, ev$wt) / wmean(p, ev$wt),
               t(calib_pts(lp, ev$case, ev$wt)), refer = wmean(ref, ev$wt), refer_bang = wmean(ev$ref_bang, ev$wt),
               sens = s(ref), sens_bang = s(ev$ref_bang),
               young = s(ref, ev$age_grp == "20-39"), old = s(ref, ev$age_grp == ">=60"),
               intercept_update = upd, row.names = NULL)
  }))
  cat("\n===== 评价期（方程、阈值均未改动；intercept_update 属于“更新”，不是验证）=====\n")
  print(eval_tab, digits = 3)
  write.csv(eval_tab, file.path(OUT_DIR, "eval_results.csv"), row.names = FALSE)
}


# =====================================================================
# 第 10 部分：保存和打印
# =====================================================================
thr_tab <- data.frame(model = names(MODELS),
                      threshold = sapply(main, `[[`, "thr"),
                      bang_target = sapply(main, `[[`, "target"), row.names = NULL)
for (m in names(MODELS)) {
  write.csv(data.frame(term = c(names(main[[m]]$beta), "THRESHOLD", "KNOT1", "KNOT2", "KNOT3", "C_L2"),
                       value = c(main[[m]]$beta, main[[m]]$thr, KNOTS, C_L2)),
            file.path(OUT_DIR, paste0("locked_", m, ".csv")), row.names = FALSE)
}
write.csv(perf,        file.path(OUT_DIR, "performance_loco.csv"), row.names = FALSE)
write.csv(coefs,       file.path(OUT_DIR, "coefficients_L2.csv"),  row.names = FALSE)
write.csv(svy_coefs,   file.path(OUT_DIR, "coefficients_svyglm.csv"), row.names = FALSE)
write.csv(pairs,       file.path(OUT_DIR, "paired_comparisons.csv"), row.names = FALSE)
write.csv(taylor_sens, file.path(OUT_DIR, "sensitivity_taylor.csv"), row.names = FALSE)
write.csv(calib,       file.path(OUT_DIR, "calibration.csv"), row.names = FALSE)
write.csv(thr_tab,     file.path(OUT_DIR, "thresholds.csv"), row.names = FALSE)

# ---- 锁定主方程：写入 locked_main.csv（评价脚本只读这个文件）----
bm <- main[[MAIN_MODEL]]
write.csv(data.frame(term = c(names(bm$beta), "THRESHOLD", "KNOT1", "KNOT2", "KNOT3", "C_L2"),
                     value = c(bm$beta, bm$thr, KNOTS, C_L2)),
          file.path(OUT_DIR, "locked_main.csv"), row.names = FALSE)
writeLines(c(paste("主方程：", MAIN_MODEL),
             paste("变量：", paste(MODELS[[MAIN_MODEL]], collapse = ", ")),
             paste("开发周期：", paste(DEV_CYCLES, collapse = ",")),
             paste("L2 C =", C_L2, "；重抽次数 =", N_BOOT),
             paste("阈值：", format(bm$thr, digits = 12), "（开发期留一周期预测，匹配 Bang 转检比例", round(bm$target, 4), "）"),
             paste("锁定时间：", Sys.time()),
             "此后不再修改；如需修改记为“继续开发”，重新存档并说明原因。"),
           file.path(OUT_DIR, "locked_info.txt"))

cat("\n===== L2 系数（原单位）及 bootstrap 区间 =====\n"); print(coefs, digits = 4)
cat("\n===== 留一周期表现 =====\n"); print(round(main_tab, 4))
cat("\n===== 配对比较 =====\n"); print(pairs, digits = 3)
cat("\n===== 校准 =====\n"); print(calib[, 1:13], digits = 3)
cat("\n===== 阈值 =====\n"); print(thr_tab, digits = 6)
cat("\n结果已保存到：", OUT_DIR, "\n")
