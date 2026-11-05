# =====================================================================
# 可及性差距：各筛查方案下高 / 低可及组的敏感度差（access_gap_v5.R）
# ---------------------------------------------------------------------
# 研究主问题：Bang 对高可及组与低可及组的敏感度差（15.7 个百分点）在新方程和
# 老年切点方案下是缩小还是扩大？差距有多少来自病例年龄构成？
#
# 比较的方案（方程和切点全部只用 2007–2012 确定，与之前脚本完全相同）：
#   Bang ≥5；V5_S0；V5_S1（老年单独切点，主方案）；V5_S4（≥60 全检）；
#   F4_S0 / F4_S1（4 项简化方程：年龄、腰高比、男性、家族史）
#
# 输出：
#   1) 每个方案在高 / 低可及组的转检、敏感度、特异度、PPV（Taylor 区间）
#   2) 敏感度差（高 − 低），及“差距相对 Bang 的变化”（配对，Taylor 区间）
#   3) 年龄标准化差距：按全部病例的年龄构成（20–39 / 40–59 / ≥60）直接标准化
#   4) 同名额曲线：总转检固定时，高 / 低可及组敏感度及差距，PSU bootstrap 95% 带
#
# 可及性只用来分组，不进方程。2013–2018 为事后再评价，不是独立验证。
# 用法：放到 scripts 文件夹，只改第 0 部分，点 Source。约 3–10 分钟。
# 输出：results\access_gap_v5_日期\
# =====================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(survey)
})
options(survey.lonely.psu = "adjust")
set.seed(20261103)


# =====================================================================
# 第 0 部分：设置（只改这里）
# =====================================================================
# Portable execution adapter. Optional positional arguments: raw_directory output_directory.
.audit_file_arg <- grep("^--file=", commandArgs(), value = TRUE)[1]
.audit_script_dir <- dirname(normalizePath(sub("^--file=", "", .audit_file_arg), winslash = "/", mustWork = TRUE))
.audit_args <- commandArgs(trailingOnly = TRUE)
PROJECT_DIR <- .audit_script_dir # optional original locked_main.csv diagnostic only
DATA_DIR <- if (length(.audit_args) >= 1) normalizePath(.audit_args[1], winslash="/", mustWork=TRUE) else file.path(.audit_script_dir, "..", "analysis", "raw")
OUT_DIR <- if (length(.audit_args) >= 2) .audit_args[2] else file.path(.audit_script_dir, "access_gap_rerun")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

DEV_CYCLES  <- c("E", "F", "G")
EVAL_CYCLES <- c("H", "I", "J")
C_L2      <- 1
KNOTS     <- c(2.5, 4.6, 7.3)
SPLIT_AGE <- 60
VARS      <- c("a", "w", "male", "famhx", "htn", "inactive")   # V5 六项（已锁定）
VARS4     <- c("a", "w", "male", "famhx")                      # 4 项简化版
B_BOOT    <- 500
BUDGETS   <- seq(0.25, 0.50, by = 0.005)


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
  # 取列；某周期没有的列补 NA（新增候选变量可能不是每个周期都有）
  sel <- function(d, v) { for (x in setdiff(v, names(d))) d[[x]] <- NA; d[, v] }
  sel(m("DEMO"), c("SEQN", "RIDSTATR", "RIDAGEYR", "RIAGENDR", "RIDEXPRG", "WTMEC2YR",
                   "SDMVSTRA", "SDMVPSU", "RIDRETH1")) %>%
    left_join(sel(m("DIQ"), c("SEQN", "DIQ010", "DIQ160", "DIQ170")), by = "SEQN") %>%
    left_join(m("BMX") %>% select(SEQN, BMXWAIST, BMXHT, BMXBMI), by = "SEQN") %>%
    left_join(sel(m("BPQ"), c("SEQN", "BPQ020", "BPQ080")), by = "SEQN") %>%
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

# ---- 新增候选特征（只做探索；缺失 / 不知道 / 拒答的自报项记为 0，不改变样本）----
full <- full %>%
  mutate(
    chol    = as.numeric(BPQ080 %in% 1),           # 自报曾被告知高胆固醇
    prediab = as.numeric(DIQ160 %in% 1),           # 自报曾被告知糖尿病前期（需要以前查过血）
    atrisk  = as.numeric(DIQ170 %in% 1),           # 自报曾被告知有糖尿病风险
    hisp    = as.numeric(RIDRETH1 %in% c(1, 2)),   # 种族/族裔（参照：非西裔白人）
    black   = as.numeric(RIDRETH1 %in% 4),
    othrace = as.numeric(RIDRETH1 %in% 5),
    bmi_c   = (bmi - 28) / 5,                      # BMI（与 w 高度相关，用来检验冗余）
    h       = (height - 168) / 10,                 # 身高（每 10 cm）
    ah      = pmax(age - 60, 0) / 10,              # 60 岁以后的年龄斜率变化（线性样条）
    ah_w = ah * w, ah_male = ah * male, ah_famhx = ah * famhx,
    ah_htn = ah * htn, ah_inactive = ah * inactive,
    male_w = male * w
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
# 第 5 部分：方程与切点（开发期，与 segmented_v5.R / denoise_features_v5.R 相同）
# =====================================================================
dd    <- full[full$dev, ]
k_old <- dd$age >= SPLIT_AGE
target_ref  <- wmean(dd$ref_bang, dd$wt)
target_sens <- wmean(dd$ref_bang[k_old & dd$case == 1], dd$wt[k_old & dd$case == 1])

loco_p <- function(d, vars) {
  lp <- rep(NA_real_, nrow(d))
  for (cy in DEV_CYCLES) {
    tr <- d$cycle != cy; te <- d$cycle == cy
    lp[te] <- lp_of(fit_l2(as.matrix(d[tr, vars]), d$case[tr], d$wt[tr]), d[te, ])
  }
  plogis(lp)
}
cut_for_sens <- function(p, y, w, target) {
  u <- sort(unique(p[y == 1]), decreasing = TRUE)
  sens <- sapply(u, function(t) sum(w[y == 1 & p >= t]) / sum(w[y == 1]))
  u[which(sens >= target - 1e-12)[1]]
}

MODELS <- list(V5 = VARS, F4 = VARS4)
old_all <- full$age >= SPLIT_AGE
THR <- list()
for (m in names(MODELS)) {
  v <- MODELS[[m]]
  p_dev <- loco_p(dd, v)
  t0 <- v5_threshold(p_dev, dd$wt, target_ref)
  t1 <- min(cut_for_sens(p_dev[k_old], dd$case[k_old], dd$wt[k_old], target_sens), t0)
  THR[[m]] <- c(single = t0, old = t1)
  fit <- fit_l2(as.matrix(dd[, v]), dd$case, dd$wt)
  pm <- paste0("p_", m); full[[pm]] <- NA_real_
  full[[pm]][full$dev] <- p_dev
  if (any(full$eval)) full[[pm]][full$eval] <- plogis(lp_of(fit, full[full$eval, ]))
  full[[paste0("ref_", m, "_S0")]] <- as.numeric(full[[pm]] >= t0)
  full[[paste0("ref_", m, "_S1")]] <- as.numeric(full[[pm]] >= ifelse(old_all, t1, t0))
}
full$ref_V5_S4 <- as.numeric(ifelse(old_all, 1, full$ref_V5_S0))
RULES <- c("Bang", "V5_S0", "V5_S1", "V5_S4", "F4_S0", "F4_S1")
refcol <- function(r) if (r == "Bang") "ref_bang" else paste0("ref_", r)
for (r in RULES[-1]) full[[paste0("dB_", r)]] <- full[[refcol(r)]] - full$ref_bang

lk <- list.files(file.path(PROJECT_DIR, "results"), pattern = "^locked_main\\.csv$",
                 recursive = TRUE, full.names = TRUE)
if (length(lk) > 0) {
  lk <- read.csv(lk[which.max(file.mtime(lk))]); vv <- setNames(lk$value, lk$term)
  cat("\nV5 切点与 locked_main.csv 的差：", signif(abs(THR$V5["single"] - vv["THRESHOLD"]), 3), "（应接近 0）\n")
}


# =====================================================================
# 第 6 部分：按可及性分组的指标与差距（Taylor 区间）
# =====================================================================
full$acc_high <- as.numeric(full$access == "high")
des <- svydesign(ids = ~psu, strata = ~strata, weights = ~wt, nest = TRUE, data = full)
periods <- c(dev = "dev", eval = "eval")
if (!any(full$eval)) periods <- periods["dev"]

ci_row <- function(est, se) 100 * c(estimate = est, lower95 = est - 1.96 * se, upper95 = est + 1.96 * se)

# 组内均值（变量 var，在子集 idx 内按可及性分组）
by_access <- function(var, idx) {
  b <- svyby(as.formula(paste0("~", var)), ~access, des[idx, ], svymean, na.rm = TRUE, covmat = TRUE)
  list(b = b, est = setNames(coef(b), rownames(b)))
}
# 高 − 低的差
gap_of <- function(var, idx) {
  ba <- by_access(var, idx)
  ct <- svycontrast(ba$b, setNames(c(1, -1), c("high", "low")))
  ci_row(unname(coef(ct)), unname(SE(ct)))
}
# 年龄标准化差距：按该时期全部病例的年龄构成加权
gap_std <- function(var, idx_cases) {
  b <- svyby(as.formula(paste0("~", var)), ~interaction(access, age_grp), des[idx_cases, ],
             svymean, na.rm = TRUE, covmat = TRUE)
  agew <- prop.table(tapply(full$wt[idx_cases], full$age_grp[idx_cases], sum))
  nm <- rownames(b); cf <- setNames(rep(0, length(nm)), nm)
  for (g in names(agew)) {
    cf[paste("high", g, sep = ".")] <- agew[g]
    cf[paste("low", g, sep = ".")]  <- -agew[g]
  }
  ct <- svycontrast(b, cf)
  ci_row(unname(coef(ct)), unname(SE(ct)))
}

grp_rows <- list(); gap_rows <- list()
for (pe in periods) {
  base <- full[[pe]] & !is.na(full$access)
  cases <- base & full$case == 1; nonc <- base & full$case == 0
  for (r in RULES) {
    v <- refcol(r)
    for (acc in c("high", "low")) {
      ka <- full$access %in% acc
      add <- function(measure, x) grp_rows[[length(grp_rows) + 1]] <<-
        data.frame(period = pe, strategy = r, access = acc, measure = measure,
                   estimate = x[1], lower95 = x[2], upper95 = x[3], row.names = NULL)
      tay <- function(var, idx) {
        m <- svymean(as.formula(paste0("~", var)), des[idx, ], na.rm = TRUE)
        ci <- confint(m); 100 * c(unname(coef(m)), ci[1], ci[2])
      }
      add("referral",    tay(v, base & ka))
      add("sensitivity", tay(v, cases & ka))
      add("specificity", 100 - tay(v, nonc & ka)[c(1, 3, 2)])
      add("PPV",         tay("case", base & ka & full[[v]] == 1))
    }
    addg <- function(measure, x) gap_rows[[length(gap_rows) + 1]] <<-
      data.frame(period = pe, strategy = r, measure = measure,
                 estimate = x[1], lower95 = x[2], upper95 = x[3], row.names = NULL)
    addg("sens_gap_high_minus_low", gap_of(v, cases))
    addg("sens_gap_age_standardised", gap_std(v, cases))
    addg("referral_gap_high_minus_low", gap_of(v, base))
    if (r != "Bang") {
      addg("change_in_sens_gap_vs_Bang", gap_of(paste0("dB_", r), cases))
      addg("change_in_std_gap_vs_Bang",  gap_std(paste0("dB_", r), cases))
    }
  }
}
grp_tab <- bind_rows(grp_rows)
gap_tab <- bind_rows(gap_rows)

# 病例年龄构成（解释差距的来源）
age_mix <- full %>% filter((dev | eval) & case == 1 & !is.na(access)) %>%
  mutate(period = ifelse(dev, "dev", "eval")) %>%
  group_by(period, access, age_grp) %>% summarise(w = sum(wt), n = n(), .groups = "drop") %>%
  group_by(period, access) %>% mutate(pct = 100 * w / sum(w)) %>% ungroup() %>% select(-w)


# =====================================================================
# 第 7 部分：同名额曲线（评价期；方程和切点固定；PSU bootstrap）
# 总转检比例固定为 b 时，高 / 低可及组各自的敏感度和差距
# =====================================================================
fam_access <- function(p, y, w, old, acc, t_old, budgets) {
  W <- sum(w)
  if (is.na(t_old)) {
    o <- order(-p); x <- c(0, cumsum(w[o])) / W
    base_ref <- rep(FALSE, length(p)); ord <- o
  } else {
    base_ref <- old & p >= t_old
    yo <- which(!old); ord <- yo[order(-p[yo])]
    x <- c(sum(w[base_ref]), sum(w[base_ref]) + cumsum(w[ord])) / W
  }
  sens_acc <- function(a) {
    k <- acc == a; Wc <- sum((w * y)[k])
    s <- c(sum((w * y)[k & base_ref]), sum((w * y)[k & base_ref]) + cumsum((w * y * k)[ord])) / Wc
    approx(x, s, xout = budgets, ties = "ordered", rule = 1)$y
  }
  cbind(high = sens_acc("high"), low = sens_acc("low"))
}
CURVES <- list(V5_single = list(m = "V5", t = NA), V5_oldfix = list(m = "V5", t = THR$V5["old"]),
               F4_single = list(m = "F4", t = NA), F4_oldfix = list(m = "F4", t = THR$F4["old"]))
curves_for <- function(d, w) {
  out <- list()
  for (cn in names(CURVES)) {
    cc <- fam_access(d[[paste0("p_", CURVES[[cn]]$m)]], d$case, w, d$age >= SPLIT_AGE,
                     d$access, CURVES[[cn]]$t, BUDGETS)
    out[[paste0(cn, ".high")]] <- cc[, "high"]; out[[paste0(cn, ".low")]] <- cc[, "low"]
    out[[paste0(cn, ".gap")]]  <- cc[, "high"] - cc[, "low"]
  }
  do.call(cbind, out)
}
boot_weights <- function(d) {
  key <- paste(d$strata, d$psu); f <- setNames(numeric(0), character(0))
  for (h in unique(d$strata)) {
    ps <- unique(d$psu[d$strata == h]); nh <- length(ps)
    if (nh < 2) { f[paste(h, ps)] <- 1; next }
    cnt <- table(factor(sample(ps, nh - 1, replace = TRUE), levels = ps))
    f[paste(h, ps)] <- as.numeric(cnt) * nh / (nh - 1)
  }
  .audit_last_f <<- f
  d$wt * f[key]
}
qna <- function(x, pr) if (mean(is.na(x)) > 0.05) NA_real_ else unname(quantile(x, pr, na.rm = TRUE))

pe_c <- if ("eval" %in% periods) "eval" else "dev"
dC <- full[full[[pe_c]] & !is.na(full$access), ]
est <- curves_for(dC, dC$wt)
bt <- array(NA_real_, c(nrow(est), ncol(est), B_BOOT), dimnames = list(NULL, colnames(est), NULL))
cat("bootstrap（B =", B_BOOT, "）...\n")
audit_status <- vector("list", B_BOOT)
audit_multipliers <- vector("list", B_BOOT)
audit_rng_states <- vector("list", B_BOOT)
for (b in seq_len(B_BOOT)) {
  audit_rng_states[[b]] <- .Random.seed
  .audit_start <- proc.time()[3]
  .audit_weights <- boot_weights(dC)
  bt[, , b] <- curves_for(dC, .audit_weights)
  audit_multipliers[[b]] <- data.frame(replicate=b, stratum_psu=names(.audit_last_f), multiplier=as.numeric(.audit_last_f))
  audit_status[[b]] <- data.frame(replicate=b, status="completed", finite_cells=sum(is.finite(bt[,,b])), missing_cells=sum(is.na(bt[,,b])), total_cells=length(bt[,,b]), elapsed_seconds=unname(proc.time()[3]-.audit_start))
  if (b %% 50 == 0) cat("Audit completed", b, "of", B_BOOT, "replicates\n")
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

curve_tab <- bind_rows(lapply(colnames(est), function(j) data.frame(
  period = pe_c, curve = j, budget = BUDGETS, value = est[, j],
  lower95 = apply(bt[, j, ], 1, qna, pr = 0.025), upper95 = apply(bt[, j, ], 1, qna, pr = 0.975))))

# 同名额下各方案差距之差（例如 V5_oldfix 的差距 − V5_single 的差距）
GCMP <- list(`V5 oldfix - V5 single` = c("V5_oldfix.gap", "V5_single.gap"),
             `F4 single - V5 single` = c("F4_single.gap", "V5_single.gap"),
             `F4 oldfix - V5 oldfix` = c("F4_oldfix.gap", "V5_oldfix.gap"))
key_budgets <- c(Bang = wmean(dC$ref_bang, dC$wt), V5_S1 = wmean(dC$ref_V5_S1, dC$wt))
ip <- function(v, b0) approx(BUDGETS, v, xout = b0)$y
key_rows <- list()
for (kb in names(key_budgets)) {
  b0 <- key_budgets[kb]
  for (cn in names(CURVES)) {
    j <- paste0(cn, ".gap")
    db <- apply(bt[, j, ], 2, ip, b0 = b0)
    key_rows[[length(key_rows) + 1]] <- data.frame(at = kb, budget = b0, item = paste(cn, "gap"),
                                                   value = ip(est[, j], b0), lower95 = qna(db, 0.025), upper95 = qna(db, 0.975), row.names = NULL)
  }
  for (nm in names(GCMP)) {
    a <- GCMP[[nm]]
    db <- apply(bt[, a[1], ] - bt[, a[2], ], 2, ip, b0 = b0)
    key_rows[[length(key_rows) + 1]] <- data.frame(at = kb, budget = b0, item = nm,
                                                   value = ip(est[, a[1]] - est[, a[2]], b0), lower95 = qna(db, 0.025), upper95 = qna(db, 0.975), row.names = NULL)
  }
}
key_tab <- bind_rows(key_rows)


# =====================================================================
# 第 8 部分：图
# =====================================================================
band <- function(x, lo, hi, col) {
  k <- !is.na(lo) & !is.na(hi)
  if (any(k)) polygon(c(x[k], rev(x[k])), c(lo[k], rev(hi[k])), col = adjustcolor(col, 0.2), border = NA)
}
png(file.path(OUT_DIR, "access_gap_curves.png"), width = 2600, height = 1100, res = 200)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
for (cn in c("V5_single", "V5_oldfix")) {
  plot(NA, xlim = 100 * range(BUDGETS), ylim = c(0, 100), xlab = "Overall referral (%)",
       ylab = "Sensitivity (%)", main = paste0(cn, ", ", pe_c, ": high vs low access"))
  for (a in c("high", "low")) {
    r <- curve_tab[curve_tab$curve == paste(cn, a, sep = "."), ]
    col <- if (a == "high") "#2166ac" else "#b2182b"
    band(100 * r$budget, 100 * r$lower95, 100 * r$upper95, col)
    lines(100 * r$budget, 100 * r$value, lwd = 2, col = col)
  }
  for (a in c("high", "low")) {
    g <- grp_tab[grp_tab$period == pe_c & grp_tab$strategy == "Bang" & grp_tab$access == a, ]
    ov <- wmean(dC$ref_bang, dC$wt) * 100
    points(ov, g$estimate[g$measure == "sensitivity"], pch = 19, col = if (a == "high") "#2166ac" else "#b2182b")
  }
  legend("bottomright", c("high access", "low access", "Bang (dots)"), col = c("#2166ac", "#b2182b", "grey30"),
         lty = c(1, 1, NA), pch = c(NA, NA, 19), lwd = 2, bty = "n", cex = 0.8)
}
invisible(dev.off())

png(file.path(OUT_DIR, "access_gap_difference.png"), width = 1400, height = 1100, res = 200)
par(mar = c(4.5, 4.5, 3, 1))
cols <- c(V5_single = "#2166ac", V5_oldfix = "#b2182b", F4_single = "#1b7837", F4_oldfix = "#d95f02")
yl <- range(100 * curve_tab$value[grepl("\\.gap$", curve_tab$curve)], na.rm = TRUE) + c(-10, 10)
plot(NA, xlim = 100 * range(BUDGETS), ylim = yl, xlab = "Overall referral (%)",
     ylab = "Sensitivity gap, high - low (pp)", main = paste0("Access gap at matched budget (", pe_c, ")"))
abline(h = 0, lty = 2)
for (cn in c("V5_single", "V5_oldfix")) {
  r <- curve_tab[curve_tab$curve == paste0(cn, ".gap"), ]
  band(100 * r$budget, 100 * r$lower95, 100 * r$upper95, cols[cn])
}
for (cn in names(cols)) {
  r <- curve_tab[curve_tab$curve == paste0(cn, ".gap"), ]
  lines(100 * r$budget, 100 * r$value, lwd = 2, col = cols[cn], lty = ifelse(grepl("^F4", cn), 2, 1))
}
gb <- gap_tab[gap_tab$period == pe_c & gap_tab$strategy == "Bang" & gap_tab$measure == "sens_gap_high_minus_low", ]
points(100 * wmean(dC$ref_bang, dC$wt), gb$estimate, pch = 19, col = "grey30")
legend("topright", c(names(cols), "Bang"), col = c(cols, "grey30"), lty = c(1, 1, 2, 2, NA),
       pch = c(NA, NA, NA, NA, 19), lwd = 2, bty = "n", cex = 0.75)
invisible(dev.off())


# =====================================================================
# 第 9 部分：保存和打印
# =====================================================================
write.csv(grp_tab,   file.path(OUT_DIR, "access_group_metrics.csv"), row.names = FALSE)
write.csv(gap_tab,   file.path(OUT_DIR, "access_gaps.csv"), row.names = FALSE)
write.csv(age_mix,   file.path(OUT_DIR, "case_age_mix_by_access.csv"), row.names = FALSE)
write.csv(curve_tab, file.path(OUT_DIR, "access_gap_curves.csv"), row.names = FALSE)
write.csv(key_tab,   file.path(OUT_DIR, "access_gap_at_key_budgets.csv"), row.names = FALSE)
write.csv(data.frame(model = names(THR), cut_single = sapply(THR, `[`, "single"),
                     cut_60plus = sapply(THR, `[`, "old")), file.path(OUT_DIR, "thresholds.csv"), row.names = FALSE)
writeLines(c(paste("bootstrap 次数：", B_BOOT),
             "方程和切点只用 2007–2012 确定；2013–2018 为事后再评价，不是独立验证。",
             paste("运行时间：", Sys.time())), file.path(OUT_DIR, "access_gap_info.txt"))
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))

fmt <- function(d) sprintf("%.1f (%.1f to %.1f)", d$estimate, d$lower95, d$upper95)
for (pe in periods) {
  cat("\n=====", ifelse(pe == "dev", "开发期 2007–2012（留一周期）", "2013–2018"), "=====\n")
  g <- grp_tab[grp_tab$period == pe & grp_tab$measure %in% c("referral", "sensitivity", "PPV"), ]
  g$value <- fmt(g); g$key <- paste(g$measure, g$access)
  w_ <- reshape(g[, c("strategy", "key", "value")], idvar = "strategy", timevar = "key", direction = "wide")
  names(w_) <- sub("^value\\.", "", names(w_)); print(w_, row.names = FALSE)
  cat("-- 差距（高 − 低，百分点）--\n")
  x <- gap_tab[gap_tab$period == pe, ]; x$value <- fmt(x)
  w_ <- reshape(x[, c("strategy", "measure", "value")], idvar = "strategy", timevar = "measure", direction = "wide")
  names(w_) <- sub("^value\\.", "", names(w_)); print(w_, row.names = FALSE)
}
cat("\n===== 病例年龄构成（%）=====\n"); print(as.data.frame(age_mix), digits = 3, row.names = FALSE)
cat("\n===== 同名额下的差距（", pe_c, "，百分点，bootstrap 95% 区间）=====\n")
key_tab$value <- sprintf("%+.1f (%+.1f to %+.1f)", 100 * key_tab$value, 100 * key_tab$lower95, 100 * key_tab$upper95)
print(key_tab[, c("at", "budget", "item", "value")], digits = 3, row.names = FALSE)
cat("\n结果已保存到：", OUT_DIR, "\n")
