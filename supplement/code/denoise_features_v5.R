# =====================================================================
# 去噪曲线 + 冗余检查 + 老年针对性特征工程（denoise_features_v5.R）
# ---------------------------------------------------------------------
# A. 去噪：不再在每个名额下重新挑“最优分法”（会在几乎等价的切点间来回跳），
#    改为画三条固定方案族曲线，在评价期按实际转检比例索引（类似 ROC，没有选择步骤）：
#      single  所有人同一个切点，改变名额            （S0 这一族）
#      oldfix  老年切点固定为该方程的 S1 老年切点，只调青年切点（S1 这一族）
#      old_all 老年全部转检，只调青年切点            （S4 这一族）
#    并按 PSU 做 Rao–Wu 重标度 bootstrap（方程和切点固定），给出 95% 置信带，
#    以及“同名额下两条曲线之差”的置信带。抽样噪声无法消除，只能这样量化。
#
# B. 冗余检查（只用开发期）：加权相关矩阵、VIF、逐个删除 V5 六项变量
#    （留一周期加权对数损失 LL 与 AUC），去掉删了也不变差的变量。
#
# C. 针对性特征（只用开发期，前向选择，规则事先写死）：
#    核心池（不涉及种族、不需要以前查过血）：
#      s 年龄样条、ah 60 岁后斜率变化、ah×（w/男性/家族史/高血压/不活动）、
#      男性×w、aw 年龄×腰高比、BMI、身高、自报高胆固醇
#    扩展池：种族/族裔、自报被告知糖尿病前期 / 有风险（后两项需要以前查过血，
#            且与“未诊断”人群定义有关，只作敏感性分析）
#    选择用平滑指标（留一周期加权对数损失 LL），而不是某个名额下的敏感度，
#    避免选择本身被噪声带偏。
#    加入规则：LL 改善 > TOL_LL，且总体 AUC 提升 ≥ TOL_AUC 或老年组 AUC 提升 ≥ TOL_AUC_OLD。
#
# D. 新方程（F_core、F_ext）沿用 S0 / S1 切点规则（开发期确定），
#    与 V5 在 2013–2018 上比较一次（Taylor 区间 + 曲线 bootstrap）。
#
# 注意：任何新特征都意味着方程不再是 10-28 锁定的 V5。2013–2018 已多次看过，
#       这里的结果只能称为事后探索，不能替代外部验证；论文里应作为探索性分析报告。
#
# 用法：放到 scripts 文件夹，只改第 0 部分，点 Source。约 5–15 分钟（bootstrap 占大部分）。
# 输出：results\denoise_features_v5_日期\
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
PROJECT_DIR <- "D:/Nottingham/Diabetes_Screening_NHANES"
DATA_DIR    <- file.path(PROJECT_DIR, "data", "raw")
OUT_DIR     <- file.path(PROJECT_DIR, "results",
                         paste0("denoise_features_v5_", format(Sys.Date(), "%Y%m%d")))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

DEV_CYCLES  <- c("E", "F", "G")   # 2007–2012 开发
EVAL_CYCLES <- c("H", "I", "J")   # 2013–2018 再评价
C_L2      <- 1
KNOTS     <- c(2.5, 4.6, 7.3)     # 年龄样条结点（与 V5 附录 B 一致）
SPLIT_AGE <- 60
VARS      <- c("a", "w", "male", "famhx", "htn", "inactive")   # V5 六项（已锁定）

CORE_POOL <- c("s", "ah", "ah_w", "ah_male", "ah_famhx", "ah_htn", "ah_inactive",
               "male_w", "aw", "bmi_c", "h", "chol")
EXT_POOL  <- c("hisp", "black", "othrace", "prediab", "atrisk")

TOL_LL      <- 0.0005   # 每人加权对数损失至少改善这么多
TOL_AUC     <- 0.002    # 总体 AUC 至少提升
TOL_AUC_OLD <- 0.005    # 或老年组 AUC 至少提升
B_BOOT      <- 500      # bootstrap 次数（先用 100 试跑）
BUDGETS     <- seq(0.25, 0.50, by = 0.005)   # 曲线横轴：总转检比例


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
# 第 5 部分：开发期工具（留一周期预测、评分、切点）
# =====================================================================
dd    <- full[full$dev, ]
k_old <- dd$age >= SPLIT_AGE
target_ref  <- wmean(dd$ref_bang, dd$wt)
target_sens <- wmean(dd$ref_bang[k_old & dd$case == 1], dd$wt[k_old & dd$case == 1])

# 候选变量可用性：开发期或评价期某个周期的原始问题整列缺失，就从候选池中去掉
src_of <- c(chol = "BPQ080", prediab = "DIQ160", atrisk = "DIQ170",
            hisp = "RIDRETH1", black = "RIDRETH1", othrace = "RIDRETH1")
avail <- sapply(c(CORE_POOL, EXT_POOL), function(v) {
  if (!v %in% names(src_of)) return(TRUE)
  x <- full[full$in_sample, ]
  all(tapply(!is.na(x[[src_of[v]]]), x$cycle, mean) > 0.5)
})
if (any(!avail)) cat(">>> 以下候选变量在某些周期缺失过多，已排除：", names(avail)[!avail], "\n")
CORE_POOL <- CORE_POOL[avail[CORE_POOL]]
EXT_POOL  <- EXT_POOL[avail[EXT_POOL]]

loco_p <- function(d, vars) {
  lp <- rep(NA_real_, nrow(d))
  for (cy in DEV_CYCLES) {
    tr <- d$cycle != cy; te <- d$cycle == cy
    lp[te] <- lp_of(fit_l2(as.matrix(d[tr, vars]), d$case[tr], d$wt[tr]), d[te, ])
  }
  plogis(lp)
}
wll <- function(p, y, w) { p <- pmin(pmax(p, 1e-9), 1 - 1e-9); -sum(w * (y * log(p) + (1 - y) * log(1 - p))) / sum(w) }
wsens <- function(ref, y, w) sum(w * ref * y) / sum(w * y)

score_cache <- new.env()
score_set <- function(vars) {
  key <- paste(sort(vars), collapse = "+")
  if (!is.null(score_cache[[key]])) return(score_cache[[key]])
  p <- loco_p(dd, vars); y <- dd$case; w <- dd$wt
  thr <- v5_threshold(p, w, target_ref)
  out <- c(LL = wll(p, y, w), AUC = w_auc(p, y, w),
           AUC_old = w_auc(p[k_old], y[k_old], w[k_old]),
           AUC_young = w_auc(p[!k_old], y[!k_old], w[!k_old]),
           sens_at_Bang_budget = wsens(as.numeric(p >= thr), y, w))
  score_cache[[key]] <- out
  out
}
cut_for_sens <- function(p, y, w, target) {
  u <- sort(unique(p[y == 1]), decreasing = TRUE)
  sens <- sapply(u, function(t) sum(w[y == 1 & p >= t]) / sum(w[y == 1]))
  u[which(sens >= target - 1e-12)[1]]
}


# =====================================================================
# 第 6 部分：冗余检查（开发期）
# =====================================================================
all_feats <- unique(c(VARS, CORE_POOL, EXT_POOL))
cor_mat <- cov.wt(as.matrix(dd[, c("age", "waist", "height", "bmi", "whtr", all_feats)]),
                  wt = dd$wt / sum(dd$wt), cor = TRUE)$cor
vif_of <- function(v) sapply(v, function(x) {
  r2 <- summary(lm(reformulate(setdiff(v, x), x), data = dd, weights = wt))$r.squared
  1 / (1 - r2)
})
vif_v5 <- vif_of(VARS)

# 逐个删除（贪婪后退）：删掉后 LL 不变差、AUC 下降 < TOL_AUC 的变量视为冗余
back_steps <- list(); cur <- VARS
repeat {
  sc <- score_set(cur)
  rows <- bind_rows(lapply(cur, function(v) {
    s2 <- score_set(setdiff(cur, v))
    data.frame(drop = v, dLL = s2["LL"] - sc["LL"], dAUC = s2["AUC"] - sc["AUC"],
               dAUC_old = s2["AUC_old"] - sc["AUC_old"],
               dSens_at_Bang = s2["sens_at_Bang_budget"] - sc["sens_at_Bang_budget"], row.names = NULL)
  }))
  rows$redundant <- rows$dLL <= 0 & rows$dAUC > -TOL_AUC
  back_steps[[length(back_steps) + 1]] <- cbind(step = length(back_steps) + 1,
                                                current = paste(cur, collapse = "+"), rows)
  if (!any(rows$redundant) || length(cur) <= 2) break
  cur <- setdiff(cur, rows$drop[rows$redundant][which.min(rows$dLL[rows$redundant])])
}
back_tab <- bind_rows(back_steps)
BASE_RED <- cur
cat("\n冗余检查后保留：", paste(BASE_RED, collapse = ", "), "\n")


# =====================================================================
# 第 7 部分：前向选择（开发期，规则事先写死）
# =====================================================================
forward <- function(start, pool, label) {
  cur <- start; steps <- list()
  repeat {
    sc <- score_set(cur); cand <- setdiff(pool, cur)
    if (length(cand) == 0) break
    rows <- bind_rows(lapply(cand, function(v) {
      s2 <- score_set(c(cur, v))
      data.frame(add = v, gain_LL = sc["LL"] - s2["LL"], dAUC = s2["AUC"] - sc["AUC"],
                 dAUC_old = s2["AUC_old"] - sc["AUC_old"],
                 dSens_at_Bang = s2["sens_at_Bang_budget"] - sc["sens_at_Bang_budget"], row.names = NULL)
    }))
    rows$qualifies <- rows$gain_LL > TOL_LL & (rows$dAUC >= TOL_AUC | rows$dAUC_old >= TOL_AUC_OLD)
    steps[[length(steps) + 1]] <- cbind(model = label, step = length(steps) + 1,
                                        current = paste(cur, collapse = "+"), rows)
    if (!any(rows$qualifies)) break
    cur <- c(cur, rows$add[rows$qualifies][which.max(rows$gain_LL[rows$qualifies])])
  }
  list(vars = cur, steps = bind_rows(steps))
}
fw_core <- forward(BASE_RED, CORE_POOL, "F_core")
fw_ext  <- forward(fw_core$vars, EXT_POOL, "F_ext")
fwd_tab <- bind_rows(fw_core$steps, fw_ext$steps)

MODELS <- list(V5 = VARS, F_core = fw_core$vars, F_ext = fw_ext$vars)
cat("F_core：", paste(MODELS$F_core, collapse = ", "), "\n")
cat("F_ext ：", paste(MODELS$F_ext, collapse = ", "), "\n")
if (setequal(MODELS$F_core, VARS)) cat(">>> F_core 与 V5 相同：核心池中没有变量达到加入标准，冗余检查也没有删掉变量。\n")

dev_scores <- bind_rows(lapply(names(MODELS), function(m)
  data.frame(model = m, vars = paste(MODELS[[m]], collapse = "+"), t(score_set(MODELS[[m]])))))
vif_models <- bind_rows(lapply(names(MODELS), function(m)
  data.frame(model = m, term = MODELS[[m]], VIF = vif_of(MODELS[[m]]), row.names = NULL)))


# =====================================================================
# 第 8 部分：各方程的切点（开发期）、最终系数、评价期预测
# =====================================================================
old_all <- full$age >= SPLIT_AGE
THR <- list(); COEF <- list()
for (m in names(MODELS)) {
  v <- MODELS[[m]]
  p_dev <- loco_p(dd, v)
  t0 <- v5_threshold(p_dev, dd$wt, target_ref)
  t1 <- min(cut_for_sens(p_dev[k_old], dd$case[k_old], dd$wt[k_old], target_sens), t0)
  THR[[m]] <- c(single = t0, old = t1)
  fit <- fit_l2(as.matrix(dd[, v]), dd$case, dd$wt)
  COEF[[m]] <- data.frame(model = m, term = names(fit), coef = fit, row.names = NULL)
  pm <- paste0("p_", m)
  full[[pm]] <- NA_real_
  full[[pm]][full$dev] <- p_dev
  if (any(full$eval)) full[[pm]][full$eval] <- plogis(lp_of(fit, full[full$eval, ]))
  full[[paste0("ref_", m, "_S0")]] <- as.numeric(full[[pm]] >= t0)
  full[[paste0("ref_", m, "_S1")]] <- as.numeric(full[[pm]] >= ifelse(old_all, t1, t0))
}
coef_tab <- bind_rows(COEF) %>% mutate(OR = ifelse(term == "(Intercept)", NA, exp(coef)))
thr_tab  <- data.frame(model = names(THR), cut_single = sapply(THR, `[`, "single"),
                       cut_60plus_S1 = sapply(THR, `[`, "old"), row.names = NULL)

lk <- list.files(file.path(PROJECT_DIR, "results"), pattern = "^locked_main\\.csv$",
                 recursive = TRUE, full.names = TRUE)
if (length(lk) > 0) {
  lk <- read.csv(lk[which.max(file.mtime(lk))]); vv <- setNames(lk$value, lk$term)
  b5 <- setNames(COEF$V5$coef, COEF$V5$term)
  cat("\nV5 与 locked_main.csv 的最大系数差：", signif(max(abs(b5 - vv[names(b5)])), 3),
      "；切点差：", signif(abs(THR$V5["single"] - vv["THRESHOLD"]), 3), "（应接近 0）\n")
}

RULES <- c("Bang" = "ref_bang", setNames(paste0("ref_", rep(names(MODELS), each = 2), c("_S0", "_S1")),
                                         paste0(rep(names(MODELS), each = 2), c("_S0", "_S1"))))
for (r in names(RULES)[-1]) {
  full[[paste0("dB_", r)]] <- full[[RULES[r]]] - full$ref_bang
  ref_v5 <- if (grepl("_S0$", r)) full$ref_V5_S0 else full$ref_V5_S1
  full[[paste0("dV_", r)]] <- full[[RULES[r]]] - ref_v5
}


# =====================================================================
# 第 9 部分：Taylor 区间（评价期与开发期）
# =====================================================================
des <- svydesign(ids = ~psu, strata = ~strata, weights = ~wt, nest = TRUE, data = full)
tay <- function(var, idx) {
  m <- svymean(as.formula(paste0("~", var)), des[idx, ], na.rm = TRUE)
  ci <- confint(m)
  100 * c(estimate = unname(coef(m)), lower95 = ci[1], upper95 = ci[2])
}
GROUPS <- list(overall = TRUE, `20-39` = quote(age_grp == "20-39"),
               `40-59` = quote(age_grp == "40-59"), `>=60` = quote(age_grp == ">=60"))
periods <- c(dev = "dev", eval = "eval")
if (!any(full$eval)) periods <- periods["dev"]
res <- list()
for (pe in periods) for (gname in names(GROUPS)) {
  g <- if (isTRUE(GROUPS[[gname]])) rep(TRUE, nrow(full)) else (eval(GROUPS[[gname]], full) %in% TRUE)
  allg <- full[[pe]] & g; cases <- allg & full$case == 1; nonc <- allg & full$case == 0
  add <- function(strategy, measure, v) res[[length(res) + 1]] <<-
    data.frame(period = pe, group = gname, strategy = strategy, measure = measure,
               estimate = v[1], lower95 = v[2], upper95 = v[3], row.names = NULL)
  for (r in names(RULES)) {
    add(r, "referral",    tay(RULES[r], allg))
    add(r, "sensitivity", tay(RULES[r], cases))
    add(r, "specificity", 100 - tay(RULES[r], nonc)[c(1, 3, 2)])
    if (r != "Bang") {
      add(r, "sens_diff_vs_Bang", tay(paste0("dB_", r), cases))
      if (!grepl("^V5_", r)) {
        add(r, "sens_diff_vs_V5_same_rule",  tay(paste0("dV_", r), cases))
        add(r, "refer_diff_vs_V5_same_rule", tay(paste0("dV_", r), allg))
      }
    }
  }
}
results <- bind_rows(res)
auc_tab <- bind_rows(lapply(periods, function(pe) {
  k <- full[[pe]]
  bind_rows(lapply(names(MODELS), function(m) data.frame(
    period = pe, model = m,
    AUC = w_auc(full[[paste0("p_", m)]][k], full$case[k], full$wt[k]),
    AUC_60plus = w_auc(full[[paste0("p_", m)]][k & old_all], full$case[k & old_all], full$wt[k & old_all]),
    AUC_under60 = w_auc(full[[paste0("p_", m)]][k & !old_all], full$case[k & !old_all], full$wt[k & !old_all]))))
}))


# =====================================================================
# 第 10 部分：固定方案族曲线 + PSU bootstrap 置信带（评价期；方程和切点固定）
# =====================================================================
# 按转检比例索引：在每个名额处线性插值（相当于在两个相邻切点之间随机化），没有任何选择步骤
fam_curve <- function(p, y, w, old, t_old, budgets) {
  W <- sum(w); Wc <- sum(w * y)
  if (is.na(t_old)) {
    o <- order(-p)
    x <- c(0, cumsum(w[o])) / W; s <- c(0, cumsum((w * y)[o])) / Wc
  } else {
    ro <- old & p >= t_old
    bw <- sum(w[ro]); bc <- sum((w * y)[ro])
    yo <- which(!old); o <- yo[order(-p[yo])]
    x <- c(bw, bw + cumsum(w[o])) / W; s <- c(bc, bc + cumsum((w * y)[o])) / Wc
  }
  approx(x, s, xout = budgets, ties = "ordered", rule = 1)$y
}
FAMS <- function(m) list(single = NA, oldfix = THR[[m]]["old"], old_all = 0)

curves_for <- function(d, w) {
  out <- list()
  for (m in names(MODELS)) for (f in names(FAMS(m))) {
    out[[paste(m, f, sep = ".")]] <- fam_curve(d[[paste0("p_", m)]], d$case, w, d$age >= SPLIT_AGE,
                                               FAMS(m)[[f]], BUDGETS)
  }
  do.call(cbind, out)
}

# Rao–Wu 重标度 bootstrap：每层抽 n_h − 1 个 PSU（有放回），权重乘 n_h/(n_h − 1) × 抽中次数
boot_weights <- function(d) {
  key <- paste(d$strata, d$psu); f <- setNames(numeric(0), character(0))
  for (h in unique(d$strata)) {
    ps <- unique(d$psu[d$strata == h]); nh <- length(ps)
    if (nh < 2) { f[paste(h, ps)] <- 1; next }
    cnt <- table(factor(sample(ps, nh - 1, replace = TRUE), levels = ps))
    f[paste(h, ps)] <- as.numeric(cnt) * nh / (nh - 1)
  }
  d$wt * f[key]
}

curve_period <- function(pe, B) {
  d <- full[full[[pe]], ]
  est <- curves_for(d, d$wt)
  bt <- array(NA_real_, c(length(BUDGETS), ncol(est), B), dimnames = list(NULL, colnames(est), NULL))
  for (b in seq_len(B)) bt[, , b] <- curves_for(d, boot_weights(d))
  list(est = est, boot = bt)
}
# 分位数：超过 5% 的重抽样在该名额处无定义（例如名额低于“老年全检”的起点）时不给区间
qna <- function(x, pr) if (mean(is.na(x)) > 0.05) NA_real_ else unname(quantile(x, pr, na.rm = TRUE))
CMP <- list(                     # 同名额下的差：A − B
  `V5: oldfix - single`      = c("V5.oldfix", "V5.single"),
  `V5: old_all - oldfix`     = c("V5.old_all", "V5.oldfix"),
  `F_core - V5 (single)`     = c("F_core.single", "V5.single"),
  `F_core - V5 (oldfix)`     = c("F_core.oldfix", "V5.oldfix"),
  `F_ext - V5 (oldfix)`      = c("F_ext.oldfix", "V5.oldfix"))

curve_rows <- list(); diff_rows <- list(); CV <- list()
for (pe in periods) {
  cat("bootstrap 曲线：", pe, "（B =", if (pe == "eval") B_BOOT else 0, "）...\n")
  cv <- curve_period(pe, if (pe == "eval") B_BOOT else 0); CV[[pe]] <- cv
  q <- function(a, pr) apply(a, c(1, 2), qna, pr = pr)
  lo <- if (pe == "eval") q(cv$boot, 0.025) else cv$est * NA
  hi <- if (pe == "eval") q(cv$boot, 0.975) else cv$est * NA
  for (j in colnames(cv$est)) curve_rows[[length(curve_rows) + 1]] <- data.frame(
    period = pe, curve = j, budget = BUDGETS, sens = cv$est[, j], lower95 = lo[, j], upper95 = hi[, j])
  if (pe == "eval") for (nm in names(CMP)) {
    a <- CMP[[nm]]
    dv <- cv$est[, a[1]] - cv$est[, a[2]]
    db <- cv$boot[, a[1], ] - cv$boot[, a[2], ]
    diff_rows[[length(diff_rows) + 1]] <- data.frame(
      comparison = nm, budget = BUDGETS, diff = dv,
      lower95 = apply(db, 1, qna, pr = 0.025),
      upper95 = apply(db, 1, qna, pr = 0.975))
  }
}
curve_tab <- bind_rows(curve_rows)
diff_tab  <- bind_rows(diff_rows)

# 关键名额处的差（评价期 Bang 与 V5-S1 的实际转检比例）
key_budgets <- c(Bang = wmean(full$ref_bang[full$eval], full$wt[full$eval]),
                 V5_S1 = wmean(full$ref_V5_S1[full$eval], full$wt[full$eval]))
key_tab <- NULL
if ("eval" %in% periods) {
  key_tab <- bind_rows(lapply(names(CMP), function(nm) bind_rows(lapply(names(key_budgets), function(kb) {
    a <- CMP[[nm]]; b0 <- key_budgets[kb]
    ip <- function(v) approx(BUDGETS, v, xout = b0)$y
    dv <- ip(CV$eval$est[, a[1]] - CV$eval$est[, a[2]])
    db <- apply(CV$eval$boot[, a[1], ] - CV$eval$boot[, a[2], ], 2, ip)
    data.frame(comparison = nm, at = kb, budget = b0, diff = dv,
               lower95 = qna(db, 0.025), upper95 = qna(db, 0.975),
               row.names = NULL)
  }))))
}


# =====================================================================
# 第 11 部分：图
# =====================================================================
band <- function(x, lo, hi, col) {
  k <- !is.na(lo) & !is.na(hi)
  if (any(k)) polygon(c(x[k], rev(x[k])), c(lo[k], rev(hi[k])), col = adjustcolor(col, 0.2), border = NA)
}
pe_plot <- if ("eval" %in% periods) "eval" else "dev"
ct <- curve_tab[curve_tab$period == pe_plot, ]
png(file.path(OUT_DIR, "family_curves.png"), width = 2600, height = 1100, res = 200)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
cols <- c(single = "#2166ac", oldfix = "#b2182b", old_all = "#d95f02")
for (m in c("V5", "F_core")) {
  plot(NA, xlim = 100 * range(BUDGETS), ylim = range(100 * ct$sens, na.rm = TRUE) + c(-3, 3),
       xlab = "Overall referral (%)", ylab = "Overall sensitivity (%)",
       main = paste0(m, ", ", ifelse(pe_plot == "eval", "2013-2018 (95% bootstrap band)", "2007-2012 LOCO")))
  for (f in names(cols)) {
    r <- ct[ct$curve == paste(m, f, sep = "."), ]
    band(100 * r$budget, 100 * r$lower95, 100 * r$upper95, cols[f])
    lines(100 * r$budget, 100 * r$sens, col = cols[f], lwd = 2)
  }
  rb <- results[results$period == pe_plot & results$group == "overall" & results$strategy == "Bang", ]
  points(rb$estimate[rb$measure == "referral"], rb$estimate[rb$measure == "sensitivity"], pch = 19, col = "grey40")
  legend("bottomright", c("single cut-off", "60+ cut-off fixed (S1 family)", "60+ all tested (S4 family)", "Bang"),
         col = c(cols, "grey40"), lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 19), lwd = 2, bty = "n", cex = 0.75)
}
invisible(dev.off())

if (nrow(diff_tab) > 0) {
  png(file.path(OUT_DIR, "family_differences.png"), width = 2600, height = 1500, res = 200)
  par(mfrow = c(2, 3), mar = c(4.5, 4.5, 3, 1))
  for (nm in names(CMP)) {
    r <- diff_tab[diff_tab$comparison == nm, ]
    plot(NA, xlim = 100 * range(BUDGETS), ylim = range(c(0, 100 * r$lower95, 100 * r$upper95), na.rm = TRUE),
         xlab = "Overall referral (%)", ylab = "Sensitivity difference (pp)", main = nm)
    band(100 * r$budget, 100 * r$lower95, 100 * r$upper95, "#2166ac")
    lines(100 * r$budget, 100 * r$diff, lwd = 2, col = "#2166ac"); abline(h = 0, lty = 2)
    abline(v = 100 * key_budgets, col = "grey60", lty = 3)
  }
  invisible(dev.off())
}


# =====================================================================
# 第 12 部分：保存和打印
# =====================================================================
write.csv(round(cor_mat, 3), file.path(OUT_DIR, "redundancy_correlation.csv"))
write.csv(data.frame(term = names(vif_v5), VIF = vif_v5), file.path(OUT_DIR, "redundancy_vif_v5.csv"), row.names = FALSE)
write.csv(back_tab,   file.path(OUT_DIR, "redundancy_backward.csv"), row.names = FALSE)
write.csv(fwd_tab,    file.path(OUT_DIR, "feature_forward_steps.csv"), row.names = FALSE)
write.csv(dev_scores, file.path(OUT_DIR, "dev_scores_by_model.csv"), row.names = FALSE)
write.csv(vif_models, file.path(OUT_DIR, "vif_by_model.csv"), row.names = FALSE)
write.csv(coef_tab,   file.path(OUT_DIR, "coefficients_by_model.csv"), row.names = FALSE)
write.csv(thr_tab,    file.path(OUT_DIR, "thresholds_by_model.csv"), row.names = FALSE)
write.csv(results,    file.path(OUT_DIR, "results_by_model.csv"), row.names = FALSE)
write.csv(auc_tab,    file.path(OUT_DIR, "auc_by_model.csv"), row.names = FALSE)
write.csv(curve_tab,  file.path(OUT_DIR, "family_curves.csv"), row.names = FALSE)
write.csv(diff_tab,   file.path(OUT_DIR, "family_differences.csv"), row.names = FALSE)
if (!is.null(key_tab)) write.csv(key_tab, file.path(OUT_DIR, "family_differences_at_key_budgets.csv"), row.names = FALSE)
writeLines(c(paste("开发期 Bang 总转检目标：", round(target_ref, 4)),
             paste("开发期 Bang 老年组敏感度目标：", round(target_sens, 4)),
             paste("冗余检查后保留：", paste(BASE_RED, collapse = "+")),
             paste("F_core：", paste(MODELS$F_core, collapse = "+")),
             paste("F_ext：", paste(MODELS$F_ext, collapse = "+")),
             paste("bootstrap 次数：", B_BOOT),
             "特征和切点只用 2007–2012 确定；2013–2018 为事后探索，不是独立验证。",
             paste("运行时间：", Sys.time())),
           file.path(OUT_DIR, "denoise_features_info.txt"))
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))

cat("\n===== V5 六项 VIF =====\n"); print(round(vif_v5, 2))
cat("\n===== 冗余检查（逐个删除，开发期留一周期）=====\n")
print(back_tab[, c("step", "drop", "dLL", "dAUC", "dAUC_old", "dSens_at_Bang", "redundant")], digits = 3, row.names = FALSE)
cat("\n===== 前向选择 =====\n")
print(fwd_tab[, c("model", "step", "add", "gain_LL", "dAUC", "dAUC_old", "dSens_at_Bang", "qualifies")], digits = 3, row.names = FALSE)
cat("\n===== 开发期各方程得分 =====\n"); print(dev_scores, digits = 4, row.names = FALSE)
cat("\n===== AUC =====\n"); print(auc_tab, digits = 3, row.names = FALSE)
for (pe in periods) {
  cat("\n=====", ifelse(pe == "dev", "开发期 2007–2012（留一周期）", "2013–2018"), ": 转检 / 敏感度 % =====\n")
  r <- results[results$period == pe & results$measure %in% c("referral", "sensitivity"), ]
  r$value <- sprintf("%.1f (%.1f-%.1f)", r$estimate, r$lower95, r$upper95)
  r$key <- paste(r$measure, r$group)
  out <- reshape(r[, c("strategy", "key", "value")], idvar = "strategy", timevar = "key", direction = "wide")
  names(out) <- sub("^value\\.", "", names(out)); print(out, row.names = FALSE)
}
if (!is.null(key_tab)) {
  cat("\n===== 同名额下的敏感度差（评价期，百分点，bootstrap 95% 区间）=====\n")
  key_tab$value <- sprintf("%+.1f (%+.1f to %+.1f)", 100 * key_tab$diff, 100 * key_tab$lower95, 100 * key_tab$upper95)
  print(key_tab[, c("comparison", "at", "budget", "value")], digits = 3, row.names = FALSE)
}
cat("\n结果已保存到：", OUT_DIR, "\n")
