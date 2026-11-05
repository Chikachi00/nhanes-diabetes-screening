# =====================================================================
# 固定总转检名额下，联合优化 <60 岁和 ≥60 岁两个切点（joint_cutoff_v5.R）
# ---------------------------------------------------------------------
# 要回答的问题：
#   S0 与 Bang 名额相同，但老年组敏感度输给 Bang；
#   S1 老年组与 Bang 持平，但总转检多了约 5.6 个百分点。
#   在“名额相同”的条件下，能不能同时做到老年组不输 Bang？
#   如果不能，老年组持平最少要多花多少名额？
#
# 方程：已锁定的 V5 单一方程（六项、线性年龄、L2 C = 1），完全不改。
# 只改切点：<60 岁一个切点 t_young，≥60 岁一个切点 t_old。
# 所有切点只在 2007–2012 开发期的留一周期预测上确定，锁定后在 2013–2018 上跑一次。
#
# 比较的策略：
#   Bang、S0、S1、S4                    （与 segmented_v5.R 定义相同，作参照）
#   J1 Bang 名额 + 总敏感度最高          （不加约束，看名额怎么分最好）
#   J2 Bang 名额 + 老年组敏感度 ≥ Bang + 总敏感度最高（主要问题；可能无解）
#   J3 S1 名额 + 总敏感度最高            （检验 S1 多出来的名额是否分配得最好）
#   J4 最小名额：老年组敏感度 ≥ Bang 且总敏感度 ≥ S0 （“老年组持平最少要多花多少名额”）
# 说明：若要求 <60 岁敏感度也不低于 S0，最小名额的解按构造就是 S1，所以这里
#       J4 只要求总敏感度不低于 S0，允许青年组让出一点。
# 另外输出“名额–最高总敏感度”曲线（有 / 无老年组约束），S0、S1、S4、Bang 标在图上。
#
# 注意：2013–2018 已多次看过，这里的结果只能称为事后再评价，不是独立验证。
#       J1–J4 在开发期上对两个切点做了搜索，开发期成绩略偏乐观，以 2013–2018 为准。
# 区间：调查设计 Taylor 线性化（权重、层、PSU），切点视为固定。
#
# 用法：放到 scripts 文件夹，只改第 0 部分，点 Source。约 1–3 分钟。
# 输出：results\joint_cutoff_v5_日期\
# =====================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(survey)
})
options(survey.lonely.psu = "adjust")
set.seed(20261102)


# =====================================================================
# 第 0 部分：设置（只改这里）
# =====================================================================
PROJECT_DIR <- "D:/Nottingham/Diabetes_Screening_NHANES"
DATA_DIR    <- file.path(PROJECT_DIR, "data", "raw")
OUT_DIR     <- file.path(PROJECT_DIR, "results",
                         paste0("joint_cutoff_v5_", format(Sys.Date(), "%Y%m%d")))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

DEV_CYCLES  <- c("E", "F", "G")   # 2007–2012 开发
EVAL_CYCLES <- c("H", "I", "J")   # 2013–2018 再评价
C_L2      <- 1
KNOTS     <- c(2.5, 4.6, 7.3)     # 只为与原脚本一致，线性方程不用样条
SPLIT_AGE <- 60                   # 两个切点的年龄界线
VARS      <- c("a", "w", "male", "famhx", "htn", "inactive")   # V5 六项，线性年龄
BUDGET_GRID <- seq(0.28, 0.46, by = 0.0025)   # 名额曲线扫描范围（总转检比例）


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
# 第 5 部分：单一方程的留一周期预测与参照切点（与 segmented_v5.R 相同）
# =====================================================================
loco_p <- function(d) {
  lp <- rep(NA_real_, nrow(d))
  for (cy in DEV_CYCLES) {
    tr <- d$cycle != cy; te <- d$cycle == cy
    lp[te] <- lp_of(fit_l2(as.matrix(d[tr, VARS]), d$case[tr], d$wt[tr]), d[te, ])
  }
  plogis(lp)
}
cut_for_sens <- function(p, y, w, target) {   # 敏感度不低于目标的最高阈值
  u <- sort(unique(p[y == 1]), decreasing = TRUE)
  sens <- sapply(u, function(t) sum(w[y == 1 & p >= t]) / sum(w[y == 1]))
  u[which(sens >= target - 1e-12)[1]]
}

dd    <- full[full$dev, ]
k_old <- dd$age >= SPLIT_AGE
target_ref  <- wmean(dd$ref_bang, dd$wt)                                                # Bang 总转检（开发期）
target_sens <- wmean(dd$ref_bang[k_old & dd$case == 1], dd$wt[k_old & dd$case == 1])   # Bang 老年敏感度（开发期）

p_dev <- loco_p(dd)
thr_single     <- v5_threshold(p_dev, dd$wt, target_ref)
thr_single_old <- min(cut_for_sens(p_dev[k_old], dd$case[k_old], dd$wt[k_old], target_sens), thr_single)

wsens <- function(ref, y, w) sum(w * ref * y) / sum(w * y)
ref_S0_dev <- as.numeric(p_dev >= thr_single)
ref_S1_dev <- as.numeric(p_dev >= ifelse(k_old, thr_single_old, thr_single))
S0_sens_dev    <- wsens(ref_S0_dev, dd$case, dd$wt)
S1_budget_dev  <- wmean(ref_S1_dev, dd$wt)


# =====================================================================
# 第 6 部分：联合切点搜索（只用开发期）
# ---------------------------------------------------------------------
# 对每一个候选老年切点（老年组所有不同的预测值，加上“老年不转检”），
# 用剩余名额确定青年切点（使总转检最接近目标；偏差相同取较高阈值），
# 记录总敏感度、老年 / 青年敏感度和转检比例。这样就枚举了该名额下
# 所有可能的“老年–青年”名额分配。
# =====================================================================
cum_tab <- function(p, y, w) {   # 阈值从高到低：每个阈值下的累计转检权重和病例权重
  u <- sort(unique(p), decreasing = TRUE)
  g <- match(p, u)
  list(thr  = c(Inf, u),
       refw = c(0, cumsum(rowsum(w, g)[, 1])),
       casw = c(0, cumsum(rowsum(w * y, g)[, 1])))
}

joint_grid <- function(p, y, w, old, budget) {
  W  <- sum(w); Wc <- sum(w * y)
  To <- cum_tab(p[old], y[old], w[old])
  Ty <- cum_tab(p[!old], y[!old], w[!old])
  rem <- budget * W - To$refw
  j <- findInterval(rem, Ty$refw)              # Ty$refw[j] <= rem
  j <- pmax(j, 1)
  j2 <- pmin(j + 1, length(Ty$refw))
  pick_hi <- abs(Ty$refw[j2] - rem) < abs(Ty$refw[j] - rem) - 1e-12
  j <- ifelse(pick_hi, j2, j)
  data.frame(t_old = To$thr, t_young = Ty$thr[j],
             ref_total  = (To$refw + Ty$refw[j]) / W,
             sens_total = (To$casw + Ty$casw[j]) / Wc,
             sens_old   = To$casw / sum(w[old] * y[old]),
             sens_young = Ty$casw[j] / sum(w[!old] * y[!old]),
             ref_old    = To$refw / sum(w[old]),
             ref_young  = Ty$refw[j] / sum(w[!old]))
}

# 在网格中按规则选一个方案：满足约束者中总敏感度最高；
# 并列时取老年切点最接近 S0 切点者（改动最小）
pick <- function(G, ok = rep(TRUE, nrow(G))) {
  G <- G[ok, ]
  if (nrow(G) == 0) return(NULL)
  G <- G[G$sens_total >= max(G$sens_total) - 1e-12, ]
  G[which.min(abs(log(pmax(G$t_old, 1e-12)) - log(thr_single))), ]
}

y_dev <- dd$case; w_dev <- dd$wt
G_bang <- joint_grid(p_dev, y_dev, w_dev, k_old, target_ref)
G_S1   <- joint_grid(p_dev, y_dev, w_dev, k_old, S1_budget_dev)

J1 <- pick(G_bang)
J2 <- pick(G_bang, G_bang$sens_old >= target_sens - 1e-12)
J3 <- pick(G_S1)

# 名额曲线：每个名额下的最高总敏感度（无约束 / 老年组 ≥ Bang）
curve <- bind_rows(lapply(BUDGET_GRID, function(b) {
  G <- joint_grid(p_dev, y_dev, w_dev, k_old, b)
  a <- pick(G); c2 <- pick(G, G$sens_old >= target_sens - 1e-12)
  bind_rows(
    data.frame(budget = b, constraint = "none", a),
    if (!is.null(c2)) data.frame(budget = b, constraint = "old_sens>=Bang", c2))
}))
rownames(curve) <- NULL

# J4：老年 ≥ Bang 且总敏感度 ≥ S0 的最小名额（在名额网格上找，再在该名额附近细化）
feas <- curve[curve$constraint == "old_sens>=Bang" & curve$sens_total >= S0_sens_dev - 1e-12, ]
J4 <- NULL; J4_budget <- NA
if (nrow(feas) > 0) {
  b_hi <- min(feas$budget)
  for (b in seq(max(b_hi - 0.0025, min(BUDGET_GRID)), b_hi, by = 0.0001)) {
    G <- joint_grid(p_dev, y_dev, w_dev, k_old, b)
    c2 <- pick(G, G$sens_old >= target_sens - 1e-12 & G$sens_total >= S0_sens_dev - 1e-12)
    if (!is.null(c2)) { J4 <- c2; J4_budget <- b; break }
  }
}

cat("\n开发期目标：Bang 总转检 =", round(100 * target_ref, 2), "%，Bang 老年敏感度 =",
    round(100 * target_sens, 2), "%，S0 总敏感度 =", round(100 * S0_sens_dev, 2), "%\n")
if (is.null(J2)) cat(">>> J2 无解：在 Bang 名额下，任何切点组合都无法让老年组敏感度不低于 Bang。\n")
if (is.null(J4)) cat(">>> J4 在扫描范围内无解，可扩大 BUDGET_GRID。\n")


# =====================================================================
# 第 7 部分：锁定方程，套用到 2013–2018
# =====================================================================
fit_single <- fit_l2(as.matrix(dd[, VARS]), dd$case, dd$wt)
lk <- list.files(file.path(PROJECT_DIR, "results"), pattern = "^locked_main\\.csv$",
                 recursive = TRUE, full.names = TRUE)
if (length(lk) > 0) {
  lk <- read.csv(lk[which.max(file.mtime(lk))]); v <- setNames(lk$value, lk$term)
  cat("\n与 locked_main.csv 的最大系数差：", signif(max(abs(fit_single - v[names(fit_single)])), 3),
      "；切点差：", signif(abs(thr_single - v["THRESHOLD"]), 3), "（应接近 0）\n")
}

full$p <- NA_real_
full$p[full$dev] <- p_dev
if (any(full$eval)) full$p[full$eval] <- plogis(lp_of(fit_single, full[full$eval, ]))

STRATS <- list(
  S0_single          = c(young = thr_single, old = thr_single),
  S1_single_agecut   = c(young = thr_single, old = thr_single_old),
  S4_test_all_60plus = c(young = thr_single, old = 0)
)
add_j <- function(name, J) if (!is.null(J)) STRATS[[name]] <<- c(young = J$t_young, old = J$t_old)
add_j("J1_bangbudget_maxsens", J1)
add_j("J2_bangbudget_oldfloor", J2)
add_j("J3_S1budget_maxsens", J3)
add_j("J4_minbudget_oldfloor", J4)

old_all <- full$age >= SPLIT_AGE
for (s in names(STRATS)) {
  thr <- ifelse(old_all, STRATS[[s]]["old"], STRATS[[s]]["young"])
  full[[paste0("ref_", s)]] <- as.numeric(full$p >= thr)
}
for (s in names(STRATS)) {
  full[[paste0("dB_", s)]] <- full[[paste0("ref_", s)]] - full$ref_bang
  full[[paste0("d0_", s)]] <- full[[paste0("ref_", s)]] - full$ref_S0_single
  full[[paste0("d1_", s)]] <- full[[paste0("ref_", s)]] - full$ref_S1_single_agecut
}

# 名额曲线上每个点的切点也套用到 2013–2018（点估计，用来画评价期曲线）
if (any(full$eval)) {
  ee <- full$eval; pe_ <- full$p[ee]; ye <- full$case[ee]; we <- full$wt[ee]; oe <- old_all[ee]
  ev <- t(sapply(seq_len(nrow(curve)), function(i) {
    r <- as.numeric(pe_ >= ifelse(oe, curve$t_old[i], curve$t_young[i]))
    c(eval_ref_total = wmean(r, we), eval_sens_total = wsens(r, ye, we),
      eval_sens_old = wsens(r[oe], ye[oe], we[oe]), eval_sens_young = wsens(r[!oe], ye[!oe], we[!oe]))
  }))
  curve <- cbind(curve, ev)
}


# =====================================================================
# 第 8 部分：调查设计下的指标（Taylor 区间）
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
for (pe in periods) {
  base <- full[[pe]]
  for (gname in names(GROUPS)) {
    g <- if (isTRUE(GROUPS[[gname]])) rep(TRUE, nrow(full)) else (eval(GROUPS[[gname]], full) %in% TRUE)
    allg <- base & g; cases <- allg & full$case == 1; nonc <- allg & full$case == 0
    add <- function(strategy, measure, v) res[[length(res) + 1]] <<-
      data.frame(period = pe, group = gname, strategy = strategy, measure = measure,
                 estimate = v[1], lower95 = v[2], upper95 = v[3], row.names = NULL)
    add("Bang", "referral",    tay("ref_bang", allg))
    add("Bang", "sensitivity", tay("ref_bang", cases))
    add("Bang", "specificity", 100 - tay("ref_bang", nonc)[c(1, 3, 2)])
    for (s in names(STRATS)) {
      add(s, "referral",    tay(paste0("ref_", s), allg))
      add(s, "sensitivity", tay(paste0("ref_", s), cases))
      add(s, "specificity", 100 - tay(paste0("ref_", s), nonc)[c(1, 3, 2)])
      add(s, "sens_diff_vs_Bang",  tay(paste0("dB_", s), cases))
      add(s, "refer_diff_vs_Bang", tay(paste0("dB_", s), allg))
      if (s != "S0_single") {
        add(s, "sens_diff_vs_S0",  tay(paste0("d0_", s), cases))
        add(s, "refer_diff_vs_S0", tay(paste0("d0_", s), allg))
      }
      if (s != "S1_single_agecut") {
        add(s, "sens_diff_vs_S1",  tay(paste0("d1_", s), cases))
        add(s, "refer_diff_vs_S1", tay(paste0("d1_", s), allg))
      }
    }
  }
}
results <- bind_rows(res)


# =====================================================================
# 第 9 部分：图（名额–最高总敏感度曲线；Bang 名额下的老年–青年取舍曲线）
# =====================================================================
pt_of <- function(pe, s, g = "overall") {
  r <- results[results$period == pe & results$group == g & results$strategy == s, ]
  c(ref = r$estimate[r$measure == "referral"], sens = r$estimate[r$measure == "sensitivity"])
}
plot_curve <- function(file) {
  png(file, width = 2400, height = 1100, res = 200)
  par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
  for (pe in periods) {
    xs <- if (pe == "dev") "ref_total" else "eval_ref_total"
    ys <- if (pe == "dev") "sens_total" else "eval_sens_total"
    a <- curve[curve$constraint == "none", ]; b <- curve[curve$constraint != "none", ]
    sts <- c("Bang", "S0_single", "S1_single_agecut", "S4_test_all_60plus")
    cols <- c("grey40", "#2166ac", "#67a9cf", "#d95f02")
    Q <- sapply(sts, function(s) pt_of(pe, s))
    plot(100 * a[[xs]], 100 * a[[ys]], type = "l", lwd = 2, col = "#2166ac",
         xlab = "Overall referral (%)", ylab = "Overall sensitivity (%)",
         main = ifelse(pe == "dev", "2007-2012 (LOCO): best split per budget", "2013-2018: same cut-offs applied"),
         xlim = range(c(100 * a[[xs]], Q["ref", ]), na.rm = TRUE),
         ylim = range(c(100 * a[[ys]], 100 * b[[ys]], Q["sens", ]), na.rm = TRUE))
    if (nrow(b) > 0) lines(100 * b[[xs]], 100 * b[[ys]], lwd = 2, lty = 2, col = "#b2182b")
    points(Q["ref", ], Q["sens", ], pch = 19, col = cols, cex = 1.3)
    legend("bottomright", c("best split, no constraint", "best split, >=60 sens >= Bang", sts),
           col = c("#2166ac", "#b2182b", cols[seq_along(sts)]), lty = c(1, 2, rep(NA, length(sts))),
           pch = c(NA, NA, rep(19, length(sts))), lwd = 2, bty = "n", cex = 0.7)
  }
  invisible(dev.off())
}
plot_curve(file.path(OUT_DIR, "budget_curve.png"))

png(file.path(OUT_DIR, "tradeoff_at_bang_budget.png"), width = 1300, height = 1100, res = 200)
par(mar = c(4.5, 4.5, 3, 1))
plot(100 * G_bang$sens_old, 100 * G_bang$sens_young, type = "l", lwd = 2, col = "#2166ac",
     xlab = ">=60 sensitivity (%)", ylab = "<60 sensitivity (%)",
     main = "Bang budget, 2007-2012 (LOCO): all age splits")
abline(v = 100 * target_sens, lty = 2, col = "grey40")
s0 <- G_bang[which.min(abs(G_bang$t_old - thr_single)), ]
points(100 * s0$sens_old, 100 * s0$sens_young, pch = 19, col = "#2166ac", cex = 1.3)
if (!is.null(J1)) points(100 * J1$sens_old, 100 * J1$sens_young, pch = 17, col = "#1b7837", cex = 1.3)
if (!is.null(J2)) points(100 * J2$sens_old, 100 * J2$sens_young, pch = 15, col = "#b2182b", cex = 1.3)
legend("bottomleft", c("S0", "J1", "J2", "Bang >=60 sensitivity"), col = c("#2166ac", "#1b7837", "#b2182b", "grey40"),
       pch = c(19, 17, 15, NA), lty = c(NA, NA, NA, 2), bty = "n", cex = 0.8)
invisible(dev.off())


# =====================================================================
# 第 10 部分：保存和打印
# =====================================================================
thr_tab <- data.frame(strategy = names(STRATS),
                      cut_under60 = sapply(STRATS, `[`, "young"),
                      cut_60plus  = sapply(STRATS, `[`, "old"), row.names = NULL)
write.csv(results, file.path(OUT_DIR, "joint_results.csv"), row.names = FALSE)
write.csv(thr_tab, file.path(OUT_DIR, "joint_thresholds.csv"), row.names = FALSE)
write.csv(curve,   file.path(OUT_DIR, "budget_curve.csv"), row.names = FALSE)
write.csv(G_bang,  file.path(OUT_DIR, "grid_at_bang_budget_dev.csv"), row.names = FALSE)
writeLines(c(paste("开发期 Bang 总转检目标：", round(target_ref, 4)),
             paste("开发期 Bang 老年组敏感度目标：", round(target_sens, 4)),
             paste("开发期 S0 总敏感度：", round(S0_sens_dev, 4)),
             paste("开发期 S1 总转检：", round(S1_budget_dev, 4)),
             paste("J2 是否有解：", !is.null(J2)),
             paste("J4 最小名额（开发期）：", J4_budget),
             "所有切点只用 2007–2012 确定；2013–2018 为事后再评价，不是独立验证。",
             paste("运行时间：", Sys.time())),
           file.path(OUT_DIR, "joint_info.txt"))
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))

cat("\n===== 切点 =====\n"); print(thr_tab, digits = 4)
for (pe in periods) {
  cat("\n=====", ifelse(pe == "dev", "开发期 2007–2012（留一周期）", "2013–2018 再评价"), "=====\n")
  r <- results[results$period == pe & results$measure %in% c("referral", "sensitivity"), ]
  r$value <- sprintf("%.1f (%.1f-%.1f)", r$estimate, r$lower95, r$upper95)
  r$key <- paste(r$measure, r$group)
  out <- reshape(r[, c("strategy", "key", "value")], idvar = "strategy", timevar = "key", direction = "wide")
  names(out) <- sub("^value\\.", "", names(out)); print(out, row.names = FALSE)
  cat("-- 与 Bang / S0 / S1 的配对差（总体与 >=60，百分点）--\n")
  d <- results[results$period == pe & results$group %in% c("overall", ">=60") & grepl("diff", results$measure), ]
  d$value <- sprintf("%+.1f (%+.1f to %+.1f)", d$estimate, d$lower95, d$upper95)
  d$key <- paste(d$measure, d$group)
  out <- reshape(d[, c("strategy", "key", "value")], idvar = "strategy", timevar = "key", direction = "wide")
  names(out) <- sub("^value\\.", "", names(out)); print(out, row.names = FALSE)
}
cat("\n结果已保存到：", OUT_DIR, "\n")
