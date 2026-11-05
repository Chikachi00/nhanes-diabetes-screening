# =====================================================================
# 老年组敏感度问题：按年龄分切点 vs 分段方程（青年 / 老年各一个方程）
# ---------------------------------------------------------------------
# 四种策略，全部只在 2007–2012 开发期内确定（留一周期预测），锁定后
# 再在 2013–2018 上跑一次：
#   S0 单一方程 + 单一切点      = 已锁定的 V5 方程（参照）
#   S1 单一方程 + 老年组单独切点 （老年组切点：开发期敏感度不低于 Bang）
#   S2 分段方程 + 单一切点      （<60 岁、≥60 岁各拟合一个六项方程）
#   S3 分段方程 + 老年组单独切点
# 主方案：S1（10-31 与导师讨论后确定）；S2/S3 作为对照。
# 补充分析：USPSTF 2021 对照（35–70 岁且 BMI≥25）；按年龄组净获益；
#           Riley 等（BMJ 2020）最低样本量标准，判断分段方程样本量是否足够。
# 方程形式、六项变量、L2 C=1、权重与 V5 附录 B / refine_v5_equation.R 相同。
#
# 注意：这次修改是在看过 2013–2018 结果后提出的，所以这里的 2013–2018 结果
#       只能称为“事后修订后的再评价”，不是独立验证。
# 区间：调查设计 Taylor 线性化（权重、层、PSU），切点视为固定；
#       开发期区间未计入重新拟合的不确定性。
#
# 用法：放到 scripts 文件夹，只改第 0 部分，点 Source。约 1–2 分钟。
# 输出：results\segmented_v5_日期\
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
                         paste0("segmented_v5_", format(Sys.Date(), "%Y%m%d")))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

DEV_CYCLES  <- c("E", "F", "G")   # 2007–2012 开发
EVAL_CYCLES <- c("H", "I", "J")   # 2013–2018 再评价
C_L2      <- 1
KNOTS     <- c(2.5, 4.6, 7.3)     # 只为与原脚本一致，线性方程不用样条
SPLIT_AGE <- 60                   # 分段年龄界线：<60 为“青年段”，≥60 为“老年段”
VARS      <- c("a", "w", "male", "famhx", "htn", "inactive")   # V5 六项，线性年龄


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
# 第 5 部分：拟合与预测（单一方程 / 分段方程）
# =====================================================================
seg_of <- function(d) ifelse(d$age >= SPLIT_AGE, "old", "young")

fit_strategy <- function(d, w, segmented) {
  if (!segmented) return(list(all = fit_l2(as.matrix(d[, VARS]), d$case, w)))
  g <- seg_of(d)
  lapply(c(young = "young", old = "old"), function(s) {
    k <- g == s
    fit_l2(as.matrix(d[k, VARS]), d$case[k], w[k])
  })
}
predict_lp <- function(fit, d) {
  if (!is.null(fit$all)) return(lp_of(fit$all, d))
  g <- seg_of(d); lp <- rep(NA_real_, nrow(d))
  for (s in c("young", "old")) lp[g == s] <- lp_of(fit[[s]], d[g == s, ])
  lp
}
loco_lp <- function(d, segmented) {   # 开发期留一周期预测
  lp <- rep(NA_real_, nrow(d))
  for (cy in DEV_CYCLES) {
    tr <- d$cycle != cy; te <- d$cycle == cy
    lp[te] <- predict_lp(fit_strategy(d[tr, ], d$wt[tr], segmented), d[te, ])
  }
  lp
}

# 老年组切点：敏感度（加权）不低于目标的最高阈值
cut_for_sens <- function(p, y, w, target) {
  u <- sort(unique(p[y == 1]), decreasing = TRUE)
  sens <- sapply(u, function(t) sum(w[y == 1 & p >= t]) / sum(w[y == 1]))
  u[which(sens >= target - 1e-12)[1]]
}

dd <- full[full$dev, ]
k_old <- dd$age >= SPLIT_AGE
target_ref  <- wmean(dd$ref_bang, dd$wt)                             # Bang 转检比例（开发期）
target_sens <- wmean(dd$ref_bang[k_old & dd$case == 1], dd$wt[k_old & dd$case == 1])  # Bang 老年敏感度

lp_single <- loco_lp(dd, FALSE)
lp_seg    <- loco_lp(dd, TRUE)
p1 <- plogis(lp_single); p2 <- plogis(lp_seg)

thr_single     <- v5_threshold(p1, dd$wt, target_ref)
thr_single_old <- cut_for_sens(p1[k_old], dd$case[k_old], dd$wt[k_old], target_sens)
thr_seg        <- v5_threshold(p2, dd$wt, target_ref)
thr_seg_old    <- cut_for_sens(p2[k_old], dd$case[k_old], dd$wt[k_old], target_sens)
# 只放宽、不收紧：老年组切点不高于原切点
thr_single_old <- min(thr_single_old, thr_single)
thr_seg_old    <- min(thr_seg_old, thr_seg)

STRATS <- list(
  S0_single        = list(seg = FALSE, young = thr_single, old = thr_single),
  S1_single_agecut = list(seg = FALSE, young = thr_single, old = thr_single_old),
  S2_segmented     = list(seg = TRUE,  young = thr_seg,    old = thr_seg),
  S3_segmented_agecut = list(seg = TRUE, young = thr_seg,  old = thr_seg_old)
)

# 最终方程：全部开发期数据拟合（S0 / S1 的方程即已锁定的 V5 方程）
fit_single <- fit_strategy(dd, dd$wt, FALSE)
fit_seg    <- fit_strategy(dd, dd$wt, TRUE)

# 核对：单一方程应与已锁定的 locked_main.csv 一致
lk <- list.files(file.path(PROJECT_DIR, "results"), pattern = "^locked_main\\.csv$",
                 recursive = TRUE, full.names = TRUE)
if (length(lk) > 0) {
  lk <- read.csv(lk[which.max(file.mtime(lk))]); v <- setNames(lk$value, lk$term)
  cat("\n与 locked_main.csv 的最大系数差：", signif(max(abs(fit_single$all - v[names(fit_single$all)])), 3),
      "；切点差：", signif(abs(thr_single - v["THRESHOLD"]), 3), "（应接近 0）\n")
}

# 写回 full：开发期用留一周期预测，评价期用锁定方程
full$p_single <- NA_real_; full$p_seg <- NA_real_
full$p_single[full$dev] <- p1
full$p_seg[full$dev]    <- p2
ee <- full$eval
if (any(ee)) {
  full$p_single[ee] <- plogis(predict_lp(fit_single, full[ee, ]))
  full$p_seg[ee]    <- plogis(predict_lp(fit_seg, full[ee, ]))
}
old_all <- full$age >= SPLIT_AGE
for (s in names(STRATS)) {
  st <- STRATS[[s]]
  p  <- if (st$seg) full$p_seg else full$p_single
  thr <- ifelse(old_all, st$old, st$young)
  full[[paste0("ref_", s)]] <- as.numeric(p >= thr)
  full[[paste0("dB_", s)]]  <- full[[paste0("ref_", s)]] - full$ref_bang      # 与 Bang 的配对差
  full[[paste0("d0_", s)]]  <- full[[paste0("ref_", s)]] - full$ref_S0_single  # 与 S0 的配对差
}


# =====================================================================
# 第 6 部分：调查设计下的指标（Taylor 区间）
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
    allg  <- base & g
    cases <- allg & full$case == 1
    nonc  <- allg & full$case == 0
    add <- function(strategy, measure, v) res[[length(res) + 1]] <<-
      data.frame(period = pe, group = gname, strategy = strategy, measure = measure,
                 estimate = v[1], lower95 = v[2], upper95 = v[3], row.names = NULL)
    add("Bang", "referral", tay("ref_bang", allg))
    add("Bang", "sensitivity", tay("ref_bang", cases))
    add("Bang", "specificity", 100 - tay("ref_bang", nonc)[c(1, 3, 2)])
    add("USPSTF2021", "referral", tay("ref_uspstf", allg))
    add("USPSTF2021", "sensitivity", tay("ref_uspstf", cases))
    add("USPSTF2021", "specificity", 100 - tay("ref_uspstf", nonc)[c(1, 3, 2)])
    for (s in names(STRATS)) {
      add(s, "referral",    tay(paste0("ref_", s), allg))
      add(s, "sensitivity", tay(paste0("ref_", s), cases))
      add(s, "specificity", 100 - tay(paste0("ref_", s), nonc)[c(1, 3, 2)])
      add(s, "sens_diff_vs_Bang", tay(paste0("dB_", s), cases))
      add(s, "refer_diff_vs_Bang", tay(paste0("dB_", s), allg))
      if (s != "S0_single") {
        add(s, "sens_diff_vs_S0",  tay(paste0("d0_", s), cases))
        add(s, "refer_diff_vs_S0", tay(paste0("d0_", s), allg))
      }
    }
  }
}
results <- bind_rows(res)

# AUC（点估计；S0 与 S1 同一个分数，S2 与 S3 同一个分数）
auc_tab <- bind_rows(lapply(periods, function(pe) {
  k <- full[[pe]]
  data.frame(period = pe,
             AUC_single = w_auc(full$p_single[k], full$case[k], full$wt[k]),
             AUC_segmented = w_auc(full$p_seg[k], full$case[k], full$wt[k]),
             AUC_bang = w_auc(full$bang[k], full$case[k], full$wt[k]))
}))

# 系数表
coef_tab <- bind_rows(
  data.frame(equation = "single (= locked V5)", term = names(fit_single$all), coef = fit_single$all),
  data.frame(equation = paste0("segment age <", SPLIT_AGE), term = names(fit_seg$young), coef = fit_seg$young),
  data.frame(equation = paste0("segment age >=", SPLIT_AGE), term = names(fit_seg$old), coef = fit_seg$old)
) %>% mutate(OR = ifelse(term == "(Intercept)", NA, exp(coef)))
rownames(coef_tab) <- NULL

thr_tab <- data.frame(strategy = names(STRATS),
                      equation = ifelse(sapply(STRATS, `[[`, "seg"), "segmented", "single"),
                      cut_under60 = sapply(STRATS, `[[`, "young"),
                      cut_60plus  = sapply(STRATS, `[[`, "old"), row.names = NULL)

n_tab <- full %>% filter(dev | eval) %>%
  mutate(period = ifelse(dev, "dev", "eval"), segment = ifelse(age >= SPLIT_AGE, "60+", "<60")) %>%
  group_by(period, segment) %>% summarise(n = n(), cases = sum(case), .groups = "drop")


# =====================================================================
# 第 6b 部分：按年龄组净获益（决策曲线；Benitez-Aurioles 等 2024 的亚组净获益思路）
# 净获益 = 真阳性/N − 假阳性/N × pt/(1−pt)，加权；pt 为“愿意为找到 1 例而多转检的人数”
# 对应的阈值概率。规则型策略（Bang、USPSTF、S0–S3）在各 pt 下转检人群不变；
# “方程按 pt 转检”为连续风险的决策曲线。点估计，不带区间。
# =====================================================================
PT <- seq(0.01, 0.15, by = 0.005)
nb_rule <- function(ref, y, w, pt) {
  N <- sum(w); tp <- sum(w * ref * y) / N; fp <- sum(w * ref * (1 - y)) / N
  tp - fp * pt / (1 - pt)
}
nb <- list()
for (pe in periods) for (gname in names(GROUPS)) {
  g <- if (isTRUE(GROUPS[[gname]])) rep(TRUE, nrow(full)) else (eval(GROUPS[[gname]], full) %in% TRUE)
  k <- full[[pe]] & g; y <- full$case[k]; w <- full$wt[k]
  rules <- c(list(refer_all = rep(1, sum(k)), Bang = full$ref_bang[k], USPSTF2021 = full$ref_uspstf[k]),
             setNames(lapply(names(STRATS), function(s) full[[paste0("ref_", s)]][k]), names(STRATS)))
  for (pt in PT) {
    v <- sapply(rules, nb_rule, y = y, w = w, pt = pt)
    v["equation_by_pt"] <- nb_rule(as.numeric(full$p_single[k] >= pt), y, w, pt)
    nb[[length(nb) + 1]] <- data.frame(period = pe, group = gname, pt = pt, strategy = names(v),
                                       net_benefit_per1000 = 1000 * unname(v))
  }
}
nb_tab <- bind_rows(nb)


# =====================================================================
# 第 6c 部分：最低样本量（Riley 等 BMJ 2020，二分类结局三条标准；与 R 包 pmsampsize 同一公式）
# 用开发期每段的“表观 Cox-Snell R²”（六项方程，未加权 glm）作为预期 R²，参数数 = 6。
# =====================================================================
riley_n <- function(n, events, r2cs, P = length(VARS), S = 0.9, delta = 0.05) {
  phi <- events / n
  r2max <- 1 - (phi^phi * (1 - phi)^(1 - phi))^2
  n1 <- P / ((S - 1) * log(1 - r2cs / S))                    # 标准 1：收缩系数 ≥0.9
  S2 <- r2cs / (r2cs + delta * r2max)
  n2 <- P / ((S2 - 1) * log(1 - r2cs / S2))                  # 标准 2：R² 乐观偏差 ≤0.05
  n3 <- (1.96 / 0.05)^2 * phi * (1 - phi)                    # 标准 3：总体风险误差 ≤0.05
  need <- ceiling(max(n1, n2, n3))
  c(available_n = n, available_events = events, prevalence = phi, R2_CS = r2cs,
    n_crit1 = ceiling(n1), n_crit2 = ceiling(n2), n_crit3 = ceiling(n3),
    required_n = need, required_events = ceiling(need * phi), required_EPP = need * phi / P,
    enough = as.numeric(n >= need))
}
ss_rows <- list()
for (sg in c("all", "<60", "60+")) {
  k <- if (sg == "all") rep(TRUE, nrow(dd)) else if (sg == "<60") dd$age < SPLIT_AGE else dd$age >= SPLIT_AGE
  dk <- dd[k, ]
  f1 <- glm(reformulate(VARS, "case"), family = binomial, data = dk)
  r2 <- 1 - exp(-(f1$null.deviance - f1$deviance) / nrow(dk))
  ss_rows[[sg]] <- data.frame(segment = sg, t(riley_n(nrow(dk), sum(dk$case), r2)), row.names = NULL)
}
ss_tab <- bind_rows(ss_rows)


# =====================================================================
# 第 7 部分：图（各年龄组敏感度和转检比例）
# =====================================================================
plot_period <- function(pe, file) {
  labs <- c("Bang", "USPSTF2021", names(STRATS)); cols <- c("grey70", "grey40", "#2166ac", "#67a9cf", "#b2182b", "#ef8a62")
  png(file, width = 2400, height = 1100, res = 200)
  par(mfrow = c(1, 2), mar = c(4, 4.5, 3, 1))
  for (ms in c("sensitivity", "referral")) {
    r <- results[results$period == pe & results$measure == ms & results$strategy %in% labs, ]
    M <- matrix(NA, length(labs), 4, dimnames = list(labs, names(GROUPS)))
    L <- U <- M
    for (i in seq_len(nrow(r))) {
      M[r$strategy[i], r$group[i]] <- r$estimate[i]
      L[r$strategy[i], r$group[i]] <- r$lower95[i]; U[r$strategy[i], r$group[i]] <- r$upper95[i]
    }
    bp <- barplot(M, beside = TRUE, col = cols, ylim = c(0, 100), border = NA,
                  ylab = paste(ms, "(%)"), main = paste0(ms, ", ", ifelse(pe == "dev", "2007-2012 (LOCO)", "2013-2018")))
    suppressWarnings(arrows(bp, L, bp, U, angle = 90, code = 3, length = 0.02, col = "grey30"))
    if (ms == "sensitivity") legend("topleft", labs, fill = cols, border = NA, bty = "n", cex = 0.75)
  }
  invisible(dev.off())
}
for (pe in periods) plot_period(pe, file.path(OUT_DIR, paste0("segmented_", pe, ".png")))

# 按年龄组净获益图（每个时期一张）
plot_nb <- function(pe, file) {
  sts <- c("refer_all", "Bang", "USPSTF2021", "S0_single", "S1_single_agecut", "equation_by_pt")
  cols <- c("grey80", "grey60", "grey30", "#2166ac", "#67a9cf", "#1b7837")
  ltys <- c(3, 1, 1, 1, 2, 1)
  png(file, width = 2600, height = 750, res = 200)
  par(mfrow = c(1, 4), mar = c(4, 4.5, 3, 1))
  for (gname in names(GROUPS)) {
    r <- nb_tab[nb_tab$period == pe & nb_tab$group == gname, ]
    yl <- range(c(0, r$net_benefit_per1000[r$strategy != "refer_all"]), na.rm = TRUE)
    yl[1] <- min(yl[1], -0.2 * yl[2])
    plot(NA, xlim = range(PT), ylim = yl, xlab = "Threshold probability",
         ylab = "Net benefit per 1,000", main = paste0(gname, ifelse(pe == "dev", " (2007-2012)", " (2013-2018)")))
    abline(h = 0, col = "grey85")
    for (i in seq_along(sts)) {
      ri <- r[r$strategy == sts[i], ]
      lines(ri$pt, ri$net_benefit_per1000, col = cols[i], lty = ltys[i], lwd = 2)
    }
    if (gname == "overall") legend("topright", sts, col = cols, lty = ltys, lwd = 2, bty = "n", cex = 0.7)
  }
  invisible(dev.off())
}
for (pe in periods) plot_nb(pe, file.path(OUT_DIR, paste0("net_benefit_by_age_", pe, ".png")))


# =====================================================================
# 第 8 部分：保存和打印
# =====================================================================
write.csv(results,  file.path(OUT_DIR, "seg_results.csv"), row.names = FALSE)
write.csv(auc_tab,  file.path(OUT_DIR, "seg_auc.csv"), row.names = FALSE)
write.csv(coef_tab, file.path(OUT_DIR, "seg_coefficients.csv"), row.names = FALSE)
write.csv(thr_tab,  file.path(OUT_DIR, "seg_thresholds.csv"), row.names = FALSE)
write.csv(n_tab,    file.path(OUT_DIR, "seg_sample_sizes.csv"), row.names = FALSE)
write.csv(nb_tab,   file.path(OUT_DIR, "net_benefit_by_age.csv"), row.names = FALSE)
write.csv(ss_tab,   file.path(OUT_DIR, "riley_sample_size.csv"), row.names = FALSE)
writeLines(c(paste("分段年龄：", SPLIT_AGE),
             paste("开发期 Bang 转检比例目标：", round(target_ref, 4)),
             paste("开发期 Bang 老年组敏感度目标：", round(target_sens, 4)),
             "所有切点和方程只用 2007–2012 确定；2013–2018 为事后修订后的再评价，不是独立验证。",
             paste("运行时间：", Sys.time())),
           file.path(OUT_DIR, "seg_info.txt"))
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))

show <- function(pe) {
  r <- results[results$period == pe & results$measure %in% c("referral", "sensitivity"), ]
  r$value <- sprintf("%.1f (%.1f-%.1f)", r$estimate, r$lower95, r$upper95)
  print(tidyr_free_wide(r))
}
tidyr_free_wide <- function(r) {
  r$key <- paste(r$measure, r$group)
  out <- reshape(r[, c("strategy", "key", "value")], idvar = "strategy", timevar = "key", direction = "wide")
  names(out) <- sub("^value\\.", "", names(out)); out
}
cat("\n===== 每段样本量和病例数 =====\n"); print(as.data.frame(n_tab))
cat("\n===== 切点 =====\n"); print(thr_tab, digits = 4)
cat("\n===== AUC =====\n"); print(auc_tab, digits = 3)
cat("\n===== 最低样本量（Riley 标准）：现有 vs 需要 =====\n")
print(ss_tab[, c("segment", "available_n", "available_events", "R2_CS", "required_n", "required_events", "enough")], digits = 3)
cat("\n===== 每千人净获益（pt = 0.03 / 0.05 / 0.10）=====\n")
nb_show <- nb_tab[abs(nb_tab$pt - 0.03) < 1e-9 | abs(nb_tab$pt - 0.05) < 1e-9 | abs(nb_tab$pt - 0.10) < 1e-9, ]
nb_show$key <- paste0(nb_show$group, "@", nb_show$pt)
for (pe in periods) {
  cat("--", pe, "--\n")
  w_ <- reshape(nb_show[nb_show$period == pe, c("strategy", "key", "net_benefit_per1000")],
                idvar = "strategy", timevar = "key", direction = "wide")
  names(w_) <- sub("^net_benefit_per1000\\.", "", names(w_)); print(w_, digits = 3, row.names = FALSE)
}
for (pe in periods) {
  cat("\n=====", ifelse(pe == "dev", "开发期 2007–2012（留一周期）", "2013–2018 再评价"), ": 转检比例 / 敏感度 % =====\n")
  show(pe)
}
cat("\n结果已保存到：", OUT_DIR, "\n")
