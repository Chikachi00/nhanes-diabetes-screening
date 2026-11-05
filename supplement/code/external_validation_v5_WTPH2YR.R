# =====================================================================
# 外部验证：NHANES 2021–2023（August 2021–August 2023，文件后缀 _L）
# external_validation_v5.R
# ---------------------------------------------------------------------
# 方程和切点全部沿用 2007–2012 开发期确定的版本（与之前所有脚本完全相同），
# 不在 2021–2023 上做任何拟合、选择或调参。2021–2023 从未参与任何分析，
# 这是 V5 的第一次真正独立验证（时间外推）。
#
# 比较的方案：Bang ≥5；USPSTF 2021；V5_R0；V5_R1（主方案）；V5_R4；F4_R0；F4_R1
# 输出（2013–2018 与 2021–2023 并排）：
#   1) 总体、各年龄组、高 / 低可及组的转检、敏感度、特异度、PPV（Taylor 区间）
#   2) 与 Bang 的配对差；可及性差距（原始、年龄标准化、相对 Bang 的变化）
#   3) AUC、校准（O:E、校准截距和斜率）
#   4) 同名额曲线（单一切点 vs 老年切点固定）及差值的 PSU bootstrap 95% 带
#
# 2021–2023 的已知差异（论文方法部分需写明）：
#   0. 2021–2023 的 MCQ 删除了家族糖尿病史（MCQ300C）。因此外部验证的主分析用开发期预先拟合的
#      “无家族史”版本：V5nf（V5 去掉家族史）、F3（年龄、腰高比、男性）、Bang_nf（Bang 去掉家族史，
#      整数切点在开发期按原 Bang 转检比例选定）。同一批无家族史版本也在 2013–2018 上评价，
#      用来量化“去掉家族史”本身的代价，使 2021–2023 的结果可以解释。
#   a. 体力活动问卷改版：旧版 GPAQ（工作、交通、休闲 5 个领域的天数）被换成
#      只问休闲时间中等 / 高强度活动的频率（PAD790Q/U、PAD810Q/U）。本脚本把
#      “每周中等 + 高强度次数合计 ≥5”视为活动。这会使“不活动”比例偏高，
#      影响 V5（不活动系数 0.19）和 Bang（活动减 1 分）；4 项简化版 F4 不受影响。
#      如果你的 PAQ_L 仍有 PAQ605 等旧变量，脚本会自动用旧定义。
#   b. 2021–2023 的 HbA1c 在采血子样本中测量，权重用采血权重 WTPH2YR（GHB_L 说明文件），
#      不用 WTMEC2YR；只有 1 个周期，不除以 3；只保留 WTPH2YR > 0 的人。抽样设计与以前不同
#      （COVID 后恢复），Taylor 区间仍按层和 PSU 计算。2013–2018 并排参照仍用 WTMEC2YR / 3，与主分析一致。
#   修订（2026-11-04）：权重改为 WTPH2YR；策略标签 S0 / S1 / S4 改为 R0 / R1 / R4，计算不变。
#   c. 可及性变量（HIQ011、HUQ030）如有缺失，只在可及性分析中排除，不影响总体分析。
#
# 用法：放到 scripts 文件夹，只改第 0 部分，点 Source。缺少的 _L 文件会自动从 CDC 下载。
# 输出：results\external_validation_v5_日期\
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
                         paste0("external_validation_v5_WTPH2YR_", format(Sys.Date(), "%Y%m%d")))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

DEV_CYCLES  <- c("E", "F", "G")   # 2007–2012 开发（只用来重建锁定方程和切点）
EVAL_CYCLES <- c("H", "I", "J")   # 2013–2018 再评价（并排参照；设为 character(0) 可跳过）
EXT_CYCLES  <- c("L")             # 2021–2023 外部验证
C_L2      <- 1
KNOTS     <- c(2.5, 4.6, 7.3)
SPLIT_AGE <- 60
VARS      <- c("a", "w", "male", "famhx", "htn", "inactive")   # V5 六项（已锁定）
VARS4     <- c("a", "w", "male", "famhx")                      # 4 项简化版
# 2021–2023 没有家族史题目：以下两个版本在开发期预先拟合（不看 2021–2023），作为外部验证的主分析
VARS_NF   <- c("a", "w", "male", "htn", "inactive")             # V5 去掉家族史
VARS3     <- c("a", "w", "male")                                # 4 项版去掉家族史
B_BOOT    <- 500
BUDGETS   <- seq(0.25, 0.55, by = 0.005)
AUTO_DOWNLOAD <- TRUE
FAMHX_VAR_EXT <- character(0)   # 若 2021–2023 的家族史题目改了名，在这里填新变量名（取值 1 = 是，2 = 否）
CDC_URL   <- "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2021/DataFiles/"   # 2021–2023 文件所在目录
EXT_MODULES <- c("DEMO", "DIQ", "BMX", "BPQ", "MCQ", "PAQ", "GHB", "HIQ", "HUQ")


# ---- 缺少的 2021–2023 文件自动下载 ----
if (AUTO_DOWNLOAD) {
  ext_dir <- file.path(DATA_DIR, "2021-2023"); dir.create(ext_dir, recursive = TRUE, showWarnings = FALSE)
  for (mod in EXT_MODULES) for (sfx in EXT_CYCLES) {
    fn <- paste0(mod, "_", sfx, ".xpt")
    if (length(list.files(DATA_DIR, pattern = paste0("^", fn, "$"), recursive = TRUE, ignore.case = TRUE)) > 0) next
    cat("下载", fn, "...\n")
    ok <- tryCatch({ download.file(paste0(CDC_URL, fn), file.path(ext_dir, fn), mode = "wb", quiet = TRUE); TRUE },
                   error = function(e) FALSE)
    if (!ok) cat(">>> 下载失败：", fn, "。请到 CDC NHANES 2021–2023 页面手动下载，放进", ext_dir, "\n")
  }
}


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
    left_join(sel(m("BMX"), c("SEQN", "BMXWAIST", "BMXHT", "BMXBMI")), by = "SEQN") %>%
    left_join(sel(m("BPQ"), c("SEQN", "BPQ020", "BPQ080")), by = "SEQN") %>%
    left_join(sel(m("MCQ"), unique(c("SEQN", "MCQ300C", FAMHX_VAR_EXT))), by = "SEQN") %>%
    left_join(sel(m("PAQ"), c("SEQN", "PAQ605", "PAQ610", "PAQ620", "PAQ625", "PAQ635",
                              "PAQ640", "PAQ650", "PAQ655", "PAQ665", "PAQ670",
                              "PAD790Q", "PAD790U", "PAD810Q", "PAD810U")), by = "SEQN") %>%
    left_join(sel(m("GHB"), c("SEQN", "LBXGH", "WTPH2YR")), by = "SEQN") %>%
    left_join(sel(m("HIQ"), c("SEQN", "HIQ011")), by = "SEQN") %>%
    left_join(sel(m("HUQ"), c("SEQN", "HUQ030")), by = "SEQN") %>%
    mutate(cycle = suffix)
}

per_week <- function(q, u) {     # 频率换算成每周次数（2021–2023 的 PAD790Q/U、PAD810Q/U）
  u <- toupper(trimws(as.character(u))); q <- suppressWarnings(as.numeric(q))
  case_when(q %in% 0 ~ 0, q >= 7777 ~ NA_real_,
            u == "D" ~ q * 7, u == "W" ~ q, u == "M" ~ q / 4.345, u == "Y" ~ q / 52.18,
            TRUE ~ NA_real_)
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

all_cycles <- c(DEV_CYCLES, EVAL_CYCLES, EXT_CYCLES)
raw <- bind_rows(lapply(all_cycles, read_cycle))

# ---- 外部验证期变量检查：哪些原始变量整列缺失（2021–2023 个别题目改名或删除）----
key_src <- c(DEMO = "RIDAGEYR", BMX = "BMXWAIST", BMX = "BMXHT", BMX = "BMXBMI", BPQ = "BPQ020",
             MCQ = "MCQ300C", DIQ = "DIQ010", GHB = "LBXGH", HIQ = "HIQ011", HUQ = "HUQ030")
rx <- raw[raw$cycle %in% EXT_CYCLES & raw$RIDSTATR %in% 2, ]
miss_src <- key_src[sapply(key_src, function(v) mean(!is.na(rx[[v]])) < 0.2)]
if (length(miss_src) > 0) {
  cat("\n>>> 2021–2023 中以下变量缺失或几乎全缺：", paste(miss_src, collapse = ", "), "\n")
  for (mod in unique(names(miss_src))) {
    cn <- names(read_nhanes(mod, EXT_CYCLES[1]))
    cat(">>>", mod, "_", EXT_CYCLES[1], " 的全部变量名：\n", sep = ""); print(cn)
  }
  if (length(FAMHX_VAR_EXT) == 0 && "MCQ300C" %in% miss_src)
    cat(">>> 2021–2023 没有家族糖尿病史：外部验证使用开发期预先拟合的“无家族史”版本（V5nf、F3、Bang_nf）。\n")
  if (any(setdiff(miss_src, "MCQ300C") %in% key_src))
    stop("除家族史外还有关键变量缺失（见上），请把打印的变量名发给 Claude。")
}

if (all(is.na(raw$WTPH2YR[raw$cycle %in% EXT_CYCLES])))
  stop("GHB_L 里没有读到 WTPH2YR：请确认 GHB_L.xpt 是最新版本（采血权重在该文件中）。")
cat("外部验证期 WTPH2YR > 0 的人数：", sum(raw$WTPH2YR[raw$cycle %in% EXT_CYCLES] > 0, na.rm = TRUE), "\n")

if (length(FAMHX_VAR_EXT) > 0) {   # 2021–2023 用新变量名的家族史
  kx <- raw$cycle %in% EXT_CYCLES
  raw$MCQ300C[kx] <- raw[[FAMHX_VAR_EXT]][kx]
}
full <- raw %>%
  filter(RIDSTATR == 2) %>%
  mutate(
    wt     = ifelse(cycle %in% EXT_CYCLES, coalesce(WTPH2YR, 0), WTMEC2YR / 3),   # 外部验证：采血权重 WTPH2YR，1 个周期，不除以 3；未采血者记 0，不进入样本但保留在设计中
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
    inactive_old = case_when(any_active %in% TRUE ~ 0,
                             rowSums(is.na(cbind(d1, d2, d3, d4, d5))) == 0 ~ 1,
                             TRUE ~ NA_real_),
    # 2021–2023 问卷改版（只问休闲时间中等 / 高强度活动的频率）：每周次数合计 ≥5 视为活动
    mod_wk = per_week(PAD790Q, PAD790U), vig_wk = per_week(PAD810Q, PAD810U),
    inactive_new = ifelse(is.na(mod_wk) | is.na(vig_wk), NA_real_, as.numeric(mod_wk + vig_wk < 5)),
    paq_version = ifelse(rowSums(!is.na(cbind(PAQ605, PAQ620, PAQ635, PAQ650, PAQ665))) > 0, "GPAQ",
                         ifelse(!is.na(PAD790Q) | !is.na(PAD810Q), "PAD790", NA)),
    inactive = ifelse(paq_version %in% "PAD790", inactive_new, inactive_old),
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
    bang_nf   = age_pts + male + htn + obese_pts - active,   # 不含家族史的 Bang 积分
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

need <- c("age", "male", "waist", "height", "bmi", "whtr", "htn",
          "inactive", "case")
# 开发期 / 2013–2018 与之前完全相同（要求可及性不缺失）；外部验证期可及性缺失者保留，只在可及性分析中排除
full <- full %>%
  mutate(in_sample = age >= 20 & DIQ010 %in% 2 & !(RIDEXPRG %in% 1) &
           complete.cases(across(all_of(need))) &
           ((!is.na(access) & !is.na(famhx)) | cycle %in% EXT_CYCLES) &
           (!(cycle %in% EXT_CYCLES) | (!is.na(wt) & wt > 0)))   # 外部验证期只保留采血权重 > 0 的人
full$in_sample[is.na(full$in_sample)] <- FALSE
full$dev  <- full$in_sample & full$cycle %in% DEV_CYCLES
full$eval <- full$in_sample & full$cycle %in% EVAL_CYCLES
full$ext  <- full$in_sample & full$cycle %in% EXT_CYCLES
if (!any(full$ext)) stop("外部验证期没有可用样本：检查 *_L.xpt 文件是否已下载，以及变量是否齐全。")
cat("外部验证期样本：", sum(full$ext), "人，病例", sum(full$case[full$ext]), "人\n")

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
# 第 5 部分：重建锁定方程和切点（只用开发期，与之前所有脚本相同）
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

MODELS <- list(V5 = VARS, F4 = VARS4, V5nf = VARS_NF, F3 = VARS3)
old_all <- full$age >= SPLIT_AGE
THR <- list(); FIT <- list()
for (m in names(MODELS)) {
  v <- MODELS[[m]]
  p_dev <- loco_p(dd, v)
  t0 <- v5_threshold(p_dev, dd$wt, target_ref)
  t1 <- min(cut_for_sens(p_dev[k_old], dd$case[k_old], dd$wt[k_old], target_sens), t0)
  THR[[m]] <- c(single = t0, old = t1)
  FIT[[m]] <- fit_l2(as.matrix(dd[, v]), dd$case, dd$wt)
  pm <- paste0("p_", m); full[[pm]] <- NA_real_
  full[[pm]][full$dev] <- p_dev
  oth <- full$eval | full$ext
  full[[pm]][oth] <- plogis(lp_of(FIT[[m]], full[oth, ]))
  full[[paste0("lp_", m)]] <- qlogis(full[[pm]])
  full[[paste0("ref_", m, "_R0")]] <- as.numeric(full[[pm]] >= t0)
  full[[paste0("ref_", m, "_R1")]] <- as.numeric(full[[pm]] >= ifelse(old_all, t1, t0))
}
full$ref_V5_R4   <- as.numeric(ifelse(old_all & !is.na(full$ref_V5_R0), 1, full$ref_V5_R0))
full$ref_V5nf_R4 <- as.numeric(ifelse(old_all & !is.na(full$ref_V5nf_R0), 1, full$ref_V5nf_R0))

# 无家族史 Bang：开发期选整数切点，使转检比例最接近原 Bang（与方程切点的规则一致）
ks <- 0:12
sh <- sapply(ks, function(k) wmean(as.numeric(dd$bang_nf >= k), dd$wt))
BANG_NF_CUT <- ks[which.min(abs(sh - target_ref))]
full$ref_bang_nf <- as.numeric(full$bang_nf >= BANG_NF_CUT)
cat("无家族史 Bang 的切点（开发期）：≥", BANG_NF_CUT, "\n")

RULES <- c("Bang", "Bang_nf", "USPSTF2021", "V5_R0", "V5_R1", "V5_R4", "F4_R0", "F4_R1",
           "V5nf_R0", "V5nf_R1", "V5nf_R4", "F3_R0", "F3_R1")
refcol <- function(r) switch(r, Bang = "ref_bang", Bang_nf = "ref_bang_nf", USPSTF2021 = "ref_uspstf", paste0("ref_", r))
for (r in RULES) {
  full[[paste0("dB_", r)]] <- full[[refcol(r)]] - full$ref_bang
  full[[paste0("dN_", r)]] <- full[[refcol(r)]] - full$ref_bang_nf
}
avail_in <- function(r, pe) !all(is.na(full[[refcol(r)]][full[[pe]]]))

lk <- list.files(file.path(PROJECT_DIR, "results"), pattern = "^locked_main\\.csv$",
                 recursive = TRUE, full.names = TRUE)
if (length(lk) > 0) {
  lk <- read.csv(lk[which.max(file.mtime(lk))]); vv <- setNames(lk$value, lk$term)
  cat("\nV5 与 locked_main.csv：最大系数差", signif(max(abs(FIT$V5 - vv[names(FIT$V5)])), 3),
      "；切点差", signif(abs(THR$V5["single"] - vv["THRESHOLD"]), 3), "（应接近 0）\n")
}

# 体力活动定义的影响：各时期“不活动”比例
paq_tab <- full %>% filter(dev | eval | ext) %>%
  mutate(period = ifelse(dev, "dev", ifelse(eval, "eval", "ext"))) %>%
  group_by(period, paq_version) %>%
  summarise(n = n(), inactive_pct = 100 * sum(wt * inactive) / sum(wt), .groups = "drop")


# =====================================================================
# 第 6 部分：指标（Taylor 区间）
# =====================================================================
des <- svydesign(ids = ~psu, strata = ~strata, weights = ~wt, nest = TRUE, data = full)
periods <- c(eval = "eval", ext = "ext")
if (!any(full$eval)) periods <- periods["ext"]
tay <- function(var, idx) {
  m <- svymean(as.formula(paste0("~", var)), des[idx, ], na.rm = TRUE)
  ci <- confint(m); 100 * c(unname(coef(m)), ci[1], ci[2])
}
ci_row <- function(est, se) 100 * c(est, est - 1.96 * se, est + 1.96 * se)
GROUPS <- list(overall = TRUE, `20-39` = quote(age_grp == "20-39"), `40-59` = quote(age_grp == "40-59"),
               `>=60` = quote(age_grp == ">=60"), high_access = quote(access == "high"),
               low_access = quote(access == "low"))

res <- list()
for (pe in periods) for (gname in names(GROUPS)) {
  g <- if (isTRUE(GROUPS[[gname]])) rep(TRUE, nrow(full)) else (eval(GROUPS[[gname]], full) %in% TRUE)
  allg <- full[[pe]] & g; cases <- allg & full$case == 1; nonc <- allg & full$case == 0
  add <- function(strategy, measure, v) res[[length(res) + 1]] <<-
    data.frame(period = pe, group = gname, strategy = strategy, measure = measure,
               estimate = v[1], lower95 = v[2], upper95 = v[3], row.names = NULL)
  for (r in RULES) {
    if (!avail_in(r, pe)) next
    v <- refcol(r)
    add(r, "referral",    tay(v, allg))
    add(r, "sensitivity", tay(v, cases))
    add(r, "specificity", 100 - tay(v, nonc)[c(1, 3, 2)])
    add(r, "PPV",         tay("case", allg & full[[v]] %in% 1))
    if (r != "Bang" && avail_in("Bang", pe)) add(r, "sens_diff_vs_Bang", tay(paste0("dB_", r), cases))
    if (r != "Bang_nf") add(r, "sens_diff_vs_Bang_nf", tay(paste0("dN_", r), cases))
  }
}
results <- bind_rows(res)

# 可及性差距
gap_of <- function(var, idx) {
  b <- svyby(as.formula(paste0("~", var)), ~access, des[idx, ], svymean, na.rm = TRUE, covmat = TRUE)
  ct <- svycontrast(b, setNames(c(1, -1), c("high", "low")))
  ci_row(unname(coef(ct)), unname(SE(ct)))
}
gap_std <- function(var, idx) {
  b <- svyby(as.formula(paste0("~", var)), ~interaction(access, age_grp), des[idx, ], svymean,
             na.rm = TRUE, covmat = TRUE)
  agew <- prop.table(tapply(full$wt[idx], full$age_grp[idx], sum))
  cf <- setNames(rep(0, nrow(b)), rownames(b))
  for (g in names(agew)) { cf[paste("high", g, sep = ".")] <- agew[g]; cf[paste("low", g, sep = ".")] <- -agew[g] }
  ct <- svycontrast(b, cf); ci_row(unname(coef(ct)), unname(SE(ct)))
}
gap_rows <- list()
for (pe in periods) {
  cases <- full[[pe]] & full$case == 1 & !is.na(full$access)
  base  <- full[[pe]] & !is.na(full$access)
  for (r in RULES) {
    if (!avail_in(r, pe)) next
    addg <- function(measure, x) gap_rows[[length(gap_rows) + 1]] <<-
        data.frame(period = pe, strategy = r, measure = measure, estimate = x[1], lower95 = x[2],
                   upper95 = x[3], row.names = NULL)
    v <- refcol(r)
    addg("sens_gap_high_minus_low", tryCatch(gap_of(v, cases), error = function(e) rep(NA, 3)))
    addg("sens_gap_age_standardised", tryCatch(gap_std(v, cases), error = function(e) rep(NA, 3)))
    addg("referral_gap_high_minus_low", tryCatch(gap_of(v, base), error = function(e) rep(NA, 3)))
    if (r != "Bang_nf") addg("change_in_sens_gap_vs_Bang_nf",
                             tryCatch(gap_of(paste0("dN_", r), cases), error = function(e) rep(NA, 3)))
  }
}
gap_tab <- bind_rows(gap_rows)

# AUC 与校准
calib <- function(lp, y, w) {
  nw <- w / mean(w)
  c(O_E = sum(w * y) / sum(w * plogis(lp)),
    cal_intercept = unname(coef(suppressWarnings(glm(y ~ offset(lp), family = binomial, weights = nw)))[1]),
    cal_slope = unname(coef(suppressWarnings(glm(y ~ lp, family = binomial, weights = nw)))[2]))
}
perf_tab <- bind_rows(lapply(periods, function(pe) {
  k <- full[[pe]]
  bind_rows(lapply(names(MODELS), function(m) {
    p <- full[[paste0("p_", m)]]
    if (all(is.na(p[k]))) return(NULL)
    data.frame(period = pe, model = m,
               AUC = w_auc(p[k], full$case[k], full$wt[k]),
               AUC_60plus = w_auc(p[k & old_all], full$case[k & old_all], full$wt[k & old_all]),
               AUC_under60 = w_auc(p[k & !old_all], full$case[k & !old_all], full$wt[k & !old_all]),
               t(calib(full[[paste0("lp_", m)]][k], full$case[k], full$wt[k])))
  }), bind_rows(lapply(c("bang", "bang_nf"), function(sc) {
    if (all(is.na(full[[sc]][k]))) return(NULL)
    x <- full[[sc]]
    data.frame(period = pe, model = paste(sc, "score"), AUC = w_auc(x[k], full$case[k], full$wt[k]),
               AUC_60plus = w_auc(x[k & old_all], full$case[k & old_all], full$wt[k & old_all]),
               AUC_under60 = w_auc(x[k & !old_all], full$case[k & !old_all], full$wt[k & !old_all]))
  })))
}))


# =====================================================================
# 第 7 部分：同名额曲线（外部验证期；方程和切点固定；PSU bootstrap）
# =====================================================================
fam_curve <- function(p, y, w, old, t_old, budgets) {
  W <- sum(w); Wc <- sum(w * y)
  if (is.na(t_old)) {
    o <- order(-p); x <- c(0, cumsum(w[o])) / W; s <- c(0, cumsum((w * y)[o])) / Wc
  } else {
    ro <- old & p >= t_old; bw <- sum(w[ro]); bc <- sum((w * y)[ro])
    yo <- which(!old); o <- yo[order(-p[yo])]
    x <- c(bw, bw + cumsum(w[o])) / W; s <- c(bc, bc + cumsum((w * y)[o])) / Wc
  }
  approx(x, s, xout = budgets, ties = "ordered", rule = 1)$y
}
CURVES <- list(V5nf.single = list(m = "V5nf", t = NA), V5nf.oldfix = list(m = "V5nf", t = THR$V5nf["old"]),
               F3.single = list(m = "F3", t = NA), F3.oldfix = list(m = "F3", t = THR$F3["old"]))
curves_for <- function(d, w) sapply(names(CURVES), function(cn)
  fam_curve(d[[paste0("p_", CURVES[[cn]]$m)]], d$case, w, d$age >= SPLIT_AGE, CURVES[[cn]]$t, BUDGETS))
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
qna <- function(x, pr) if (mean(is.na(x)) > 0.05) NA_real_ else unname(quantile(x, pr, na.rm = TRUE))

dX <- full[full$ext, ]
est <- curves_for(dX, dX$wt)
bt <- array(NA_real_, c(nrow(est), ncol(est), B_BOOT), dimnames = list(NULL, colnames(est), NULL))
cat("bootstrap（B =", B_BOOT, "）...\n")
for (b in seq_len(B_BOOT)) bt[, , b] <- curves_for(dX, boot_weights(dX))
curve_tab <- bind_rows(lapply(colnames(est), function(j) data.frame(
  curve = j, budget = BUDGETS, sens = est[, j],
  lower95 = apply(bt[, j, ], 1, qna, pr = 0.025), upper95 = apply(bt[, j, ], 1, qna, pr = 0.975))))
CMP <- list(`V5nf: oldfix - single` = c("V5nf.oldfix", "V5nf.single"),
            `F3 - V5nf (single)` = c("F3.single", "V5nf.single"),
            `F3 - V5nf (oldfix)` = c("F3.oldfix", "V5nf.oldfix"))
key_budgets <- c(Bang_nf = wmean(dX$ref_bang_nf, dX$wt), V5nf_R0 = wmean(dX$ref_V5nf_R0, dX$wt),
                 V5nf_R1 = wmean(dX$ref_V5nf_R1, dX$wt))
ip <- function(v, b0) approx(BUDGETS, v, xout = b0)$y
diff_tab <- bind_rows(lapply(names(CMP), function(nm) {
  a <- CMP[[nm]]; db <- bt[, a[1], ] - bt[, a[2], ]
  data.frame(comparison = nm, budget = BUDGETS, diff = est[, a[1]] - est[, a[2]],
             lower95 = apply(db, 1, qna, pr = 0.025), upper95 = apply(db, 1, qna, pr = 0.975))
}))
key_tab <- bind_rows(lapply(names(CMP), function(nm) bind_rows(lapply(names(key_budgets), function(kb) {
  a <- CMP[[nm]]; b0 <- key_budgets[kb]
  db <- apply(bt[, a[1], ] - bt[, a[2], ], 2, ip, b0 = b0)
  data.frame(comparison = nm, at = kb, budget = b0, diff = ip(est[, a[1]] - est[, a[2]], b0),
             lower95 = qna(db, 0.025), upper95 = qna(db, 0.975), row.names = NULL)
}))))


# =====================================================================
# 第 8 部分：图
# =====================================================================
band <- function(x, lo, hi, col) {
  k <- !is.na(lo) & !is.na(hi)
  if (any(k)) polygon(c(x[k], rev(x[k])), c(lo[k], rev(hi[k])), col = adjustcolor(col, 0.2), border = NA)
}
png(file.path(OUT_DIR, "ext_curves.png"), width = 2600, height = 1100, res = 200)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
pt_of <- function(r) { x <- results[results$period == "ext" & results$group == "overall" & results$strategy == r, ]
c(x$estimate[x$measure == "referral"], x$estimate[x$measure == "sensitivity"]) }
for (m in c("V5nf", "F3")) {
  plot(NA, xlim = 100 * range(BUDGETS), ylim = range(100 * curve_tab$sens, na.rm = TRUE) + c(-5, 3),
       xlab = "Overall referral (%)", ylab = "Overall sensitivity (%)", main = paste0(m, ", 2021-2023 (external)"))
  cols <- c(single = "#2166ac", oldfix = "#b2182b")
  for (f in names(cols)) {
    r <- curve_tab[curve_tab$curve == paste(m, f, sep = "."), ]
    band(100 * r$budget, 100 * r$lower95, 100 * r$upper95, cols[f]); lines(100 * r$budget, 100 * r$sens, lwd = 2, col = cols[f])
  }
  pts <- c("Bang_nf", "USPSTF2021", paste0(m, c("_R0", "_R1")), if (m == "V5nf") "V5nf_R4")
  pc <- c("grey40", "black", "#2166ac", "#b2182b", "#d95f02")[seq_along(pts)]
  for (i in seq_along(pts)) { q <- pt_of(pts[i]); points(q[1], q[2], pch = 19, col = pc[i], cex = 1.2) }
  legend("bottomright", c("single cut-off", "60+ cut-off fixed", pts), col = c(cols, pc),
         lty = c(1, 1, rep(NA, length(pts))), pch = c(NA, NA, rep(19, length(pts))), lwd = 2, bty = "n", cex = 0.7)
}
invisible(dev.off())

png(file.path(OUT_DIR, "ext_age_sensitivity.png"), width = 1800, height = 1100, res = 200)
par(mar = c(4.5, 4.5, 3, 1))
gs <- c("overall", "20-39", "40-59", ">=60", "high_access", "low_access")
sts <- c("Bang_nf", "V5nf_R0", "V5nf_R1", "F3_R1")
M <- L <- U <- matrix(NA, length(sts), length(gs), dimnames = list(sts, gs))
for (s in sts) for (g in gs) {
  x <- results[results$period == "ext" & results$group == g & results$strategy == s & results$measure == "sensitivity", ]
  if (nrow(x)) { M[s, g] <- x$estimate; L[s, g] <- x$lower95; U[s, g] <- x$upper95 }
}
bp <- barplot(M, beside = TRUE, ylim = c(0, 105), col = c("grey70", "#2166ac", "#b2182b", "#d95f02"), border = NA,
              ylab = "Sensitivity (%)", main = "2021-2023 external validation: sensitivity by group")
suppressWarnings(arrows(bp, L, bp, U, angle = 90, code = 3, length = 0.02, col = "grey30"))
legend("topleft", sts, fill = c("grey70", "#2166ac", "#b2182b", "#d95f02"), border = NA, bty = "n", cex = 0.8)
invisible(dev.off())


# =====================================================================
# 第 9 部分：保存和打印
# =====================================================================
write.csv(results,   file.path(OUT_DIR, "ext_results.csv"), row.names = FALSE)
write.csv(gap_tab,   file.path(OUT_DIR, "ext_access_gaps.csv"), row.names = FALSE)
write.csv(perf_tab,  file.path(OUT_DIR, "ext_auc_calibration.csv"), row.names = FALSE)
write.csv(curve_tab, file.path(OUT_DIR, "ext_family_curves.csv"), row.names = FALSE)
write.csv(diff_tab,  file.path(OUT_DIR, "ext_family_differences.csv"), row.names = FALSE)
write.csv(key_tab,   file.path(OUT_DIR, "ext_differences_at_key_budgets.csv"), row.names = FALSE)
write.csv(paq_tab,   file.path(OUT_DIR, "ext_inactive_definition_check.csv"), row.names = FALSE)
write.csv(data.frame(model = names(THR), cut_single = sapply(THR, `[`, "single"), cut_60plus = sapply(THR, `[`, "old")),
          file.path(OUT_DIR, "thresholds_locked.csv"), row.names = FALSE)
writeLines(c("方程和切点只用 2007–2012 确定；2021–2023 未参与任何拟合或选择，是独立的时间外推验证。",
             paste("外部验证样本：", sum(full$ext), "人，病例", sum(full$case[full$ext])),
             paste("bootstrap 次数：", B_BOOT), paste("运行时间：", Sys.time())),
           file.path(OUT_DIR, "ext_info.txt"))
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))

fmt <- function(d) sprintf("%.1f (%.1f-%.1f)", d$estimate, d$lower95, d$upper95)
cat("\n===== 不活动定义检查（各时期“不活动”加权比例）=====\n"); print(as.data.frame(paq_tab), digits = 3, row.names = FALSE)
cat("\n===== AUC 与校准 =====\n"); print(perf_tab, digits = 3, row.names = FALSE)
for (pe in periods) {
  cat("\n=====", ifelse(pe == "eval", "2013–2018（参照）", "2021–2023 外部验证"), ": 转检 / 敏感度 % =====\n")
  r <- results[results$period == pe & results$measure %in% c("referral", "sensitivity"), ]
  r$value <- fmt(r); r$key <- paste(r$measure, r$group)
  out <- reshape(r[, c("strategy", "key", "value")], idvar = "strategy", timevar = "key", direction = "wide")
  names(out) <- sub("^value\\.", "", names(out)); print(out, row.names = FALSE)
  cat("-- 与 Bang 的敏感度差（总体 / >=60）--\n")
  d <- results[results$period == pe & results$measure %in% c("sens_diff_vs_Bang", "sens_diff_vs_Bang_nf") &
                 results$group %in% c("overall", ">=60"), ]
  d$value <- fmt(d); print(d[, c("strategy", "measure", "group", "value")], row.names = FALSE)
  cat("-- 可及性差距 --\n")
  x <- gap_tab[gap_tab$period == pe, ]; x$value <- fmt(x)
  out <- reshape(x[, c("strategy", "measure", "value")], idvar = "strategy", timevar = "measure", direction = "wide")
  names(out) <- sub("^value\\.", "", names(out)); print(out, row.names = FALSE)
}
cat("\n===== 2021–2023 同名额下的敏感度差（百分点，bootstrap 95% 区间）=====\n")
key_tab$value <- sprintf("%+.1f (%+.1f to %+.1f)", 100 * key_tab$diff, 100 * key_tab$lower95, 100 * key_tab$upper95)
print(key_tab[, c("comparison", "at", "budget", "value")], digits = 3, row.names = FALSE)
cat("\n结果已保存到：", OUT_DIR, "\n")
