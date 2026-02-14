# functions.R - 辅助函数定义
# 项目：围术期 ACEI/ARB 停药策略因果森林分析

# ── 数据读取 ──
read_data <- function(paths) {
  for (path in paths) {
    if (file.exists(path)) {
      cat("Data file found:", path, "\n")
      return(readxl::read_excel(path))
    }
  }
  stop("Cannot find data file. Tried paths:\n", paste(paths, collapse = "\n"))
}

# ── 列名重命名 ──
rename_columns <- function(data, name_mapping) {
  for (old_name in names(name_mapping)) {
    if (old_name %in% names(data)) {
      names(data)[names(data) == old_name] <- name_mapping[old_name]
    } else {
      cat("Warning: Column '", old_name, "' not found in data.\n")
    }
  }
  return(data)
}

# ── 变量类型转换 ──
convert_variable_types <- function(data, continuous_vars, categorical_vars) {
  for (var in continuous_vars) {
    if (var %in% names(data)) {
      data[[var]] <- as.numeric(as.character(data[[var]]))
    }
  }
  for (var in categorical_vars) {
    if (var %in% names(data)) {
      data[[var]] <- as.factor(data[[var]])
    }
  }
  return(data)
}

# ── RR / RD / OR 计算（含 95% CI）──
calculate_rr_rd <- function(data, outcome_var, exposure_var) {
  n_exp <- sum(data[[exposure_var]] == 1, na.rm = TRUE)
  n_ctl <- sum(data[[exposure_var]] == 0, na.rm = TRUE)
  ev_exp <- sum(data[[exposure_var]] == 1 & data[[outcome_var]] == 1, na.rm = TRUE)
  ev_ctl <- sum(data[[exposure_var]] == 0 & data[[outcome_var]] == 1, na.rm = TRUE)
  ne_exp <- n_exp - ev_exp
  ne_ctl <- n_ctl - ev_ctl

  r_exp <- ev_exp / n_exp
  r_ctl <- ev_ctl / n_ctl

  # Relative Risk
  rr <- r_exp / r_ctl
  rr_se <- sqrt((1 - r_exp) / ev_exp + (1 - r_ctl) / ev_ctl)
  rr_lo <- exp(log(rr) - 1.96 * rr_se)
  rr_hi <- exp(log(rr) + 1.96 * rr_se)

  # Risk Difference (per 1000)
  rd <- (r_exp - r_ctl) * 1000
  rd_se <- sqrt(r_exp * (1 - r_exp) / n_exp + r_ctl * (1 - r_ctl) / n_ctl)
  rd_lo <- rd - 1.96 * rd_se * 1000
  rd_hi <- rd + 1.96 * rd_se * 1000

  # Odds Ratio
  or_val <- (ev_exp / ne_exp) / (ev_ctl / ne_ctl)
  or_se <- sqrt(1 / ev_exp + 1 / ne_exp + 1 / ev_ctl + 1 / ne_ctl)
  or_lo <- exp(log(or_val) - 1.96 * or_se)
  or_hi <- exp(log(or_val) + 1.96 * or_se)

  list(
    n_exposed = n_exp, n_control = n_ctl,
    events_exposed = ev_exp, events_control = ev_ctl,
    risk_exposed = r_exp, risk_control = r_ctl,
    rr = rr, rr_ci_lower = rr_lo, rr_ci_upper = rr_hi,
    rd = rd, rd_ci_lower = rd_lo, rd_ci_upper = rd_hi,
    or = or_val, or_ci_lower = or_lo, or_ci_upper = or_hi
  )
}

# ── 因果森林 ATE 提取（含 CI 和 p-value）──
extract_ate_results <- function(cf_model, outcome_name = "adverse_event") {
  ate <- average_treatment_effect(cf_model)
  est <- ate["estimate"]
  se  <- ate["std.err"]
  z   <- est / se
  p   <- 2 * pnorm(-abs(z))
  data.frame(
    outcome  = outcome_name,
    ATE      = est,
    SE       = se,
    CI_lower = est - 1.96 * se,
    CI_upper = est + 1.96 * se,
    z_value  = z,
    p_value  = p,
    row.names = NULL
  )
}

# ── 准备因果森林输入矩阵 ──
prepare_cf_matrices <- function(data,
                                outcome_col = "adverse_event",
                                treatment_col = "preop_discontinuation",
                                covariate_cols = NULL) {
  if (is.null(covariate_cols)) {
    covariate_cols <- c("age", "sex", "center", "heart_failure", "hypertension",
                        "diabetes", "ckd", "respiratory_disease", "cpb", "egfr",
                        "bmi", "ccb", "statin", "antiplatelet", "beta_blocker")
  }
  # 只保留数据中实际存在的协变量
  covariate_cols <- intersect(covariate_cols, names(data))

  subset_data <- data[, c(covariate_cols, outcome_col, treatment_col)]
  complete_idx <- complete.cases(subset_data)
  subset_data <- subset_data[complete_idx, ]

  Y <- as.numeric(as.character(subset_data[[outcome_col]]))
  W <- as.numeric(as.character(subset_data[[treatment_col]]))
  X <- model.matrix(~ . - 1, data = subset_data[, covariate_cols, drop = FALSE])

  list(Y = Y, W = W, X = X, data = subset_data, complete_idx = complete_idx,
       covariate_names = colnames(X))
}

cat("Functions loaded.\n")
