# sensitivity_analysis.R - 敏感性分析模块
# 项目：围术期 ACEI/ARB 停药策略因果森林分析
#
# 内容：
#   1. IPTW 分析（替代 PSM 的验证）
#   2. E-value 计算（未测量混杂）
#   3. 不同 Caliper 敏感性分析
#   4. 亚组交互检验

cat("\n=== Sensitivity Analyses ===\n")

# ── 1. IPTW 分析 ──
run_iptw_analysis <- function(data,
                              outcome_col   = "adverse_event",
                              treatment_col = "preop_discontinuation",
                              ps_formula    = NULL) {
  cat("  Running IPTW analysis...\n")

  if (is.null(ps_formula)) {
    ps_formula <- as.formula(paste(
      treatment_col, "~ age + center + egfr + sex + surgery_type +",
      "hypertension + diabetes + heart_failure + ckd +",
      "preop_cr + bmi + preop_hb +",
      "ccb + statin + antiplatelet + beta_blocker +",
      "cpb + los + preop_bun + preop_plt + preop_alb"
    ))
  }

  # 使用 WeightIt 计算 IPTW 权重
  W_num <- as.numeric(as.character(data[[treatment_col]]))
  Y_num <- as.numeric(as.character(data[[outcome_col]]))

  # 倾向性评分
  ps_model <- glm(ps_formula, data = data, family = binomial())
  ps <- ps_model$fitted.values

  # IPTW 权重（stabilized）
  p_treatment <- mean(W_num, na.rm = TRUE)
  weights <- ifelse(W_num == 1,
                    p_treatment / ps,
                    (1 - p_treatment) / (1 - ps))

  # 截断极端权重（1st and 99th percentile）
  w_lo <- quantile(weights, 0.01, na.rm = TRUE)
  w_hi <- quantile(weights, 0.99, na.rm = TRUE)
  weights_trimmed <- pmin(pmax(weights, w_lo), w_hi)

  # IPTW 加权 logistic 回归
  iptw_model <- glm(Y_num ~ W_num, weights = weights_trimmed, family = quasibinomial())

  coef_summary <- summary(iptw_model)$coefficients
  or_iptw    <- exp(coef_summary["W_num", "Estimate"])
  or_se      <- coef_summary["W_num", "Std. Error"]
  or_lo      <- exp(coef_summary["W_num", "Estimate"] - 1.96 * or_se)
  or_hi      <- exp(coef_summary["W_num", "Estimate"] + 1.96 * or_se)
  p_iptw     <- coef_summary["W_num", "Pr(>|t|)"]

  # Robust SE (sandwich estimator)
  robust_se <- tryCatch({
    sqrt(sandwich::vcovHC(iptw_model, type = "HC1")["W_num", "W_num"])
  }, error = function(e) or_se)

  or_robust_lo <- exp(coef_summary["W_num", "Estimate"] - 1.96 * robust_se)
  or_robust_hi <- exp(coef_summary["W_num", "Estimate"] + 1.96 * robust_se)

  # Risk difference via weighted means
  rd_iptw <- weighted.mean(Y_num[W_num == 1], weights_trimmed[W_num == 1]) -
             weighted.mean(Y_num[W_num == 0], weights_trimmed[W_num == 0])

  result <- data.frame(
    method       = "IPTW (stabilized, trimmed)",
    outcome      = outcome_col,
    OR           = or_iptw,
    OR_CI_lower  = or_robust_lo,
    OR_CI_upper  = or_robust_hi,
    RD           = rd_iptw,
    p_value      = p_iptw,
    row.names    = NULL
  )

  cat("  IPTW OR:", round(or_iptw, 3),
      "[", round(or_robust_lo, 3), "-", round(or_robust_hi, 3), "]",
      "p =", round(p_iptw, 4), "\n")

  return(result)
}

# ── 2. E-value 计算 ──
calculate_evalue <- function(or_estimate, or_ci_lower, or_ci_upper, rare_outcome = FALSE) {
  cat("  Calculating E-value...\n")

  # 如果结局不罕见，将 OR 转换为近似 RR
  if (!rare_outcome) {
    # 使用 Zhang and Yu (1998) 转换: RR = OR / (1 - p0 + p0 * OR)
    # 简化：对于中等 OR，可以直接用 EValue 包
    # 这里使用 EValue 包的方法
  }

  tryCatch({
    ev <- EValue::evalues.OR(est = or_estimate, lo = or_ci_lower, hi = or_ci_upper,
                             rare = rare_outcome)
    result <- data.frame(
      measure     = c("point_estimate", "CI_lower"),
      OR          = c(or_estimate, or_ci_lower),
      E_value     = c(ev[2, "E-values"], ev[3, "E-values"]),
      row.names   = NULL
    )
    cat("  E-value for point estimate:", round(ev[2, "E-values"], 2), "\n")
    cat("  E-value for CI bound:", round(ev[3, "E-values"], 2), "\n")
    return(result)
  }, error = function(e) {
    cat("  E-value calculation failed:", e$message, "\n")
    # 手动计算 E-value: E = OR + sqrt(OR * (OR - 1))
    if (or_estimate >= 1) {
      ev_point <- or_estimate + sqrt(or_estimate * (or_estimate - 1))
    } else {
      rr_inv <- 1 / or_estimate
      ev_point <- rr_inv + sqrt(rr_inv * (rr_inv - 1))
    }
    data.frame(measure = "point_estimate", OR = or_estimate,
               E_value = ev_point, row.names = NULL)
  })
}

# ── 3. 不同 Caliper 敏感性分析 ──
caliper_sensitivity <- function(data,
                                outcome_col   = "adverse_event",
                                treatment_col = "preop_discontinuation",
                                ps_formula    = NULL,
                                calipers      = c(0.05, 0.1, 0.15, 0.2)) {
  cat("  Running caliper sensitivity analysis...\n")

  if (is.null(ps_formula)) {
    ps_formula <- as.formula(paste(
      treatment_col, "~ age + center + egfr + sex + surgery_type +",
      "hypertension + diabetes + heart_failure + ckd +",
      "preop_cr + bmi + preop_hb +",
      "ccb + statin + antiplatelet + beta_blocker +",
      "cpb + los + preop_bun + preop_plt + preop_alb"
    ))
  }

  results <- list()
  for (cal in calipers) {
    cat("    Caliper =", cal, "... ")
    m <- tryCatch({
      matchit(ps_formula, data = data, method = "nearest",
              ratio = 1, caliper = cal, replace = FALSE)
    }, error = function(e) NULL)

    if (!is.null(m)) {
      md <- match.data(m)
      n_matched <- nrow(md)

      # 计算效应
      W_m <- as.numeric(as.character(md[[treatment_col]]))
      Y_m <- as.numeric(as.character(md[[outcome_col]]))
      rd <- mean(Y_m[W_m == 1]) - mean(Y_m[W_m == 0])

      mod <- glm(Y_m ~ W_m, family = binomial())
      or_val <- exp(coef(mod)["W_m"])
      or_se  <- summary(mod)$coefficients["W_m", "Std. Error"]

      results[[as.character(cal)]] <- data.frame(
        caliper    = cal,
        n_matched  = n_matched,
        RD         = rd,
        OR         = or_val,
        OR_CI_lower = exp(coef(mod)["W_m"] - 1.96 * or_se),
        OR_CI_upper = exp(coef(mod)["W_m"] + 1.96 * or_se),
        p_value    = summary(mod)$coefficients["W_m", "Pr(>|z|)"],
        row.names  = NULL
      )
      cat("n =", n_matched, ", OR =", round(or_val, 3), "\n")
    } else {
      cat("FAILED\n")
    }
  }

  do.call(rbind, results)
}

# ── 4. 亚组交互检验 ──
subgroup_interaction_tests <- function(matched_data,
                                       outcome_col   = "adverse_event",
                                       treatment_col = "preop_discontinuation",
                                       subgroup_vars = c("age_group", "egfr_group",
                                                         "diabetes", "hypertension",
                                                         "heart_failure", "ckd", "sex")) {
  cat("  Running subgroup interaction tests...\n")

  W <- as.numeric(as.character(matched_data[[treatment_col]]))
  Y <- as.numeric(as.character(matched_data[[outcome_col]]))

  results <- list()
  for (sv in subgroup_vars) {
    if (!sv %in% names(matched_data)) {
      cat("    Skipping", sv, "(not found)\n")
      next
    }

    subgroup <- matched_data[[sv]]

    # 主效应模型
    mod_main <- glm(Y ~ W, family = binomial())

    # 交互模型
    mod_int <- tryCatch({
      glm(Y ~ W * subgroup, family = binomial())
    }, error = function(e) NULL)

    if (is.null(mod_int)) {
      cat("    Skipping", sv, "(model failed)\n")
      next
    }

    # LRT for interaction
    lrt <- anova(mod_main, mod_int, test = "Chisq")
    p_interaction <- lrt$`Pr(>Chi)`[2]

    # 每个亚组内的效应
    levels_sv <- unique(subgroup)
    subgroup_effects <- list()
    for (lev in levels_sv) {
      idx <- subgroup == lev & !is.na(subgroup)
      if (sum(idx) < 20) next
      mod_sub <- tryCatch({
        glm(Y[idx] ~ W[idx], family = binomial())
      }, error = function(e) NULL)

      if (!is.null(mod_sub)) {
        coefs <- summary(mod_sub)$coefficients
        if ("W[idx]" %in% rownames(coefs)) {
          subgroup_effects[[as.character(lev)]] <- data.frame(
            subgroup_var   = sv,
            subgroup_level = as.character(lev),
            n              = sum(idx),
            OR             = exp(coefs["W[idx]", "Estimate"]),
            OR_CI_lower    = exp(coefs["W[idx]", "Estimate"] - 1.96 * coefs["W[idx]", "Std. Error"]),
            OR_CI_upper    = exp(coefs["W[idx]", "Estimate"] + 1.96 * coefs["W[idx]", "Std. Error"]),
            p_value        = coefs["W[idx]", "Pr(>|z|)"],
            p_interaction  = p_interaction,
            row.names      = NULL
          )
        }
      }
    }

    if (length(subgroup_effects) > 0) {
      results[[sv]] <- do.call(rbind, subgroup_effects)
    }

    cat("    ", sv, ": p_interaction =", round(p_interaction, 4), "\n")
  }

  if (length(results) > 0) {
    combined <- do.call(rbind, results)
    # FDR 校正
    combined$p_interaction_fdr <- p.adjust(
      combined$p_interaction[!duplicated(combined$subgroup_var)],
      method = "fdr"
    )[match(combined$subgroup_var,
            unique(combined$subgroup_var))]
    return(combined)
  }

  return(NULL)
}

cat("Sensitivity analysis functions loaded.\n")
