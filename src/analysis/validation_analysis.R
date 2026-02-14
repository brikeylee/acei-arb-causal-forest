# validation_analysis.R - 验证和稳健性分析模块
# 项目：围术期 ACEI/ARB 停药策略因果森林分析
#
# 改动说明：
#   1. Bootstrap 函数使用实际列名参数化
#   2. 交叉验证改用 RATE (Rank-Weighted ATE) 替代循环论证
#   3. 超参数敏感性分析保留
#   4. 校准分析修正：使用 grf 内置方法

cat("\n=== Phase 4: Validation Analysis ===\n")

# ── 1. Bootstrap 验证（修复版）──
bootstrap_validation <- function(matched_data,
                                 outcome_col   = "adverse_event",
                                 treatment_col = "preop_discontinuation",
                                 cate_col      = "tau_hat",
                                 covariates    = NULL,
                                 n_bootstrap   = 500,
                                 seed          = 123) {
  cat("  Running bootstrap validation (n =", n_bootstrap, ")...\n")
  set.seed(seed)

  if (is.null(covariates)) {
    covariates <- c("age", "sex", "center", "heart_failure", "hypertension",
                    "diabetes", "ckd", "respiratory_disease", "cpb", "egfr",
                    "bmi", "ccb", "statin", "antiplatelet", "beta_blocker")
  }
  covariates <- intersect(covariates, names(matched_data))

  # 准备完整数据
  all_cols <- c(covariates, outcome_col, treatment_col)
  if (cate_col %in% names(matched_data)) {
    work_data <- matched_data[complete.cases(matched_data[, c(all_cols, cate_col)]), ]
  } else {
    work_data <- matched_data[complete.cases(matched_data[, all_cols]), ]
  }

  Y <- as.numeric(as.character(work_data[[outcome_col]]))
  W <- as.numeric(as.character(work_data[[treatment_col]]))
  X <- model.matrix(~ . - 1, data = work_data[, covariates, drop = FALSE])

  # 原始 ATE
  cf_original <- causal_forest(X, Y, W, num.trees = 2000,
                                honesty = TRUE, tune.parameters = "all", seed = seed)
  original_ate <- average_treatment_effect(cf_original)
  original_estimate <- original_ate["estimate"]

  # Bootstrap
  boot_ates <- numeric(n_bootstrap)
  for (i in seq_len(n_bootstrap)) {
    idx <- sample(nrow(work_data), replace = TRUE)
    cf_b <- causal_forest(X[idx, ], Y[idx], W[idx],
                          num.trees = 2000, honesty = TRUE,
                          tune.parameters = "all")
    boot_ates[i] <- average_treatment_effect(cf_b)["estimate"]
    if (i %% 100 == 0) cat("    Bootstrap iteration", i, "/", n_bootstrap, "\n")
  }

  boot_summary <- data.frame(
    metric = c("original_estimate", "bootstrap_mean", "bootstrap_se",
               "ci_lower", "ci_upper", "bias"),
    value = c(
      original_estimate,
      mean(boot_ates),
      sd(boot_ates),
      quantile(boot_ates, 0.025),
      quantile(boot_ates, 0.975),
      mean(boot_ates) - original_estimate
    ),
    row.names = NULL
  )

  cat("  Bootstrap ATE: mean =", round(mean(boot_ates), 4),
      ", 95% CI [", round(quantile(boot_ates, 0.025), 4), ",",
      round(quantile(boot_ates, 0.975), 4), "]\n")

  return(list(
    summary     = boot_summary,
    boot_ates   = boot_ates,
    original    = original_estimate,
    cf_original = cf_original
  ))
}

# ── 2. RATE 异质性验证（替代循环交叉验证）──
rate_heterogeneity_test <- function(cf_model, matched_data, cate_col = "tau_hat") {
  cat("  Running RATE (Rank-weighted ATE) test...\n")

  tau_hat <- matched_data[[cate_col]]
  valid_idx <- !is.na(tau_hat)

  # 使用 grf 内置的 RATE
  rate_result <- tryCatch({
    rank_average_treatment_effect(cf_model, tau_hat[valid_idx])
  }, error = function(e) {
    cat("  RATE not available, using test_calibration instead.\n")
    NULL
  })

  # test_calibration 作为补充
  cal_test <- tryCatch(test_calibration(cf_model), error = function(e) NULL)

  list(rate = rate_result, calibration_test = cal_test)
}

# ── 3. 超参数敏感性分析 ──
hyperparameter_sensitivity <- function(matched_data,
                                       outcome_col   = "adverse_event",
                                       treatment_col = "preop_discontinuation",
                                       covariates    = NULL,
                                       num_trees_range     = c(1000, 2000, 4000),
                                       min_node_size_range = c(5, 10, 20)) {
  cat("  Running hyperparameter sensitivity...\n")

  if (is.null(covariates)) {
    covariates <- c("age", "sex", "center", "heart_failure", "hypertension",
                    "diabetes", "ckd", "respiratory_disease", "cpb", "egfr",
                    "bmi", "ccb", "statin", "antiplatelet", "beta_blocker")
  }
  covariates <- intersect(covariates, names(matched_data))

  mats <- prepare_cf_matrices(matched_data, outcome_col, treatment_col, covariates)

  params <- expand.grid(num_trees = num_trees_range, min_node_size = min_node_size_range)
  results <- list()

  for (i in seq_len(nrow(params))) {
    cf_i <- causal_forest(
      mats$X, mats$Y, mats$W,
      num.trees     = params$num_trees[i],
      min.node.size = params$min_node_size[i],
      honesty       = TRUE,
      tune.parameters = "all",
      seed          = 123
    )
    ate_i <- average_treatment_effect(cf_i)
    results[[i]] <- data.frame(
      num_trees     = params$num_trees[i],
      min_node_size = params$min_node_size[i],
      ATE           = ate_i["estimate"],
      SE            = ate_i["std.err"],
      row.names     = NULL
    )
  }

  results_df <- do.call(rbind, results)
  results_df$CI_lower <- results_df$ATE - 1.96 * results_df$SE
  results_df$CI_upper <- results_df$ATE + 1.96 * results_df$SE

  cat("  Hyperparameter sensitivity results:\n")
  print(results_df)

  return(results_df)
}

# ── 4. 变量重要性稳定性 ──
variable_importance_stability <- function(matched_data,
                                          outcome_col   = "adverse_event",
                                          treatment_col = "preop_discontinuation",
                                          covariates    = NULL,
                                          n_iterations  = 100,
                                          sample_prop   = 0.8) {
  cat("  Running variable importance stability (n =", n_iterations, ")...\n")

  if (is.null(covariates)) {
    covariates <- c("age", "sex", "center", "heart_failure", "hypertension",
                    "diabetes", "ckd", "respiratory_disease", "cpb", "egfr",
                    "bmi", "ccb", "statin", "antiplatelet", "beta_blocker")
  }
  covariates <- intersect(covariates, names(matched_data))

  mats <- prepare_cf_matrices(matched_data, outcome_col, treatment_col, covariates)
  imp_matrix <- matrix(NA, nrow = n_iterations, ncol = ncol(mats$X))
  colnames(imp_matrix) <- mats$covariate_names

  set.seed(123)
  for (i in seq_len(n_iterations)) {
    idx <- sample(length(mats$Y), size = floor(length(mats$Y) * sample_prop))
    cf_i <- causal_forest(mats$X[idx, ], mats$Y[idx], mats$W[idx],
                          num.trees = 1000, seed = i)
    imp_matrix[i, ] <- variable_importance(cf_i)
  }

  stability <- data.frame(
    variable        = mats$covariate_names,
    mean_importance = colMeans(imp_matrix, na.rm = TRUE),
    sd_importance   = apply(imp_matrix, 2, sd, na.rm = TRUE),
    row.names       = NULL
  )
  stability$cv <- stability$sd_importance / pmax(stability$mean_importance, 1e-10)
  stability <- stability[order(-stability$mean_importance), ]

  return(stability)
}

cat("Validation analysis functions loaded.\n")
