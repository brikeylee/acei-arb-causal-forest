# clinical_decision_support.R - 临床决策支持工具模块
# 项目：围术期 ACEI/ARB 停药策略因果森林分析
#
# 改动说明：
#   1. 修正变量名（无 smoking/previous_mi/gender）
#   2. 临床评分添加交叉验证
#   3. 校准分析修正

cat("Clinical decision support functions loaded.\n")

# ── 1. 临床评分系统（含 CV 验证）──
create_clinical_score <- function(matched_data,
                                  cate_col = "tau_hat") {
  cat("  Building clinical decision score...\n")

  scored <- matched_data %>%
    filter(!is.na(.data[[cate_col]])) %>%
    mutate(
      age_score  = case_when(age < 60 ~ 0L, age < 70 ~ 1L, age < 80 ~ 2L, TRUE ~ 3L),
      egfr_score = case_when(egfr >= 90 ~ 0L, egfr >= 60 ~ 1L, egfr >= 30 ~ 2L, TRUE ~ 3L),
      dm_score   = ifelse(as.numeric(as.character(diabetes)) == 1, 2L, 0L),
      htn_score  = ifelse(as.numeric(as.character(hypertension)) == 1, 1L, 0L),
      bmi_score  = case_when(bmi < 25 ~ 0L, bmi < 30 ~ 1L, TRUE ~ 2L),
      total_score = age_score + egfr_score + dm_score + htn_score + bmi_score
    )

  # 评分性能
  score_perf <- scored %>%
    group_by(total_score) %>%
    summarise(
      n = n(),
      mean_cate  = mean(.data[[cate_col]], na.rm = TRUE),
      se_cate    = sd(.data[[cate_col]], na.rm = TRUE) / sqrt(n()),
      ci_lower   = mean_cate - 1.96 * se_cate,
      ci_upper   = mean_cate + 1.96 * se_cate,
      .groups = "drop"
    ) %>%
    mutate(
      confidence = case_when(
        total_score >= 6 ~ "High",
        total_score >= 3 ~ "Moderate",
        TRUE ~ "Low"
      )
    )

  # 5-fold CV 验证评分预测能力
  cat("  Cross-validating score system (5-fold)...\n")
  set.seed(123)
  folds <- sample(rep(1:5, length.out = nrow(scored)))
  cv_cors <- numeric(5)

  for (f in 1:5) {
    train_idx <- folds != f
    test_idx  <- folds == f
    # 在训练集上拟合 score ~ cate 关系
    mod <- lm(as.formula(paste(cate_col, "~ total_score")), data = scored[train_idx, ])
    pred <- predict(mod, newdata = scored[test_idx, ])
    cv_cors[f] <- cor(pred, scored[[cate_col]][test_idx], use = "complete.obs")
  }

  cat("  CV correlation (score vs CATE):", round(mean(cv_cors), 3),
      "+/-", round(sd(cv_cors), 3), "\n")

  cv_result <- data.frame(
    fold = 1:5,
    correlation = cv_cors,
    row.names = NULL
  )

  list(
    scored_data      = scored,
    score_performance = score_perf,
    cv_validation     = cv_result,
    cv_mean_cor       = mean(cv_cors)
  )
}

# ── 2. 决策树 ──
build_decision_tree <- function(matched_data,
                                 cate_col    = "tau_hat",
                                 min_split   = 100,
                                 tree_vars   = c("age", "egfr", "bmi", "sex",
                                                 "diabetes", "hypertension",
                                                 "heart_failure", "ckd")) {
  cat("  Building decision tree...\n")

  tree_vars <- intersect(tree_vars, names(matched_data))
  tree_data <- matched_data[!is.na(matched_data[[cate_col]]), c(cate_col, tree_vars)]
  tree_data <- tree_data[complete.cases(tree_data), ]

  # 将 CATE 分为三类
  q33 <- quantile(tree_data[[cate_col]], 1/3, na.rm = TRUE)
  q67 <- quantile(tree_data[[cate_col]], 2/3, na.rm = TRUE)
  tree_data$benefit_group <- factor(
    case_when(
      tree_data[[cate_col]] >= q67 ~ "High_Benefit",
      tree_data[[cate_col]] >= q33 ~ "Moderate_Benefit",
      TRUE ~ "Lower_Benefit"
    ),
    levels = c("Lower_Benefit", "Moderate_Benefit", "High_Benefit")
  )

  formula_str <- paste("benefit_group ~", paste(tree_vars, collapse = " + "))
  tree_model <- rpart(
    as.formula(formula_str),
    data    = tree_data,
    method  = "class",
    control = rpart.control(minsplit = min_split, cp = 0.01, maxdepth = 4)
  )

  list(tree_model = tree_model, tree_data = tree_data)
}

# ── 3. 校准分析（修正版）──
calibration_analysis <- function(cf_model, matched_data,
                                  cate_col      = "tau_hat",
                                  outcome_col   = "adverse_event",
                                  treatment_col = "preop_discontinuation",
                                  n_bins        = 10) {
  cat("  Running calibration analysis...\n")

  tau <- matched_data[[cate_col]]
  Y   <- as.numeric(as.character(matched_data[[outcome_col]]))
  W   <- as.numeric(as.character(matched_data[[treatment_col]]))

  valid <- !is.na(tau) & !is.na(Y) & !is.na(W)
  tau <- tau[valid]; Y <- Y[valid]; W <- W[valid]

  # 按 CATE 分位数分组
  breaks <- quantile(tau, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE)
  breaks <- unique(breaks)  # 避免重复断点
  bins <- cut(tau, breaks = breaks, include.lowest = TRUE, labels = FALSE)

  cal_data <- data.frame(tau = tau, Y = Y, W = W, bin = bins)

  cal_summary <- cal_data %>%
    group_by(bin) %>%
    summarise(
      n = n(),
      predicted_cate = mean(tau),
      # 组内观察 ATE：treated 均值 - control 均值
      n_treated  = sum(W == 1),
      n_control  = sum(W == 0),
      rate_treated = mean(Y[W == 1]),
      rate_control = mean(Y[W == 0]),
      observed_cate = rate_treated - rate_control,
      calibration_error = predicted_cate - observed_cate,
      .groups = "drop"
    )

  # 校准斜率
  cal_slope <- tryCatch({
    mod <- lm(observed_cate ~ predicted_cate, data = cal_summary, weights = n)
    coef(mod)["predicted_cate"]
  }, error = function(e) NA)

  # 校准图
  p_cal <- ggplot(cal_summary, aes(x = predicted_cate, y = observed_cate)) +
    geom_point(aes(size = n), color = "#00468B", alpha = 0.8) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
    geom_smooth(method = "lm", se = TRUE, color = "#ED0000", alpha = 0.2) +
    labs(title = "Calibration: Predicted vs Observed CATE",
         subtitle = paste("Calibration slope =", round(cal_slope, 3)),
         x = "Predicted CATE", y = "Observed CATE", size = "N patients") +
    theme_publication()

  list(
    calibration_data  = cal_summary,
    calibration_slope = cal_slope,
    calibration_plot  = p_cal
  )
}

# ── 4. NNT/NNH 计算 ──
calculate_nnt_nnh <- function(matched_data,
                               outcome_col   = "adverse_event",
                               treatment_col = "preop_discontinuation",
                               subgroup_var  = NULL) {
  cat("  Calculating NNT/NNH...\n")

  calc_one <- function(data, label = "Overall") {
    W <- as.numeric(as.character(data[[treatment_col]]))
    Y <- as.numeric(as.character(data[[outcome_col]]))
    r1 <- mean(Y[W == 1], na.rm = TRUE)
    r0 <- mean(Y[W == 0], na.rm = TRUE)
    rd <- r1 - r0
    nnt <- ifelse(rd > 0, 1 / rd, NA)  # 停用增加风险 → NNT to harm
    nnh <- ifelse(rd < 0, 1 / abs(rd), NA)  # 停用减少风险 → NNH
    data.frame(subgroup = label, n = length(Y),
               risk_discont = r1, risk_contin = r0,
               risk_diff = rd, NNT = nnt, NNH = nnh,
               row.names = NULL)
  }

  results <- list(calc_one(matched_data, "Overall"))

  if (!is.null(subgroup_var) && subgroup_var %in% names(matched_data)) {
    for (lev in unique(matched_data[[subgroup_var]])) {
      if (is.na(lev)) next
      sub <- matched_data[matched_data[[subgroup_var]] == lev & !is.na(matched_data[[subgroup_var]]), ]
      if (nrow(sub) >= 50) {
        results[[length(results) + 1]] <- calc_one(sub, paste0(subgroup_var, "=", lev))
      }
    }
  }

  do.call(rbind, results)
}
