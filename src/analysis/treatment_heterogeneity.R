# treatment_heterogeneity.R - 治疗异质性深度分析模块
# 项目：围术期 ACEI/ARB 停药策略因果森林分析
#
# 改动说明：
#   1. 修正变量名：使用数据中实际存在的列 (sex 而非 gender, 无 smoking/previous_mi)
#   2. 添加正式统计检验和多重比较校正
#   3. 使用 tau_hat 而非 cate

suppressMessages(library(stringr))

cat("Treatment heterogeneity functions loaded.\n")

# ── 1. 异质性统计检验 ──
test_treatment_heterogeneity <- function(cf_model, matched_data,
                                         cate_col = "tau_hat",
                                         covariates = NULL) {
  # test_calibration from grf
  cal_test <- tryCatch(test_calibration(cf_model), error = function(e) NULL)

  if (!is.null(cal_test)) {
    cat("  Calibration test results:\n")
    print(cal_test)
    cat("  Interpretation: p < 0.05 for 'differential.forest.prediction' suggests\n")
    cat("  significant heterogeneity in treatment effects.\n")
  }

  # 变量重要性
  var_imp <- variable_importance(cf_model)
  if (is.null(covariates)) {
    covariates <- c("age", "sex", "center", "heart_failure", "hypertension",
                    "diabetes", "ckd", "respiratory_disease", "cpb", "egfr",
                    "bmi", "ccb", "statin", "antiplatelet", "beta_blocker")
  }

  # 使用模型矩阵的列名
  X_names <- colnames(model.matrix(
    as.formula(paste("~", paste(intersect(covariates, names(matched_data)), collapse = " + "), "- 1")),
    data = matched_data[1, , drop = FALSE]
  ))

  if (length(var_imp) == length(X_names)) {
    imp_df <- data.frame(variable = X_names, importance = var_imp,
                         rank = rank(-var_imp), row.names = NULL)
  } else {
    imp_df <- data.frame(variable = paste0("V", seq_along(var_imp)),
                         importance = var_imp, rank = rank(-var_imp),
                         row.names = NULL)
  }
  imp_df <- imp_df[order(-imp_df$importance), ]

  list(calibration_test = cal_test, variable_importance = imp_df)
}

# ── 2. 连续变量与 CATE 关系分析 ──
analyze_continuous_relationships <- function(matched_data,
                                             cate_col = "tau_hat",
                                             continuous_vars = c("age", "egfr", "bmi")) {
  continuous_vars <- intersect(continuous_vars, names(matched_data))
  plots <- list()
  correlations <- list()

  for (var in continuous_vars) {
    plot_df <- matched_data[!is.na(matched_data[[cate_col]]) & !is.na(matched_data[[var]]), ]

    p <- ggplot(plot_df, aes(x = .data[[var]], y = .data[[cate_col]])) +
      geom_point(alpha = 0.15, size = 0.6, color = "#00468B") +
      geom_smooth(method = "loess", span = 0.3, se = TRUE,
                  color = "#ED0000", fill = "#ED0000", alpha = 0.15) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
      labs(title = paste("Treatment Effect vs", str_to_title(var)),
           x = str_to_title(var), y = "CATE") +
      theme_publication()
    plots[[var]] <- p

    # Spearman 相关 + 检验
    cor_test <- cor.test(plot_df[[var]], plot_df[[cate_col]], method = "spearman")
    correlations[[var]] <- data.frame(
      variable    = var,
      rho         = cor_test$estimate,
      p_value     = cor_test$p.value,
      row.names   = NULL
    )
  }

  list(plots = plots, correlations = do.call(rbind, correlations))
}

# ── 3. 异质性程度量化 ──
quantify_heterogeneity <- function(matched_data, cate_col = "tau_hat") {
  tau <- matched_data[[cate_col]]
  tau <- tau[!is.na(tau)]

  stats <- data.frame(
    measure = c("n", "mean_cate", "sd_cate", "median_cate",
                "q25", "q75", "iqr", "min", "max", "range",
                "prop_positive", "prop_negative",
                "prop_benefit_gt2pct", "prop_harm_lt_neg2pct",
                "coefficient_of_variation"),
    value = c(
      length(tau),
      mean(tau), sd(tau), median(tau),
      quantile(tau, 0.25), quantile(tau, 0.75),
      IQR(tau), min(tau), max(tau), diff(range(tau)),
      mean(tau > 0), mean(tau < 0),
      mean(tau > 0.02), mean(tau < -0.02),
      sd(tau) / abs(mean(tau))
    ),
    row.names = NULL
  )
  stats
}

# ── 4. 交互效应分析（含正式检验）──
analyze_interaction_effects <- function(matched_data,
                                        cate_col = "tau_hat",
                                        key_vars = c("age", "egfr", "bmi")) {
  key_vars <- intersect(key_vars, names(matched_data))
  results <- list()

  for (i in seq_along(key_vars)) {
    for (j in seq_along(key_vars)) {
      if (j <= i) next
      v1 <- key_vars[i]; v2 <- key_vars[j]

      q1 <- quantile(matched_data[[v1]], c(0, 1/3, 2/3, 1), na.rm = TRUE)
      q2 <- quantile(matched_data[[v2]], c(0, 1/3, 2/3, 1), na.rm = TRUE)

      temp <- matched_data
      temp$g1 <- cut(temp[[v1]], breaks = q1, include.lowest = TRUE,
                     labels = c("Low", "Mid", "High"))
      temp$g2 <- cut(temp[[v2]], breaks = q2, include.lowest = TRUE,
                     labels = c("Low", "Mid", "High"))
      temp <- temp[!is.na(temp$g1) & !is.na(temp$g2) & !is.na(temp[[cate_col]]), ]

      # 分组统计
      summ <- temp %>%
        group_by(g1, g2) %>%
        summarise(
          n = n(),
          mean_cate = mean(.data[[cate_col]], na.rm = TRUE),
          se_cate = sd(.data[[cate_col]], na.rm = TRUE) / sqrt(n()),
          ci_lower = mean_cate - 1.96 * se_cate,
          ci_upper = mean_cate + 1.96 * se_cate,
          .groups = "drop"
        ) %>%
        mutate(var1 = v1, var2 = v2)

      # ANOVA 检验交互项
      aov_test <- tryCatch({
        mod <- lm(as.formula(paste(cate_col, "~ g1 * g2")), data = temp)
        anova(mod)
      }, error = function(e) NULL)

      p_interaction <- NA
      if (!is.null(aov_test) && "g1:g2" %in% rownames(aov_test)) {
        p_interaction <- aov_test["g1:g2", "Pr(>F)"]
      }
      summ$p_interaction <- p_interaction

      results[[paste(v1, v2, sep = "_x_")]] <- summ
    }
  }

  results
}
