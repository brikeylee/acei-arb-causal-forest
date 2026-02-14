# causal_forest.R - 因果森林模型分析（主结局 + 各独立结局）
# 项目：围术期 ACEI/ARB 停药策略因果森林分析
#
# Treatment coding: W=1 = 术前停用 ACEI/ARB (discontinuation)
#                   W=0 = 术前继续 (continuation)
# ATE > 0 → 停用增加不良事件风险 → 继续用药有保护作用
# ATE < 0 → 停用减少不良事件风险

cat("\n=== Phase 3: Causal Forest Analysis ===\n")

# ── 定义协变量集 ──
cf_covariates <- c("age", "sex", "center", "heart_failure", "hypertension",
                   "diabetes", "ckd", "respiratory_disease", "cpb", "egfr",
                   "bmi", "ccb", "statin", "antiplatelet", "beta_blocker")

treatment_col <- "preop_discontinuation"

# ── 定义所有需要分析的结局 ──
outcomes_to_analyze <- c(
  "adverse_event",  # 主要复合结局
  "death",          # 次要结局
  "aki1", "aki2", "aki3", "aki_all"
)

# ── 函数：对单个结局运行因果森林 ──
run_single_cf <- function(data, outcome_col, treatment_col, covariates, num_trees = 2000) {
  cat("  Running CF for outcome:", outcome_col, "...\n")

  mats <- prepare_cf_matrices(data, outcome_col, treatment_col, covariates)

  cf <- causal_forest(
    mats$X, mats$Y, mats$W,
    num.trees        = num_trees,
    honesty          = TRUE,
    tune.parameters  = "all",
    seed             = 123
  )

  # 提取结果
  ate_result <- extract_ate_results(cf, outcome_col)
  tau_hat    <- predict(cf)$predictions
  var_imp    <- variable_importance(cf)

  # test_calibration
  cal_test <- tryCatch(test_calibration(cf), error = function(e) NULL)

  list(
    model       = cf,
    ate         = ate_result,
    tau_hat     = tau_hat,
    var_imp     = data.frame(variable = mats$covariate_names,
                             importance = var_imp,
                             rank = rank(-var_imp)),
    calibration = cal_test,
    matrices    = mats
  )
}

# ── 1. 主要分析：复合结局 adverse_event ──
cat("Running primary analysis (adverse_event)...\n")
primary_result <- run_single_cf(matchdat, "adverse_event", treatment_col, cf_covariates)

# 保存主模型
saveRDS(primary_result$model,
        here("results", "main_analysis", "models", "causal_forest_model.rds"))

# 添加 CATE 到匹配数据
matchdat$tau_hat <- NA_real_
matchdat$tau_hat[primary_result$matrices$complete_idx] <- primary_result$tau_hat

# ATE 结果
cat("\n--- Primary Outcome: Composite Adverse Event ---\n")
print(primary_result$ate)

# ── 2. 次要分析：各独立结局 ──
cat("\nRunning secondary analyses...\n")
all_ate_results <- list(primary_result$ate)
secondary_models <- list()

for (outcome in setdiff(outcomes_to_analyze, "adverse_event")) {
  res <- run_single_cf(matchdat, outcome, treatment_col, cf_covariates)
  all_ate_results[[length(all_ate_results) + 1]] <- res$ate
  secondary_models[[outcome]] <- res
}

# 汇总所有结局的 ATE
outcome_summary <- do.call(rbind, all_ate_results)
cat("\n--- All Outcome ATEs ---\n")
print(outcome_summary)

write.csv(outcome_summary,
          here("results", "outcome_summary.csv"),
          row.names = FALSE)

# ── 3. 保存主结局详细结果 ──
# ATE
write.csv(primary_result$ate,
          here("results", "main_analysis", "causal_forest", "average_treatment_effect.csv"),
          row.names = FALSE)

# 带 CATE 的匹配数据
write.csv(matchdat,
          here("results", "main_analysis", "causal_forest", "matched_data_with_cate.csv"),
          row.names = FALSE)

# 变量重要性
write.csv(primary_result$var_imp[order(-primary_result$var_imp$importance), ],
          here("results", "main_analysis", "causal_forest", "variable_importance.csv"),
          row.names = FALSE)

# Calibration test
if (!is.null(primary_result$calibration)) {
  cal_df <- data.frame(
    test = c("mean.forest.prediction", "differential.forest.prediction"),
    estimate = primary_result$calibration[, 1],
    std.err  = primary_result$calibration[, 2],
    t.value  = primary_result$calibration[, 3],
    p.value  = primary_result$calibration[, 4]
  )
  write.csv(cal_df,
            here("results", "main_analysis", "causal_forest", "calibration_test.csv"),
            row.names = FALSE)
}

# ── 4. 亚组特征比较（高获益 vs 高风险）──
cat("Subgroup characterization...\n")
q25 <- quantile(matchdat$tau_hat, 0.25, na.rm = TRUE)
q75 <- quantile(matchdat$tau_hat, 0.75, na.rm = TRUE)

high_benefit <- matchdat[!is.na(matchdat$tau_hat) & matchdat$tau_hat >= q75, ]
high_risk    <- matchdat[!is.na(matchdat$tau_hat) & matchdat$tau_hat <= q25, ]

# 注意：tau_hat > 0 意味着停药增加风险 → 继续用药获益更大
# 所以 high_benefit = tau_hat 最高的 25%（继续用药获益最大）
write.csv(high_benefit,
          here("results", "main_analysis", "causal_forest", "high_benefit_subgroup.csv"),
          row.names = FALSE)
write.csv(high_risk,
          here("results", "main_analysis", "causal_forest", "high_risk_subgroup.csv"),
          row.names = FALSE)

# 亚组对比表
compare_vars <- cf_covariates
subgroup_comparison <- data.frame(
  variable = compare_vars,
  high_benefit_mean = sapply(compare_vars, function(v) {
    if (is.factor(high_benefit[[v]])) {
      mean(as.numeric(as.character(high_benefit[[v]])), na.rm = TRUE)
    } else {
      mean(high_benefit[[v]], na.rm = TRUE)
    }
  }),
  high_risk_mean = sapply(compare_vars, function(v) {
    if (is.factor(high_risk[[v]])) {
      mean(as.numeric(as.character(high_risk[[v]])), na.rm = TRUE)
    } else {
      mean(high_risk[[v]], na.rm = TRUE)
    }
  }),
  row.names = NULL
)
write.csv(subgroup_comparison,
          here("results", "main_analysis", "causal_forest", "subgroup_comparison.csv"),
          row.names = FALSE)

# ── 5. 可视化 ──
cat("Generating causal forest visualizations...\n")

cf_vis_dir <- here("results", "main_analysis", "causal_forest")

# CATE 分布图
p_ite <- ggplot(matchdat[!is.na(matchdat$tau_hat), ], aes(x = tau_hat)) +
  geom_histogram(bins = 50, fill = "#00468B", alpha = 0.7, color = "white") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#ED0000", linewidth = 0.8) +
  labs(title = "Distribution of Individual Treatment Effects",
       subtitle = "Positive values indicate higher risk with discontinuation",
       x = "Conditional Average Treatment Effect (CATE)",
       y = "Number of Patients") +
  theme_publication()
ggsave(file.path(cf_vis_dir, "individual_treatment_effect_distribution.pdf"),
       p_ite, width = 8, height = 5, dpi = 300)

# CATE vs Age
p_age <- ggplot(matchdat[!is.na(matchdat$tau_hat), ],
                aes(x = age, y = tau_hat)) +
  geom_point(alpha = 0.2, size = 0.8, color = "#00468B") +
  geom_smooth(method = "loess", span = 0.3, color = "#ED0000", fill = "#ED0000", alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  labs(title = "Treatment Effect by Age",
       x = "Age (years)", y = "CATE") +
  theme_publication()
ggsave(file.path(cf_vis_dir, "age_treatment_effect.pdf"),
       p_age, width = 7, height = 5, dpi = 300)

# CATE vs eGFR
p_egfr <- ggplot(matchdat[!is.na(matchdat$tau_hat), ],
                 aes(x = egfr, y = tau_hat)) +
  geom_point(alpha = 0.2, size = 0.8, color = "#00468B") +
  geom_smooth(method = "loess", span = 0.3, color = "#42B540", fill = "#42B540", alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  labs(title = "Treatment Effect by eGFR",
       x = "eGFR (mL/min/1.73m²)", y = "CATE") +
  theme_publication()
ggsave(file.path(cf_vis_dir, "egfr_treatment_effect.pdf"),
       p_egfr, width = 7, height = 5, dpi = 300)

# 变量重要性
vimp <- primary_result$var_imp[order(-primary_result$var_imp$importance), ]
vimp_top <- head(vimp, 15)
vimp_top$variable <- factor(vimp_top$variable, levels = rev(vimp_top$variable))

p_vimp <- ggplot(vimp_top, aes(x = variable, y = importance)) +
  geom_col(fill = "#00468B", alpha = 0.85) +
  coord_flip() +
  labs(title = "Variable Importance in Causal Forest",
       x = "", y = "Importance") +
  theme_publication()
ggsave(file.path(cf_vis_dir, "variable_importance.pdf"),
       p_vimp, width = 8, height = 6, dpi = 300)

cat("=== Causal forest analysis complete ===\n")
