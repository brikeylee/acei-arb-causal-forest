# main.R - 主程序入口：因果森林分析完整流程
# 项目：围术期 ACEI/ARB 停药策略对心脏手术患者预后的影响
# 目标期刊：Annals of Intensive Care
#
# 使用方法：
#   setwd("/path/to/CF")   # 设置工作目录为项目根目录
#   source("src/main.R")   # 运行完整流程
#
# Treatment coding: W=1 = 术前停用 (discontinuation), W=0 = 继续 (continuation)
# ATE > 0 → 停用增加风险 → 继续用药有保护作用

cat("============================================================\n")
cat("  Causal Forest Analysis: Perioperative ACEI/ARB Management\n")
cat("============================================================\n\n")

# ── 加载核心模块（函数定义 + 环境设置）──
source("src/core/setup.R")
source("src/core/functions.R")

# ── 加载分析模块（函数定义）──
source("src/analysis/validation_analysis.R")
source("src/analysis/sensitivity_analysis.R")
source("src/analysis/treatment_heterogeneity.R")
source("src/analysis/individualized_treatment.R")
source("src/analysis/clinical_decision_support.R")
source("src/analysis/supplementary_analysis.R")

# ══════════════════════════════════════════════════════════════
# 主分析流程函数
# ══════════════════════════════════════════════════════════════
run_full_analysis <- function() {

  # ── Phase 1: 数据预处理 ──
  source("src/core/data_prep.R")

  # ── Phase 2: 倾向性评分匹配 ──
  source("src/core/matching.R")

  # ── Phase 3: 因果森林分析 ──
  source("src/core/causal_forest.R")

  # ── Phase 4: 验证分析 ──
  cat("\n=== Phase 4: Validation & Robustness ===\n")

  # 4a. Bootstrap 验证
  boot_result <- bootstrap_validation(matchdat, n_bootstrap = 500)
  write.csv(boot_result$summary,
            here("results", "enhanced_analysis", "validation_analysis", "bootstrap_summary.csv"),
            row.names = FALSE)

  # 4b. RATE 异质性验证
  rate_result <- rate_heterogeneity_test(primary_result$model, matchdat, "tau_hat")

  # 4c. 超参数敏感性
  hp_sens <- hyperparameter_sensitivity(matchdat)
  write.csv(hp_sens,
            here("results", "enhanced_analysis", "validation_analysis", "hyperparameter_sensitivity.csv"),
            row.names = FALSE)

  # 4d. 变量重要性稳定性
  vimp_stability <- variable_importance_stability(matchdat, n_iterations = 100)
  write.csv(vimp_stability,
            here("results", "enhanced_analysis", "validation_analysis", "variable_importance_stability.csv"),
            row.names = FALSE)

  # ── Phase 5: 敏感性分析 ──
  cat("\n=== Phase 5: Sensitivity Analyses ===\n")

  sens_dir <- here("results", "sensitivity_analysis")

  # 5a. IPTW
  iptw_result <- run_iptw_analysis(imputed_data)
  write.csv(iptw_result, file.path(sens_dir, "iptw_result.csv"), row.names = FALSE)

  # 5b. E-value
  # 使用主分析的 OR (从匹配数据计算)
  W_matched <- as.numeric(as.character(matchdat$preop_discontinuation))
  Y_matched <- as.numeric(as.character(matchdat$adverse_event))
  or_model <- glm(Y_matched ~ W_matched, family = binomial())
  or_est <- exp(coef(or_model)["W_matched"])
  or_se  <- summary(or_model)$coefficients["W_matched", "Std. Error"]
  or_lo  <- exp(coef(or_model)["W_matched"] - 1.96 * or_se)
  or_hi  <- exp(coef(or_model)["W_matched"] + 1.96 * or_se)

  evalue_result <- calculate_evalue(or_est, or_lo, or_hi, rare_outcome = FALSE)
  write.csv(evalue_result, file.path(sens_dir, "evalue_result.csv"), row.names = FALSE)

  # 5c. Caliper 敏感性
  caliper_result <- caliper_sensitivity(imputed_data)
  write.csv(caliper_result, file.path(sens_dir, "caliper_sensitivity.csv"), row.names = FALSE)

  # 5d. 亚组交互检验
  interaction_result <- subgroup_interaction_tests(matchdat)
  if (!is.null(interaction_result)) {
    write.csv(interaction_result, file.path(sens_dir, "subgroup_interactions.csv"),
              row.names = FALSE)
  }

  # ── Phase 6: 治疗异质性深度分析 ──
  cat("\n=== Phase 6: Treatment Heterogeneity ===\n")

  het_dir <- here("results", "enhanced_analysis", "heterogeneity_analysis")

  # 6a. 异质性统计检验
  het_test <- test_treatment_heterogeneity(primary_result$model, matchdat)
  write.csv(het_test$variable_importance, file.path(het_dir, "variable_importance_detailed.csv"),
            row.names = FALSE)

  # 6b. 连续变量关系
  cont_rel <- analyze_continuous_relationships(matchdat, "tau_hat", c("age", "egfr", "bmi"))
  write.csv(cont_rel$correlations, file.path(het_dir, "continuous_correlations.csv"),
            row.names = FALSE)
  for (var_name in names(cont_rel$plots)) {
    ggsave(here("results", "enhanced_analysis", "visualizations",
                paste0(var_name, "_cate_relationship.pdf")),
           cont_rel$plots[[var_name]], width = 7, height = 5, dpi = 300)
  }

  # 6c. 异质性量化
  het_quant <- quantify_heterogeneity(matchdat, "tau_hat")
  write.csv(het_quant, file.path(het_dir, "heterogeneity_quantification.csv"),
            row.names = FALSE)

  # 6d. 交互效应
  int_effects <- analyze_interaction_effects(matchdat, "tau_hat", c("age", "egfr", "bmi"))
  int_combined <- do.call(rbind, int_effects)
  write.csv(int_combined, file.path(het_dir, "interaction_effects.csv"),
            row.names = FALSE)

  # ── Phase 7: 个体化治疗 ──
  cat("\n=== Phase 7: Individualized Treatment ===\n")

  ind_dir <- here("results", "enhanced_analysis", "individualized_treatment")

  # 7a. 患者聚类
  cluster_result <- perform_patient_clustering(matchdat, "tau_hat")
  write.csv(cluster_result$cluster_profiles, file.path(ind_dir, "cluster_profiles.csv"),
            row.names = FALSE)
  write.csv(cluster_result$silhouette_scores, file.path(ind_dir, "silhouette_scores.csv"),
            row.names = FALSE)

  # 7b. 治疗建议
  rec_result <- generate_treatment_recommendations(matchdat, "tau_hat")
  write.csv(rec_result$recommendation_summary,
            file.path(ind_dir, "recommendation_summary.csv"), row.names = FALSE)
  write.csv(rec_result$individual_recommendations,
            file.path(ind_dir, "individual_recommendations.csv"), row.names = FALSE)

  # 7c. CATE 分布可视化
  cate_plots <- plot_cate_distribution(matchdat, "tau_hat", by_group = "cluster")
  if (is.list(cate_plots) && "overall" %in% names(cate_plots)) {
    ggsave(here("results", "enhanced_analysis", "visualizations", "cate_distribution.pdf"),
           cate_plots$overall, width = 8, height = 5, dpi = 300)
    ggsave(here("results", "enhanced_analysis", "visualizations", "cate_by_cluster.pdf"),
           cate_plots$by_group, width = 8, height = 5, dpi = 300)
  }

  # ── Phase 8: 临床决策支持 ──
  cat("\n=== Phase 8: Clinical Decision Support ===\n")

  cds_dir <- here("results", "enhanced_analysis", "clinical_decision_support")

  # 8a. 临床评分
  score_result <- create_clinical_score(matchdat, "tau_hat")
  write.csv(score_result$score_performance, file.path(cds_dir, "clinical_score_performance.csv"),
            row.names = FALSE)
  write.csv(score_result$cv_validation, file.path(cds_dir, "score_cv_validation.csv"),
            row.names = FALSE)

  # 8b. 决策树
  tree_result <- build_decision_tree(matchdat, "tau_hat")
  saveRDS(tree_result$tree_model, file.path(cds_dir, "decision_tree_model.rds"))
  pdf(here("results", "enhanced_analysis", "visualizations", "decision_tree.pdf"),
      width = 12, height = 8)
  rpart.plot(tree_result$tree_model,
             main = "Clinical Decision Tree for ACEI/ARB Management",
             extra = 104, under = TRUE, faclen = 3)
  dev.off()

  # 8c. 校准分析
  cal_result <- calibration_analysis(primary_result$model, matchdat, "tau_hat")
  write.csv(cal_result$calibration_data, file.path(cds_dir, "calibration_analysis.csv"),
            row.names = FALSE)
  ggsave(here("results", "enhanced_analysis", "visualizations", "calibration_plot.pdf"),
         cal_result$calibration_plot, width = 7, height = 6, dpi = 300)

  # 8d. NNT/NNH
  nnt_result <- calculate_nnt_nnh(matchdat, subgroup_var = "age_group")
  nnt_egfr   <- calculate_nnt_nnh(matchdat, subgroup_var = "egfr_group")
  nnt_all    <- rbind(nnt_result, nnt_egfr[-1, ])  # 避免重复 Overall
  write.csv(nnt_all,
            here("results", "main_analysis", "benefit_risk", "nnt_nnh_analysis.csv"),
            row.names = FALSE)

  # ── Phase 9: 补充分析（AKI 悖论 + 竞争风险 + 传统对比 + 综合图）──
  cat("\n=== Phase 9: Supplementary Analyses ===\n")

  supp_dir <- here("results", "sensitivity_analysis")

  # 9a. 竞争风险分析
  cr_result <- run_competing_risk_analysis(matchdat)
  write.csv(cr_result, file.path(supp_dir, "competing_risk_analysis.csv"),
            row.names = FALSE)

  # 9b. 传统 vs CF 亚组对比
  trad_vs_cf <- run_traditional_subgroup_comparison(matchdat, cate_col = "tau_hat")
  write.csv(trad_vs_cf, file.path(supp_dir, "traditional_vs_cf_comparison.csv"),
            row.names = FALSE)

  # 9c. 综合 Forest Plot
  forest_result <- create_comprehensive_forest_plot(outcome_summary, matchdat)
  write.csv(forest_result$data, file.path(supp_dir, "comprehensive_forest_data.csv"),
            row.names = FALSE)
  ggsave(here("results", "enhanced_analysis", "visualizations", "comprehensive_forest_plot.pdf"),
         forest_result$plot, width = 11, height = 10, dpi = 300)

  # 9d. AKI 悖论深入分析
  aki_paradox <- analyze_aki_paradox(matchdat)
  write.csv(aki_paradox$paradox_table,
            file.path(supp_dir, "aki_paradox_detail.csv"), row.names = FALSE)
  write.csv(aki_paradox$severity_spectrum,
            file.path(supp_dir, "aki_severity_spectrum.csv"), row.names = FALSE)
  ggsave(here("results", "enhanced_analysis", "visualizations", "aki_paradox_plot.pdf"),
         aki_paradox$paradox_plot, width = 9, height = 6, dpi = 300)

  # ── 完成 ──
  cat("\n============================================================\n")
  cat("  Analysis complete! All results saved to results/\n")
  cat("============================================================\n")

  invisible(list(
    primary_result   = primary_result,
    boot_result      = boot_result,
    iptw_result      = iptw_result,
    evalue_result    = evalue_result,
    cluster_result   = cluster_result,
    score_result     = score_result,
    cal_result       = cal_result,
    cr_result        = cr_result,
    aki_paradox      = aki_paradox,
    forest_result    = forest_result
  ))
}

# ── 执行 ──
all_results <- run_full_analysis()
