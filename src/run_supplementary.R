# run_supplementary.R - 仅运行补充分析（Phase 9）
# 前提：主分析已完成，匹配数据和 CATE 已存在

cat("Loading modules...\n")
source("src/core/setup.R")
source("src/core/functions.R")
source("src/analysis/supplementary_analysis.R")

cat("Loading matched data with CATE...\n")
matchdat <- read.csv(here("results", "main_analysis", "causal_forest", "matched_data_with_cate.csv"))
outcome_summary <- read.csv(here("results", "outcome_summary.csv"))

cat("Dataset loaded: n =", nrow(matchdat), "\n")

supp_dir <- here("results", "sensitivity_analysis")

# 1. 竞争风险分析
cr_result <- run_competing_risk_analysis(matchdat)
write.csv(cr_result, file.path(supp_dir, "competing_risk_analysis.csv"), row.names = FALSE)

# 2. 传统 vs CF 亚组对比
trad_vs_cf <- run_traditional_subgroup_comparison(matchdat, cate_col = "tau_hat")
write.csv(trad_vs_cf, file.path(supp_dir, "traditional_vs_cf_comparison.csv"), row.names = FALSE)

# 3. 综合 Forest Plot
forest_result <- create_comprehensive_forest_plot(outcome_summary, matchdat)
write.csv(forest_result$data, file.path(supp_dir, "comprehensive_forest_data.csv"), row.names = FALSE)
ggsave(here("results", "enhanced_analysis", "visualizations", "comprehensive_forest_plot.pdf"),
       forest_result$plot, width = 11, height = 10, dpi = 300)

# 4. AKI 悖论深入分析
aki_paradox <- analyze_aki_paradox(matchdat)
write.csv(aki_paradox$paradox_table, file.path(supp_dir, "aki_paradox_detail.csv"), row.names = FALSE)
write.csv(aki_paradox$severity_spectrum, file.path(supp_dir, "aki_severity_spectrum.csv"), row.names = FALSE)
ggsave(here("results", "enhanced_analysis", "visualizations", "aki_paradox_plot.pdf"),
       aki_paradox$paradox_plot, width = 9, height = 6, dpi = 300)

cat("\nSupplementary analyses complete!\n")
