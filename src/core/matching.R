# matching.R - 倾向性评分匹配
# 项目：围术期 ACEI/ARB 停药策略因果森林分析
#
# 改动说明：
#   1. 使用 here::here() 路径
#   2. 使用 BMI（已在 data_prep.R 中计算）
#   3. 输出完整 SMD 表

cat("\n=== Phase 2: Propensity Score Matching ===\n")

# ── 1. 定义 PSM 公式 ──
match_formula <- as.formula(
  "preop_discontinuation ~ age + center + egfr + sex + surgery_type +
   hypertension + diabetes + heart_failure + ckd +
   preop_cr + bmi + preop_hb +
   ccb + statin + antiplatelet + beta_blocker +
   cpb + los + preop_bun + preop_plt + preop_alb"
)

# ── 2. 执行匹配 ──
cat("Running 1:1 nearest-neighbor matching (caliper = 0.1)...\n")
matching <- matchit(
  match_formula,
  data    = imputed_data,
  method  = "nearest",
  ratio   = 1,
  caliper = 0.1,
  replace = FALSE
)

matchdat <- match.data(matching)
cat("Matched sample: n =", nrow(matchdat), "\n")
cat("  Discontinued (W=1):", sum(as.numeric(as.character(matchdat$preop_discontinuation)) == 1), "\n")
cat("  Continued    (W=0):", sum(as.numeric(as.character(matchdat$preop_discontinuation)) == 0), "\n")

# ── 3. 创建分组变量 ──
matchdat$age_group <- cut(matchdat$age, breaks = c(0, 65, Inf),
                          labels = c("<65", ">=65"), right = FALSE)
matchdat$egfr_group <- cut(matchdat$egfr, breaks = c(0, 60, 90, Inf),
                           labels = c("<60", "60-89", ">=90"), right = FALSE)

# ── 4. 匹配平衡诊断 ──
cat("Generating balance diagnostics...\n")

# Love plot (SMD)
love_plot <- love.plot(
  matching,
  stats     = "mean.diffs",
  stars     = "std",
  threshold = c(m = 0.1),
  var.order = "unadjusted",
  abs       = TRUE,
  colors    = c("#00468B", "#ED0000")
) +
  labs(title = "Covariate Balance Before and After Matching",
       x = "Absolute Standardized Mean Difference") +
  theme_publication()

ggsave(here("results", "main_analysis", "matching", "balance_plot.pdf"),
       love_plot, width = 10, height = 8, dpi = 300)

# PS 分布图
ps_df <- data.frame(
  ps    = matching$distance,
  group = factor(imputed_data$preop_discontinuation,
                 levels = c("0", "1"),
                 labels = c("Continuation", "Discontinuation"))
)
ps_plot <- ggplot(ps_df, aes(x = ps, fill = group)) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("#00468B", "#ED0000"), name = "Group") +
  labs(x = "Propensity Score", y = "Density",
       title = "Propensity Score Distribution") +
  theme_publication()

ggsave(here("results", "main_analysis", "matching", "ps_distribution.pdf"),
       ps_plot, width = 8, height = 5, dpi = 300)

# ── 5. 基线特征表 ──
cat("Creating baseline characteristic tables...\n")

table_vars <- c(
  # 连续变量
  "age", "bmi", "egfr", "preop_cr", "preop_hb", "preop_bun",
  "preop_plt", "preop_alb", "preop_k", "los",
  # 分类变量
  "sex", "surgery_type", "heart_failure", "hypertension", "diabetes",
  "ckd", "respiratory_disease", "cpb",
  "ccb", "statin", "antiplatelet", "beta_blocker",
  # 结局
  "adverse_event", "aki_all", "aki1", "aki2", "aki3",
  "death", "ecmo", "iabp", "dialysis"
)

# 匹配前
tab_before <- CreateTableOne(
  vars       = table_vars,
  strata     = "preop_discontinuation",
  data       = imputed_data,
  factorVars = intersect(table_vars, categorical_vars),
  smd        = TRUE
)

# 匹配后
tab_after <- CreateTableOne(
  vars       = table_vars,
  strata     = "preop_discontinuation",
  data       = matchdat,
  factorVars = intersect(table_vars, categorical_vars),
  smd        = TRUE
)

write.csv(print(tab_before, smd = TRUE),
          here("results", "main_analysis", "matching", "baseline_before_matching.csv"))
write.csv(print(tab_after, smd = TRUE),
          here("results", "main_analysis", "matching", "baseline_after_matching.csv"))

# ── 6. 保存匹配数据 ──
write.csv(matchdat,
          here("results", "main_analysis", "matching", "matched_data.csv"),
          row.names = FALSE)

cat("=== Matching complete ===\n")
