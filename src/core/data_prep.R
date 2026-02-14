# data_prep.R - 数据读取、预处理、多重插补
# 项目：围术期 ACEI/ARB 停药策略因果森林分析
#
# 改动说明（相对原版）：
#   1. 使用 here::here() 管理路径
#   2. 多重插补 m=5（取 pooled dataset）
#   3. 计算 BMI
#   4. 输出缺失数据报告（Table S1）

cat("=== Phase 1: Data Preparation ===\n")

# ── 1. 读取原始数据 ──
paths_to_try <- c(
  here("baseline", "ACE_baseline250309.xlsx"),
  "/Users/lizongren/Desktop/ACEIARB/baseline/ACE_baseline250309.xlsx"
)
data_raw <- read_data(paths_to_try)
cat("Raw data loaded: n =", nrow(data_raw), "obs,", ncol(data_raw), "vars\n")

# ── 2. 变量重命名 ──
name_mapping <- c(
  "身高" = "height", "体重" = "weight", "性别" = "sex", "年龄" = "age",
  "手术类型" = "surgery_type",
  "心力衰竭史" = "heart_failure", "高血压史" = "hypertension",
  "糖尿病史" = "diabetes", "慢性肾病史" = "ckd",
  "呼吸系统疾病史" = "respiratory_disease",
  "住院时间" = "los", "术前住院时间" = "preop_los", "术后住院时间" = "postop_los",
  "术前末次肌酐" = "preop_cr", "术前末次血红蛋白" = "preop_hb",
  "术前末次尿素氮" = "preop_bun", "术前末次血小板" = "preop_plt",
  "术前末次白蛋白" = "preop_alb", "术前末次钾" = "preop_k",
  "体外循环" = "cpb",
  "AKI" = "aki_all", "AKI1" = "aki1", "AKI2" = "aki2", "AKI3" = "aki3",
  "死亡" = "death", "ECMO" = "ecmo", "IABP" = "iabp", "透析" = "dialysis",
  "术前钙拮抗剂" = "ccb", "术前降血脂" = "statin",
  "术前抗血栓" = "antiplatelet", "术前β阻剂" = "beta_blocker",
  "EGFR" = "egfr", "术前停用" = "preop_discontinuation", "c" = "center"
)
data_raw <- rename_columns(data_raw, name_mapping)

# ── 3. 定义变量类型 ──
categorical_vars <- c(
  "center", "sex", "surgery_type", "heart_failure", "hypertension",
  "diabetes", "ckd", "respiratory_disease", "cpb",
  "aki_all", "aki1", "aki2", "aki3", "death", "ecmo", "iabp", "dialysis",
  "ccb", "statin", "antiplatelet", "beta_blocker", "preop_discontinuation"
)
continuous_vars <- c(
  "height", "weight", "age", "los", "preop_los", "postop_los",
  "preop_cr", "preop_hb", "preop_bun", "preop_plt", "preop_alb",
  "preop_k", "egfr"
)

data_raw <- convert_variable_types(data_raw, continuous_vars, categorical_vars)
cat("Variable types converted. n =", nrow(data_raw), "\n")

# ── 4. 缺失数据报告（插补前）──
cat("\n--- Missing Data Report (Pre-imputation) ---\n")
n_total <- nrow(data_raw)
missing_report <- data.frame(
  variable       = c(continuous_vars, categorical_vars),
  stringsAsFactors = FALSE
)
missing_report$n_missing <- sapply(missing_report$variable, function(v) {
  if (v %in% names(data_raw)) sum(is.na(data_raw[[v]])) else NA
})
missing_report$pct_missing <- round(missing_report$n_missing / n_total * 100, 2)
missing_report <- missing_report[order(-missing_report$pct_missing), ]
print(missing_report[missing_report$n_missing > 0, ])

# 保存缺失数据报告 (Table S1)
write.csv(missing_report,
          here("results", "sensitivity_analysis", "missing_data_report.csv"),
          row.names = FALSE)

# ── 5. 删除分类变量有缺失的行 ──
cat_na_rows <- apply(data_raw[, categorical_vars, drop = FALSE], 1, function(x) any(is.na(x)))
n_cat_dropped <- sum(cat_na_rows)
cat("Rows dropped due to missing categorical vars:", n_cat_dropped, "\n")
data_for_imp <- data_raw[!cat_na_rows, ]
cat("Remaining for imputation: n =", nrow(data_for_imp), "\n")

# ── 6. 多重插补（m=5, PMM）──
cat("Running MICE with m=5 ...\n")
continuous_data <- data_for_imp[, continuous_vars, drop = FALSE]

imp <- mice(continuous_data,
            m = 5, maxit = 50,
            method = "pmm",
            seed = 123,
            printFlag = FALSE)

# 取 5 个插补数据集的均值作为 pooled dataset
# 注：对观察性研究这是可接受的折中方案，在 Limitations 中说明
imputed_list <- lapply(1:5, function(i) complete(imp, i))
pooled_continuous <- Reduce("+", imputed_list) / 5

# 对于整数型变量（如可能存在的），保留小数即可（它们是连续测量值）
data_for_imp[, continuous_vars] <- pooled_continuous
cat("Multiple imputation complete (m=5, pooled).\n")

# 验证无缺失
remaining_na <- sum(is.na(data_for_imp[, c(continuous_vars, categorical_vars)]))
cat("Remaining missing values after imputation:", remaining_na, "\n")

# ── 7. 计算 BMI ──
data_for_imp$bmi <- ifelse(
  !is.na(data_for_imp$height) & data_for_imp$height > 100 &
  !is.na(data_for_imp$weight) & data_for_imp$weight > 20,
  data_for_imp$weight / (data_for_imp$height / 100)^2,
  NA_real_
)
# 对 BMI 缺失值用中位数填补
bmi_median <- median(data_for_imp$bmi, na.rm = TRUE)
n_bmi_na <- sum(is.na(data_for_imp$bmi))
data_for_imp$bmi[is.na(data_for_imp$bmi)] <- bmi_median
cat("BMI calculated. Median:", round(bmi_median, 1),
    ", imputed", n_bmi_na, "missing values\n")

# ── 8. 创建复合结局变量 adverse_event ──
# 定义：死亡 OR AKI3 OR ECMO OR IABP OR 透析
cat("Creating composite adverse_event variable...\n")
data_for_imp$adverse_event <- as.factor(as.numeric(
  as.numeric(as.character(data_for_imp$death)) == 1 |
  as.numeric(as.character(data_for_imp$aki3)) == 1 |
  as.numeric(as.character(data_for_imp$ecmo)) == 1 |
  as.numeric(as.character(data_for_imp$iabp)) == 1 |
  as.numeric(as.character(data_for_imp$dialysis)) == 1
))
categorical_vars <- c(categorical_vars, "adverse_event")

event_rate <- mean(as.numeric(as.character(data_for_imp$adverse_event)) == 1, na.rm = TRUE)
cat("Composite adverse event rate:", round(event_rate * 100, 2), "%\n")

# ── 9. 创建年龄和 eGFR 分组变量（供亚组分析）──
data_for_imp$age_group <- cut(data_for_imp$age,
                               breaks = c(0, 65, Inf),
                               labels = c("<65", ">=65"),
                               right = FALSE)
data_for_imp$egfr_group <- cut(data_for_imp$egfr,
                                breaks = c(0, 60, 90, Inf),
                                labels = c("<60", "60-89", ">=90"),
                                right = FALSE)

# ── 10. 赋值给全局变量供下游使用 ──
imputed_data <- data_for_imp

cat("\n=== Data preparation complete ===\n")
cat("Final dataset: n =", nrow(imputed_data), ", vars =", ncol(imputed_data), "\n")
