# export_model.R - 拟合 CATE 近似模型并导出系数供前端使用
library(here)
library(jsonlite)

cat("Loading matched data with CATE...\n")
d <- read.csv(here("results", "main_analysis", "causal_forest", "matched_data_with_cate.csv"))
cat("n =", nrow(d), "\n")

# 确保数值型
d$sex_num      <- ifelse(d$sex == "M", 1, 0)  # 1=Male, 0=Female
d$diabetes_num <- as.integer(d$diabetes)
d$htn_num      <- as.integer(d$hypertension)
d$hf_num       <- as.integer(d$heart_failure)
d$ckd_num      <- as.integer(d$ckd)

# 拟合近似模型（含交互项）
model <- lm(tau_hat ~ age + egfr + bmi + sex_num +
              diabetes_num + htn_num + hf_num + ckd_num +
              age:egfr + age:diabetes_num + age:hf_num,
            data = d)

cat("\nModel summary:\n")
print(summary(model))
cat("\nR-squared:", summary(model)$r.squared, "\n")

# 提取系数
coefs <- as.list(coef(model))
names(coefs) <- gsub(":", "_x_", names(coefs))

# CATE 分布统计（供前端绘制直方图）
tau <- d$tau_hat[!is.na(d$tau_hat)]
hist_breaks <- seq(min(tau) - 0.001, max(tau) + 0.001, length.out = 51)
hist_obj <- hist(tau, breaks = hist_breaks, plot = FALSE)

cate_distribution <- list(
  counts = hist_obj$counts,
  breaks = hist_obj$breaks,
  mids   = hist_obj$mids
)

# 分布统计
cate_stats <- list(
  n        = length(tau),
  mean     = mean(tau),
  sd       = sd(tau),
  median   = median(tau),
  q25      = quantile(tau, 0.25, names = FALSE),
  q75      = quantile(tau, 0.75, names = FALSE),
  min      = min(tau),
  max      = max(tau),
  pct_positive = mean(tau > 0)
)

# 组装 JSON
export_data <- list(
  coefficients     = coefs,
  r_squared        = summary(model)$r.squared,
  cate_distribution = cate_distribution,
  cate_stats       = cate_stats,
  variable_means   = list(
    age     = mean(d$age, na.rm = TRUE),
    egfr    = mean(d$egfr, na.rm = TRUE),
    bmi     = mean(d$bmi, na.rm = TRUE)
  )
)

out_path <- here("tools", "clinical_calculator", "model_data.json")
write_json(export_data, out_path, auto_unbox = TRUE, digits = 8, pretty = TRUE)
cat("\nExported to:", out_path, "\n")
