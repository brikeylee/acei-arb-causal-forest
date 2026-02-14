# setup.R - 环境设置：安装包、加载库、创建目录、定义主题
# 项目：围术期 ACEI/ARB 停药策略因果森林分析

# ── Locale（确保中文列名匹配）──
invisible(Sys.setlocale("LC_ALL", "en_US.UTF-8"))

# ── CRAN 镜像 & 用户库 ──
options(repos = c(CRAN = "https://cloud.r-project.org"))
user_lib <- Sys.getenv("R_LIBS_USER")
if (nchar(user_lib) > 0) {
  dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(user_lib, .libPaths()))
}

# ── 必要的 R 包 ──
required_packages <- c(
  # 基础数据处理
  "dplyr", "tidyr", "readxl", "here",
  # 可视化

  "ggplot2", "ggpubr", "gridExtra", "grid", "scales", "viridis",
  # 统计分析
  "tableone", "mice", "survival", "sandwich", "lmtest",
  # 匹配与因果推断
  "MatchIt", "cobalt", "WeightIt", "survey",
  # 因果森林
  "grf",
  # 敏感性分析
  "EValue",
  # 决策树与机器学习
  "rpart", "rpart.plot", "pROC", "boot",
  # 聚类

  "cluster", "factoextra"
)

# 检查并安装缺失的包
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages)) {
  cat("Installing missing packages:", paste(new_packages, collapse = ", "), "\n")
  install.packages(new_packages, lib = .libPaths()[1])
}

# 加载所有包
suppressMessages({
  for (pkg in required_packages) {
    library(pkg, character.only = TRUE)
  }
})

# ── 使用 here::here() 创建目录结构 ──
dirs_to_create <- c(

  here("results"),
  here("results", "main_analysis", "matching"),
  here("results", "main_analysis", "causal_forest"),
  here("results", "main_analysis", "models"),
  here("results", "main_analysis", "benefit_risk"),
  here("results", "enhanced_analysis", "individualized_treatment"),
  here("results", "enhanced_analysis", "heterogeneity_analysis"),
  here("results", "enhanced_analysis", "clinical_decision_support"),
  here("results", "enhanced_analysis", "validation_analysis"),
  here("results", "enhanced_analysis", "visualizations"),
  here("results", "sensitivity_analysis")
)

for (d in dirs_to_create) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ── 随机数种子 ──
set.seed(123)

# ── 出版级别主题 (Lancet / Nature style) ──
theme_publication <- function(base_size = 11) {
  theme_bw(base_size = base_size) %+replace%
    theme(
      panel.background  = element_rect(fill = "white"),
      panel.grid.major  = element_line(color = "grey93", linewidth = 0.3),
      panel.grid.minor  = element_blank(),
      panel.border      = element_rect(color = "grey70", fill = NA, linewidth = 0.5),
      axis.line         = element_blank(),
      axis.ticks        = element_line(color = "grey70", linewidth = 0.3),
      axis.text         = element_text(size = rel(0.9), color = "black"),
      axis.title        = element_text(size = rel(1), face = "bold"),
      legend.background = element_rect(fill = "white", color = NA),
      legend.key        = element_rect(fill = "white", color = NA),
      legend.title      = element_text(face = "bold", size = rel(0.9)),
      legend.text       = element_text(size = rel(0.85)),
      plot.title        = element_text(face = "bold", size = rel(1.15), hjust = 0),
      plot.subtitle     = element_text(size = rel(0.95), hjust = 0),
      strip.background  = element_rect(fill = "grey95", color = "grey70"),
      strip.text        = element_text(face = "bold", size = rel(0.9)),
      plot.margin       = margin(10, 10, 10, 10)
    )
}

theme_set(theme_publication())

# Lancet-style 配色
lancet_colors <- c("#00468B", "#ED0000", "#42B540", "#0099B4", "#925E9F",
                   "#FDAF91", "#AD002A", "#ADB6B6")

# ── 美化映射表 ──
outcome_labels <- c(
  "adverse_event" = "Major Adverse Events",
  "aki_all"       = "AKI (All Stages)",
  "aki1"          = "AKI Stage 1",
  "aki2"          = "AKI Stage 2",
  "aki3"          = "AKI Stage 3",
  "death"         = "30-day Mortality",
  "dialysis"      = "Dialysis"
)

cat("Setup complete: all packages loaded, directories created.\n")
