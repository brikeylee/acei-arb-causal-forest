# individualized_treatment.R - 个体化治疗策略分析模块
# 项目：围术期 ACEI/ARB 停药策略因果森林分析
#
# 改动说明：
#   1. 修正变量名：使用 sex/tau_hat，去掉 smoking/previous_mi/gender
#   2. 聚类数 k 添加 silhouette 验证
#   3. CATE 分类基于分位数而非任意阈值

cat("Individualized treatment functions loaded.\n")

# ── 1. 患者聚类分析（含 silhouette 验证）──
perform_patient_clustering <- function(matched_data,
                                       cate_col = "tau_hat",
                                       cluster_features = c("age", "egfr", "bmi",
                                                            "hypertension", "diabetes"),
                                       k_range = 2:6) {
  cat("  Running patient clustering...\n")

  cluster_features <- intersect(cluster_features, names(matched_data))
  features_and_cate <- c(cate_col, cluster_features)

  work <- matched_data[complete.cases(matched_data[, features_and_cate]), features_and_cate]

  # 标准化
  work_scaled <- as.data.frame(scale(sapply(work, function(x) as.numeric(as.character(x)))))

  # Silhouette 最优 k
  cat("  Evaluating optimal k...\n")
  sil_scores <- numeric(length(k_range))
  for (i in seq_along(k_range)) {
    k <- k_range[i]
    set.seed(123)
    km <- kmeans(work_scaled, centers = k, nstart = 25)
    ss <- cluster::silhouette(km$cluster, dist(work_scaled))
    sil_scores[i] <- mean(ss[, 3])
    cat("    k =", k, " silhouette =", round(sil_scores[i], 3), "\n")
  }
  optimal_k <- k_range[which.max(sil_scores)]
  cat("  Optimal k =", optimal_k, "(silhouette =", round(max(sil_scores), 3), ")\n")

  # 最终聚类
  set.seed(123)
  km_final <- kmeans(work_scaled, centers = optimal_k, nstart = 25)
  work$cluster <- as.factor(km_final$cluster)

  # 聚类画像
  cluster_profiles <- work %>%
    group_by(cluster) %>%
    summarise(
      n_patients = n(),
      mean_cate  = mean(.data[[cate_col]], na.rm = TRUE),
      sd_cate    = sd(.data[[cate_col]], na.rm = TRUE),
      mean_age   = mean(age, na.rm = TRUE),
      mean_egfr  = mean(egfr, na.rm = TRUE),
      mean_bmi   = mean(bmi, na.rm = TRUE),
      prop_htn   = mean(as.numeric(as.character(hypertension)) == 1, na.rm = TRUE),
      prop_dm    = mean(as.numeric(as.character(diabetes)) == 1, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      pct = round(n_patients / sum(n_patients) * 100, 1),
      recommendation = case_when(
        mean_cate > quantile(work[[cate_col]], 0.75, na.rm = TRUE) ~
          "Strong recommendation for continuation",
        mean_cate > median(work[[cate_col]], na.rm = TRUE) ~
          "Moderate recommendation for continuation",
        TRUE ~ "Individualized decision"
      )
    )

  # 将聚类标签加回原数据
  matched_data$cluster <- NA_character_
  complete_idx <- complete.cases(matched_data[, features_and_cate])
  matched_data$cluster[complete_idx] <- as.character(km_final$cluster)

  list(
    clustered_data    = matched_data,
    cluster_profiles  = cluster_profiles,
    optimal_k         = optimal_k,
    silhouette_scores = data.frame(k = k_range, silhouette = sil_scores)
  )
}

# ── 2. 个体化治疗建议（基于分位数阈值）──
generate_treatment_recommendations <- function(matched_data,
                                                cate_col = "tau_hat") {
  cat("  Generating treatment recommendations...\n")

  tau <- matched_data[[cate_col]]
  tau <- tau[!is.na(tau)]

  # 基于分位数定义阈值（数据驱动）
  q_lo  <- quantile(tau, 0.25)
  q_hi  <- quantile(tau, 0.75)

  recommendations <- matched_data %>%
    filter(!is.na(.data[[cate_col]])) %>%
    mutate(
      cate_category = case_when(
        .data[[cate_col]] >= q_hi ~ "High",
        .data[[cate_col]] >= q_lo ~ "Moderate",
        TRUE ~ "Low"
      ),
      final_recommendation = case_when(
        cate_category == "High" ~ "Continue ACE/ARB (strong evidence of benefit)",
        cate_category == "Moderate" ~ "Continue ACE/ARB (moderate evidence)",
        cate_category == "Low" ~ "Individualized decision"
      ),
      confidence_level = cate_category
    )

  # 汇总
  rec_summary <- recommendations %>%
    group_by(final_recommendation, confidence_level) %>%
    summarise(
      n_patients = n(),
      percentage = round(n() / nrow(recommendations) * 100, 1),
      mean_cate  = mean(.data[[cate_col]], na.rm = TRUE),
      sd_cate    = sd(.data[[cate_col]], na.rm = TRUE),
      mean_age   = mean(age, na.rm = TRUE),
      mean_egfr  = mean(egfr, na.rm = TRUE),
      .groups = "drop"
    )

  list(
    individual_recommendations = recommendations,
    recommendation_summary     = rec_summary,
    thresholds                 = c(q_low = q_lo, q_high = q_hi)
  )
}

# ── 3. CATE 分布可视化 ──
plot_cate_distribution <- function(matched_data, cate_col = "tau_hat",
                                    by_group = NULL) {
  plot_df <- matched_data[!is.na(matched_data[[cate_col]]), ]

  p1 <- ggplot(plot_df, aes(x = .data[[cate_col]])) +
    geom_histogram(bins = 50, fill = "#00468B", alpha = 0.7, color = "white") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "#ED0000", linewidth = 0.7) +
    labs(title = "Distribution of Conditional Average Treatment Effects",
         x = "CATE (Risk Difference)", y = "Number of Patients") +
    theme_publication()

  if (!is.null(by_group) && by_group %in% names(plot_df)) {
    p2 <- ggplot(plot_df, aes(x = .data[[by_group]], y = .data[[cate_col]],
                               fill = .data[[by_group]])) +
      geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "#ED0000") +
      scale_fill_manual(values = lancet_colors) +
      labs(title = paste("CATE Distribution by", str_to_title(by_group)),
           x = str_to_title(by_group), y = "CATE") +
      theme_publication() + theme(legend.position = "none")
    return(list(overall = p1, by_group = p2))
  }

  p1
}
