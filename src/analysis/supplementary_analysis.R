# supplementary_analysis.R - 补充分析模块
# 内容：
#   1. 竞争风险分析（死亡作为 AKI 的竞争事件）
#   2. 传统 logistic 回归亚组分析（与 CF 对比）
#   3. 综合 Forest Plot（全部结局一图展示）
#   4. AKI 悖论深入分析

suppressMessages({
  library(dplyr)
  library(ggplot2)
  library(survival)
  library(here)
})

cat("Supplementary analysis functions loaded.\n")

# ══════════════════════════════════════════════════════════════
# 1. 竞争风险分析
# ══════════════════════════════════════════════════════════════
run_competing_risk_analysis <- function(matched_data,
                                        treatment_col = "preop_discontinuation") {
  cat("  Running competing risk analysis (death as competing event for AKI)...\n")

  W <- as.numeric(as.character(matched_data[[treatment_col]]))
  death   <- as.numeric(as.character(matched_data$death))
  aki_all <- as.numeric(as.character(matched_data$aki_all))
  aki1    <- as.numeric(as.character(matched_data$aki1))
  aki3    <- as.numeric(as.character(matched_data$aki3))

  # 竞争风险状态编码：0=无事件, 1=AKI, 2=死亡(竞争事件)
  # 对于 AKI 结局，死亡是竞争事件
  make_cr_status <- function(aki_vec, death_vec) {
    status <- rep(0L, length(aki_vec))
    status[aki_vec == 1] <- 1L      # AKI 发生
    status[death_vec == 1 & aki_vec == 0] <- 2L  # 死亡但无 AKI（竞争事件）
    status
  }

  # 使用 logistic 回归近似竞争风险分析
  # （因为这是短期二元结局而非时间-事件数据）
  results <- list()

  aki_outcomes <- list(
    "AKI_all"   = aki_all,
    "AKI_Stage1" = aki1,
    "AKI_Stage3" = aki3
  )

  for (outcome_name in names(aki_outcomes)) {
    aki_vec <- aki_outcomes[[outcome_name]]
    cr_status <- make_cr_status(aki_vec, death)

    # 标准分析（不考虑竞争风险）
    mod_standard <- glm(aki_vec ~ W, family = binomial())
    or_std <- exp(coef(mod_standard)["W"])
    se_std <- summary(mod_standard)$coefficients["W", "Std. Error"]
    p_std  <- summary(mod_standard)$coefficients["W", "Pr(>|z|)"]

    # 排除死亡后的分析（仅在存活者中）
    alive_idx <- death == 0
    if (sum(alive_idx) > 100) {
      mod_alive <- glm(aki_vec[alive_idx] ~ W[alive_idx], family = binomial())
      or_alive <- exp(coef(mod_alive)["W[alive_idx]"])
      se_alive <- summary(mod_alive)$coefficients["W[alive_idx]", "Std. Error"]
      p_alive  <- summary(mod_alive)$coefficients["W[alive_idx]", "Pr(>|z|)"]
    } else {
      or_alive <- NA; se_alive <- NA; p_alive <- NA
    }

    # 多项式 logistic 回归处理竞争事件
    # 简化为：分别计算 cause-specific hazard
    # Cause 1: AKI
    mod_cs1 <- glm((cr_status == 1) ~ W, family = binomial())
    or_cs1 <- exp(coef(mod_cs1)["W"])
    p_cs1  <- summary(mod_cs1)$coefficients["W", "Pr(>|z|)"]
    se_cs1 <- summary(mod_cs1)$coefficients["W", "Std. Error"]

    # Cause 2: Death (competing)
    mod_cs2 <- glm((cr_status == 2) ~ W, family = binomial())
    or_cs2 <- exp(coef(mod_cs2)["W"])
    p_cs2  <- summary(mod_cs2)$coefficients["W", "Pr(>|z|)"]
    se_cs2 <- summary(mod_cs2)$coefficients["W", "Std. Error"]

    results[[outcome_name]] <- data.frame(
      outcome = outcome_name,
      # 标准分析
      OR_standard = or_std,
      OR_std_CI_lo = exp(coef(mod_standard)["W"] - 1.96 * se_std),
      OR_std_CI_hi = exp(coef(mod_standard)["W"] + 1.96 * se_std),
      p_standard = p_std,
      # 排除死亡
      OR_survivors_only = or_alive,
      p_survivors = p_alive,
      # Cause-specific: AKI
      OR_cause_specific_AKI = or_cs1,
      OR_csAKI_CI_lo = exp(coef(mod_cs1)["W"] - 1.96 * se_cs1),
      OR_csAKI_CI_hi = exp(coef(mod_cs1)["W"] + 1.96 * se_cs1),
      p_cause_specific_AKI = p_cs1,
      # Cause-specific: Death
      OR_cause_specific_death = or_cs2,
      p_cause_specific_death = p_cs2,
      row.names = NULL
    )
  }

  combined <- do.call(rbind, results)
  cat("  Competing risk analysis complete.\n")
  combined
}

# ══════════════════════════════════════════════════════════════
# 2. 传统 logistic 回归亚组分析（对比 CF）
# ══════════════════════════════════════════════════════════════
run_traditional_subgroup_comparison <- function(matched_data,
                                                 treatment_col = "preop_discontinuation",
                                                 cate_col = "tau_hat") {
  cat("  Running traditional vs CF subgroup comparison...\n")

  W <- as.numeric(as.character(matched_data[[treatment_col]]))
  Y <- as.numeric(as.character(matched_data$adverse_event))

  subgroups <- list(
    age_lt65  = matched_data$age < 65,
    age_ge65  = matched_data$age >= 65,
    egfr_lt60 = matched_data$egfr < 60,
    egfr_60_89 = matched_data$egfr >= 60 & matched_data$egfr < 90,
    egfr_ge90 = matched_data$egfr >= 90,
    dm_yes    = as.numeric(as.character(matched_data$diabetes)) == 1,
    dm_no     = as.numeric(as.character(matched_data$diabetes)) == 0,
    hf_yes    = as.numeric(as.character(matched_data$heart_failure)) == 1,
    hf_no     = as.numeric(as.character(matched_data$heart_failure)) == 0,
    male      = as.numeric(as.character(matched_data$sex)) == 1,
    female    = as.numeric(as.character(matched_data$sex)) == 0
  )

  results <- list()
  for (sg_name in names(subgroups)) {
    idx <- subgroups[[sg_name]] & !is.na(subgroups[[sg_name]])
    if (sum(idx) < 50) next

    # 传统 logistic
    mod <- glm(Y[idx] ~ W[idx], family = binomial())
    coefs <- summary(mod)$coefficients
    if (!"W[idx]" %in% rownames(coefs)) next

    or_trad <- exp(coefs["W[idx]", "Estimate"])
    se_trad <- coefs["W[idx]", "Std. Error"]
    rd_trad <- mean(Y[idx & W == 1], na.rm = TRUE) - mean(Y[idx & W == 0], na.rm = TRUE)

    # CF CATE 均值
    tau_sub <- matched_data[[cate_col]][idx]
    mean_cate_cf <- mean(tau_sub, na.rm = TRUE)

    results[[sg_name]] <- data.frame(
      subgroup = sg_name,
      n = sum(idx),
      # 传统方法
      OR_traditional = or_trad,
      OR_trad_CI_lo = exp(coefs["W[idx]", "Estimate"] - 1.96 * se_trad),
      OR_trad_CI_hi = exp(coefs["W[idx]", "Estimate"] + 1.96 * se_trad),
      RD_traditional = rd_trad,
      p_traditional = coefs["W[idx]", "Pr(>|z|)"],
      # CF 方法
      mean_CATE_cf = mean_cate_cf,
      # 方向一致性
      direction_consistent = (rd_trad > 0) == (mean_cate_cf > 0),
      row.names = NULL
    )
  }

  combined <- do.call(rbind, results)
  cat("  Traditional vs CF comparison complete.\n")
  combined
}

# ══════════════════════════════════════════════════════════════
# 3. 综合 Forest Plot（全部结局 + 亚组一图展示）
# ══════════════════════════════════════════════════════════════
create_comprehensive_forest_plot <- function(outcome_summary, matched_data,
                                              treatment_col = "preop_discontinuation") {
  cat("  Creating comprehensive forest plot...\n")

  W <- as.numeric(as.character(matched_data[[treatment_col]]))

  # A. 各结局的 OR（从匹配数据直接计算）
  outcome_vars <- c("adverse_event", "death", "aki1", "aki2", "aki3", "aki_all",
                     "ecmo", "iabp", "dialysis")
  outcome_labels <- c("Composite Adverse Events", "30-day Mortality",
                       "AKI Stage 1", "AKI Stage 2", "AKI Stage 3", "AKI (All Stages)",
                       "ECMO", "IABP", "Dialysis")

  or_results <- list()
  for (i in seq_along(outcome_vars)) {
    ov <- outcome_vars[i]
    if (!ov %in% names(matched_data)) next
    Y_i <- as.numeric(as.character(matched_data[[ov]]))
    valid <- !is.na(Y_i) & !is.na(W)
    mod <- tryCatch(glm(Y_i[valid] ~ W[valid], family = binomial()), error = function(e) NULL)
    if (is.null(mod)) next
    coefs <- summary(mod)$coefficients
    if (!"W[valid]" %in% rownames(coefs)) next

    est <- coefs["W[valid]", "Estimate"]
    se  <- coefs["W[valid]", "Std. Error"]

    or_results[[length(or_results) + 1]] <- data.frame(
      label    = outcome_labels[i],
      OR       = exp(est),
      CI_lower = exp(est - 1.96 * se),
      CI_upper = exp(est + 1.96 * se),
      p_value  = coefs["W[valid]", "Pr(>|z|)"],
      category = "Outcomes",
      row.names = NULL
    )
  }

  # B. 主要亚组的复合结局 OR
  subgroup_defs <- list(
    "Age < 65"  = matched_data$age < 65,
    "Age >= 65" = matched_data$age >= 65,
    "eGFR < 60" = matched_data$egfr < 60,
    "eGFR 60-89" = matched_data$egfr >= 60 & matched_data$egfr < 90,
    "eGFR >= 90" = matched_data$egfr >= 90,
    "Diabetes"    = as.numeric(as.character(matched_data$diabetes)) == 1,
    "No Diabetes" = as.numeric(as.character(matched_data$diabetes)) == 0,
    "Heart Failure" = as.numeric(as.character(matched_data$heart_failure)) == 1,
    "No Heart Failure" = as.numeric(as.character(matched_data$heart_failure)) == 0
  )

  Y_ae <- as.numeric(as.character(matched_data$adverse_event))
  for (sg_name in names(subgroup_defs)) {
    idx <- subgroup_defs[[sg_name]] & !is.na(subgroup_defs[[sg_name]])
    if (sum(idx) < 50) next
    mod <- tryCatch(glm(Y_ae[idx] ~ W[idx], family = binomial()), error = function(e) NULL)
    if (is.null(mod)) next
    coefs <- summary(mod)$coefficients
    if (!"W[idx]" %in% rownames(coefs)) next
    est <- coefs["W[idx]", "Estimate"]
    se  <- coefs["W[idx]", "Std. Error"]

    or_results[[length(or_results) + 1]] <- data.frame(
      label    = sg_name,
      OR       = exp(est),
      CI_lower = exp(est - 1.96 * se),
      CI_upper = exp(est + 1.96 * se),
      p_value  = coefs["W[idx]", "Pr(>|z|)"],
      category = "Subgroups",
      row.names = NULL
    )
  }

  forest_data <- do.call(rbind, or_results)

  # 整理顺序
  forest_data$category <- factor(forest_data$category, levels = c("Outcomes", "Subgroups"))
  forest_data$label <- factor(forest_data$label, levels = rev(forest_data$label))

  # 构建 forest plot
  forest_data$sig <- ifelse(forest_data$p_value < 0.05, "Significant", "Not significant")
  forest_data$or_text <- sprintf("%.2f (%.2f-%.2f)", forest_data$OR,
                                  forest_data$CI_lower, forest_data$CI_upper)

  p <- ggplot(forest_data, aes(x = OR, y = label, color = sig)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
    geom_point(size = 3) +
    geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.2) +
    geom_text(aes(x = max(forest_data$CI_upper, na.rm = TRUE) * 1.1, label = or_text),
              hjust = 0, size = 3, color = "black") +
    scale_color_manual(values = c("Significant" = "#ED0000", "Not significant" = "#00468B"),
                       name = "") +
    scale_x_log10() +
    facet_grid(category ~ ., scales = "free_y", space = "free_y") +
    labs(title = "Odds Ratios for Discontinuation vs Continuation of ACEI/ARB",
         subtitle = "Values > 1 favor continuation (discontinuation increases risk)",
         x = "Odds Ratio (log scale)", y = "") +
    theme_publication() +
    theme(strip.text = element_text(face = "bold", size = 11),
          legend.position = "bottom")

  list(plot = p, data = forest_data)
}

# ══════════════════════════════════════════════════════════════
# 4. AKI 悖论深入分析
# ══════════════════════════════════════════════════════════════
analyze_aki_paradox <- function(matched_data,
                                 treatment_col = "preop_discontinuation") {
  cat("  Analyzing AKI paradox in detail...\n")

  W <- as.numeric(as.character(matched_data[[treatment_col]]))
  aki1 <- as.numeric(as.character(matched_data$aki1))
  aki3 <- as.numeric(as.character(matched_data$aki3))
  death <- as.numeric(as.character(matched_data$death))
  ae   <- as.numeric(as.character(matched_data$adverse_event))
  ecmo <- as.numeric(as.character(matched_data$ecmo))
  iabp <- as.numeric(as.character(matched_data$iabp))
  dialysis <- as.numeric(as.character(matched_data$dialysis))

  # AKI1 患者中后续严重事件的发生率
  aki1_patients <- aki1 == 1

  paradox_table <- data.frame(
    group = c("AKI1 with Continuation", "AKI1 with Discontinuation",
              "No AKI1 with Continuation", "No AKI1 with Discontinuation"),
    n = c(sum(aki1_patients & W == 0, na.rm = TRUE),
          sum(aki1_patients & W == 1, na.rm = TRUE),
          sum(!aki1_patients & W == 0, na.rm = TRUE),
          sum(!aki1_patients & W == 1, na.rm = TRUE)),
    death_rate = c(
      mean(death[aki1_patients & W == 0], na.rm = TRUE),
      mean(death[aki1_patients & W == 1], na.rm = TRUE),
      mean(death[!aki1_patients & W == 0], na.rm = TRUE),
      mean(death[!aki1_patients & W == 1], na.rm = TRUE)
    ),
    aki3_rate = c(
      mean(aki3[aki1_patients & W == 0], na.rm = TRUE),
      mean(aki3[aki1_patients & W == 1], na.rm = TRUE),
      mean(aki3[!aki1_patients & W == 0], na.rm = TRUE),
      mean(aki3[!aki1_patients & W == 1], na.rm = TRUE)
    ),
    composite_ae_rate = c(
      mean(ae[aki1_patients & W == 0], na.rm = TRUE),
      mean(ae[aki1_patients & W == 1], na.rm = TRUE),
      mean(ae[!aki1_patients & W == 0], na.rm = TRUE),
      mean(ae[!aki1_patients & W == 1], na.rm = TRUE)
    ),
    ecmo_iabp_rate = c(
      mean((ecmo | iabp)[aki1_patients & W == 0], na.rm = TRUE),
      mean((ecmo | iabp)[aki1_patients & W == 1], na.rm = TRUE),
      mean((ecmo | iabp)[!aki1_patients & W == 0], na.rm = TRUE),
      mean((ecmo | iabp)[!aki1_patients & W == 1], na.rm = TRUE)
    ),
    dialysis_rate = c(
      mean(dialysis[aki1_patients & W == 0], na.rm = TRUE),
      mean(dialysis[aki1_patients & W == 1], na.rm = TRUE),
      mean(dialysis[!aki1_patients & W == 0], na.rm = TRUE),
      mean(dialysis[!aki1_patients & W == 1], na.rm = TRUE)
    ),
    row.names = NULL
  )

  # AKI 严重度分层的 OR
  aki_severity <- data.frame(
    outcome = c("AKI Stage 1 (mild)", "AKI Stage 2 (moderate)",
                "AKI Stage 3 (severe)", "Death", "ECMO/IABP", "Dialysis",
                "Composite Adverse Events"),
    stringsAsFactors = FALSE
  )

  outcome_vecs <- list(aki1, as.numeric(as.character(matched_data$aki2)),
                       aki3, death, as.numeric(ecmo | iabp), dialysis, ae)

  for (i in seq_along(outcome_vecs)) {
    y <- outcome_vecs[[i]]
    valid <- !is.na(y) & !is.na(W)
    mod <- tryCatch(glm(y[valid] ~ W[valid], family = binomial()), error = function(e) NULL)
    if (!is.null(mod)) {
      coefs <- summary(mod)$coefficients
      key <- "W[valid]"
      if (key %in% rownames(coefs)) {
        aki_severity$OR[i] <- exp(coefs[key, "Estimate"])
        aki_severity$CI_lower[i] <- exp(coefs[key, "Estimate"] - 1.96 * coefs[key, "Std. Error"])
        aki_severity$CI_upper[i] <- exp(coefs[key, "Estimate"] + 1.96 * coefs[key, "Std. Error"])
        aki_severity$p_value[i] <- coefs[key, "Pr(>|z|)"]
      }
    }
  }

  # AKI 悖论图：severity spectrum
  aki_severity$outcome <- factor(aki_severity$outcome,
                                  levels = rev(aki_severity$outcome))
  aki_severity$direction <- ifelse(aki_severity$OR < 1, "Favors Discontinuation",
                                    "Favors Continuation")

  p_paradox <- ggplot(aki_severity, aes(x = OR, y = outcome, color = direction)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
    geom_point(size = 4) +
    geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.25, linewidth = 0.8) +
    scale_color_manual(values = c("Favors Continuation" = "#ED0000",
                                   "Favors Discontinuation" = "#00468B"),
                       name = "Direction of Effect") +
    scale_x_log10(breaks = c(0.5, 0.7, 1, 1.3, 1.5, 2)) +
    labs(title = "The AKI Paradox: Severity-Dependent Treatment Effects",
         subtitle = "OR > 1: discontinuation increases risk; OR < 1: discontinuation decreases risk",
         x = "Odds Ratio (log scale)", y = "") +
    theme_publication() +
    theme(legend.position = "bottom")

  list(
    paradox_table = paradox_table,
    severity_spectrum = aki_severity,
    paradox_plot = p_paradox
  )
}
