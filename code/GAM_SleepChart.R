# ==============================================================================
# analysis:    Sleep duration vs. Biological Age Clocks (GAM)
# paper:       Sleep chart of biological aging clocks in middle and late life (The MULTI Consortium et al., 2026)
# author:      Junhao (Hao) Wen
# maintainer:  Junhao (Hao) Wen
# email:       junhao.wen89@gmail.com
# version:     0.1.0
# status:      Release
# license:     MIT (see LICENSE)
# copyright:   (c) 2025–2026 Junhao (Hao) Wen
# ==============================================================================

# ---- Packages ----
library(mgcv)
library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyverse)

# ---- I/O: set project-relative paths ----
PROJECT_ROOT <- "path/to/project"
IN_DIR       <- file.path(PROJECT_ROOT, "data")
OUT_DIR      <- file.path(PROJECT_ROOT, "outputs", "sleep_chart")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

sleep_file <- file.path(IN_DIR, "sleep.csv")
bag_file   <- file.path(IN_DIR, "BAG.tsv")
cov_file   <- file.path(IN_DIR, "covariate.csv")

# ---- Load data ----
sleep <- read.csv(sleep_file, stringsAsFactors = FALSE)
bag <- read.delim(bag_file, header = TRUE,
                  na.strings = c("NA", "", ".", "-9999"),
                  stringsAsFactors = FALSE)
covs <- read.csv(cov_file, stringsAsFactors = FALSE)

# ---- Merge ----
df <- sleep %>%
  select(participant_id, sleep_duration = sleep_duration_f1160) %>%
  filter(!sleep_duration %in% c(-1, -3)) %>%   # remove DK/Prefer-not-to-answer, etc.
  full_join(bag,  by = "participant_id") %>%
  full_join(covs, by = "participant_id")

# ---- 23 BAGs & covariates ----
BAG_list <- c(
  "Reproductive_female_ProtBAG", "Pulmonary_ProtBAG", "Heart_ProtBAG",
  "Brain_ProtBAG", "Eye_ProtBAG", "Hepatic_ProtBAG", "Renal_ProtBAG",
  "Reproductive_male_ProtBAG", "Endocrine_ProtBAG", "Immune_ProtBAG", "Skin_ProtBAG",
  "Endocrine_MetBAG", "Digestive_MetBAG", "Hepatic_MetBAG", "Immune_MetBAG",
  "Metabolic_MetBAG", "Brain_MRIBAG", "Adipose_MRIBAG", "Kidney_MRIBAG", "Heart_MRIBAG",
  "Liver_MRIBAG", "Pancreas_MRIBAG", "Spleen_MRIBAG"
)

covariates <- c(
  "Sex", "Age",
  "Weight", "Standing_height",
  "Waist_circumference", "BMI",
  "Diastolic_blood_pressure", "Systolic_blood_pressure",
  "Disease_status", "Assessment_center",
  "Time_diff")

# ---- Helper: robust extractor for summary tables ----
safe_extract <- function(obj, keys) {
  if (is.null(obj)) return(NA_real_)
  for (k in keys) if (k %in% names(obj)) return(suppressWarnings(as.numeric(obj[[k]])))
  return(NA_real_)
}

# ---- Core routine: model selection + inference + plotting ----
fit_and_test_effects <- function(outcome) {
  vars_to_keep <- c("sleep_duration", outcome, covariates)
  df_model <- df %>%
    select(all_of(vars_to_keep)) %>%
    tidyr::drop_na() %>%
    filter(dplyr::between(sleep_duration, 4, 10)) %>% # GAM modeling is prone to tail data sparcity, so we restricted the analyses for sleep hours from 4-10 hours per day.
    mutate(
      sex = factor(Sex, levels = c(0, 1), labels = c("female", "male")),
      # Normalize outcome to ~[1e-5, 1 + 1e-5] to stabilize Gamma/log fits if needed
      !!outcome := {
        x <- get(outcome)
        rng <- range(x, na.rm = TRUE)
        (x - rng[1]) / (rng[2] - rng[1]) + 1e-5
      }
    )

  # Candidate families & k
  families <- list(
    gaussian = gaussian(),
    tdist    = scat(),
    gamma    = Gamma(link = "log")
  )
  k_values <- c(3, 5, 10, 15, 20) # Data-driven selection for the optimal k, which allows both linear and non-linear relationships between sleep duration and BAG

  best_model  <- NULL
  best_aic    <- Inf
  best_family <- NA_character_
  best_k      <- NA_integer_

  # Grid search over (family, k)
  for (k in k_values) {
    model_formula <- as.formula(paste0(
      outcome, " ~ s(sleep_duration, k=", k, ", bs='cr') + ",
      "sex + s(sleep_duration, by=sex, k=", k, ", bs='cr') + ",
      paste(covariates, collapse = " + ")
    ))
    for (fam_name in names(families)) {
      fam <- families[[fam_name]]
      fit <- tryCatch(
        gam(model_formula, data = df_model, method = "REML", family = fam),
        error = function(e) NULL
      )
      if (!is.null(fit)) {
        aic <- AIC(fit)
        if (is.finite(aic) && aic < best_aic) {
          best_model  <- fit
          best_aic    <- aic
          best_family <- fam_name
          best_k      <- k
        }
      }
    }
  }

  if (is.null(best_model)) {
    warning(sprintf("No converged model for %s; returning NA.", outcome))
    return(list(plot = ggplot() + theme_void(), stats = data.frame(Outcome = outcome)))
  }

  # ---- Summaries ----
  sum_gam <- summary(best_model)

  # s() main sleep smooth (by default keyed as "s(sleep_duration)")
  main_sleep <- tryCatch(sum_gam$s.table["s(sleep_duration)", ], error = function(e) NULL)
  # sex main effect (label differs by link/family; match on name)
  ptab <- tryCatch(sum_gam$p.table, error = function(e) NULL)
  sex_row <- if (!is.null(ptab)) {
    rn <- rownames(ptab)
    hit <- which(rn %in% c("sexmale", "sexmale (Intercept)"))
    if (length(hit) == 0) which(grepl("^sexmale$", rn)) else hit
  } else integer(0)
  sex_diff <- if (length(sex_row) == 1) ptab[sex_row, ] else NULL

  # interaction smooth (by=sex)
  inter_row <- tryCatch(rownames(sum_gam$s.table), error = function(e) character(0))
  inter_hit <- which(inter_row == "s(sleep_duration):sexmale")
  sex_interaction_term <- if (length(inter_hit) == 1) sum_gam$s.table[inter_hit, ] else NULL

  # ---- Predictions ----
  pred_data <- expand.grid(
    sleep_duration = seq(4, 10, length.out = 100),
    sex = factor(c("female", "male"), levels = c("female", "male"))
  )
  # Hold covariates at mean/mode
  for (cov in covariates) {
    if (is.numeric(df_model[[cov]])) {
      pred_data[[cov]] <- mean(df_model[[cov]], na.rm = TRUE)
    } else {
      pred_data[[cov]] <- names(sort(table(df_model[[cov]]), decreasing = TRUE))[1]
    }
  }

  pred <- predict(best_model, newdata = pred_data, se.fit = TRUE)
  pred_data$fit <- as.numeric(pred$fit)

  # Empirical minima of the estimated BAG for sleep (min predicted BAG) per sex
  optimals <- pred_data %>%
    group_by(sex) %>%
    slice(which.min(fit)) %>%
    ungroup()

  # 95% CI on link scale (approx)
  ci_data <- pred_data %>%
    mutate(
      lower = fit - 1.96 * pred$se.fit,
      upper = fit + 1.96 * pred$se.fit
    )

  # Bonferroni gate for showing “optimal” BAG across 23 outcomes
  add_optimal_lines <- (!is.null(main_sleep)) &&
    (safe_extract(main_sleep, c("p-value", "pvalue", "p.val")) < 0.05 / length(BAG_list))
  optimals_filtered <- if (isTRUE(add_optimal_lines)) optimals else optimals[0, ]

  # ---- Plot: points + fit ----
  p_points <- ggplot(df_model, aes(x = sleep_duration, y = !!sym(outcome), color = sex)) +
    geom_point(alpha = 0.3, size = 1) +
    scale_color_manual(values = c(female = "#EB5F2C", male = "#0072B5")) +
    labs(title = paste(outcome), x = NULL, y = "BAG") +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 10),
      legend.position = "none",
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text   = element_text(color = "black")
    )

  p_fit <- ggplot(ci_data, aes(x = sleep_duration, y = fit, color = sex)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = sex), alpha = 0.2, color = NA) +
    geom_line(linewidth = 0.8) +
    geom_vline(
      data = optimals_filtered,
      aes(xintercept = sleep_duration, color = sex),
      linetype = "dashed", linewidth = 0.6
    ) +
    scale_color_manual(values = c(female = "#EB5F2C", male = "#0072B5")) +
    scale_fill_manual(values = c(female = "#EB5F2C", male = "#0072B5")) +
    labs(x = "Sleep Duration (hours)", y = "BAG") +
    theme_minimal(base_size = 10) +
    theme(legend.position = "none", axis.text = element_text(color = "black"))

  p <- (p_points / p_fit) +
    plot_layout(heights = c(1, 1.2)) +
    plot_annotation(subtitle = paste("Best-fitting family:", best_family))

  # ---- Stats table ----
  stats <- data.frame(
    Outcome                 = outcome,
    Family                  = best_family,
    Optimal_k               = best_k,
    Sleep_edf               = safe_extract(main_sleep, c("edf", "Ref.df")),
    Sleep_pvalue            = safe_extract(main_sleep, c("p-value", "pvalue", "p.val")),
    Sex_coef                = safe_extract(sex_diff,  c("Estimate")),
    Sex_pvalue              = safe_extract(sex_diff,  c("Pr(>|t|)", "Pr(>|z|)", "p-value")),
    Sex_interaction_pvalue  = safe_extract(sex_interaction_term, c("p-value", "pvalue", "p.val")),
    Female_optimal          = if ("female" %in% optimals$sex) optimals$sleep_duration[optimals$sex == "female"] else NA_real_,
    Male_optimal            = if ("male"   %in% optimals$sex)   optimals$sleep_duration[optimals$sex   == "male"]   else NA_real_
  )

  # ---- Export per-outcome data (optional) ----
  ci_export <- ci_data %>%
    select(sleep_duration, sex, fit, lower, upper) %>%
    rename(
      BAG_predict       = fit,
      BAG_predict_lower = lower,
      BAG_predict_upper = upper
    )

  readr::write_tsv(ci_export,
                   file.path(OUT_DIR, paste0("Sleepchart_figure_data_", outcome, ".tsv"))) # this is used to generate the sleep-BAG visualization at our SleepChart portal: e.g., https://labs-laboratory.com/sleepchart/bag/
  readr::write_tsv(stats,
                   file.path(OUT_DIR, paste0("Sleepchart_figure_stats_", outcome, ".tsv")))

  list(plot = p, stats = stats)
}

# ---- Run across all outcomes ----
all_results <- lapply(BAG_list, fit_and_test_effects)
names(all_results) <- BAG_list

# ---- Combine plots ----
plot_list    <- lapply(all_results, function(x) x$plot)
combined_plot <- wrap_plots(plot_list, ncol = 6)

# ---- Export pooled stats ----
stats_table <- bind_rows(lapply(all_results, function(x) x$stats))
readr::write_tsv(stats_table, file.path(OUT_DIR, "BAG_sleep_stats_GAM_CI.tsv"))
