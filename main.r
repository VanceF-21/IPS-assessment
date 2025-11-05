### Introductory Probability and Statistics: Biomarkers and Pain Analysis


## 1. SETUP
library(readxl)
library(dplyr)
library(tidyr)
library(janitor)
library(ggplot2)
library(broom)
library(car)

set.seed(123)

cat("Setup complete.\n\n")


# 2. DATA LOADING AND PREPROCESSING
suppressMessages({
  biomarkers_raw <- read_excel("data/biomarkers.xlsx")
  covariates_raw <- read_excel("data/covariates.xlsx")
})

biomarkers_cleaned <- clean_names(biomarkers_raw)
covariates_cleaned <- clean_names(covariates_raw)

# Process biomarkers (extract patient ID & timepoint)
biomarkers_processed <- biomarkers_cleaned %>%
  separate(col = 1, into = c("patient_id", "timepoint"), sep = "-", remove = FALSE) %>%
  mutate(patient_id = as.integer(patient_id))

biomarkers_inclusion <- biomarkers_processed %>%
  filter(timepoint == "0weeks")

biomarker_names <- c(
  "il_8", "vegf_a", "opg", "tgf_beta_1", "il_6",
  "cxcl9", "cxcl1", "il_18", "csf_1"
)

# Process covariates
covariates_processed <- covariates_cleaned %>%
  mutate(sex = factor(sex, levels = c(1, 2), labels = c("Male", "Female")))

# Merge datasets
merged_data <- inner_join(biomarkers_inclusion, covariates_processed, by = "patient_id")

cat("Data loading complete.\n")
cat(paste("Observations after merging:", nrow(merged_data), "\n"))

# Handle missing 12-month VAS values
missing_count <- sum(is.na(merged_data$vas_12months))
cat(paste("Missing values in vas_12months:", missing_count, "\n"))

# Task-specific subsets
data_task1 <- merged_data                                  # 117 obs
data_task2 <- merged_data %>% filter(!is.na(vas_12months)) # 115 obs

cat(paste("Task 1 sample size:", nrow(data_task1), "\n"))
cat(paste("Task 2 sample size:", nrow(data_task2), "\n\n"))

cat(paste("Sex distribution (Task 1):",
          sum(data_task1$sex == "Male"), "males,",
          sum(data_task1$sex == "Female"), "females\n"))
cat(paste("Sex distribution (Task 2):",
          sum(data_task2$sex == "Male"), "males,",
          sum(data_task2$sex == "Female"), "females\n\n"))

# Descriptive statistics by sex (Task 1)
cat("Descriptive Statistics by Sex (Task 1 dataset):\n")
summary_by_sex <- data_task1 %>%
  group_by(sex) %>%
  summarise(
    n = n(),
    across(all_of(biomarker_names),
           list(mean = ~mean(., na.rm = TRUE),
                sd   = ~sd(.,   na.rm = TRUE)),
           .names = "{.col}_{.fn}")
  )
print(summary_by_sex, width = Inf)
cat("\n")


## 3. TASK 1: HYPOTHESIS TESTING
cat("TASK 1: HYPOTHESIS TESTING\n")
cat("Research question: Do biomarker levels at inclusion differ between males and females?\n\n")

# Assumption checks (example IL-8)
male_data   <- data_task1 %>% filter(sex == "Male")   %>% pull(il_8)
female_data <- data_task1 %>% filter(sex == "Female") %>% pull(il_8)

shapiro_male   <- shapiro.test(male_data)
shapiro_female <- shapiro.test(female_data)

cat("Shapiro–Wilk test for normality (IL-8):\n")
cat(paste("  Males:   W =", round(shapiro_male$statistic, 4),
          "p =", round(shapiro_male$p.value, 4), "\n"))
cat(paste("  Females: W =", round(shapiro_female$statistic, 4),
          "p =", round(shapiro_female$p.value, 4), "\n\n"))

levene_test <- leveneTest(il_8 ~ sex, data = data_task1)
cat("Levene’s test for homogeneity of variance (IL-8):\n")
print(levene_test)

# t-tests for all biomarkers
test_results <- data.frame(
  biomarker = biomarker_names,
  p_value = NA, t_statistic = NA, df = NA,
  mean_male = NA, mean_female = NA
)

for (i in seq_along(biomarker_names)) {
  f <- as.formula(paste(biomarker_names[i], "~ sex"))
  res <- t.test(f, data = data_task1)
  test_results[i, 2:6] <- c(res$p.value, res$statistic, res$parameter,
                            res$estimate[1], res$estimate[2])
}

cat("\nResults of t-tests (uncorrected):\n")
print(test_results, digits = 4)

# Bonferroni correction
test_results_corrected <- test_results %>%
  mutate(
    p_adj_bonferroni = p.adjust(p_value, method = "bonferroni"),
    significant_uncorrected = ifelse(p_value < 0.05, "Yes", "No"),
    significant_bonferroni  = ifelse(p_adj_bonferroni < 0.05, "Yes", "No")
  )

cat("\nResults with Bonferroni correction:\n")
print(test_results_corrected, digits = 4)

sig_biomarkers <- test_results_corrected %>%
  filter(significant_uncorrected == "Yes") %>%
  mutate(direction = ifelse(mean_male > mean_female, "Male > Female", "Female > Male")) %>%
  select(biomarker, mean_male, mean_female, direction, p_value,
         p_adj_bonferroni, significant_bonferroni)

if (nrow(sig_biomarkers) > 0) {
  cat("\nBiomarkers with significant differences (p < 0.05 before correction):\n")
  print(sig_biomarkers, digits = 4, row.names = FALSE)
} else {
  cat("\nNo significant differences found.\n")
}
cat("\nTask 1 complete.\n\n")


## 4. TASK 2: REGRESSION MODELLING
cat("TASK 2: REGRESSION MODELLING\n\n")

train_index <- sample(1:nrow(data_task2), size = 0.8 * nrow(data_task2))
train_data  <- data_task2[train_index, ]
test_data   <- data_task2[-train_index, ]

cat(paste("Training set:", nrow(train_data), "observations\n"))
cat(paste("Testing set:",  nrow(test_data),  "observations\n\n"))

# FULL MODEL
cat("--- FULL MODEL ---\n")
covariate_names <- c("age", "vas_at_inclusion")
all_predictors  <- c(biomarker_names, covariate_names)

model_formula_full <- as.formula(
  paste("vas_12months ~", paste(all_predictors, collapse = " + "))
)

model_full <- lm(model_formula_full, data = train_data)
model_coefficients_full <- tidy(model_full)

cat("Parameter estimates:\n")
print(model_coefficients_full, digits = 4)

sig_predictors_from_full <- model_coefficients_full %>%
  filter(p.value < 0.05, term != "(Intercept)") %>%
  arrange(p.value)

if (nrow(sig_predictors_from_full) > 0) {
  cat("\nSignificant predictors (p < 0.05):\n")
  print(sig_predictors_from_full %>% select(term, estimate, p.value),
        digits = 4, row.names = FALSE)
}

model_summary_full <- summary(model_full)
cat("\nModel fit (Training):\n")
cat(paste("R²:", round(model_summary_full$r.squared, 4), "\n"))
cat(paste("Adjusted R²:", round(model_summary_full$adj.r.squared, 4), "\n"))
cat(paste("Residual SE:", round(model_summary_full$sigma, 4), "\n"))

cat("\nVariance Inflation Factors:\n")
print(round(vif(model_full), 2))

cooks_d <- cooks.distance(model_full)
cat(paste("\nObservations with Cook’s D > 4/n:",
          sum(cooks_d > 4 / nrow(train_data)), "\n"))

pred_full <- predict(model_full, newdata = test_data)
rmse_full <- sqrt(mean((test_data$vas_12months - pred_full)^2))
mae_full  <- mean(abs(test_data$vas_12months - pred_full))
cor_full  <- cor(test_data$vas_12months, pred_full)

cat("\nOut-of-sample evaluation (Test):\n")
cat(paste("RMSE:", round(rmse_full, 4),
          "MAE:", round(mae_full, 4),
          "Correlation:", round(cor_full, 4), "\n\n"))

# Reduced model
cat("--- REDUCED MODEL - VARIABLE SELECTION ---\n\n")

# Rank ALL predictors by p-value
all_predictors_ranked <- model_coefficients_full %>%
  filter(term != "(Intercept)") %>%
  arrange(p.value)

if (nrow(all_predictors_ranked) >= 1) {
  max_predictors <- nrow(all_predictors_ranked) 
  
  # Store results for each k
  selection_results <- data.frame(
    k = integer(),
    adj_r2_train = numeric(),
    rmse_train = numeric(),
    adj_r2_test = numeric(),
    rmse_test = numeric(),
    mae_test = numeric(),
    cor_test = numeric(),
    aic = numeric(),
    bic = numeric()
  )
  
  # Try models with k = 1, 2, 3, ... predictors
  for (k in 1:max_predictors) {
    selected_predictors <- all_predictors_ranked$term[1:k]
    formula_k <- as.formula(
      paste("vas_12months ~", paste(selected_predictors, collapse = " + "))
    )
    model_k <- lm(formula_k, data = train_data)
    
    # Training performance
    summary_k <- summary(model_k)
    adj_r2_train <- summary_k$adj.r.squared
    rmse_train <- sqrt(mean(residuals(model_k)^2))
    
    # Test performance
    pred_k <- predict(model_k, newdata = test_data)
    rmse_test <- sqrt(mean((test_data$vas_12months - pred_k)^2))
    mae_test <- mean(abs(test_data$vas_12months - pred_k))
    cor_test <- cor(test_data$vas_12months, pred_k)
    
    ss_total_test <- sum((test_data$vas_12months - mean(test_data$vas_12months))^2)
    ss_residual_test <- sum((test_data$vas_12months - pred_k)^2)
    r2_test <- 1 - (ss_residual_test / ss_total_test)
    n_test <- nrow(test_data)
    adj_r2_test <- 1 - ((1 - r2_test) * (n_test - 1) / (n_test - k - 1))
    
    # Store results
    selection_results <- rbind(selection_results, data.frame(
      k = k,
      adj_r2_train = adj_r2_train,
      rmse_train = rmse_train,
      adj_r2_test = adj_r2_test,
      rmse_test = rmse_test,
      mae_test = mae_test,
      cor_test = cor_test,
      aic = AIC(model_k),
      bic = BIC(model_k)
    ))
  }
  
  cat("Model selection results (by number of predictors 'k'):\n")
  print(selection_results, digits = 4, row.names = FALSE)
  
  # Select best model
  best_k_rmse <- selection_results$k[which.min(selection_results$rmse_test)]
  best_k_bic <- selection_results$k[which.min(selection_results$bic)]
  
  cat("\nBest k based on test RMSE:", best_k_rmse, "\n")
  cat("Best k based on BIC:", best_k_bic, "\n")
  
  final_k <- best_k_rmse
  cat(paste("\nSelected k =", final_k, "for final reduced model\n\n"))
  
  # Build final reduced model
  selected_final_predictors <- all_predictors_ranked$term[1:final_k]
  
  cat("Final reduced model predictors:\n")
  for (i in 1:final_k) {
    p_val <- all_predictors_ranked$p.value[i]
    cat(paste("  ", i, ".", selected_final_predictors[i], 
              "(p =", sprintf("%.4f", p_val), ")\n"))
  }
  cat("\n")
  
  model_formula_reduced <- as.formula(
    paste("vas_12months ~", paste(selected_final_predictors, collapse = " + "))
  )
  model_reduced <- lm(model_formula_reduced, data = train_data)
  
  cat("Parameter estimates:\n")
  model_coefficients_reduced <- tidy(model_reduced)
  print(model_coefficients_reduced, digits = 4)
  
  # Model fit
  cat("\nModel fit (Training):\n")
  model_summary_reduced <- summary(model_reduced)
  adj_r_squared_reduced <- model_summary_reduced$adj.r.squared
  cat(paste("R-squared:", round(model_summary_reduced$r.squared, 4), "\n"))
  cat(paste("Adjusted R-squared:", round(adj_r_squared_reduced, 4), "\n"))
  cat(paste("Residual standard error:", round(model_summary_reduced$sigma, 4), "\n"))
  cat(paste("AIC:", round(AIC(model_reduced), 2), "\n"))
  cat(paste("BIC:", round(BIC(model_reduced), 2), "\n\n"))
  
  # Out-of-sample evaluation
  cat("Out-of-sample evaluation (Test):\n")
  predictions_reduced <- predict(model_reduced, newdata = test_data)
  rmse_reduced <- sqrt(mean((test_data$vas_12months - predictions_reduced)^2))
  mae_reduced <- mean(abs(test_data$vas_12months - predictions_reduced))
  cor_reduced <- cor(test_data$vas_12months, predictions_reduced)
  
  cat(paste("RMSE:", round(rmse_reduced, 4), "\n"))
  cat(paste("MAE:", round(mae_reduced, 4), "\n"))
  cat(paste("Correlation:", round(cor_reduced, 4), "\n\n"))
  
  # Final comparison
  cat("--- FINAL MODEL COMPARISON ---\n")
  comparison_df <- data.frame(
    Model = c("Full", "Reduced"),
    N_Predictors = c(length(all_predictors), final_k),
    Adj_R2_Train = c(round(adj_r_squared_full, 4), round(adj_r_squared_reduced, 4)),
    RMSE_Test = c(round(rmse_full, 4), round(rmse_reduced, 4)),
    MAE_Test = c(round(mae_full, 4), round(mae_reduced, 4)),
    Correlation = c(round(cor_full, 4), round(cor_reduced, 4)),
    AIC = c(round(AIC(model_full), 2), round(AIC(model_reduced), 2)),
    BIC = c(round(BIC(model_full), 2), round(BIC(model_reduced), 2))
  )
  print(comparison_df, row.names = FALSE)
  
  # Visual comparison of model selection
  # Main plot: Test RMSE
  selection_plot <- ggplot(selection_results, aes(x = k)) +
    geom_line(aes(y = rmse_test, color = "Test RMSE"), linewidth = 1) +
    geom_point(aes(y = rmse_test, color = "Test RMSE"), size = 3) +
    geom_vline(xintercept = final_k, linetype = "dashed", color = "red") +
    labs(
      x = "Number of Predictors (k)",
      y = "RMSE on Test Set",
      color = ""
    ) +
    theme_minimal(base_size = 14) +
    scale_x_continuous(breaks = seq(1, max_predictors, by = ifelse(max_predictors > 8, 2, 1))) +
    annotate("text", x = final_k, y = max(selection_results$rmse_test),
             label = paste("Selected k =", final_k),
             hjust = -0.1, color = "red", size = 4)
  
  ggsave("figs/model_selection.pdf", selection_plot, width = 8, height = 5)
  print(selection_plot)
  
  # Additional plot: Multiple metrics comparison
  selection_long <- selection_results %>%
    select(k, rmse_test, aic, bic, adj_r2_train) %>%
    pivot_longer(cols = c(rmse_test, aic, bic, adj_r2_train),
                 names_to = "metric", values_to = "value")
  
  metric_labels <- c(
    rmse_test = "Test RMSE",
    aic = "AIC",
    bic = "BIC",
    adj_r2_train = "Adjusted R2 (Train)"
  )
  
  multi_plot <- ggplot(selection_long, aes(x = k, y = value)) +
    geom_line(linewidth = 1, color = "steelblue") +
    geom_point(size = 2, color = "steelblue") +
    geom_vline(xintercept = final_k, linetype = "dashed", color = "red", alpha = 0.5) +
    facet_wrap(~ metric, scales = "free_y", ncol = 2,
               labeller = labeller(metric = metric_labels)) +
    labs(x = "Number of Predictors (k)", y = "Value") +
    theme_minimal(base_size = 12) +
    scale_x_continuous(breaks = seq(1, max_predictors, by = ifelse(max_predictors > 8, 2, 1))) +
    theme(strip.text = element_text(face = "bold"))
  
  ggsave("figs/model_selection_detailed.pdf", multi_plot, width = 10, height = 8)
  print(multi_plot)
  
} else {
  cat("\nNo significant predictors found.\n")
  cat("Only full model is presented.\n")
  model_reduced <- NULL
  predictions_reduced <- NULL
}

cat("\nTask 2 complete.\n\n")


## Generate diagnostic plots
pdf("figs/diagnostic_full.pdf", width = 10, height = 10)
par(mfrow = c(2, 2), cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.8, 
    mar = c(5, 5, 3, 2))
plot(model_full)
dev.off()

if (!is.null(model_reduced)) {
  pdf("figs/diagnostic_reduced.pdf", width = 10, height = 10)
  par(mfrow = c(2, 2), cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.8, 
      mar = c(5, 5, 3, 2))
  plot(model_reduced)
  dev.off()
}

par(mfrow = c(1, 1))


## Prediction plots
plot_full <- ggplot(
  data.frame(actual = test_data$vas_12months, predicted = predictions_full),
  aes(x = actual, y = predicted)
) +
  geom_point(alpha = 0.6, size = 3, color = "blue") +
  geom_abline(slope = 1, intercept = 0, color = "red", 
              linetype = "dashed", linewidth = 1) +
  labs(x = "Actual 12-month VAS", y = "Predicted 12-month VAS") +
  theme_minimal(base_size = 16) +
  coord_fixed(xlim = c(0, 10), ylim = c(0, 10)) +
  annotate("text", x = 2, y = 9,
           label = paste0("r = ", round(cor_full, 3), 
                          "\nRMSE = ", round(rmse_full, 2)),
           size = 5, color = "blue")

ggsave("figs/prediction_full.pdf", plot_full, width = 7, height = 7)
print(plot_full)

if (!is.null(predictions_reduced)) {
  plot_reduced <- ggplot(
    data.frame(actual = test_data$vas_12months, predicted = predictions_reduced),
    aes(x = actual, y = predicted)
  ) +
    geom_point(alpha = 0.6, size = 3, color = "darkgreen") +
    geom_abline(slope = 1, intercept = 0, color = "red", 
                linetype = "dashed", linewidth = 1) +
    labs(x = "Actual 12-month VAS", y = "Predicted 12-month VAS") +
    theme_minimal(base_size = 16) +
    coord_fixed(xlim = c(0, 10), ylim = c(0, 10)) +
    annotate("text", x = 2, y = 9,
             label = paste0("r = ", round(cor_reduced, 3), 
                            "\nRMSE = ", round(rmse_reduced, 2)),
             size = 5, color = "darkgreen")
  
  ggsave("figs/prediction_reduced.pdf", plot_reduced, width = 7, height = 7)
  print(plot_reduced)
}

cat("\nAnalysis complete.\n")