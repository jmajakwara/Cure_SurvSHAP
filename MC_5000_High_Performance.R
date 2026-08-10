rm(list=ls())
options(encoding = "UTF-8")

library(ggplot2)
library(MASS)
library(caret)
library(data.table)
library(survival)
library(rsample)
library(flexsurv)
library(flexsurvcure)
library(tidyverse)
library(censored)
library(tidymodels)
library(survex)
library(patchwork)

##### Housekeeping ######
cat("R Version:", R.version.string, "\n")
cat("Start time:", date(), "\n")

# ============================================================================
# DATA GENERATION FUNCTION
# ============================================================================

generate_logn_cure_data <- function(n = 5000, seed = NULL) {
  
  if(!is.null(seed)) set.seed(seed)
  
  # Covariates
  X1 <- sample(c("A", "B", "C", "D"), n, replace = TRUE, prob = c(0.04, 0.2, 0.62, 0.14)) 
  X2 <- pmin(pmax(rnorm(n, 0.18, 0.05), 0.05), 0.35)   # BorrowerRate - strong latency
  X3 <- pmin(rlnorm(n, log(0.15), 0.85), 1.0)
  X4 <- rpois(n, 8)
  X5 <- runif(n, 0, 1)
  X6 <- rpois(n, 1.2)
  X7 <- rbinom(n, 1, 0.55)
  X8 <- sample(c("1", "2", "3"), n, replace = TRUE, prob = c(0.13, 0.68, 0.19))
  
  X_inc <- model.matrix(~ X1)[, -1]
  X_lat <- cbind(X2, X3, X4, X5, X6, X7)
  
  # Stronger incidence separation
  b_true <- c(1.2, -0.6, -0.7, -0.9)   # Higher intercept → higher cure rate
  eta_inc <- b_true[1] + X_inc %*% b_true[-1]
  pi_susceptible <- 1 / (1 + exp(-eta_inc))
  susceptible <- rbinom(n, 1, pi_susceptible)
  
  # Latency with stronger effects + slight non-linearity
  g_true <- c(-5.5, -0.65, 0.04, 0.45, -0.08, -0.09)
  sigma <- 0.77
  mu0 <- 7.5
  
  log_T <- rep(NA_real_, n)
  idx <- which(susceptible == 1)
  # Add mild time interaction simulation by scaling with a factor
  log_T[idx] <- mu0 + X_lat[idx, ] %*% g_true +  rnorm(length(idx), 0, sigma) * 
                (2 + 1.5* X2[idx])  # mild interaction
  
  Y <- rep(Inf, n)
  Y[idx] <- exp(log_T[idx])
  
  # More realistic censoring
  C <- rexp(n, rate = 0.0018)
  time_obs <- pmin(Y, C)
  status <- as.numeric(Y <= C)
  
  data.frame(time = time_obs/2, status = status, X1 = X1, X2 = X2, 
             cure = 1 - susceptible, X3 = X3, X4 = X4, X5 = X5, 
             X6 = X6, X7 = X7)
}

# ============================================================================
# PREDICTION FUNCTION FOR CURE MODEL
# ============================================================================

predict_survival_cure <- function(object, newdata, times, train_data = NULL) {
  
  n <- nrow(newdata)
  nt <- length(times)
  
  s <- summary(
    object,
    newdata = newdata,
    t = times,
    type = "survival",
    tidy = TRUE,
    ci = FALSE
  )
  
  s$.id <- rep(seq_len(n), each = nt)
  
  mat <- s %>%
    dplyr::select(.id, time, est) %>%
    tidyr::pivot_wider(names_from = time, values_from = est) %>%
    dplyr::arrange(.id)
  
  result_matrix <- as.matrix(mat[, -1])
  colnames(result_matrix) <- NULL
  
  return(result_matrix)
}

# ============================================================================
# PREDICTION FUNCTION FOR STANDARD MODEL
# ============================================================================

predict_survival_standard <- function(object, newdata, times = NULL) {
  v <- augment(object, newdata, eval_time = times)
  return(v$.pred)
}

predict_risk <- function(object, newdata) {
  
  horizon <- median(sim_train$time[sim_train$status == 1])
  
  surv_mat <- predict_survival_cure(
    object,
    newdata,
    times = horizon
  )
  
  risk_score <- 1 - as.vector(surv_mat)
  
  return(risk_score)
}

# ============================================================================
# SAFE MODEL FITTING WRAPPER
# ============================================================================

safe_fit <- function(expr) {
  
  tryCatch(
    expr,
    error = function(e) {
      
      msg <- conditionMessage(e)
      
      if(grepl("system is exactly singular", msg) ||
         grepl("Lapack routine dgesv", msg)) {
        
        cat("  -> Singular Hessian encountered. Skipping model.\n")
        return(NULL)
      }
      
      cat("  -> Error:", msg, "\n")
      return(NULL)
    }
  )
}

# ============================================================================
# SAFE PERFORMANCE WRAPPER
# ============================================================================

safe_performance <- function(expr) {
  
  tryCatch(
    expr,
    error = function(e) {
      cat("  -> Performance calculation failed:", conditionMessage(e), "\n")
      return(NULL)
    }
  )
}

# ============================================================================
# SINGLE REPLICATION FUNCTION
# ============================================================================

run_single_replication <- function(rep_id,
                                   time_points = NULL,
                                   verbose = FALSE) {
  
  if(verbose) {
    cat("  Replication", rep_id, "starting...\n")
  }
  
  dat <- generate_logn_cure_data(
    n = 5000,
    seed = rep_id + 2000
  )
  
  # Split data
  data_split <- initial_split(
    dat,
    prop = 0.7,
    strata = status
  )
  
  sim_train <- training(data_split)
  sim_test  <- testing(data_split)
  
  X_train <- sim_train %>%
    dplyr::select(X1, X2, X3, X4, X5, X6, X7)
  
  X_test <- sim_test %>%
    dplyr::select(X1, X2, X3, X4, X5, X6, X7)
  
  time_points <- seq(
    0,
    max(sim_test$time),
    length.out = 100
  )
  
  # =========================================================================
  # STANDARD LOGNORMAL MODEL
  # =========================================================================
  
  log_spec <- survival_reg(dist = "lnorm") |>
    set_engine("flexsurv") |>
    set_mode("censored regression")
  
  log_norm <- safe_fit(
    log_spec |> fit(
      Surv(time, status) ~ X1 + X2 + X3 + X4 + X5 + X6 + X7,
      data = sim_train
    )
  )
  
  # =========================================================================
  # CURE MODELS
  # =========================================================================
  
  fit_lnorm <- safe_fit(
    flexsurvcure(
      Surv(time, status) ~ X1,
      anc = list(meanlog = ~ X2 + X3 + X4 + X5 + X6 + X7),
      data = sim_train,
      dist = "lnorm",
      link = "logistic",
      mixture = TRUE
    )
  )
  
  fit_weibull <- safe_fit(
    flexsurvcure(
      Surv(time, status) ~ X1,
      anc = list(scale = ~ X2 + X3 + X4 + X5 + X6 + X7),
      data = sim_train,
      dist = "weibull",
      link = "logistic",
      mixture = TRUE
    )
  )
  
  fit_exponential <- safe_fit(
    flexsurvcure(
      Surv(time, status) ~ X1,
      anc = list(rate = ~ X2 + X3 + X4 + X5 + X6 + X7),
      data = sim_train,
      dist = "exp",
      link = "logistic",
      mixture = TRUE
    )
  )
  
  fit_loglogistic <- safe_fit(
    flexsurvcure(
      Surv(time, status) ~ X1,
      anc = list(scale = ~ X2 + X3 + X4 + X5 + X6 + X7),
      data = sim_train,
      dist = "llogis",
      link = "logistic",
      mixture = TRUE
    )
  )
  
  # =========================================================================
  # INITIALIZE RESULTS
  # =========================================================================
  
  results <- list(
    
    rep_id = rep_id,
    
    standard_c_index = NA_real_,
    standard_ibs = NA_real_,
    
    lnorm_cure_c_index = NA_real_,
    lnorm_cure_ibs = NA_real_,
    
    weibull_cure_c_index = NA_real_,
    weibull_cure_ibs = NA_real_,
    
    exponential_cure_c_index = NA_real_,
    exponential_cure_ibs = NA_real_,
    
    loglogistic_cure_c_index = NA_real_,
    loglogistic_cure_ibs = NA_real_
  )
  
  # =========================================================================
  # STANDARD MODEL PERFORMANCE
  # =========================================================================
  
  if(!is.null(log_norm)) {
    
    explainer_st <- survex::explain(
      model = log_norm,
      data = X_train,
      y = Surv(sim_train$time, sim_train$status),
      predict_survival = predict_survival_standard,
      times = time_points,
      label = "Standard Model",
      verbose = FALSE
    )
    
    perf_st <- safe_performance(
      model_performance(
        explainer = explainer_st,
        data = X_test,
        y = Surv(sim_test$time, sim_test$status),
        type = "metrics",
        metrics = c(c_index, integrated_brier_score),
        times = time_points
      )
    )
    
    if(!is.null(perf_st)) {
      results$standard_c_index <- as.numeric(perf_st$result[[1]])
      results$standard_ibs <- as.numeric(perf_st$result[[2]])
    }
  }
  
  # =========================================================================
  # FUNCTION TO EVALUATE CURE MODELS
  # =========================================================================
  
  evaluate_cure_model <- function(fit_object, model_name) {
    
    if(is.null(fit_object)) {
      return(c(NA_real_, NA_real_))
    }
    
    explainer_cure <- survex::explain_survival(
      model = fit_object,
      data = X_train,
      y = Surv(sim_train$time, sim_train$status),
      predict_survival_function = predict_survival_cure,
      predict_risk_function = predict_risk,
      times = time_points,
      label = model_name,
      verbose = FALSE
    )
    
    perf <- safe_performance(
      model_performance(
        explainer = explainer_cure,
        data = X_test,
        y = Surv(sim_test$time, sim_test$status),
        type = "metrics",
        metrics = c(c_index, integrated_brier_score),
        times = time_points
      )
    )
    
    if(is.null(perf)) {
      return(c(NA_real_, NA_real_))
    }
    
    return(c(
      as.numeric(perf$result[[1]]),
      as.numeric(perf$result[[2]])
    ))
  }
  
  # =========================================================================
  # EVALUATE ALL CURE MODELS
  # =========================================================================
  
  tmp <- evaluate_cure_model(fit_lnorm, "Lognormal Cure")
  results$lnorm_cure_c_index <- tmp[1]
  results$lnorm_cure_ibs <- tmp[2]
  
  tmp <- evaluate_cure_model(fit_weibull, "Weibull Cure")
  results$weibull_cure_c_index <- tmp[1]
  results$weibull_cure_ibs <- tmp[2]
  
  tmp <- evaluate_cure_model(fit_exponential, "Exponential Cure")
  results$exponential_cure_c_index <- tmp[1]
  results$exponential_cure_ibs <- tmp[2]
  
  tmp <- evaluate_cure_model(fit_loglogistic, "Loglogistic Cure")
  results$loglogistic_cure_c_index <- tmp[1]
  results$loglogistic_cure_ibs <- tmp[2]
  
  return(results)
}

# ============================================================================
# MONTE CARLO SIMULATION
# ============================================================================

run_monte_carlo <- function(n_sim = 1000,
                            verbose = TRUE,
                            save_intermediate = TRUE) {
  
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("MONTE CARLO SIMULATION\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  cat("Number of simulations:", n_sim, "\n")
  cat("Start time:", date(), "\n\n")
  
  all_results <- vector("list", n_sim)
  
  for(i in 1:n_sim) {
    
    if(verbose && i %% 50 == 0) {
      cat(sprintf(
        "Progress: %d/%d simulations completed (%.1f%%)\n",
        i, n_sim, 100 * i / n_sim
      ))
    }
    
    all_results[[i]] <- run_single_replication(
      rep_id = i + 1000,
      verbose = FALSE
    )
  }
  
  return(list(
    raw_results = all_results,
    n_sim = n_sim,
    timestamp = date()
  ))
}

# ============================================================================
# SUMMARY TABLE
# ============================================================================

create_summary_table <- function(full_results) {
  
  raw <- bind_rows(full_results$raw_results)
  
  summary_table <- tibble(
    Model = c(
      "Standard Lognormal",
      "Lognormal Cure",
      "Weibull Cure",
      "Exponential Cure",
      "Loglogistic Cure"
    ),
    
    C_Index_Mean = c(
      mean(raw$standard_c_index, na.rm = TRUE),
      mean(raw$lnorm_cure_c_index, na.rm = TRUE),
      mean(raw$weibull_cure_c_index, na.rm = TRUE),
      mean(raw$exponential_cure_c_index, na.rm = TRUE),
      mean(raw$loglogistic_cure_c_index, na.rm = TRUE)
    ),
    
    C_Index_SD = c(
      sd(raw$standard_c_index, na.rm = TRUE),
      sd(raw$lnorm_cure_c_index, na.rm = TRUE),
      sd(raw$weibull_cure_c_index, na.rm = TRUE),
      sd(raw$exponential_cure_c_index, na.rm = TRUE),
      sd(raw$loglogistic_cure_c_index, na.rm = TRUE)
    ),
    
    IBS_Mean = c(
      mean(raw$standard_ibs, na.rm = TRUE),
      mean(raw$lnorm_cure_ibs, na.rm = TRUE),
      mean(raw$weibull_cure_ibs, na.rm = TRUE),
      mean(raw$exponential_cure_ibs, na.rm = TRUE),
      mean(raw$loglogistic_cure_ibs, na.rm = TRUE)
    ),
    
    IBS_SD = c(
      sd(raw$standard_ibs, na.rm = TRUE),
      sd(raw$lnorm_cure_ibs, na.rm = TRUE),
      sd(raw$weibull_cure_ibs, na.rm = TRUE),
      sd(raw$exponential_cure_ibs, na.rm = TRUE),
      sd(raw$loglogistic_cure_ibs, na.rm = TRUE)
    )
  )
  
  numeric_cols <- names(summary_table)[-1]
  
  summary_table[numeric_cols] <- lapply(
    summary_table[numeric_cols],
    round,
    4
  )
  
  return(summary_table)
}

# ============================================================================
# FULL SIMULATION
# ============================================================================

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("FULL MONTE CARLO SIMULATION (1000 iterations)\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

full_results <- run_monte_carlo(
  n_sim = 1000,
  verbose = TRUE,
  save_intermediate = TRUE
)

# ============================================================================
# CREATE FINAL SUMMARY TABLE
# ============================================================================

summary_table <- create_summary_table(full_results)

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("FINAL SUMMARY TABLE\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

print(summary_table)



cat("\nEnd time:", date(), "\n")