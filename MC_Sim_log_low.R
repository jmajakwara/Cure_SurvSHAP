rm(list=ls())
options(encoding = "UTF-8")

library(ggplot2)
library(MASS)
library(caret)
library(reticulate)
library(data.table)
library(survival)
library(rsample)
library(flexsurv)
library(flexsurvcure)
library(tidyverse)
library(censored)
library(tidymodels)
library(survex)


logit <- function(p) log(p / (1 - p))

generate_logn_cure_data <- function(n = 1000, seed = 2000, scenario = "low_cure") {
  set.seed(seed)
  
  # Generate Mildly Correlated Latent Variables ===
  rho <- 0.15
  cor_matrix <- matrix(rho, nrow = 4, ncol = 4)
  diag(cor_matrix) <- 1
  
  means <- rep(0, 4)
  sds <- rep(1, 4)
  cov_matrix <- diag(sds) %*% cor_matrix %*% diag(sds)
  
  Y <- MASS::mvrnorm(n = n, mu = means, Sigma = cov_matrix)
  Y <- as.data.frame(Y)
  colnames(Y) <- paste0("Y", 1:4)
  
  # =Transform to Target Distributions ===
  X <- Y |>
    mutate(
      X1 = case_when(
        Y1 <= qnorm(0.04) ~ "A",
        Y1 <= qnorm(0.24) ~ "B",
        Y1 <= qnorm(0.86) ~ "C",
        TRUE ~ "D"
      ),
      X2 = pmin(pmax(exp(Y2 * 0.45 + log(0.18)), 0.05), 0.35),
      X3 = pmin(plogis(Y3 * 0.75 + logit(0.15)), 1.0),
      X4 = pmin(pmax(Y4 * 1/6 + 0.5, 0), 1)
    ) |>
    mutate(X1 = factor(X1, levels = c("A","B","C","D")))
  
  X5 = rpois(n, lambda = 1.5)
  
  # Incidence Submodel - Different Scenarios ===
  X_inc <- model.matrix(~ X1, data = X)[, -1]
  
  b_true <- switch(scenario,
    high_cure = c(1.2, -1.1, -1.6, -2.45),   
    medium_cure = c(0.95, -1.15, -1.35, -2.15),
    low_cure = c(0.3, -1.25, -1.85, -2.35)   
  )
  
  eta_inc <- b_true[1] + X_inc %*% b_true[-1]
  pi_cure <- 1 / (1 + exp(-eta_inc))            # Probability of being cured
  cure <- rbinom(n, 1, pi_cure)                 # 1 = cured, 0 = susceptible
  
  # Latency Submodel ===
  X_lat <- cbind(X$X2, X$X3, X$X4)
  
  g_true <- c(-5.65, 0.75, 0.3)
  sigma <- 1.1
  mu0 <- 7.58
  

  
  log_T <- rep(NA_real_, n)
  idx <- which(cure == 0)   # Only susceptible individuals
  log_T[idx] <- mu0 + X_lat[idx, ] %*% g_true + 
                rnorm(length(idx), 0, sigma)
  
  Y_time <- rep(Inf, n)
  Y_time[idx] <- exp(log_T[idx])
  
  C <- rexp(n, rate = 0.0001)
  time_obs <- pmin(Y_time, C)
  status <- as.numeric(Y_time <= C)
  
  data.frame(
    time = time_obs / 25,
    status = status,
    X1 = X$X1,
    X2 = X$X2,
    X3 = X$X3,
    X4 = X$X4,
    X5 = X5,
    cure = cure
  )
}
# ====================== MONTE CARLO SIMULATION ======================

monte_carlo_cure <- function(B = 1000, n = 1000, scenario = "low_cure", seed = 2000) {
  
  set.seed(seed)
  results_list <- list()
  convergence_status <- data.frame(replication = 1:B, converged = FALSE, 
                                   loglik = NA_real_, stringsAsFactors = FALSE)
  

  
  true_values <- switch(scenario,
    high_cure = list(X1B = -1.1, X1C = -1.6, X1D = -2.45,`meanlog(X2)` = -5.65,
            `meanlog(X3)` = 0.75,`meanlog(X4)` = 0.3), 
    medium_cure = list(X1B = -1.15,X1C = -1.35,X1D = -2.15,`meanlog(X2)` = -5.65,
      `meanlog(X3)` = 0.75,`meanlog(X4)` = 0.3),
    low_cure = list(X1B = -1.25,X1C = -1.85,X1D = -2.35,`meanlog(X2)` = -5.65,
      `meanlog(X3)` = 0.75,`meanlog(X4)` = 0.3)
  )

   
 
  pb <- txtProgressBar(min = 0, max = B, style = 3)
  successful_runs <- 0
  
  for(b in 1:B) {
    
    # Generate data with unique seed per replication
    dat <- generate_logn_cure_data(n = n, seed = b + 10000, scenario = scenario)
    
    # Fit model with tryCatch for both errors and warnings
    fit <- tryCatch({
      flexsurvcure::flexsurvcure(
        Surv(time, status) ~ X1,
        anc = list(meanlog = ~ X2 + X3 + X4),
        data = dat,
        dist = "lnorm",
        link = "logistic",
        mixture = TRUE
      )
    }, error = function(e) {
      # Uncomment to debug: cat("\nError in rep", b, ":", e$message, "\n")
      return(NULL)
    })
    
    # Check if fit is valid
    if(!is.null(fit) && is.null(fit$fail) && is.finite(fit$loglik)) {
      successful_runs <- successful_runs + 1
      convergence_status$converged[b] <- TRUE
      convergence_status$loglik[b] <- fit$loglik
      
      # Extract tidy results
      est_table <- broom::tidy(fit)
      
      # Filter for covariate parameters (exclude distributional parameters)
      # Note: Column names in tidy() are typically: term, estimate, std.error, statistic, p.value
      cov_params <- est_table %>%
        filter(!term %in% c("meanlog", "sdlog", "theta")) %>%
        mutate(
          replication = b,
          scenario = scenario,
          std.error = ifelse("std.error" %in% names(.), std.error, NA_real_)
        )
      
      results_list[[b]] <- cov_params
    } else {
      convergence_status$converged[b] <- FALSE
      convergence_status$loglik[b] <- NA
    }
    
    setTxtProgressBar(pb, b)
  }
  
  close(pb)
  
  # Combine successful results only
  sim_results <- bind_rows(results_list[!sapply(results_list, is.null)])
  
  if(nrow(sim_results) == 0) {
    cat("\nERROR: No successful model convergences!\n")
    return(list(raw = NULL, metrics = NULL, convergence = convergence_status))
  }
  
  # ====================== ADD TRUE VALUES TO RESULTS ======================
  sim_results <- sim_results %>%
    mutate(true_value = case_when(
      term == "X1B" ~ true_values$X1B,
      term == "X1C" ~ true_values$X1C,
      term == "X1D" ~ true_values$X1D,
      term == "meanlog(X2)" ~ true_values$`meanlog(X2)`,
      term == "meanlog(X3)" ~ true_values$`meanlog(X3)`,
      term == "meanlog(X4)" ~ true_values$`meanlog(X4)`,
      # Alternative naming from broom::tidy()
      term == "X2" ~ true_values$`meanlog(X2)`,
      term == "X3" ~ true_values$`meanlog(X3)`,
      term == "X4" ~ true_values$`meanlog(X4)`,
      TRUE ~ NA_real_
    ))
  
  # ====================== PERFORMANCE METRICS ======================
  metrics <- sim_results %>%
    group_by(term) %>%
    summarise(
      true_value = first(true_value),
      n_successful = n(),
      mean_est = mean(estimate, na.rm = TRUE),
      bias = mean_est - true_value,
      percent_bias = ifelse(true_value != 0, bias / abs(true_value) * 100, NA_real_),
      MSE = mean((estimate - true_value)^2, na.rm = TRUE),
      RMSE = sqrt(MSE),
      sd_est = sd(estimate, na.rm = TRUE),
      # Only compute coverage if std.error available
      coverage_90 = if(all(!is.na(std.error))) {
        mean(abs(estimate - true_value) <= qnorm(0.95) * std.error, na.rm = TRUE)
      } else NA_real_,
      coverage_95 = if(all(!is.na(std.error))) {
        mean(abs(estimate - true_value) <= qnorm(0.975) * std.error, na.rm = TRUE)
      } else NA_real_,
      .groups = "drop"
    ) %>%
    arrange(term)
  
  # ====================== PRINT SUMMARY ======================
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("MONTE CARLO SIMULATION RESULTS\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  cat("Scenario:", scenario, "\n")
  cat("Attempted replications:", B, "\n")
  cat("Successful convergences:", successful_runs, 
      sprintf("(%.1f%%)\n", 100 * successful_runs / B))
  cat("\n")
  print(metrics, n = Inf)
  
  # ====================== CONVERGENCE REPORT ======================
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("CONVERGENCE DIAGNOSTICS\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  cat("Log-likelihood summary for converged models:\n")
  print(summary(convergence_status$loglik[convergence_status$converged]))
  
  # Return full results
  return(list(
    raw = sim_results, 
    metrics = metrics, 
    convergence = convergence_status,
    successful_runs = successful_runs,
    total_runs = B
  ))
}

# ====================== RUN CORRECTED SIMULATION ======================
cat("\nRunning Monte Carlo Simulation for High Cure Scenario...\n")
set.seed(2000)
mc_results_high <- monte_carlo_cure(B = 1000, n = 3000, scenario = "low_cure", seed = 2000)
print(mc_results_high)



# ====================== PLOT RESULTS ======================
if(!is.null(mc_results_high$metrics)) {
  # Bias plot
  bias_plot <- mc_results_high$metrics %>%
    filter(!is.na(percent_bias)) %>%
    ggplot(aes(x = term, y = percent_bias, fill = term)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(title = paste("Parameter Bias -", "High Cure Scenario"),
         x = "Parameter", y = "Percent Bias (%)") +
    theme_minimal() +
    theme(legend.position = "none")
  
  print(bias_plot)
  
  # Convergence plot
  conv_df <- data.frame(
    replication = 1:length(mc_results_high$convergence$converged),
    converged = mc_results_high$convergence$converged
  )
  
  conv_plot <- ggplot(conv_df, aes(x = replication, y = as.numeric(converged))) +
    geom_point(alpha = 0.5, size = 0.8) +
    geom_smooth(se = FALSE, method = "loess") +
    labs(title = "Convergence by Replication",
         x = "Replication Number", y = "Converged (1 = Yes)") +
    theme_minimal()
  
  print(conv_plot)
}
