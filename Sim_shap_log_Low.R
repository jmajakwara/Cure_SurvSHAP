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

generate_logn_cure_data <- function(n = 1000, seed = 2026) {
  set.seed(seed)
 
  # === 1. Generate Mildly Correlated Latent Variables ===
  rho <- 0.15
  cor_matrix <- matrix(rho, nrow = 4, ncol = 4)
  diag(cor_matrix) <- 1
 
  means <- rep(0, 4)
  sds <- rep(1, 4)
  cov_matrix <- diag(sds) %*% cor_matrix %*% diag(sds)
 
  Y <- MASS::mvrnorm(n = n, mu = means, Sigma = cov_matrix)
  Y <- as.data.frame(Y)
  colnames(Y) <- paste0("Y", 1:4)
 
  # === 2. Transform to Target Distributions ===
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
  
  X5 = rpois(n, lambda = 2.2)
 
  # === 3. STRONGER Incidence (X1) - This is the key change ===
  X_inc <- model.matrix(~ X1, data = X)[, -1]
 
  b_true <- c(0.3, -1.25, -1.85, -2.65) 
  eta_inc <- b_true[1] + X_inc %*% b_true[-1]
  pi_cure <- 1 / (1 + exp(-eta_inc))
  cure <- rbinom(n, 1, pi_cure)             # 1 = cured
 
  # === 4. Latency with Smoother Interaction ===
  X_lat <- cbind(X$X2, X$X3, X$X4)
 
  g_true <- c(-5.65, 0.75, 0.3)
  sigma <- 1.1
  mu0 <- 7.58
 
 
  log_T <- rep(NA_real_, n)
  idx <- which(cure == 0)
  log_T[idx] <- mu0 + X_lat[idx, ] %*% g_true + rnorm(length(idx), 0, sigma)
 
  Y_time <- rep(Inf, n)
  Y_time[idx] <- exp(log_T[idx])
 
  # Realistic censoring (adjusted for better time scale)
  C <- rexp(n, rate = 0.0001)
  time_obs <- pmin(Y_time, C)
  status <- as.numeric(Y_time <= C)
 
  data.frame(
    time = time_obs / 25,          # Better time scaling
    status = status,
    X1 = X$X1,
    X2 = X$X2,
    X3 = X$X3,
    X4 = X$X4,
    X5 = X5,
    cure = cure
  )
}


# ========================= RUN IT =========================
seed=2000
set.seed(seed)
dat <- generate_logn_cure_data(n = 1000, seed = seed)

summary(dat$time)
prop.table(table(dat$status))
prop.table(table(dat$cure))

data_split <- initial_split(dat, prop = 0.7, strata = status)     
sim_train <- training(data_split)
sim_test  <- testing(data_split)

X_train <- sim_train %>% dplyr::select(X1, X2, X3, X4, X5)
X_test  <- sim_test  %>% dplyr::select(X1, X2, X3, X4, X5)

sim_train <- as_tibble(sim_train)
sim_test <- as.data.frame(sim_test)

###### Lognormal Modeling ########

cat('The Lognormal standard Model using flexsurv package  \n')
log_spec <- survival_reg(dist = "lnorm") |>
  set_engine("flexsurv") |>
  set_mode("censored regression") 

log_norm <- log_spec |> fit(Surv(time, status) ~ X1 + X2 + X3 + X4 + X5, 
                            data = sim_train)
tidy(log_norm)
glance(log_norm)


cat('The Lognormal cure Model using flexsurvcure package  \n')

log_fit <- flexsurvcure::flexsurvcure(Surv(time, status) ~ X1, 
                                      anc = list(meanlog = ~ X2 + X3 + X4 + X5),
                                      data = sim_train,
                                      dist = 'lnorm',
                                      link = 'logistic',  
                                      mixture = TRUE)
tidy(log_fit)
glance(log_fit)


time_points <- seq(0, max(sim_test$time), length.out = 100)

# --------------------- PREDICT SURVIVAL FUNCTION (SURV_MATRIX FORMAT) ---------------------
# According to survex documentation, predict_survival_function must return a matrix
# where rows correspond to observations and columns to times

predict_survival <- function(object, newdata, times = NULL) {
  if (is.null(times)) times <- time_points
  
  n <- nrow(newdata)
  nt <- length(times)
  
  # Get predictions
  s <- summary(object, 
               newdata = newdata, 
               t = times, 
               type = "survival", 
               tidy = TRUE, 
               ci = FALSE)
  
  # CRITICAL FIX: Create id column manually
  # Each observation has nt rows, in sequence
  s$.id <- rep(seq_len(n), each = nt)
  
  # Reshape to wide format
  mat <- s %>%
    dplyr::select(.id, time, est) %>%
    tidyr::pivot_wider(names_from = time, values_from = est) %>%
    dplyr::arrange(.id)
  
  # Convert to matrix and remove id column
  result_matrix <- as.matrix(mat[, -1])
  colnames(result_matrix) <- NULL
  
  return(result_matrix)
}

# --------------------- OPTIONAL: RISK PREDICTION FUNCTION ---------------------
predict_risk <- function(object, newdata) {
  # Use median event time from training data as horizon
  horizon <- median(sim_train$time[sim_train$status == 1])
  surv_mat <- predict_survival(object, newdata, times = horizon)
  risk_score <- 1 - as.vector(surv_mat)
  return(risk_score)
}
# --------------------- EXPLAINER FOR CURE MODEL ---------------------
# For custom models, use explain_survival() or use explain() with predict_survival_function argument

# Method 1: Using explain_survival() - Recommended for custom models
explainer <- survex::explain_survival(
  model = log_fit,
  data = X_train,
  y = Surv(sim_train$time, sim_train$status),
  predict_survival_function = predict_survival,
  predict_risk_function = predict_risk,
  times = time_points,
  label = "Lognormal Cure Model"
)

# Method 2: Using explain() with predict_survival_function argument
# (Alternative approach if explain_survival() not available)
explainer_alt <- survex::explain(
  model = log_fit,
  data = X_train,
  y = Surv(sim_train$time, sim_train$status),
  predict_survival_function = predict_survival,
  predict_risk_function = predict_risk,
  times = time_points,
  label = "Lognormal Cure Model"
)

# --------------------- EXPLAINER FOR STANDARD MODEL ---------------------
# For standard flexsurv model, we can use the automated interface
# or create a custom prediction function

predict_standard_survival <- function(object, newdata, times = NULL) {
   v = augment(object, newdata, eval_time = times)
   return(v$.pred)
}

explainer_st <- survex::explain(
  model = log_norm,
  data = X_train,
  y = Surv(sim_train$time, sim_train$status),
  predict_survival = predict_standard_survival,
  times = time_points,
  label = "Standard Lognormal Model"
)

# --------------------- MODEL DIAGNOSTICS ---------------------
cure_residuals <- model_diagnostics(explainer)
standard_residuals <- model_diagnostics(explainer_st)

# Plot diagnostics - note: plot_type can be "deviance", "martingale", or "Cox-Snell"
# The xvariable parameter can be any column name in the result
plot(cure_residuals, xvariable = "X2")  # Plot against BorrowerRate
plot(cure_residuals, plot_type = "Cox-Snell")  # Cox-Snell residual diagnostic plot
plot(cure_residuals, plot_type = "martingale")
plot(standard_residuals, plot_type = "martingale")

# --------------------- MODEL PERFORMANCE ---------------------
# Use unquoted function names for metrics
perf_cure <- model_performance(
  explainer = explainer, 
  type = "metrics",
  metrics = c(c_index, integrated_brier_score, brier_score),
  times = time_points
)

perf_standard <- model_performance(
  explainer = explainer_st, 
  type = "metrics",
  metrics = c(c_index, integrated_brier_score, brier_score),
  times = time_points
)

# Print results
cat("The C-index of standard model is", perf_standard$result[[1]])
cat("The Integrated Brier Score of standard model is", perf_standard$result[[2]])
cat("The C-index of cure model is", perf_cure$result[[1]])
cat("The Integrated Brier Score of cure model is", perf_cure$result[[2]])
plot(perf_cure$result[[3]])
lines(perf_standard$result[[3]])


# --------------------- ADDITIONAL: SURVSHAP(t) ANALYSIS ---------------------
# Calculate SurvSHAP(t) values using model_parts() or model_survshap()
# Note: This may be computationally intensive

# For global feature importance with SurvSHAP(t)
shap_result <- predict_parts(
  explainer,
  new_observation = X_test,
  type = "survshap",        
  N = 10,
  output_type = "survival",
  calculation_method = "kernelshap",
  aggregation_method = "integral" #"mean_absolute"  #
)

plot(shap_result, geom = "importance",max_vars = 8) + 
  ggtitle("Global Feature Importance over Time") +
  theme_minimal()

plot(shap_result, geom = "beeswarm",max_vars = 8) + 
#  ggtitle("SurvSHAP(t) - Lognormal Mixture Cure Model") +
  theme_minimal()