rm(list=ls())
options(encoding = "UTF-8")

library(ggplot2)
library(MASS)
library(caret)
library(reticulate)
library(data.table)
library(survival)
library(survminer)
library(rsample)
library(flexsurv)
library(dplyr)
library(flexsurvcure)
library(tidyverse)
library(censored)
library(tidymodels)
library(survex)

##### Housekeeping ######
R.version.string
date()
Sys.time()




# ===============================
# DATA PREPARATION 
# ===============================
mortgage_data  <- read.csv("C:\\Users\\a0025724\\OneDrive - University of Witwatersrand\\RESEARCH\\Credit Risk Modeling\\SHAP\\Code\\mortgage.csv")
str(mortgage_data)

mortgage <- mortgage_data %>%
  dplyr::group_by(id) %>%
  dplyr::summarise(
    # Extract static (time-invariant) covariates at origination
    FICO = first(FICO_orig_time),
    LTV = first(LTV_orig_time),
    HPI = first(hpi_orig_time),
    GDP = first(gdp_time),
    UER = first(uer_time),
    Interest_rate = first(Interest_Rate_orig_time),
    balance = first(balance_orig_time),
    investor = first(investor_orig_time),
    REtype_SF = first(REtype_SF_orig_time),
    REtype_CO = first(REtype_CO_orig_time),
    REtype_PU = first(REtype_PU_orig_time),
    
    first_observed = min(time),
    last_observed = max(time),
    
    event_default = as.numeric(any(default_time == 1)),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    time_to_event = last_observed - first_observed,
    status = dplyr::case_when(
      event_default == 1 ~ 1,
      TRUE ~ 0  # Censored (includes loans that paid off or are still active at study end)
    )
  ) |>
  # Select the final columns for modeling
  dplyr::select(id, time_to_event, status, 
                FICO, LTV, GDP, UER, HPI, Interest_rate, balance, 
                investor, REtype_SF, REtype_CO, REtype_PU)

mortgage <- mortgage |> dplyr::filter(time_to_event > 0) |>
  dplyr::mutate(
    FICO_score = dplyr::case_when(
      FICO >= 400 & FICO <= 659 ~ "Fair",
      FICO >= 660 & FICO <= 739 ~ "Good",
      FICO >= 740 & FICO <= 799 ~ "Very good",
      FICO >= 800 & FICO <= 899 ~ "Excellent",
      TRUE ~ NA_character_  # Values outside these ranges become NA
    )
  )

mortgage <- mortgage |> dplyr::mutate(
    dplyr::across(c(FICO,LTV,UER,HPI,GDP, Interest_rate, balance), as.numeric),
    dplyr::across(c(FICO_score, REtype_SF, REtype_CO, REtype_PU, investor), as.factor)
  )
summary(mortgage)

cat("\n--- Survival Data Summary ---\n")
cat("Number of loans:", nrow(mortgage), "\n")
cat("Default rate:", mean(mortgage$status), "\n")
cat("Censoring rate:", 1 - mean(mortgage$status), "\n")

head(mortgage, n=10)
str(mortgage)

# Checking correlation of variables
X <- mortgage |> dplyr::select(FICO,LTV,UER,HPI,GDP,Interest_rate,balance)
cor(as.data.frame(X))
# Note that HPI and GDP has a high correlation of 0.46, therefore one of the variable will be dropped (HPI)
# Also FICO and balance is high of 0.32, therefore balance will be dropped

# Distributions of variables
hist(mortgage$time_to_event)
hist(mortgage$time_to_event[mortgage$status == 0],
     main = "Histogram of Time to Event (Censored)",
     xlab = "Time to Event",
     col = "lightblue")
hist(mortgage$FICO)
hist(mortgage$LTV)
hist(mortgage$GDP)
hist(mortgage$UER)


# ===============================
# DATA SPLITTING FIRST
# ===============================

set.seed(2000)
data_split <- initial_split(mortgage, prop = 0.7, strata = status)  
mortgage_train <- training(data_split)
mortgage_test  <- testing(data_split)

X_train <- mortgage_train %>% dplyr::select("FICO_score", "LTV", "UER", "GDP")
X_test  <- mortgage_test  %>% dplyr::select("FICO_score", "LTV", "UER", "GDP")

mortgage_train <- as_tibble(mortgage_train)
mortgage_test <- as.data.frame(mortgage_test)

prop.table(table(mortgage$FICO_score))

fit <- survfit(Surv(time_to_event, status) ~ FICO_score, data = mortgage)
p <- ggsurvplot(fit, data = mortgage, pval = TRUE, pval.method = TRUE, pval.coord = c(0, 0.15),
           pval.method.coord = c(0, 0.3),  risk.table = TRUE, linetype = "strata",
           legend = "right", risk.table.fontsize = 3,ggtheme = theme_minimal(),
           title = "Kaplan-Meier Survival Curves by FICO_score",
           xlab = "Time (months)", ylab = "Survival Probability")
p

ggsave("mort_km_plot.pdf", plot = p, width = 10, height = 8, dpi = 300)



# VIF greater than 5 indicates significant multicollinearity.
cat("calculate VIF for each predictor variable, might need to consider those VIF > 2.5.......\n")
cox_model <- coxph(Surv(time_to_event, status) ~ FICO + LTV + GDP + UER + HPI, data = mortgage, y=TRUE, x = TRUE)
vif(cox_model)


cat("Fitting the Standard survival AFT models...........\n")

logl_spec <- survival_reg(dist = "llogis") |>
  set_engine("flexsurv") |>
  set_mode("censored regression") 


logl_fit <- logl_spec |> fit(Surv(time_to_event, status) ~ FICO_score + LTV + UER + GDP, data = mortgage)



cat("Results from complete data \n")
tidy(logl_fit)
glance(logl_fit)

log_spec <- survival_reg(dist = "lnorm") |>
  set_engine("flexsurv") |>
  set_mode("censored regression") 


log_st <- log_spec |> fit(Surv(time_to_event, status) ~ FICO_score + LTV + UER + GDP, data = mortgage_train)



cat("Results from complete data \n")
tidy(log_st)
glance(log_st)

exp_spec <-  survival_reg(dist = "exp") |>
   set_engine("flexsurv") |>
   set_mode("censored regression") 

exp_fit <- exp_spec |> fit(Surv(time_to_event, status) ~ FICO_score + LTV + UER + GDP + Interest_rate, data = mortgage)
tidy(exp_fit)
glance(exp_fit)

wei_spec <-  survival_reg(dist = "weibull") |>
   set_engine("flexsurv") |>
   set_mode("censored regression") 

wei_fit <- wei_spec |> fit(Surv(time_to_event, status) ~ FICO_score + LTV + UER + GDP + Interest_rate, data = mortgage)
tidy(wei_fit)
glance(wei_fit)





log_fit <- flexsurvcure(Surv(time_to_event,status) ~ FICO_score,
                           anc=list(meanlog = ~ LTV + UER + GDP),  
			   data = mortgage_train,
			   link = "logistic",   # incidence (cure probability) on logit scale
  			 dist = "lnorm",     # latency = lognormal AFT
  			 mixture = TRUE       
)

print(tidy(log_fit))
print(glance(log_fit))



cure_par <- log_fit$res["theta", "est"]
cat("Cure fraction estimate:", exp(cure_par)/(1+exp(cure_par)), "\n")

time_points <- sort(unique(mortgage_test$time_to_event))

time_points <- seq(0, max(mortgage_test$time), length.out = 100)

# --------------------- PREDICT SURVIVAL FUNCTION (SURV_MATRIX FORMAT) ---------------------
# predict_survival_function must return a matrix
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
  
  # Create id column manually each observation has nt rows, in sequence
  s$.id <- rep(seq_len(n), each = nt)
  
  # Reshape to wide format
  mat <- s |>
    dplyr::select(.id, time, est) |>
    tidyr::pivot_wider(names_from = time, values_from = est) |>
    dplyr::arrange(.id)
  
  # Convert to matrix and remove id column
  result_matrix <- as.matrix(mat[, -1])
  colnames(result_matrix) <- NULL
  return(result_matrix)
}

# --------------------- RISK PREDICTION FUNCTION ---------------------
predict_risk <- function(object, newdata) {
  # Using median event time from training data as horizon
  horizon <- median(sim_train$time[sim_train$status == 1])
  surv_mat <- predict_survival(object, newdata, times = horizon)
  risk_score <- 1 - as.vector(surv_mat)
  return(risk_score)
}
# --------------------- EXPLAINER FOR CURE MODEL ---------------------
explainer <- survex::explain_survival(
  model = log_fit,
  data = X_train,
  y = Surv(mortgage_train$time_to_event, mortgage_train$status),
  predict_survival_function = predict_survival,
  predict_risk_function = predict_risk,
  times = time_points,
  label = "Lognormal Cure Model"
)



# --------------------- SURVSHAP(t) ANALYSIS ---------------------
# Calculate SurvSHAP(t) values using model_parts() or model_survshap()
# For global feature importance with SurvSHAP(t)
shap_result <- predict_parts(
  explainer,
  new_observation = X_test,
  type = "survshap",        
  N = 100,
  output_type = "survival",
  calculation_method = "kernelshap",
  aggregation_method = "integral" 
)

plot(shap_result, geom = "importance",max_vars = 8) + 
  ggtitle("Global Feature Importance over Time") +
  theme_minimal()

plot(shap_result, geom = "beeswarm",max_vars = 8) + 
  ggtitle("SurvSHAP(t) - Lognormal Mixture Cure Model") +
  theme_minimal()

