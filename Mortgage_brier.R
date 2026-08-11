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

##### Housekeeping ######
R.version.string
date()
Sys.time()




# ===============================
# DATA PREPARATION 
# ===============================
mortgage_data  <- read.csv("mortgage.csv")


mortgage <- mortgage_data %>%
  dplyr::group_by(id) %>%
  dplyr::summarise(
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
    status = case_when(
      event_default == 1 ~ 1,
      TRUE ~ 0 
    )
  ) %>%

  dplyr::select(id, time_to_event, status, 
                FICO, LTV, GDP, UER, HPI, Interest_rate, balance, 
                investor, REtype_SF, REtype_CO, REtype_PU)

mortgage <- mortgage |> dplyr::filter(time_to_event > 0) |>
  dplyr::mutate(
    Credit_score = dplyr::case_when(
      FICO >= 400 & FICO <= 659 ~ "Fair",
      FICO >= 660 & FICO <= 739 ~ "Good",
      FICO >= 740 & FICO <= 799 ~ "Very good",
      FICO >= 800 & FICO <= 899 ~ "Excellent",
      TRUE ~ NA_character_  # Values outside these ranges become NA
    )
  )

mortgage <- mortgage %>%
dplyr::mutate(
    dplyr::across(c(FICO,LTV,UER,HPI,GDP, Interest_rate, balance), as.numeric),
    dplyr::across(c(Credit_score, REtype_SF, REtype_CO, REtype_PU, investor), as.factor)
  )
summary(mortgage)

cat("\n--- Survival Data Summary ---\n")
cat("Number of loans:", nrow(mortgage), "\n")
cat("Default rate:", mean(mortgage$status), "\n")
cat("Censoring rate:", 1 - mean(mortgage$status), "\n")

head(mortgage,n=10)



# ===============================
# PROPER DATA SPLITTING FIRST
# ===============================

set.seed(2000)
mortgage <- mortgage |>
  dplyr::mutate(event_time = Surv(time_to_event, status))

data_split <- initial_split(mortgage, prop = 0.7, strata = status) 
mortgage_train <- training(data_split)
mortgage_test  <- testing(data_split)



###### Cox Modeling ########


cat('The Exp Model on all variables output \n')

time_points <- sort(unique(mortgage_test$time_to_event))

exp_spec <- 
  survival_reg(dist = "exp") |>
  set_engine("flexsurv") |>
  set_mode("censored regression") 
exp_spec


exp_fit <- exp_spec |> fit(event_time ~ Credit_score + HPI + LTV + GDP + UER,data = mortgage_train)

exp_pred <- augment(exp_fit, mortgage_test, eval_time = time_points)
brier_exp <-
  exp_pred |> 
  brier_survival(truth = event_time, .pred)

quantile(brier_exp$.estimate)



log_spec <- 
  survival_reg(dist = "lnorm") |>
  set_engine("flexsurv") |>
  set_mode("censored regression") 
log_spec


log_fit <- log_spec |> fit(event_time ~ Credit_score + HPI + LTV + GDP + UER,data = mortgage_train)


log_pred <- augment(log_fit, mortgage_test, eval_time = time_points)
brier_log <-
  log_pred |> 
  brier_survival(truth = event_time, .pred)

quantile(brier_log$.estimate)


logl_spec <- 
  survival_reg(dist = "llogis") |>
  set_engine("flexsurv") |>
  set_mode("censored regression") 
logl_spec


logl_fit <- logl_spec |> fit(event_time ~ Credit_score + HPI + LTV + GDP + UER,data = mortgage_train)

logl_pred <- augment(logl_fit, mortgage_test, eval_time = time_points)
brier_logl <-
  logl_pred |> 
  brier_survival(truth = event_time, .pred)

quantile(brier_logl$.estimate)

wei_spec <- 
  survival_reg(dist = "weibull") |>
  set_engine("flexsurv") |>
  set_mode("censored regression") 
wei_spec


wei_fit <- wei_spec %>% fit(event_time ~ Credit_score + HPI + LTV + GDP + UER,data = mortgage_train)

wei_pred <- augment(wei_fit, mortgage_test, eval_time = time_points)

brier_wei <-
  wei_pred |> 
  brier_survival(truth = event_time, .pred)

quantile(brier_wei$.estimate)


Model = c("Exponential",  "Lognormal", "Loglogistic", "Weibull")

# Create a data frame to store Brier scores and time points
brier_data <- data.frame(
  Time = time_points,
  Exponential = brier_exp$.estimate,
  Lognormal = brier_log$.estimate,
  Loglogistic = brier_logl$.estimate,
  Weibull = brier_wei$.estimate
)

# Pivot_longer the data frame to long format
brier_data_long <- pivot_longer(brier_data, cols = -Time, names_to = "Model", values_to = "Brier_Score")

# Create the plot using ggplot2
pdf(file="Brier score for distributions.pdf")
ggplot(brier_data_long, aes(x = Time, y = Brier_Score, color = Model)) +
  geom_line() +
  geom_hline(yintercept = 1 / 4, col = "red", lty = 2) +
  labs(x = "Time (in months)",
       y = "Brier score") +
  theme_minimal() +
  theme(legend.position = "right",
        axis.text = element_text(size = 12),  
        axis.title = element_text(size = 14)) 
dev.off()

