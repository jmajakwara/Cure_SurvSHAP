rm(list=ls())

library(ggplot2)
library(MASS)
library(caret)
library(reticulate)
library(data.table)
library(survival)
library(rsample)
library(Hmisc)
library(BART)
library(rsample)
library(SurvMetrics)
library(prodlim)
library(pec)
library(rms)
library(flexsurv)
library(flexsurvcure)
library(tidyverse)
library(censored)
library(tidymodels)
library(survex)
library(survminer)

##### Housekeeping ######
options(expressions=10000)
R.version.string
date()
Sys.time()
options(digits=6)
options(scipen=6)


set.seed(2025)
loans <- read.csv("loans.csv",header=TRUE,sep=",",na.strings = "NA")

loans$time <- as.numeric(loans$time)
loans$OpenRevolvingMonthlyPayment <- as.numeric(loans$OpenRevolvingMonthlyPayment)
loans$EmploymentStatusDuration <- as.numeric(loans$EmploymentStatusDuration)
loans$AmountDelinquent <- as.numeric(loans$AmountDelinquent)
loans$OpenRevolvingMonthlyPayment <- as.numeric(loans$OpenRevolvingMonthlyPayment)
loans$RevolvingCreditBalance <- as.numeric(loans$RevolvingCreditBalance)
loans$AvailableBankcardCredit <- as.numeric(loans$AvailableBankcardCredit)
loans$LoanOriginalAmount <- as.numeric(loans$LoanOriginalAmount)
loans$ListingCategory <- as.factor(loans$ListingCategory)
loans$EmploymentStatus <- as.factor(loans$EmploymentStatus)
loans$CreditScoreRange <- as.factor(loans$CreditScoreRange)
loans$Occupation <- as.factor(loans$Occupation)
#loans <- loans  %>% mutate_if(is.character, as.factor)
loans$CurrentlyInGroup <- as.factor(loans$CurrentlyInGroup)
loans$ProsperScore <- ifelse(loans$ProsperScore == 11,10,loans$ProsperScore)
loans$ProsperScore <- as.factor(loans$ProsperScore)
loans$Term <- as.factor(loans$Term)

# convert days to months


loans <- loans %>% dplyr::filter(DebtToIncomeRatio < 1 | is.na(DebtToIncomeRatio))
loans <- loans %>% dplyr::select(-2,-Occupation,-ProsperScore)

loans_cc <- loans  %>% na.omit()

loans_cc <- loans_cc %>% dplyr::select(time, status, OpenCreditLines, InquiriesLast6Months,   BankcardUtilization, DebtToIncomeRatio,  BorrowerRate,
				LoanCurrentDaysDelinquent, CurrentlyInGroup, CreditScoreRange, Term, TradesNeverDelinquent,IsBorrowerHomeowner)
loans <- loans %>% dplyr::select(-ProsperRating,-LenderYield,-CurrentCreditLines)

data_split <- initial_split(loans_cc, prop = 0.7, strata = status)     
#Create data frames for the two sets:
loans_train <- training(data_split)
loans_test  <- testing(data_split)

pdf(file="KM_Plot.pdf", width = 10, height = 7)
fit <- survfit(Surv(time, status) ~ CreditScoreRange, data = loans_cc)
ggsurvplot(fit, data = loans_cc, 
           pval = TRUE, 
           pval.method = TRUE, 
           pval.coord = c(0, 0.15),
           pval.method.coord = c(0, 0.3), 
           conf.int = FALSE, 
           risk.table = TRUE, 
           linetype = "strata",
           title = "Kaplan-Meier Survival Curves stratified by Credit Score",
           xlab = "Time (Days)", 
           ylab = "Survival Probability",
           legend = "right",  # Position legend on the right side
           ggtheme = theme_minimal())  # Apply minimal theme
dev.off()

X_train <- loans_train %>% select(OpenCreditLines, InquiriesLast6Months,   BankcardUtilization, DebtToIncomeRatio,  BorrowerRate,
				 CurrentlyInGroup, CreditScoreRange, Term, IsBorrowerHomeowner)
X_test  <- loans_test  %>% select(OpenCreditLines, InquiriesLast6Months,   BankcardUtilization, DebtToIncomeRatio,  BorrowerRate,
				 CurrentlyInGroup, CreditScoreRange, Term, IsBorrowerHomeowner)



loans_train_all <- loans[!(loans$LoanKey %in% c(loans_test$LoanKey)),]

loans_train <- as_tibble(loans_train)
loans_test <- as.data.frame(loans_test)


###### Exponential Modeling ########

cat('The Exp Model using selected variables output \n')

loans_cc <- loans_cc |>
  mutate(event_time = Surv(time, status))


split   <- initial_split(loans_cc, prop = 0.75, strata = status)
loans_tr  <- training(split)
loans_val <- testing(split)

exp_spec <- 
  survival_reg(dist = "exp") |>
  set_engine("flexsurv") |>
  set_mode("censored regression") 
exp_spec


exp_fit <- exp_spec |> fit(event_time ~ BorrowerRate + OpenCreditLines + DebtToIncomeRatio + InquiriesLast6Months + BankcardUtilization + TradesNeverDelinquent +
                                           	Term,
						data = loans_tr)
exp_fit1 <- exp_spec |> fit(Surv(time,status) ~ BorrowerRate + OpenCreditLines + DebtToIncomeRatio + InquiriesLast6Months + BankcardUtilization + TradesNeverDelinquent +
                                           	Term,
						data = loans_cc)


exp_fit
tidy(exp_fit)
cat("Results from training set \n")
glance(exp_fit)
cat("Results from CC set \n")
glance(exp_fit1)


time_points <- sort(unique(loans_val$time))
exp_pred <- augment(exp_fit, loans_val, eval_time = time_points)
brier_scores <-
  exp_pred |> 
  brier_survival(truth = event_time, .pred)

quantile(brier_scores$.estimate)


IBS_exp <- exp_pred |> brier_survival_integrated(truth = event_time, .pred)

IBS_exp

roc_scores <-
  exp_pred |> 
  roc_auc_survival(truth = event_time, .pred)



exp_cindex <- exp_pred |>  concordance_survival(truth = event_time, estimate = .pred_time)
exp_cindex

### Cross-validation lognormal
loans_folds <- vfold_cv(loans_cc, v = 5, strata = status)

loans_folds$splits[[1]] |> analysis() |> dim()

cat("Cross validation \n ")
exp_cross_val <- exp_spec |> fit_resamples(event_time ~ BorrowerRate + OpenCreditLines + DebtToIncomeRatio + InquiriesLast6Months + BankcardUtilization + TradesNeverDelinquent +
                                           	Term,
						resamples = loans_folds, 
						eval_time = sort(unique(loans_cc$time)), 
						metrics = metric_set(brier_survival,concordance_survival, brier_survival_integrated))
print(show_notes(.Last.tune.result))
exp_metric = collect_metrics(exp_cross_val)


cat("The IBS is \n ")
collect_metrics(exp_cross_val) |> 
  filter(.metric == "brier_survival_integrated")

cat("The c-index is \n ")
collect_metrics(exp_cross_val) |> 
  filter(.metric == "concordance_survival")


####### Loglogistic Modelling ###########
logl_spec <- 
  survival_reg(dist = "llogis") |>
  set_engine("flexsurv") |>
  set_mode("censored regression") 
logl_spec


logl_fit <- logl_spec |> fit(event_time ~ OpenCreditLines + InquiriesLast6Months + Term +
                                        BankcardUtilization + DebtToIncomeRatio + CurrentlyInGroup + CreditScoreRange, 
					data = loans_tr)

logl_fit1 <- logl_spec |> fit(event_time ~ OpenCreditLines + InquiriesLast6Months + Term +
                                        BankcardUtilization + DebtToIncomeRatio + CurrentlyInGroup + CreditScoreRange, 
					data = loans_cc)


logl_fit
tidy(logl_fit)
cat("Results from training set \n")
glance(logl_fit)
cat("Results from CC set \n")
glance(logl_fit1)


time_points <- sort(unique(loans_val$time))
logl_pred <- augment(logl_fit, loans_val, eval_time = time_points)
brier_scores <-
  logl_pred |> 
  brier_survival(truth = event_time, .pred)

quantile(brier_scores$.estimate)


IBS_logl <- logl_pred |> brier_survival_integrated(truth = event_time, .pred)

IBS_logl

roc_scores <-
  logl_pred |> 
  roc_auc_survival(truth = event_time, .pred)



logl_cindex <- logl_pred |>  concordance_survival(truth = event_time, estimate = .pred_time)
logl_cindex

### Cross-validation lognormal
loans_folds <- vfold_cv(loans_cc, v = 5, strata = status)


cat("Cross validation \n ")
logl_cross_val <- logl_spec |> fit_resamples(event_time ~ OpenCreditLines + InquiriesLast6Months + Term +
                                        		BankcardUtilization + DebtToIncomeRatio + CurrentlyInGroup + CreditScoreRange,
							# method = "Nelder-Mead", 
							resamples = loans_folds, 
							eval_time = sort(unique(loans_cc$time)), 
							metrics = metric_set(concordance_survival,brier_survival,brier_survival_integrated))

logl_metric = collect_metrics(logl_cross_val)


cat("The IBS is = ")
collect_metrics(logl_cross_val) |> 
  filter(.metric == "brier_survival_integrated")
cat("The c-index is = ")
collect_metrics(logl_cross_val) |> 
  filter(.metric == "concordance_survival")


######### Weibull Modelling ########################
wei_spec <- 
  survival_reg(dist = "weibull") |>
  set_engine("flexsurv") |>
  set_mode("censored regression") 
wei_spec


wei_fit <- wei_spec |> fit(event_time ~  OpenCreditLines + InquiriesLast6Months + BankcardUtilization + DebtToIncomeRatio + BorrowerRate  +
                                         CreditScoreRange + Term, 
					data = loans_tr)

wei_fit1 <- wei_spec |> fit(event_time ~  OpenCreditLines + InquiriesLast6Months + BankcardUtilization + DebtToIncomeRatio + BorrowerRate  +
                                         CreditScoreRange + Term, 
					data = loans_cc)


wei_fit
tidy(wei_fit)
cat("Results from training set \n")
tidy(wei_fit)
glance(wei_fit)
cat("Results from CC set \n")
tidy(wei_fit1)
glance(wei_fit1)


time_points <- sort(unique(loans_val$time))
wei_pred <- parsnip::augment(wei_fit, loans_val, eval_time = time_points)

brier_scores <-
  wei_pred |> 
  brier_survival(truth = event_time, .pred)

quantile(brier_scores$.estimate)


IBS_wei <- wei_pred |> brier_survival_integrated(truth = event_time, .pred)


roc_scores <-
  wei_pred |> 
  roc_auc_survival(truth = event_time, .pred)



wei_cindex <- wei_pred |>  concordance_survival(truth = event_time, estimate = .pred_time)
wei_cindex

### Cross-validation lognormal
loans_folds <- vfold_cv(loans_cc, v = 5, strata = status)

loans_folds$splits[[1]] |> analysis() |> dim()

cat("Cross validation \n ")
wei_cross_val <- wei_spec |> fit_resamples(event_time ~ OpenCreditLines + InquiriesLast6Months + BankcardUtilization + DebtToIncomeRatio + BorrowerRate  +
                                         	CreditScoreRange + Term,
						resamples = loans_folds, 
						eval_time = sort(unique(loans_cc$time)), 
						metrics = metric_set(concordance_survival, brier_survival, brier_survival_integrated))

wei_metric = collect_metrics(wei_cross_val)

cat("The IBS is  \n")
collect_metrics(wei_cross_val) |> 
  filter(.metric == "brier_survival_integrated")
cat("The c-index is  \n")
collect_metrics(wei_cross_val) |> 
  filter(.metric == "concordance_survival")



###### Lognormal Modeling ########
log_spec <- 
  survival_reg(dist = "lnorm") |>
  set_engine("flexsurv") |>
  set_mode("censored regression") 
log_spec


log_fit <- log_spec |> fit(event_time ~ OpenCreditLines + InquiriesLast6Months + BorrowerRate +  BankcardUtilization + 
                                          	DebtToIncomeRatio + IsBorrowerHomeowner + CreditScoreRange + Term, 
						data = loans_tr)

log_fit1 <- log_spec |> fit(event_time ~ OpenCreditLines + InquiriesLast6Months + BorrowerRate + BankcardUtilization + 
                                          	DebtToIncomeRatio +  IsBorrowerHomeowner + CreditScoreRange + Term, 
						data = loans_cc)


log_fit
tidy(log_fit)
cat("Results from training set \n")
glance(log_fit)
cat("Results from CC set \n")
tidy(log_fit1)
glance(log_fit1)


time_points <- sort(unique(loans_val$time))
log_pred <- augment(log_fit, loans_val, eval_time = time_points)
brier_scores <-
  log_pred |> 
  brier_survival(truth = event_time, .pred)

quantile(brier_scores$.estimate)


IBS_log <- log_pred |> brier_survival_integrated(truth = event_time, .pred)

IBS_log

roc_scores <-
  log_pred |> 
  roc_auc_survival(truth = event_time, .pred)



log_cindex <- log_pred |>  concordance_survival(truth = event_time, estimate = .pred_time)
log_cindex

### Cross-validation lognormal
loans_folds <- vfold_cv(loans_cc, v = 5, strata = status)

loans_folds$splits[[1]] |> analysis() |> dim()

cat("Cross validation \n ")
log_cross_val <- log_spec |> fit_resamples(event_time ~ OpenCreditLines + InquiriesLast6Months + BorrowerRate +  BankcardUtilization + 
                                          	DebtToIncomeRatio +  IsBorrowerHomeowner + CreditScoreRange + Term,
						resamples = loans_folds, 
						eval_time = sort(unique(loans_cc$time)), 
						metrics = metric_set(concordance_survival,brier_survival,brier_survival_integrated))

log_metric = collect_metrics(log_cross_val)

cat("The IBS is = ")
collect_metrics(log_cross_val) |> 
  filter(.metric == "brier_survival_integrated")
cat("The c-index is = ")
collect_metrics(log_cross_val) |> 
  filter(.metric == "concordance_survival")


cat('The Lognormal cure Model with CurrentlyInGroup  \n')


log_fit <- flexsurvcure::flexsurvcure(Surv(time,status) ~ CreditScoreRange, anc = list(meanlog = ~  DebtToIncomeRatio +  BorrowerRate + OpenCreditLines + 
							InquiriesLast6Months +  IsBorrowerHomeowner + Term + BankcardUtilization),
                                   			data = loans_cc,
                                   			dist = 'lnorm',
                                   			link = 'logistic',  
                                   			mixture = TRUE)

tidy(log_fit)
glance(log_fit)
glance(log_fit)[["logLik"]]



time_points <- seq(0,max(loans_test$time),length.out = 100) 



# --------------------- PREDICT FUNCTION ---------------------
predict_survival <- function(object, newdata, times = NULL) {
  if (is.null(times)) times <- time_points
  n <- nrow(newdata)
  nt <- length(times)
  s <- summary(object, newdata = newdata, t = times, type = "survival", tidy = TRUE, ci = FALSE)
  s$.id <- rep(seq_len(n), each = nt)
  mat <- s |> dplyr::select(.id,time, est) |> pivot_wider(names_from = time, values_from = est)
  as.matrix(mat[,-1])  # remove .id column
}


# --------------------- EXPLAINER ---------------------
explainer <- survex::explain(
  model = log_fit,
  data = X_train,
  y = Surv(loans_train$time, loans_train$status),
  predict_function = predict_survival,
  times = time_points,
  label = "Lognormal Cure Model"
)


# --------------------- SurvSHAP(t) ---------------------
cat("Computing SurvSHAP(t) on test observations...\n")

shap_result <- predict_parts(
  explainer,
  new_observation = X_test,
  type = "survshap",        
  N = 100,
  output_type = "survival",
  calculation_method = "kernelshap",
  aggregation_method = "integral"
)


# --------------------- PLOTS ---------------------
pdf("SHAP_beeswarm.pdf", width = 11, height = 7)
plot(shap_result, geom = "beeswarm",max_vars = 8) + 
  ggtitle("SurvSHAP(t) - Lognormal Mixture Cure Model") +
  theme_minimal() + 
  theme(legend.position = "bottom")
dev.off()

pdf("SHAP_importance.pdf", width = 10, height = 6)
plot(shap_result, type = "importance",max_vars = 8) + 
  ggtitle("Global Feature Importance over Time") +
  theme_minimal() + 
  theme(legend.position = "bottom")
dev.off()
