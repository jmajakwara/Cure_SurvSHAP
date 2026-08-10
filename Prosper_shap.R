rm(list=ls())
options(encoding = "UTF-8")
#Sys.setlocale("LC_CTYPE", "C")

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
				LoanCurrentDaysDelinquent, IsBorrowerHomeowner, CreditScoreRange, Term, TradesNeverDelinquent)
loans <- loans %>% dplyr::select(-ProsperRating,-LenderYield,-CurrentCreditLines)

data_split <- initial_split(loans_cc, prop = 0.7, strata = status)     
#Create data frames for the two sets:
loans_train <- training(data_split)
loans_test  <- testing(data_split)

X_train <- loans_train %>% select(OpenCreditLines, InquiriesLast6Months,   BankcardUtilization, DebtToIncomeRatio,  BorrowerRate,
				 IsBorrowerHomeowner, CreditScoreRange, Term)
X_test  <- loans_test  %>% select(OpenCreditLines, InquiriesLast6Months,   BankcardUtilization, DebtToIncomeRatio,  BorrowerRate,
				 IsBorrowerHomeowner, CreditScoreRange, Term)



loans_train_all <- loans[!(loans$LoanKey %in% c(loans_test$LoanKey)),]

loans_train <- as_tibble(loans_train)
loans_test <- as.data.frame(loans_test)

###### Lognormal Modeling ########


cat('The Lognormal cure Model using flexsurvcure package  \n')

(log_fit <- flexsurvcure::flexsurvcure(Surv(time,status) ~ CreditScoreRange, anc = list(meanlog = ~ IsBorrowerHomeowner + DebtToIncomeRatio + 
							OpenCreditLines + InquiriesLast6Months + BorrowerRate + Term + BankcardUtilization),
                                   			data = loans_train,
                                   			dist = 'lnorm',
                                   			link = 'logistic',  
                                   			mixture = TRUE))
tidy(log_fit)
glance(log_fit)
glance(log_fit)[["logLik"]]


time_points <- seq(0,max(loans_test$time),length.out = 100) #time_points <- sort(unique(loans_test$time))

# --------------------- PREDICT FUNCTION ---------------------
predict_survival <- function(object, newdata, times = NULL) {
  if (is.null(times)) times <- time_points
  n <- nrow(newdata)
  nt <- length(times)
  s <- summary(object, newdata = newdata, t = times, type = "survival", tidy = TRUE, ci = FALSE)
  s$.id <- rep(seq_len(n), each = nt)
  mat <- s %>% dplyr::select(.id,time, est) %>% pivot_wider(names_from = time, values_from = est)
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

library(future)
future::plan(multisession, workers = 16)
# --------------------- SurvSHAP(t) ---------------------
cat("Computing SurvSHAP(t) on test observations...\n")
shap_result <- predict_parts(
  explainer,
  new_observation = X_test[1:1000,],
  type = "survshap",
  N = 50,
  calculation_method = "kernelshap"
)

# --------------------- PLOTS ---------------------
pdf("SHAP_beeswarm.pdf", width = 11, height = 7)
plot(shap_result, geom = "beeswarm",max_vars = 8) + 
  ggtitle("SurvSHAP(t) - Lognormal Mixture Cure Model") +
  theme_minimal() + 
  theme(legend.position = "bottom")
dev.off()

pdf("SHAP_importance.pdf", width = 10, height = 6)
plot(shap_result, geom = "importance",max_vars = 8) + 
  ggtitle("Global Feature Importance over Time") +
  theme_minimal() + 
  theme(legend.position = "bottom")

dev.off()




cat("Done! All plots saved.\n")
