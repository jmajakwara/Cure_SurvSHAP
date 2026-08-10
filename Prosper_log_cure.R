rm(list=ls())

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

###### Lognormal Modeling ########


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


