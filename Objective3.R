# FILENAME: FinalProjectObjective3.Rmd
# AUTHOR: Devon Rojas
# Summary: Performs forecasting analysis for Project Objective 3

##########################################
##### SET SEED AND IMPORT LIBRARIES ######
##########################################
set.seed(2)

library(tidyverse)
library(ggplot2)
library(glmnet)

##########################################
###### PREPROCESS DATA FOR ANALYSIS ######
##########################################
CleanedStudyData = read.csv("data/processed/cleaned_data_obj3.csv")
data = CleanedStudyData %>% mutate(across(PatientID:AssessmentID, as.factor))
# -------------------------------------
# Records Baseline (Day 0) data into a separate data frame for later use.
# -------------------------------------
Baseline = data %>% 
  filter(VisitDay == 0) %>% 
  select(PatientID, P1:G16) %>% 
  mutate(BaselinePANSS_Total = rowSums(across(P1:G16))) %>%
  rename_with(~ paste("Baseline", .x, sep=""), P1:G16)
# -------------------------------------
# Create additional variables in data for model
# -------------------------------------
StudyData = data  %>% 
  select(-c(X, SiteID, RaterID)) %>% 
  arrange(PatientID, VisitDay) %>%
  # -------------------------------------
# Record Lag and Lead variables for PANSS scores
# -------------------------------------
mutate(across(P1:PANSS_Total, ~ lag(.x, 1), .names="Lag{.col}"),
       across(P1:PANSS_Total, ~ lead(.x, 1), .names="Lead{.col}")) %>%
  fill(LagP1:LagPANSS_Total, .direction="up") %>%
  fill(LeadP1:LeadPANSS_Total, .direction="down") %>%
  # -------------------------------------
# Assign a TreatmentWeek variable to each patient's subsequent visit day
# (POTENTIAL ERROR: Duplicate VisitDays for a Patient getting incorrect
# TreatmentWeeks assigned to them.)
# -------------------------------------
group_by(PatientID) %>%
  arrange(VisitDay) %>%
  mutate(TreatmentWeek = row_number() - 1) %>%
  ungroup() %>%
  # -------------------------------------
# Join Baseline data into entire table and populate 
# through each PatientID
# -------------------------------------
left_join(Baseline) %>%
  arrange(PatientID, VisitDay) %>%
  fill(BaselineP1:BaselinePANSS_Total)

##########################################
####### DESIGN MODEL FOR ANALYSIS ########
##########################################
P = paste("P", 1:7, sep="")
N = paste("N", 1:7, sep="")
G = paste("G", 1:15, sep="")
a = c(P, N, G)
# -------------------------------------
# Dynamically generates all Lag, Lead, and Baseline variables
# -------------------------------------
LagVars = paste("Lag", a, sep="", collapse="+")
LeadVars = paste("Lead", a, sep="", collapse="+")
BaselineVars = paste("Baseline", a, sep="", collapse="+")
# -------------------------------------
# Build model from user-specified Predictors and dynamically 
# generated variables from above
# -------------------------------------
model = as.formula(paste("PANSS_Total ~ 
                         PatientID + 
                         Country + 
                         TreatmentWeek*TxGroup + 
                         BaselinePANSS_Total + 
                         LagPANSS_Total + 
                         LeadPANSS_Total", 
                         LagVars, 
                         LeadVars, 
                         BaselineVars, 
                         sep="+", collapse="+"))

# -------------------------------------
# Pull Study E data out since we will be predicting on it later
# -------------------------------------
# ModelData = StudyData %>% filter(Study != "E")
StudyEData = StudyData %>% filter(Study == "E")
ModelData = StudyEData # * Actually, lets look at what just using StudyE on model will do.

##########################################
######## LASSO REGRESSION METHOD #########
##########################################
x = model.matrix(model, ModelData)[, -1]
y = ModelData$PANSS_Total
grid = 10^seq(10, -2, length=100)

train = sample(1:nrow(x), 0.7*nrow(x))
test = (-train)

lasso.mod = glmnet(x[train,], y[train], alpha=1, lambda=grid)
# plot(lasso.mod)

cv.out = cv.glmnet(x[train, ], y[train], alpha=1)
plot(cv.out)
bestlam.lasso = cv.out$lambda.min
# bestlam.lasso

##########################################
######## RIDGE REGRESSION METHOD #########
##########################################
# ridge.mod = glmnet(x[train,], y[train], alpha=0, lambda=grid)
# plot(ridge.mod)
# 
# cv.out = cv.glmnet(x[train, ], y[train], alpha=0)
# plot(cv.out)
# bestlam.ridge = cv.out$lambda.min
# bestlam.ridge

##########################################
###### ANALYZE TEST RMSE FOR MODEL #######
##########################################
lasso.mod = glmnet(x[train, ], y[train], alpha=1, lambda=bestlam.lasso)
# ridge.mod = glmnet(x[train, ], y[train], alpha=0, lambda=bestlam)

preds.lasso = predict(lasso.mod, s=bestlam.lasso, newx=x[test,])
# preds.ridge = predict(ridge.mod, s=bestlam, newx=x[test,])

caret::RMSE(preds.lasso, y[test])
# RMSE(preds.ridge, y[test])

##########################################
### PERFORM REGRESSION ON ACTUAL DATA ####
##########################################
# -------------------------------------
# Calculate the last TreatmentWeek value for a Patient
# -------------------------------------
# StudyEData.MaxTreatmentWeek = StudyEData %>% 
#   group_by(PatientID) %>% 
#   summarise(maxTreatmentWeek = max(TreatmentWeek))

StudyEData.Wk18 = StudyEData %>% 
  group_by(PatientID) %>% 
  filter(abs(VisitDay - 126) == min(abs(VisitDay - 126))) %>% 
  distinct(PatientID, .keep_all=T) %>%
  mutate(VisitDay = 126)

# -------------------------------------
# Calculate next TreatmentWeek in sequence for Patient
# -------------------------------------
# StudyEData.Wk18 = StudyEData.Wk18 %>% 
#   left_join(StudyEData.MaxTreatmentWeek) %>% 
#   mutate(TreatmentWeek = maxTreatmentWeek + 1) %>% 
#   select(-maxTreatmentWeek)
# -------------------------------------
# Transform data for Lasso/Ridge Regression
# -------------------------------------
x_test = model.matrix(model, StudyEData.Wk18)[, -1]
# -------------------------------------
# Capture existing PANSS_Total scores for Patients 
# who have Week 18 data
# -------------------------------------
Week18Actual = StudyEData %>% 
  filter(Week == 18) %>% 
  arrange(PatientID) %>% 
  select(PatientID, PANSS_Total)
# -------------------------------------
# Create dataframe for model PANSS predictions
# -------------------------------------
# res = predict(ridge.mod, s=bestlam.ridge, newx=x_test)
res = predict(lasso.mod, s=bestlam.lasso, newx=x_test)
preds = tibble(
  PatientID = StudyEData.Wk18$PatientID, 
  PANSS_Total = res
)
# -------------------------------------
# Merge model predictions with existing PANSS_Total
# scores for Patients at Week 18 and export to CSV.
# -------------------------------------
Predictions = Week18Actual %>%
  full_join(preds, by="PatientID") %>%
  arrange(PatientID) %>%
  mutate(PANSS_Total = ifelse(is.na(PANSS_Total.x), PANSS_Total.y, PANSS_Total.x)) %>%
  select(-PANSS_Total.x, -PANSS_Total.y)
write.csv(Predictions, "output/PANSS_Predictions_Lasso_v3.csv", row.names = F)

##########################################
##########################################
##### UNUSED CODE AND OTHER ANALYSIS #####
##########################################
##########################################

# -------------------------------------
# As part of my development process, I considered using 
# weights in my model to influence how observations were
# treated. One of my approaches was to see how far away
# each patient observation was from 'Week 18' (which I 
# set as VisitDay 126), and then assign more weight to
# observations that were closer to that day. I thought
# this would add a sort of a KNN algorithm to the model.
# Unfortunately, after testing the resulting weights on
# various model forms, it appeared that it did not help
# to increase my model performance.
#
# **Note: This was done before I developed my latest model
# so it may be that the weights had no effect because my
# model was poor quality. **
# -------------------------------------
# StudyEData %>% mutate(weight = 1 -abs(VisitDay - 126)/126) %>%
#   ggplot(aes(Week, PANSS_Total)) +
#   geom_point(aes(group=PatientID, color=weight))

# -------------------------------------
# I also looked at the relationship between Treatment Week
# and PANSS_Total across different polynomial degrees to 
# get a better idea of which to use in my model design.
# -------------------------------------
# ggplot(StudyEData, aes(TreatmentWeek, PANSS_Total)) +
#   geom_point() +
#   geom_smooth(method="lm", formula=y~x, se=F, color="red") +
#   geom_smooth(method="lm", formula=y~I(x^2), se=F, color="blue") +
#   geom_smooth(method="lm", formula=y~I(x^3), se=F, color="green") +
#   geom_smooth(method="lm", formula=y~I(x^4), se=F, color="purple") +
#   facet_wrap(~TxGroup)

# -------------------------------------
# To further explore higher order relationships between
# the data, I plotted the Train and Test MSE of some
# models to see which order gave the best performance.
# -------------------------------------
# dfs = 10
# dfs = 5
# train_mse = rep(0, dfs)
# test_mse = rep(0, dfs)
# 
# model.1 = function(i) {
#   return (PANSS_Total ~ 
#             Country + 
#             TxGroup + 
#             (P1_0:G16_0) +
#             poly(TreatmentWeek, i, raw=T) + 
#             poly(TreatmentWeek, i, raw=T)*(PANSS_Total_0) +
#             poly(TreatmentWeek, i, raw=T)*(P1_0:G16_0) +
#             poly(TreatmentWeek, i, raw=T)*TxGroup +
#             TxGroup*PANSS_Total_0 +
#             TxGroup*(P1_0:G16_0))
# }
# model.2 = function(i) {
#   PANSS_Total ~
#   poly(TreatmentWeek, i, raw=T) +
#   poly(TreatmentWeek, i, raw=T)*(PANSS_Total_0) +
#   poly(TreatmentWeek, i, raw=T)*TxGroup +
#   TxGroup*PANSS_Total_0
# }
# 
# data = StudyEData.Wk1_17
# 
# for(i in 1:dfs) {
#   train = sample(1:nrow(data), nrow(data)/2)
#   test = (-train)
#   
#   model = lm(model.1(i),
#              data,
#              subset=train)
#   
#   pred_train = predict(model)
#   pred_test = predict(model, data[test,])
# 
#   train_mse[i] = RMSE(pred_train, data$PANSS_Total[train])
#   test_mse[i] = RMSE(pred_test, data$PANSS_Total[test])
# } 
# 
# ggplot(tibble(test_mse, 
#               train_mse, 
#               avg_mse = (test_mse + train_mse) /2,
#               boot_mse = boot_mse,
#               df = 1:dfs)) + 
#   geom_line(aes(df, test_mse), color="red") +
#   geom_line(aes(df, train_mse), color="blue") +
#   geom_line(aes(df, avg_mse), color="green") +
#   scale_x_continuous(breaks =c(1:dfs)) +
#   ylab("MSE")