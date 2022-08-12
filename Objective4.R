# FILENAME: FinalProjectObjective4.Rmd
# AUTHOR: Devon Rojas
# Summary: Performs classification analysis for Project Objective 4

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
CleanedStudyData = read.csv("data/processed/cleaned_data_obj4.csv")
data = CleanedStudyData %>%mutate(across(PatientID:AssessmentID, as.factor))
# -------------------------------------
# Create additional variables in data for model
# -------------------------------------
StudyData = data  %>% 
  select(-c(X, SiteID, RaterID)) %>% 
  arrange(PatientID, VisitDay) %>%
  # -------------------------------------
# Build Lag, Lead, Change, and Flag variables
# for PANSS scores
# -------------------------------------
mutate(
  # -------------------------------------
  # Lag Variables
  # -------------------------------------
  across(P1:G16, ~ lag(.x, 1), .names="Lag{.col}"),
  across(LagP1:LagG16, ~replace_na(.x, 0)),
  # -------------------------------------
  # Lead Variables
  # -------------------------------------
  across(P1:G16, ~ lead(.x, 1), .names="Lead{.col}"),
  across(LeadP1:LeadG16, ~replace_na(.x, 0)),
  # -------------------------------------
  # Change Variables (current - lag)
  # -------------------------------------
  across(P1: G16, ~ .x - lag(.x, 1), .names="Change{.col}"),
  across(ChangeP1:ChangeG16, ~replace_na(.x, 0)),
  # -------------------------------------
  # Flag Variables (change > 2)
  # -------------------------------------
  across(ChangeP1:ChangeG16,
         ~ case_when(abs(.x) > 2 ~ 1, .x <=2 ~ 0), 
         .names="Flag{.col}")
) %>%
  # -------------------------------------
# Build Lag and Lead Variables for Composite Score
# -------------------------------------
mutate(
  CompositeScore = rowSums(across(N1:N7)) - rowSums(across(P1:P7)),
  LagCompositeScore = lag(CompositeScore, 1),
  LeadCompositeScore = lead(CompositeScore, 1)
) %>%
  replace_na(list(
    LagCompositeScore = 0, 
    LeadCompositeScore = 0)
  ) %>%
  # -------------------------------------
# Build Flags for various conditions
# -------------------------------------
mutate(
  # -------------------------------------
  # Flag if P score, N score, or Composite Score 
  # ranges pass a certain threshold.
  # -------------------------------------
  FlagInconsistentP = max(c_across(P1:P7)) - min(c_across(P1:P7)) > 2,
  FlagInconsistentN = max(c_across(N1:N7)) - min(c_across(N1:N7)) > 2,
  FlagInconsistentComposite = CompositeScore - LagCompositeScore > 5,
  # -------------------------------------
  # Flag if certain PANSS items are measured above
  # a certain value, indicating remission.
  # 
  # Reference:
  # Santor, D.A., Ascher-Svanum, H., Lindenmayer, JP. et al. 
  # Item response analysis of the Positive and Negative Syndrome Scale. 
  # BMC Psychiatry 7, 66 (2007). https://doi.org/10.1186/1471-244X-7-66
  # -------------------------------------
  FlagDelusion = P1 >= 3,
  FlagHallucination = P3 >=3,
  FlagGrandiosity = P5 >=3, 
  FlagMannerisms = G5 >= 2,
  # -------------------------------------
  # Flag if Composite score for current observation is
  # spike in patient's study duration.
  # -------------------------------------
  FlagUnusualOutcome = 
    CompositeScore > LagCompositeScore & 
    CompositeScore < LeadCompositeScore,
  # -------------------------------------
  # Flag if more than 3 PANSS items have a large
  # change from previous observation values.
  # -------------------------------------
  FlagInconsistentMeasure = rowSums(across(FlagChangeP1:FlagChangeG16)) > 3
) %>%
  # -------------------------------------
# Reclassify LeadStatus into patient's who passed
# assessment and those who did not
# -------------------------------------
mutate(LeadStatus = as.factor(ifelse(LeadStatus=="Passed", "P", "NP"))) %>%
  # -------------------------------------
# Assign a TreatmentWeek variable to each patient's subsequent visit day
# (POTENTIAL ERROR: Duplicate VisitDays for a Patient getting incorrect
# TreatmentWeeks assigned to them.)
# -------------------------------------
group_by(PatientID) %>%
  arrange(VisitDay) %>%
  mutate(TreatmentWeek = row_number() -1) %>%
  ungroup() %>%
  select(-c(PatientID))

##########################################
####### DESIGN MODEL FOR ANALYSIS ########
##########################################
P = paste("P", 1:7, sep="")
N = paste("N", 1:7, sep="")
a = c(P,N)
# -------------------------------------
# Dynamically generates all Lag, Lead, and PANSS variables
# as well as interactions between Lag and Lead variables
# -------------------------------------
PANSSVars = paste(a, collapse="+")
LagVars = paste("Lag", a, sep="")
LeadVars = paste("Lead", a, sep="")
Interactions = paste(LagVars, LeadVars, sep="*", collapse = "+")
LagVars = paste(LagVars, collapse="+")
LeadVars = paste(LeadVars, collapse="+")
# -------------------------------------
# Dynamically generates Change and FlagChange variables
# -------------------------------------
ChangeVars = paste("Change", a, sep="", collapse="+")
PANSSFlags = paste("FlagChange", a, sep="", collapse="+")
# -------------------------------------
# List of other Flags to include in model
# -------------------------------------
Flags = paste("FlagInconsistentComposite", 
              "FlagDelusion", 
              "FlagHallucination", 
              "FlagGrandiosity", 
              "FlagMannerisms", 
              "FlagUnusualOutcome",
              "FlagInconsistentMeasure", 
              sep="+")
# -------------------------------------
# Build model from user-specified Predictors and dynamically 
# generated variables from above
# -------------------------------------
model = as.formula(paste("LeadStatus ~ 
                         TxGroup*TreatmentWeek+
                         P1*P3*P5*G5*LagP1*LagP3*LagP5*LagG5", 
                         Flags, 
                         PANSSFlags, 
                         sep="+"))

# -------------------------------------
# Pull Study E data out since we will be predicting on it later
# -------------------------------------
StudyE = StudyData %>% filter(Study == "E")
ModelData = StudyData %>% filter(Study != "E")

test = sample(1:nrow(ModelData), nrow(ModelData)/2)
train = (-test)

##########################################
####### LOGISTIC REGRESSION METHOD #######
##########################################
fit = glm(model, ModelData, subset=train, family="binomial")
# summary(fit)

train_probs = predict(fit, type="response")
train_preds = rep("NP", nrow(ModelData[train,]))
train_preds[train_probs>0.8] = "P"
train_response = ModelData[train,"LeadStatus"]
TRAINMSE = mean(train_preds == train_response)

test_probs = predict(fit, ModelData[test,], type="response")
test_preds = rep("NP", nrow(ModelData[test,]))
test_preds[test_probs>0.8] = "P"
test_response = ModelData[test, "LeadStatus"]
TESTMSE = mean(test_preds == test_response)

# Generate MSE for coin flip
unif_probs = runif(length(train_response), 0, 1)
unif_preds = rep("NP", length(unif_probs))
unif_preds[unif_probs>0.8] = "P"
RANDOMMSE = mean(unif_preds == train_response)

print(paste(TRAINMSE, TESTMSE, RANDOMMSE, sep=" "))

fit = glm(model, ModelData, family="binomial")
preds = predict(fit, StudyE, type="response")
StudyEPredictions = tibble(
  AssessmentID = StudyE$AssessmentID, 
  LeadStatus = preds
)

write.csv(StudyEPredictions, "output/LeadStatusPredictionsLogRegression_REVISED_v2.csv", row.names=F)

##########################################
######## LASSO REGRESSION METHOD #########
##########################################
# -------------------------------------
# Transform data for Lasso/Ridge Regression
# -------------------------------------
# x = model.matrix(model, ModelData)[, -1]
# y = ModelData$LeadStatus
# grid = 10^seq(10, -2, length=100)
# 
# lasso.mod = glmnet(x[train,], y[train], alpha=1, lambda=grid, family="binomial")
# plot(lasso.mod)
# cv.out = cv.glmnet(x[train,], y[train], alpha=1, family="binomial")
# plot(cv.out)
# bestlam = cv.out$lambda.min
# bestlam
# 
# lasso.mod = glmnet(x[train,], y[train], alpha=1, family="binomial", lambda=bestlam)
# preds = predict(lasso.mod, x[test,], type="class", se=bestlam)
# table(preds, y[test])
# mean(preds == y[test])
# 
# x = model.matrix(model, StudyE)[, -1]
# 
# Predictions = predict(lasso.mod, StudyE, type="response", se=bestlam)
# 
# LassoPredictions = tibble(AssessmentID = StudyE$AssessmentID, LeadStatus = Predictions)
# write.csv(LassoPredictions, "output/LeadStatusPredctions_Lasso.csv")

##########################################
######## RIDGE REGRESSION METHOD #########
##########################################
# -------------------------------------
# Transform data for Lasso/Ridge Regression
# -------------------------------------
# x = model.matrix(model, ModelData)[, -1]
# y = ModelData$LeadStatus
# grid = 10^seq(10, -2, length=100)
# 
# ridge.mod = glmnet(x[train,], y[train], alpha=0, lambda=grid, family="binomial")
# plot(ridge.mod)
# cv.out = cv.glmnet(x[train,], y[train], alpha=0, family="binomial")
# plot(cv.out)
# bestlam = cv.out$lambda.min
# bestlam
# 
# ridge.mod = glmnet(x[train,], y[train], alpha=1, family="binomial", lambda=bestlam)
# preds = predict(ridge.mod, x[test,], type="class")
# table(preds, y[test])
# mean(preds == y[test])
# 
# x = model.matrix(model, StudyE)[, -1]
# 
# Predictions = predict(ridge.mod, StudyE, type="response", se=bestlam)
# 
# RidgePredictions = tibble(AssessmentID = StudyE$AssessmentID, LeadStatus = Predictions)
# write.csv(RidgePredictions, "output/LeadStatusPredctions_Ridge.csv")

##########################################
##########################################
##### UNUSED CODE AND OTHER ANALYSIS #####
##########################################
##########################################
# -------------------------------------
# Look at distribution of LeadStatus across studies
# -------------------------------------
# StudyData %>%
#   filter(Study != "E") %>% 
#   group_by(Study, LeadStatus) %>% 
#   summarise(n = n()) %>% 
#   ungroup(LeadStatus) %>% 
#   mutate(pct = n/sum(n)) %>%
#   ggplot(aes(LeadStatus, pct, fill=LeadStatus)) +
#     geom_bar(stat="identity", position="dodge") +
#     facet_wrap(~Study, scales="free_x")

# -------------------------------------
# Look at PANSS score across studies by TxWeek
# -------------------------------------
# StudyData %>% 
#   group_by(Study, TreatmentWeek, TxGroup) %>% 
#   summarise(PANSS_Total = mean(PANSS_Total), .groups="keep") %>% 
#   ungroup() %>%
#   ggplot(aes(TreatmentWeek, PANSS_Total, color=TxGroup)) +
#     geom_line() + facet_wrap(~Study, scales="free")