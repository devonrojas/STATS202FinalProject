# FILENAME: FinalProjectObjective1.Rmd
# AUTHOR: Devon Rojas
# Summary: Analyzes study data for treatment effect in patients

##########################################
##### SET SEED AND IMPORT LIBRARIES ######
##########################################
set.seed(2)

library(tidyverse)
library(ggplot2)

##########################################
###### PREPROCESS DATA FOR ANALYSIS ######
##########################################
CleanedStudyData = read.csv("data/processed/cleaned_data_obj3.csv")
CleanedStudyData = CleanedStudyData %>% 
  select(-c(X, LeadStatus)) %>%
  arrange(PatientID, VisitDay) %>%
  group_by(PatientID, Week) %>%
  mutate(PANSS_Total = ceiling(mean(PANSS_Total))) %>%
  slice(which.min(PANSS_Total)) %>%
  ungroup(Week) %>%
  mutate(TreatmentWeek = row_number(Week) - 1) %>%
  ungroup() %>%
  select(-c(P1:G16))

# -------------------------------------
# Records Baseline (Day 0) data into a separate data frame for later use.
# -------------------------------------
Baseline = CleanedStudyData %>% 
  filter(VisitDay == 0) %>% 
  select(PatientID, PANSS_Total) %>%
  rename(BaselinePANSS = PANSS_Total)

# -------------------------------------
# Create additional variables in data for model
# -------------------------------------
Data = CleanedStudyData %>%
  arrange(PatientID, VisitDay) %>%
  # -------------------------------------
# Join Baseline data into entire table and populate 
# through each PatientID
# -------------------------------------
left_join(Baseline, by="PatientID") %>%
  # -------------------------------------
# Look at only Treatment Weeks 1, 2, 4, and 6
# -------------------------------------
filter(TreatmentWeek <=2 | TreatmentWeek == 4 | TreatmentWeek == 6) %>%
  # -------------------------------------
# Classify patient PANSS scores to CGI-S levels
# -------------------------------------
mutate(
  CGI = as.factor(
    case_when(
      PANSS_Total < 75 ~ "Mild",
      PANSS_Total >= 75 & PANSS_Total < 95 ~ "Moderate",
      PANSS_Total >= 95 & PANSS_Total < 116 ~ "Marked",
      PANSS_Total >= 116 ~ "Severe"
    )
  )
) %>%
  # -------------------------------------
# Record change in PANSS score from baseline observation
# -------------------------------------
mutate(
  PctChange = round((PANSS_Total - BaselinePANSS) / BaselinePANSS, 2)
) %>%
  # -------------------------------------
# Classify patient PANSS changes to CGI-I levels
# -------------------------------------
mutate(
  ImprovementLevel = as.factor(
    case_when(
      PctChange > 0.1 ~ "Worsen",
      TreatmentWeek == 1 & PctChange > -0.19 ~ "None",
      TreatmentWeek == 1 & PctChange <= -0.19 & PctChange > -0.4 ~ "Minimal",
      TreatmentWeek == 1 & PctChange <= -0.4 ~ "Much",
      TreatmentWeek == 2 & PctChange > -0.23 ~ "None",
      TreatmentWeek == 2 & PctChange <= -0.23 & PctChange > -0.45 ~ "Minimal",
      TreatmentWeek == 2 & PctChange <= -0.45 ~ "Much",
      TreatmentWeek == 4 & PctChange > -0.26 ~ "None",
      TreatmentWeek == 4 & PctChange <= -0.26 & PctChange > -0.51 ~ "Minimal",
      TreatmentWeek == 4 & PctChange <= -0.51 ~ "Much",
      TreatmentWeek == 6 & PctChange > -0.28 ~ "None",
      TreatmentWeek == 6 & PctChange <= -0.28 & PctChange > -0.53 ~ "Minimal",
      TreatmentWeek == 6 & PctChange <= -0.53 ~ "Much",
      TRUE ~ "Baseline"
    )
  )
)

##########################################
############ VISUALIZATIONS ##############
##########################################
# -------------------------------------
# PANSS scores across treatment group by treatment week
# -------------------------------------
p1 = Data %>%
  group_by(Study, TxGroup, TreatmentWeek) %>%
  summarise(avgPANSS = mean(PANSS_Total)) %>%
  ggplot(aes(TreatmentWeek, avgPANSS, color=TxGroup)) +
  geom_line() +
  facet_wrap(~Study, scales="free", ncol=1) +
  labs(title = "PANSS Scores by Week across Study and TxGroup")
ggsave("temp/Objective1/PANSSTotalByTxGroup.png", p1, width=7, height=7)
# -------------------------------------
# PANSS scores and CGI levels across treatment groups 
# by treatment week
# -------------------------------------
p2 = Data %>% group_by(Study, TxGroup, TreatmentWeek) %>% 
  summarise(PANSS_Total = mean(PANSS_Total)) %>%
  ungroup() %>%
  ggplot(aes(TreatmentWeek, PANSS_Total)) +
  geom_hline(yintercept=58, linetype="dashed", lwd=0.5, color="red", alpha=0.7) +
  geom_hline(yintercept=75, linetype="dashed", lwd=0.5, color="red", alpha=0.7) +
  geom_hline(yintercept=95, linetype="dashed", lwd=0.5, color="red", alpha=0.7) +
  geom_hline(yintercept=116, linetype="dashed", lwd=0.5, color="red", alpha=0.7) +
  geom_line(aes(color=Study), lwd=0.8) +
  # geom_smooth(method="lm", formula = y ~ poly(x, 4, raw=T), se=T, aes(color=Study), lwd=0.8) +
  scale_y_continuous(breaks=c(58, 75, 95, 116), labels=c("Mild", "Moderate", "Marked", "Severe")) +
  facet_wrap(~TxGroup, scales="free_y", ncol=1) +
  labs(title = "PANSS Scores by Week across TxGroup with CGI Markers")
ggsave("temp/Objective1/PANSSTotalWithCGIByTxWeek.png", p2, width=7, height=7)
# -------------------------------------
# CGI-I levels across studies and treatment group 
# by treatment week
# -------------------------------------
p3 = Data %>%
  group_by(Study, TxGroup, TreatmentWeek, ImprovementLevel) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  complete(Study, TxGroup, TreatmentWeek, ImprovementLevel, fill=list(n=0)) %>%
  group_by(Study, TxGroup, TreatmentWeek) %>%
  mutate(pct = n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(factor(TreatmentWeek), pct, fill=factor(ImprovementLevel, levels=c("None", "Much", "Minimal", "Worsen", "Baseline")))) +
  geom_bar(position="stack", stat="identity") +
  facet_grid(TxGroup~Study, scales="free_x") +
  scale_fill_manual(
    values=c(
      "Baseline" = "black",
      "Worsen" = "#FF0000",
      "None" = "#CCCCCC",
      "Minimal" = "#FFCC00",
      "Much" = "#339900"
    ),
    name="Improvement Level"
  ) +
  labs(title = "Improvent Level Makeup by Treatment Group and Study", x = "TreatmentWeek", y="Percent of Observations")
ggsave("temp/Objective1/ImprovementLevelsByStudyandTxGroup.png", p3, width=15, height=7)
# -------------------------------------
# Change in CGI-I levels across studies between treatment 
# groups by treatment week
# -------------------------------------
p4 = Data %>% 
  group_by(Study, TxGroup, TreatmentWeek, ImprovementLevel) %>% 
  summarise(n = n()) %>% 
  ungroup() %>%
  complete(Study, TxGroup, TreatmentWeek, ImprovementLevel, fill=list(n=0)) %>% 
  group_by(Study, TxGroup, TreatmentWeek) %>% 
  mutate(pct = n/sum(n)) %>% 
  ungroup() %>% 
  group_by(Study, TreatmentWeek, ImprovementLevel) %>%
  arrange(Study, TreatmentWeek, ImprovementLevel) %>%
  # + = bigger treatment | - = bigger control
  summarise(diff = pct - lag(pct)) %>%
  drop_na() %>%
  ggplot(aes(factor(TreatmentWeek), diff, fill=factor(ImprovementLevel, levels=c("Much", "Minimal", "None", "Worsen", "Baseline")))) +
  geom_bar(position="dodge", stat="identity") +
  facet_wrap(~Study, scales="free_x", nrow=1) +
  scale_fill_manual(
    values=c(
      "Baseline" = "black",
      "Worsen" = "#FF0000",
      "None" = "#CCCCCC",
      "Minimal" = "#FFCC00",
      "Much" = "#339900"
    ),
    name="Improvement Level"
  ) +
  labs(title = "Change in Improvement Level Percentage by Study", x = "TreatmentWeek", y="Control vs. Treatment")
ggsave("temp/Objective1/ImprovementLevelsDiffByStudy.png", p4, width=15, height=7)

##########################################
########### HYPOTHESIS TESTING ###########
##########################################
# Model 1
fit1 = lm(PANSS_Total ~ Study*TxGroup*BaselinePANSS, Data)
summary(fit1)
# confint(fit1)

# Model 1
fit2 = lm(PANSS_Total ~Study*TxGroup, Data)
summary(fit2)
# confint(fit2)

# Model 3
fit3 = lm(PANSS_Total ~ BaselinePANSS*TxGroup, Data)
summary(fit3)
# confint(fit3)
# -------------------------------------
# Perform ANOVA testing to check model robustness
# -------------------------------------
anova(fit3, fit2, fit1) # Fit 1 best