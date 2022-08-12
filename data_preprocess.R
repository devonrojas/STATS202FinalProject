# FILENAME: data_preprocess.R
# AUTHOR: Devon Rojas
# Summary: Combines the five study datasets into one
#          main dataset and adds a Week column. Saves
#          to CSV file in data/processed folder.

################################################
### COMBINE DATA FILES AND PREP FOR ANALYSIS ###
################################################

studyA = read.csv("data/raw/Study_A.csv")
studyB = read.csv("data/raw/Study_B.csv")
studyC = read.csv("data/raw/Study_C.csv")
studyD = read.csv("data/raw/Study_D.csv")
studyE = read.csv("data/raw/Study_E.csv")

study_agg = merge(studyA, studyB, all=T)
study_agg = merge(study_agg, studyC, all=T)
study_agg = merge(study_agg, studyD, all=T)
study_agg = merge(study_agg, studyE, all=T)

# Add a Week variable to data  indicating running Week of measurements
# (Note: This may not be a continuous variable due to inconsistency in
#        patient visits.)
CleanedStudyData = study_agg %>%
  mutate(Week = as.integer(ceiling((VisitDay + 1) / 7)))

write.csv(CleanedStudyData, "data/processed/cleaned_data_obj4.csv")

# Scan for additional measurements for a given patient and visit day.
# If more than one measurement exists, eliminate all but one.
CleanedStudyData = CleanedStudyData %>%
  group_by(PatientID, VisitDay) %>%
  distinct(VisitDay, .keep_all=T) %>%
  ungroup(PatientID)

write.csv(CleanedStudyData, "data/processed/cleaned_data_obj3.csv")