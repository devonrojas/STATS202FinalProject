# FILENAME: FinalProjectObjective2.Rmd
# AUTHOR: Devon Rojas
# Summary: Performs clustering analysis on study datasets

##########################################
##### SET SEED AND IMPORT LIBRARIES ######
##########################################
set.seed(2)

library(tidyverse)
library(ggplot2)
library(cluster)
library(patchwork)

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
  ungroup()

# -------------------------------------
# Use only baseline measurements for analysis
# -------------------------------------
Baseline = CleanedStudyData %>% 
  filter(VisitDay == 0) %>%
  mutate(
    # -------------------------------------
    # Calculate patient composite score (bi-polar index)
    # -------------------------------------
    CompositeScore = rowSums(across(N1:N7)) - rowSums(across(P1:P7)),
    # -------------------------------------
    # Classify CGI-S levels numerically for PANSS scores
    # -------------------------------------
    CGIScore = case_when(
      PANSS_Total < 45 ~ 1,
      PANSS_Total >= 45 & PANSS_Total < 58 ~ 2,
      PANSS_Total >= 58 & PANSS_Total < 75 ~ 3,
      PANSS_Total >= 75 & PANSS_Total < 95 ~ 4,
      PANSS_Total >= 95 & PANSS_Total < 116 ~ 5,
      PANSS_Total >= 116 & PANSS_Total < 145 ~ 6,
      PANSS_Total >= 145 ~ 7
    )
  )
# select(PatientID, P1:PANSS_Total) %>%
# mutate(across(P1:PANSS_Total, .names="Baseline{.col}")) %>%
# select(-c(P1:PANSS_Total))

##########################################
## DETERMINING OPTIMAL K FOR CLUSTERING ##
##########################################
# -------------------------------------
# Remove non-numeric attributes and scale all remaining
# -------------------------------------
kdata = Baseline %>% 
  select(where(is.numeric)) %>%
  select(-c(PatientID, SiteID, RaterID, AssessmentID, VisitDay, Week, TreatmentWeek)) %>%
  select(-c(P1:G16)) %>%
  mutate(across(everything(), scale))
# -------------------------------------
# Compute total withinss score for various k
# -------------------------------------
krange = 2:12
kmodel = list()
kwithinss = rep(0, length(krange))
for(i in 1:length(krange)) {
  print(paste("Now computing kmeans for k =", krange[i], sep=" "))
  k = kmeans(kdata, centers=krange[i], nstart=25, iter.max=300)
  kwithinss[i] = k$tot.withinss
  kmodel[[i]] = k
}
# -------------------------------------
# Compute silhouette scores for various k and plot results
# -------------------------------------
colors = brewer_pal(type="qual", palette="Set3")
png("temp/KmeansSilhouettes.png", width=1920, height=1080)
ob = par(mfrow=c(2, 6))
sils = rep(0, length(krange))
for(i in 1:length(krange)) {
  ss = cluster::silhouette(kmodel[[i]]$cluster, dist(kdata))
  plot(ss, col=colors(i+1), border=NA)
  s = mean(ss[,3])
  sils[i] = s
}
dev.off()
par(ob)

# -------------------------------------
# Plot total withinss score across k
# -------------------------------------
WithinsPlot = tibble(k = krange, withinss = kwithinss) %>%
  ggplot(aes(x=k, y=withinss)) +
  geom_line() +
  scale_x_continuous(breaks=krange)
# -------------------------------------
# Plot average silhouette score across k
# -------------------------------------
MeanSilhouettePlot = tibble(k = krange, sil_vals = sils) %>%
  ggplot(aes(k, sil_vals)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(breaks=krange)
KSelectionPlot = WithinsPlot / MeanSilhouettePlot
ggsave("temp/Objective2/KSelection.png", KSelectionPlot, width=11, height=8)

# -------------------------------------
# Plot specific silhouette plots for k's of interest
# -------------------------------------
png("temp/Objective2/KmeansSilhouettes4.png", width=7, height=8, units="in", res=200)
ss = silhouette(kmodel[[3]]$cluster, dist(kdata))
plot(ss, col=colors(4), border=NA, main="Silhouette Plot for K=4")
abline(v=mean(ss[,3]), lty=2, col="red")
dev.off()
png("temp/Objective2/KmeansSilhouettes10.png", width=7, height=8, units="in", res=200)
ss = silhouette(kmodel[[9]]$cluster, dist(kdata))
plot(ss, col=colors(10), border=NA, main="Silhouette Plot for K=10")
abline(v=mean(ss[,3]), lty=2, col="red")
dev.off()

##########################################
############ CLUSTERING STEP #############
##########################################
CHOSENK = 4 # Chosen via elbow
# CHOSENK = which.max(sils) + 1 # Chosen via silhouette. Offset because K started at 2

# -------------------------------------
# Run clustering on final k value and append results to original data
# -------------------------------------
km.out = kmeans(kdata, centers=CHOSENK, nstart=50, iter.max=200)
ClusteredData = km.out %>% augment(Baseline)
# -------------------------------------
# Compute cluster means for CGI-S score, Composite score, and PANSS score
# -------------------------------------
ClusteredDataMeans = ClusteredData %>% group_by(.cluster) %>% summarise(across(c(CGIScore, CompositeScore, PANSS_Total), mean))

##########################################
######### CLUSTER VISUALIZATIONS #########
##########################################
# -------------------------------------
# Cluster proportions
# -------------------------------------
ClusterProportionPlot = ClusteredData %>% 
  ggplot(aes(.cluster, fill=.cluster)) +
  geom_bar(position="stack", stat="count")
ggsave("temp/Objective2/ClusterProportions.png", ClusterProportionPlot, width=11, height=8)

# -------------------------------------
# CGI-S x PANSS score by cluster
# -------------------------------------
CGIxPANSSCountwMeansPlot = ClusteredData %>%
  ggplot(aes(CGIScore, PANSS_Total, color=.cluster)) +
  geom_count() +
  geom_point(data=ClusteredDataMeans, color="black") +
  geom_vline(data=ClusteredDataMeans, aes(xintercept=CGIScore), color="red", linetype="dashed") +
  geom_hline(data=ClusteredDataMeans, aes(yintercept=PANSS_Total), color="red", linetype="dashed") +
  geom_text(data=ClusteredDataMeans, label="Cluster Avg", vjust=-1.25, color="black") +
  facet_wrap(~.cluster, ncol=2) +
  scale_x_continuous(breaks=1:6, labels=c("Normal", "Borderline", "Mild", "Moderate", "Marked", "Severe"))
ggsave("temp/Objective2/CGIxPANSSCountswClusterMeans.png", CGIxPANSSCountwMeansPlot, width=11, height=8)

# -------------------------------------
# Cluster proportions by study
# -------------------------------------
ClusterProportionsByStudyPlot = ClusteredData %>% group_by(Study, .cluster) %>% summarise(n = n()) %>%
  ggplot(aes(.cluster, n, fill=Study)) +
  geom_bar(position="fill", stat="identity")
ggsave("temp/Objective2/ClusterProportionsByStudy.png", ClusterProportionsByStudyPlot, width=11, height=8)

# CGI scores per cluster
# ClusteredData %>% group_by(.cluster, CGIScore) %>% summarise(n=n()) %>%
#   ggplot(aes(CGIScore, n, fill=.cluster)) +
#   geom_col(position=position_dodge2(width=0.9, preserve="single")) +
#   scale_x_continuous(breaks=c(1:6))

# -------------------------------------
# Composite and PANSS score densities by cluster
# -------------------------------------
CompositeViolin = ClusteredData %>%
  ggplot(aes(.cluster, CompositeScore, fill=.cluster)) +
  geom_violin()
PANSSViolin = ClusteredData %>%
  ggplot(aes(.cluster, PANSS_Total, fill=.cluster)) +
  geom_violin()
ViolinPlot = PANSSViolin | CompositeViolin

# -------------------------------------
# CGI-S density plot by cluster
# -------------------------------------
CGIDensity = ClusteredData %>%
  ggplot(aes(CGIScore)) +
  geom_density(aes(fill=.cluster), alpha=0.5) +
  facet_wrap(~.cluster) +
  scale_x_continuous(breaks=1:6, labels=c("Normal", "Borderline", "Mild", "Moderate", "Marked", "Severe"))

CombinedPlot = ViolinPlot / CGIDensity
ggsave("temp/Objective2/Violin_and_CGIDensityPlot.png", CombinedPlot, width=11, height=8)

# -------------------------------------
# Composite by PANSS score cluster plot
# -------------------------------------
ClusterPlot = ClusteredData %>% 
  ggplot(aes(CompositeScore, PANSS_Total, color=.cluster)) +
  geom_point() +
  geom_text(data=ClusteredDataMeans, aes(label=.cluster), color="black", size=8, fontface="bold")
ggsave("temp/Objective2/ClusterPlot.png", ClusterPlot, width=11, height=8)