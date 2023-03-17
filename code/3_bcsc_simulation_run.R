## Supporting code for "Identifying sources of disparities in surveillance mammography performance and personalized 
## recommendations for supplemental breast imaging: A simulation study" by Hubbard, Pujol, Alhajjar, Edoh, and Martin
##
## This file runs surveillance mammography simulations and outputs csv files of simulation results and pdf figures
## appearing in the manuscript
##
## Last updated: 3/16/23

library(truncnorm)
library(dplyr)
library(ggplot2)
library(tidyr)
library(here)
library(xtable)
library(pROC)

source(here("code","0_bcsc_simulation_prepare_data.R"))
source(here("code","1_bcsc_simulation_inputs.R"))
source(here("code","2_bcsc_simulation_funcn.R"))

## Description of simulations:
## model1: No race/ethnicity effects
## model2: Race/ethnicity effects on primary cancer and treatment characteristics
## model4: Race/ethnicity effects on SDOH
## model4: Race/ethnicity effects on primary cancer and treatment characteristics and SDOH
## model5: Race/ethnicity effects on primary cancer and treatment characteristics,  SDOH and completeness of 
##         ascertainment of second cancers

## Run simulations

NSIM <- 10000 # number of simulations
set.seed(727)

model1 <- sim1_funcn(nrace = race_freq, Nsim = NSIM, sdoh_flag = 0, cancer_flag = 0)
model2 <- sim1_funcn(nrace = race_freq, Nsim = NSIM, sdoh_flag = 0, cancer_flag = 1)
model3 <- sim1_funcn(nrace = race_freq, Nsim = NSIM, sdoh_flag = 1, cancer_flag = 0)
model4 <- sim1_funcn(nrace = race_freq, Nsim = NSIM, sdoh_flag = 1, cancer_flag = 1)
model5 <- sim1_funcn(nrace = race_freq, Nsim = NSIM, sdoh_flag = 1, cancer_flag = 1,
                     cancer_misclass_prob = c(0.936,0.918,0.94,0.956,0.96)) 
                     # misclassification probabilities based on
                     # "timely follow-up" rates from McCarthy 2016 AJPM

## Output numerical results

sens.out     <- rbind(model1$SENS,model2$SENS,model3$SENS,model4$SENS,model5$SENS)
cdr.out      <- rbind(model1$CDR,model2$CDR,model3$CDR,model4$CDR,model5$CDR)
ppv.out      <- rbind(model1$PPV,model2$PPV,model3$PPV,model4$PPV,model5$PPV)
recall.out   <- rbind(model1$RECALL,model2$RECALL,model3$RECALL,model4$RECALL,model5$RECALL)
canc.out     <- rbind(model1$CANC,model2$CANC,model3$CANC,model4$CANC,model5$CANC)
fn.out       <- rbind(model1$FN,model2$FN,model3$FN,model4$FN,model5$FN)
fp.out       <- rbind(model1$FP,model2$FP,model3$FP,model4$FP,model5$FP)
highrisk.out <- rbind(model1$HIGHRISK,model2$HIGHRISK,model3$HIGHRISK,model4$HIGHRISK,model5$HIGHRISK)
lowsens.out  <- rbind(model1$LOWSENS,model2$LOWSENS,model3$LOWSENS,model4$LOWSENS,model5$LOWSENS)

write.csv(sens.out,here("results","sens.csv"),row.names = F)
write.csv(cdr.out,here("results","cdr.csv"),row.names = F)
write.csv(ppv.out,here("results","ppv.csv"),row.names = F)
write.csv(recall.out,here("results","recall.csv"),row.names = F)
write.csv(canc.out,here("results","canc.csv"),row.names = F)
write.csv(fn.out,here("results","fn.csv"),row.names = F)
write.csv(fp.out,here("results","fp.csv"),row.names = F)
write.csv(highrisk.out,here("results","highrisk.csv"),row.names = F)
write.csv(lowsens.out,here("results","lowsens.csv"),row.names = F)

# Generate numerical summaries 
meanres <- function(x){
  apply(x,2,mean, na.rm = T)
}

mean.sens     <- rbind(meanres(model1$SENS),meanres(model2$SENS),meanres(model3$SENS),meanres(model4$SENS),meanres(model5$SENS))
mean.cdr      <- rbind(meanres(model1$CDR),meanres(model2$CDR),meanres(model3$CDR),meanres(model4$CDR),meanres(model5$CDR))
mean.ppv      <- rbind(meanres(model1$PPV),meanres(model2$PPV),meanres(model3$PPV),meanres(model4$PPV),meanres(model5$PPV))
mean.recall   <- rbind(meanres(model1$RECALL),meanres(model2$RECALL),meanres(model3$RECALL),meanres(model4$RECALL),meanres(model5$RECALL))
mean.canc     <- rbind(meanres(model1$CANC),meanres(model2$CANC),meanres(model3$CANC),meanres(model4$CANC),meanres(model5$CANC))
mean.fn       <- rbind(meanres(model1$FN),meanres(model2$FN),meanres(model3$FN),meanres(model4$FN),meanres(model5$FN))
mean.fp       <- rbind(meanres(model1$FP),meanres(model2$FP),meanres(model3$FP),meanres(model4$FP),meanres(model5$FP))
mean.highrisk <- rbind(meanres(model1$HIGHRISK),meanres(model2$HIGHRISK),meanres(model3$HIGHRISK),meanres(model4$HIGHRISK),meanres(model5$HIGHRISK))
mean.lowsens  <- rbind(meanres(model1$LOWSENS),meanres(model2$LOWSENS),meanres(model3$LOWSENS),meanres(model4$LOWSENS),meanres(model5$LOWSENS))

## Output manuscript figures

pdf(here("results","sensitivity.pdf"), height = 10, width = 10)
par(mfrow = c(3,2), mai = c(0.6, 0.4, 0.2, 0.1))
boxplot(model1$SENS, main = "(A) No Race Differences",ylim=c(0,1))
boxplot(model2$SENS, main = "(B) Cancer Varies by Race",ylim=c(0,1))
boxplot(model3$SENS, main = "(C) SDOH Varies by Race",ylim=c(0,1))
boxplot(model4$SENS, main = "(D) SDOH and Cancer Vary by Race",ylim=c(0,1))
boxplot(model5$SENS, main = "(E) Differential Outcome Misclassification by Race",ylim=c(0,1))
dev.off()

pdf(here("results","cdr.pdf"), height = 10, width = 10)
par(mfrow = c(3,2), mai = c(0.6, 0.4, 0.2, 0.1))
boxplot(model1$CDR, main = "(A) No Race Differences",ylim=c(0,0.03))
boxplot(model2$CDR, main = "(B) Cancer Varies by Race",ylim=c(0,0.03))
boxplot(model3$CDR, main = "(C) SDOH Varies by Race",ylim=c(0,0.03))
boxplot(model4$CDR, main = "(D) SDOH and Cancer Vary by Race",ylim=c(0,0.03))
boxplot(model5$CDR, main = "(E) Differential Outcome Misclassification by Race",ylim=c(0,0.03))
dev.off()

pdf(here("results","ppv.pdf"), height = 10, width = 10)
par(mfrow = c(3,2), mai = c(0.6, 0.4, 0.2, 0.1))
boxplot(model1$PPV, main = "(A) No Race Differences",ylim=c(0,0.6))
boxplot(model2$PPV, main = "(B) Cancer Varies by Race",ylim=c(0,0.6))
boxplot(model3$PPV, main = "(C) SDOH Varies by Race",ylim=c(0,0.6))
boxplot(model4$PPV, main = "(D) SDOH and Cancer Vary by Race",ylim=c(0,0.6))
boxplot(model5$PPV, main = "(E) Differential Outcome Misclassification by Race",ylim=c(0,0.6))
dev.off()

pdf(here("results","recall.pdf"), height = 10, width = 10)
par(mfrow = c(3,2), mai = c(0.6, 0.4, 0.2, 0.1))
boxplot(model1$RECALL, main = "(A) No Race Differences",ylim=c(0,0.15))
boxplot(model2$RECALL, main = "(B) Cancer Varies by Race",ylim=c(0,0.15))
boxplot(model3$RECALL, main = "(C) SDOH Varies by Race",ylim=c(0,0.15))
boxplot(model4$RECALL, main = "(D) SDOH and Cancer Vary by Race",ylim=c(0,0.15))
boxplot(model5$RECALL, main = "(E) Differential Outcome Misclassification by Race",ylim=c(0,0.15))
dev.off()

pdf(here("results","cancer.pdf"), height = 10, width = 10)
par(mfrow = c(3,2), mai = c(0.6, 0.4, 0.2, 0.1))
boxplot(model1$CANC, main = "(A) No Race Differences",ylim=c(0,0.05))
boxplot(model2$CANC, main = "(B) Cancer Varies by Race",ylim=c(0,0.05))
boxplot(model3$CANC, main = "(C) SDOH Varies by Race",ylim=c(0,0.05))
boxplot(model4$CANC, main = "(D) SDOH and Cancer Vary by Race",ylim=c(0,0.05))
boxplot(model5$CANC, main = "(E) Differential Outcome Misclassification by Race",ylim=c(0,0.05))
dev.off()

pdf(here("results","fn.pdf"), height = 10, width = 10)
par(mfrow = c(3,2), mai = c(0.6, 0.4, 0.2, 0.1))
boxplot(model1$FN, main = "(A) No Race Differences",ylim=c(0,0.02))
boxplot(model2$FN, main = "(B) Cancer Varies by Race",ylim=c(0,0.02))
boxplot(model3$FN, main = "(C) SDOH Varies by Race",ylim=c(0,0.02))
boxplot(model4$FN, main = "(D) SDOH and Cancer Vary by Race",ylim=c(0,0.02))
boxplot(model5$FN, main = "(E) Differential Outcome Misclassification by Race",ylim=c(0,0.02))
dev.off()

pdf(here("results","fp.pdf"), height = 10, width = 10)
par(mfrow = c(3,2), mai = c(0.6, 0.4, 0.2, 0.1))
boxplot(model1$FP, main = "(A) No Race Differences",ylim=c(0,0.15))
boxplot(model2$FP, main = "(B) Cancer Varies by Race",ylim=c(0,0.15))
boxplot(model3$FP, main = "(C) SDOH Varies by Race",ylim=c(0,0.15))
boxplot(model4$FP, main = "(D) SDOH and Cancer Vary by Race",ylim=c(0,0.15))
boxplot(model5$FP, main = "(E) Differential Outcome Misclassification by Race",ylim=c(0,0.15))
dev.off()

pdf(here("results","sensitivity threshold.pdf"), height = 10, width = 10)
par(mfrow = c(3,2), mai = c(0.6, 0.4, 0.2, 0.1))
boxplot(model1$LOWSENS, main = "(A) No Race Differences",ylim=c(0,0.5))
boxplot(model2$LOWSENS, main = "(B) Cancer Varies by Race",ylim=c(0,0.5))
boxplot(model3$LOWSENS, main = "(C) SDOH Varies by Race",ylim=c(0,0.5))
boxplot(model4$LOWSENS, main = "(D) SDOH and Cancer Vary by Race",ylim=c(0,0.5))
boxplot(model5$LOWSENS, main = "(E) Differential Outcome Misclassification by Race",ylim=c(0,0.5))
dev.off()

pdf(here("results","risk threshold.pdf"), height = 10, width = 10)
par(mfrow = c(3,2), mai = c(0.6, 0.4, 0.2, 0.1))
boxplot(model1$HIGHRISK, main = "(A) No Race Differences",ylim=c(0,0.2))
boxplot(model2$HIGHRISK, main = "(B) Cancer Varies by Race",ylim=c(0,0.2))
boxplot(model3$HIGHRISK, main = "(C) SDOH Varies by Race",ylim=c(0,0.2))
boxplot(model4$HIGHRISK, main = "(D) SDOH and Cancer Vary by Race",ylim=c(0,0.2))
boxplot(model5$HIGHRISK, main = "(E) Differential Outcome Misclassification by Race",ylim=c(0,0.2))
dev.off()

