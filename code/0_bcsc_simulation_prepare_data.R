## Supporting code for "Identifying sources of disparities in surveillance mammography performance and personalized 
## recommendations for supplemental breast imaging: A simulation study" by Hubbard, Pujol, Alhajjar, Edoh, and Martin
##
## This file calculates summary statistics and outputs values for simulation inputs based on analysis of BCSC data
## Note that this file requires raw BCSC data (available by request only) to run
##
## Last updated: 5/19/23


## ----------------------------------------------------------------------------
## Read in data and recode variables
## ----------------------------------------------------------------------------

## read in csv data
ab505 <- read.csv("/Users/rhubb/Library/CloudStorage/Box-Box/to do/ADJOINT/data/ADJOINT_BCSC_data.csv")
colnames(ab505)[1] <- "catype1"
ab505 <- ab505[!is.na(ab505$caraceeth),] # exclude subjects with missing race
colnames(ab505)[c(5,12,14,15,21)] <-c('her1', 'adjthrpy','bmi','density',
                                        'yrsincdx') 

## recode variables to match simulation specifications
ab505$posmam <- ifelse(ab505$assminit_eod_c %in% c(0,3,4,5),1,0)
ab505$cdr <- ifelse(ab505$posmam == 1 & ab505$seccancfu1yr_cutoff_c == 1,1,0)
ab505$dxage  <- as.double(substr(ab505$dxage,1,2)) # truncate age at 90
ab505$geo_income  <- ab505$geo_income_med_2010_ACS/10000 # rescale to stabilize OR
ab505$regscreen <- 1*(ab505$lmamintmo <15)
ab505$invasive <- ab505$catype1-1 # recode to 1 = invasive, 0 = DCIS
ab505$erplus1 <- ifelse(ab505$erplus1 %in% c(1,2),0,1) #1 = ER+ or unknown, 0 = ER-
ab505$srgrt1 <- ifelse(ab505$srgrt1 %in% c(3),1,ab505$srgrt1) 
ab505$srgrt1 <- ifelse(ab505$srgrt1 %in% c(4,5),3,ab505$srgrt1) 
                      #1 = (BCS wo or unknown RT), 2 = (BCS w RT), 3 = (mastectomy)
ab505$caraceeth <- recode(ab505$caraceeth,"1" = "White","2" = "Black", 
                        "3" = "Asian/PI", "4" = "Asian/PI", "5" = "Other",
                        "6" = "Other", "7" = "Other", "8" = "Hispanic")
ab505$adjthrpy <- ifelse(ab505$adjthrpy == 0,0,1) 
ab505$age50 <- ifelse(ab505$dxage <= 50, 1, 0)
ab505$agecat <- cut(ab505$dxage, 
                    breaks=c(-Inf, 39,49,59,69,79, Inf), 
                    labels=c("<40","40-49","50-59","60-69","70-79","80+"))
ab505$yrsinccat <- cut(ab505$yrsincdx, 
                    breaks=c(-Inf, 5,10,15, Inf), 
                    labels=c("<6","6-10","11-15",">15"))
ab505$bmicat <- cut(ab505$bmi, 
                       breaks=c(-Inf, 18.499,20,25,30,40, Inf), 
                       labels=c("<18.5","18.5-20","21-25","26-30","31-40",">40"))
ab505$incomecat <- cut(ab505$geo_income_med_2010_ACS, 
                       breaks=c(-Inf, 35000,50000,75000,100000, Inf), 
                       labels=c("<35,000","35,000-50,000","50,000-75,000",
                                "75,000-100,000",">90,000"))
ab505$lmamintcat <- cut(ab505$lmamintmo, right = FALSE,
                       breaks=c(-Inf, 12,24,36, Inf), 
                       labels=c("<12","12-23","24-35",">35"))
ab505$mod1cat <- ifelse(ab505$mod1_c %in% c(2,4),2,ifelse(ab505$mod1_c %in% c(3,5),3,ab505$mod1_c))
ab505$invasiveY <- ifelse(ab505$invasive == 1, "Yes","No")
ab505$erplus1Y <- ifelse(ab505$erplus1 == 1, "Positive","Negative")
ab505$her1Y <- ifelse(ab505$her1 == 1, "Positive","Negative")
ab505$chemo1Y <- ifelse(ab505$chemo1_r == 1, "Yes","No")
ab505$srgcat <- ifelse(ab505$srg1_c%in% c(3,4),"Mastectomy","Breast Conserving Surgery") 
ab505$rad1Y <- ifelse(ab505$rad1_r == 1, "Yes","No")
ab505$stage <- floor(ab505$ajcc1_v8_c)

## ----------------------------------------------------------------------------
## Create table of summary statistics stratified by race
## ----------------------------------------------------------------------------
# function to generate summary statistics table
summ.data <- function(x, strat){
  if (x[2] == 1) LEV <- levels(factor(ab505[,x[1]])) # if var is categorical save all levels of variable before stratifying
  
  if (!is.null(strat)) dat <- ab505[strat,]
  else dat <- ab505
  
  d <- dat[,x[1]]
  if (x[2] == 1){ # summarize categorical var with N (%)
    d <- factor(d, levels = LEV)
    n <- table(d, useNA = "always")
    tot <- sum(table(d)) # compute total non-missing
    p <-round(100*n/tot,1)
    out <- vector(length = length(n))
    out[!is.na(names(n))] <- paste(n[!is.na(names(n))]," (",sprintf("%.1f",p[!is.na(names(n))]),")",sep="")
    out[is.na(names(n))] <- n[is.na(names(n))] # for missings just give total without %
    out <- c("",out)
    names(out) <- c(x[3],names(table(d)))
  } else{ # summarize continuous variable with Median (IQR)
    med <- round(median(d,na.rm = T),1)
    iqr <- round(quantile(d, probs = c(0.25,0.75), na.rm = T),1)
    out <- paste(sprintf("%.1f",med), " (",sprintf("%.1f",iqr[1]),", ",sprintf("%.1f",iqr[2]),")",sep = "")
    out <- c(out,sum(is.na(d)))
    names(out) <- c(x[3],"Missing")
  }  
  return(out)
}

tabvars <- c("agecat","yrsinccat","bmicat","density","incomecat","educat","lmamintcat",
             "mod1cat","stage","erplus1Y","her1Y","chemo1Y", "srgcat","rad1Y")
tabtype <- rep(1,length(tabvars))
tablab  <- c("Age (years), N (%)","Year since primary diagnosis, N (%)","BMI, N (%)",
             "Breast density, N (%)","Income, N (%)","Education, N (%)", "Time since last mammogram (months), N (%)",
             "Mode of detection, N (%)","Stage, N (%)",
             "ER status, N (%)","HER2 status, N (%)", 
             "Adjuvant systemic, N (%)","Surgery, N (%)","Radiation, N (%)")

all     <- unlist(apply(cbind(tabvars,tabtype,tablab),1,summ.data, strat = NULL))
asian   <- unlist(apply(cbind(tabvars,tabtype,tablab),1,summ.data, strat = (ab505$caraceeth == "Asian/PI")))
black   <- unlist(apply(cbind(tabvars,tabtype,tablab),1,summ.data, strat = (ab505$caraceeth == "Black")))
hispanic<- unlist(apply(cbind(tabvars,tabtype,tablab),1,summ.data, strat = (ab505$caraceeth == "Hispanic")))
white   <- unlist(apply(cbind(tabvars,tabtype,tablab),1,summ.data, strat = (ab505$caraceeth == "White")))
other   <- unlist(apply(cbind(tabvars,tabtype,tablab),1,summ.data, strat = (ab505$caraceeth == "Other")))

tab1 <- cbind(all,asian,black,hispanic,white,other)
print.xtable(tab1)

## ----------------------------------------------------------------------------
## Patient characteristics
## ----------------------------------------------------------------------------

# Stratified age distribution parameters
write(by(ab505$dxage,ab505$caraceeth,min),here("inputs","dxage_min_race"))
write(by(ab505$dxage,ab505$caraceeth,mean),here("inputs","dxage_mean_race"))
write(by(ab505$dxage,ab505$caraceeth,max),here("inputs","dxage_max_race"))
write(by(ab505$dxage,ab505$caraceeth,sd),here("inputs","dxage_sd_race"))

# Unstratified breast density 
write.csv(count(ab505, vars = density),here("inputs","density_freq.csv"), row.names = F)

# Stratified breast density
write.csv(do.call(rbind,by(ab505$density, ab505$caraceeth, table)),here("inputs","density_freq_race.csv"), row.names = F)

## ----------------------------------------------------------------------------
## Patient social determinants of health
## ----------------------------------------------------------------------------

# Statified BMI
write(by(ab505$bmi,ab505$caraceeth,mean,na.rm = T),here("inputs","bmi_mean_race"))
write(by(ab505$bmi,ab505$caraceeth,sd,na.rm = T),here("inputs","bmi_sd_race"))

# Stratified income
write(tapply(ab505$geo_income,ab505$caraceeth,mean,na.rm = T),here("inputs","geo_income_mean_race"))
write(tapply(ab505$geo_income,ab505$caraceeth,sd,na.rm = T),here("inputs","geo_income_sd_race"))

# Unstratified regular mammography screening 
write.csv(count(ab505, vars = regscreen),here("inputs","regscreen_freq.csv"), row.names = F)

# Stratified regular mammography screening 
write.csv(do.call(rbind,tapply(ab505$regscreen,ab505$caraceeth,table)),here("inputs","regscreen_freq_race.csv"), row.names = F)

# Unstratified education
write.csv(count(ab505, vars = educat),here("inputs","educat_freq.csv"), row.names = F)

# Stratified education
write.csv(do.call(rbind,by(ab505$educat, ab505$caraceeth, table)),here("inputs","educat_freq_race.csv"), row.names = F)

## ----------------------------------------------------------------------------
## Primary cancer and treatment characteristics
## ----------------------------------------------------------------------------

# Unstratified mode of detection
write.csv(count(ab505, vars = mod1cat),here("inputs","mod1cat_freq.csv"), row.names = F)

# Stratified  mode of detection
write.csv(do.call(rbind,by(ab505$mod1cat, ab505$caraceeth, table)),here("inputs","mod1cat_freq_race.csv"), row.names = F)

# Unstratified ER/PR status
write.csv(count(ab505, vars = erplus1),here("inputs","erplus1_freq.csv"), row.names = F)

# Stratified ER/PR status
write.csv(do.call(rbind,by(ab505$erplus1, ab505$caraceeth, table)),here("inputs","erplus1_freq_race.csv"), row.names = F)

# Unstratified HER2 status
write.csv(table(ab505$her1),here("inputs","her1c_freq.csv"), row.names = F)

# Stratified HER2 status
write.csv(do.call(rbind,tapply(ab505$her1,ab505$caraceeth,table)),here("inputs","her1_freq_race.csv"), row.names = F)

# Unstratified stage
write.csv(count(ab505, vars = stage),here("inputs","stage_freq.csv"), row.names = F)

# Stratified stage
write.csv(do.call(rbind,by(ab505$stage, ab505$caraceeth, table)),here("inputs","stage_freq_race.csv"), row.names = F)

# Unstratified surgery and RT
write.csv(count(ab505, vars = srgrt1),here("inputs","srgrt1_freq.csv"), row.names = F)

# Stratified surgery and RT
write.csv(do.call(rbind,by(ab505$srgrt1, ab505$caraceeth, table)),here("inputs","srgrt1_freq_race.csv"), row.names = F)

# Unstratified adjuvant therapy
write.csv(table(ab505$adjthrpy),here("inputs","adjthrpy_freq.csv"), row.names = F)

# Stratified HER2 status
write.csv(do.call(rbind,tapply(ab505$adjthrpy,ab505$caraceeth,table)),here("inputs","adjthrpy_freq_race.csv"), row.names = F)

## ----------------------------------------------------------------------------
## Logistic regression relating patient characteristics, cancer characteristics,
## and SDOH to outcome measures
## ----------------------------------------------------------------------------

# Cancer risk model
ctrl <- trainControl(method = "cv",  # Cross-validation method (e.g., "cv", "repeatedcv")
                     number = 10,    # Number of folds
                     classProbs = TRUE,  # Preserve class probabilities for AUC calculation
                     summaryFunction = twoClassSummary)  # AUC summary function

# Train the logistic regression model
ab505canc.dat <- ab505[,c("seccancfu1yr_cutoff_c","dxage", "yrsincdx", "bmi", 
                          "density", "geo_income", "regscreen", "stage", 
                          "erplus1", "her1", "srgrt1", "adjthrpy", "educat", 
                          "mod1cat")]
ab505canc.dat <- na.omit(ab505canc.dat)
ab505canc.dat$seccancfu1yr_cutoff_c <- factor(ab505canc.dat$seccancfu1yr_cutoff_c, labels = c("No","Yes"))

cancmod <- train(seccancfu1yr_cutoff_c ~  dxage + yrsincdx + bmi + 
                 factor(density) + geo_income + regscreen + factor(stage) + 
                 erplus1 + her1 + factor(srgrt1) + adjthrpy + factor(educat) +
                 factor(mod1cat),           
               data = ab505canc.dat,      
               method = "glm",      # Use "glm" for logistic regression
               trControl = ctrl,    # Cross-validation control
               metric = "ROC")      # Performance metric (e.g., "Accuracy", "ROC")

# Compute the cross-validated AUC
auc_cancmodcv <- cancmod$results$ROC
#> auc_cancmodcv
#[1] 0.6242052

cancmod <- glm(seccancfu1yr_cutoff_c ~  dxage + yrsincdx + bmi + 
                 factor(density) + geo_income + regscreen + factor(stage) + 
                 erplus1 + her1 + factor(srgrt1) + adjthrpy + factor(educat) +
                 factor(mod1cat),
               family = "binomial", data = ab505)

# evaluate model performance
roc_cancmod <- roc(cancmod$y, cancmod$fitted.values)
auc_cancmod <- auc(roc_cancmod)
#auc(roc_cancmod)
#Area under the curve: 0.7009

cancmod$coef[1] <- -6.92 # manually adjust intercept to calibrate cancer rate
write.csv(cancmod$coef,here("inputs","cancmodcoef.csv"), row.names = F)

# Sensitivity - positive mammogram & cancer
# Cancer risk model
ctrl <- trainControl(method = "cv",  # Cross-validation method (e.g., "cv", "repeatedcv")
                     number = 10,    # Number of folds
                     classProbs = TRUE,  # Preserve class probabilities for AUC calculation
                     summaryFunction = twoClassSummary)  # AUC summary function

# Train the logistic regression model
ab505sens.dat <- ab505[ab505$seccancfu1yr_cutoff_c == 1,c("posmam","dxage", "yrsincdx", 
                          "density", "geo_income", "regscreen", "bmi",
                          "erplus1", "her1","srgrt1", "adjthrpy", "educat", 
                          "mod1cat")]
ab505sens.dat <- na.omit(ab505sens.dat)
ab505sens.dat$posmam <- factor(ab505sens.dat$posmam, labels = c("No","Yes"))

sensmod <- train(posmam ~  dxage + yrsincdx + bmi +
                   factor(density) + geo_income + regscreen + 
                   erplus1 + her1 + factor(srgrt1) + adjthrpy + factor(educat) +
                   factor(mod1cat),           
                 data = ab505sens.dat,        
                 method = "glm",      # Use "glm" for logistic regression
                 trControl = ctrl,    # Cross-validation control
                 metric = "ROC")      # Performance metric (e.g., "Accuracy", "ROC")

# Compute the cross-validated AUC
auc_sensmodcv <- sensmod$results$ROC
#> auc_sensmodcv
#[1] 0.6591667

sensmod <- glm(posmam ~  dxage + yrsincdx + bmi +
                 factor(density) + geo_income + regscreen  + 
                 erplus1 + her1 + factor(srgrt1) + adjthrpy + factor(educat)+
                 factor(mod1cat), family = "binomial",
               data = ab505[ab505$seccancfu1yr_cutoff_c == 1,])
# evaluate model performance
roc_sensmod <- roc(sensmod$y, sensmod$fitted.values)
auc_sensmod <- auc(roc_sensmod)
#auc_sensmod
#Area under the curve: 0.8361

sensmod$coef[1] <- -5.1 # manually adjust intercept to calibrate sensitivity
write.csv(sensmod$coef,here("inputs","sensmodcoef.csv"), row.names = F)

# Specificity - positive mammogram & no cancer 
specmod <- glm(posmam ~  dxage + yrsincdx + bmi + 
                 factor(density) + geo_income + regscreen + factor(stage) + 
                 erplus1 + her1 + factor(srgrt1) + adjthrpy + factor(educat) +
                 factor(mod1cat), family = "binomial",
               data = ab505[ab505$seccancfu1yr_cutoff_c == 0,])

specmod$coef[1] <- -3.3 # manually adjust intercept to calibrate specificity
write.csv(specmod$coef,here("inputs","specmodcoef.csv"), row.names = F)

## Output table of regression coefficients for paper supplement
stage_ind <- which(regexpr("stage",names(coef(cancmod)))>0)[1]
bmi_ind   <- which(regexpr("bmi",names(coef(cancmod)))>0)[1]
her1_ind   <- which(regexpr("her1",names(coef(cancmod)))>0)[1]
sens_coef_wstage <- c(round(coef(sensmod)[1:(stage_ind-1)],3),rep("NA",3),
                      round(coef(sensmod)[(stage_ind):length(coef(sensmod))],3))
coef.tab <- cbind(round(coef(cancmod),3),sens_coef_wstage)
row.names(coef.tab) <- c("Intercept","Age at diagnosis (years)", "Time since diagnosis (years)",
                         "BMI (kg/m$^2$)", "B: scattered areas of fibroglandular density",
                         "C: heterogeneously dense", "D: extremely dense",
                         "Median neighborhood income (per \\$10,000)",
                         "Prior mammogram within 15 months", "Stage 1","Stage 2","Stage 3",
                         "ER+ primary cancer","HER2+ primary cancer","Breast conserving surgery with radiation",
                         "Mastectomy","Adjuvant therapy", "High school grad or GED",
                         "Some college/tech","College grad or postgrad",
                         "Interval","Clinical")
colnames(coef.tab) <- c("Cancer risk model","Mammography sensitivity model")
xtable(coef.tab)

