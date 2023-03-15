## Supporting code for "Identifying sources of disparities in surveillance mammography performance and personalized 
## recommendations for supplemental breast imaging: A simulation study" by Hubbard, Pujol, Alhajjar, Edoh, and Martin
##
## This file provides the simulation functions that generate simulated data for surveillance mammography cohort
##
## Last updated: 3/15/23

## ----------------------------------------------------------------------------
# Function to iterate Nsim replicates of simulation
#   nrace =vector with number of patients in each race/ethnicity group ordered as Asian, Black, Hispanic, Other, White
#   sdoh_flag = binary indicator of whether to simulate SDOH effects
#   cancer_flag = binary indicator of whether to simulate cancer/treatment effects
#   riskthresh = threshhold for classigying "high second cancer risk"
#   sensethresh = threshhold for classifying "low sensitivity"
#   cancer_misclass_prob = vector of probability of capturing second cancers by race

sim1_funcn <- function(nrace, Nsim = 10000, sdoh_flag = 1, cancer_flag = 1, 
                       riskthresh = 0.025, sensthresh = 0.5, cancer_misclass_prob = rep(1,length(nrace))){
  
  race_cat   <- c('Asian/PI','Black', 'Hispanic' ,'Other','White')
  cdr.hat    <- matrix(nrow = Nsim, ncol = (length(race_cat)+1))
  ppv.hat    <- matrix(nrow = Nsim, ncol = (length(race_cat)+1))
  sens.hat   <- matrix(nrow = Nsim, ncol = (length(race_cat)+1))
  recall.hat <- matrix(nrow = Nsim, ncol = (length(race_cat)+1))
  fn.hat     <- matrix(nrow = Nsim, ncol = (length(race_cat)+1))
  fp.hat     <- matrix(nrow = Nsim, ncol = (length(race_cat)+1))
  cancer.hat <- matrix(nrow = Nsim, ncol = (length(race_cat)+1))
  highrisk   <- matrix(nrow = Nsim, ncol = (length(race_cat)+1))
  lowsens    <- matrix(nrow = Nsim, ncol = (length(race_cat)+1))
  
  for (j in 1:Nsim){
    sim.data <- NULL
    for(i in 1:length(race_cat)){
      temp <- do.call(cbind,bc_simulation(nrace[i], i, sdoh_flag, cancer_flag, obscancprob = cancer_misclass_prob[i]))
      
      # calculate race specific performance measures
      sens.hat[j,(i+1)]    <- mean(temp[temp$seccancfu1yr_obs == 1,"posmam_hat"])
      cdr.hat[j,(i+1)]     <- mean(temp[,"posmam_hat"]*temp[,"seccancfu1yr_obs"])
      ppv.hat[j,(i+1)]     <- mean(temp[temp$posmam_hat == 1,"seccancfu1yr_obs"])
      recall.hat[j,(i+1)]  <- mean(temp[,"posmam_hat"])
      cancer.hat[j,(i+1)]  <- mean(temp[,"seccancfu1yr_obs"])
      fn.hat[j,(i+1)]      <- mean((temp[,"posmam_hat"]==0)*temp[,"seccancfu1yr_obs"])
      fp.hat[j,(i+1)]      <- mean(temp[temp$seccancfu1yr_obs == 0,"posmam_hat"])
      highrisk[j,(i+1)]    <- mean(temp[,"cancmod_prob"] > riskthresh)
      lowsens[j,(i+1)]     <- mean(temp[,"pos_prob"] < sensthresh)
      
      if (i > 1) colnames(sim.data) <- colnames(temp)
      sim.data <- rbind(sim.data,temp)
      
    }
    
    # calculate overall performance measures
    sens.hat[j,1]    <- mean(sim.data[sim.data$seccancfu1yr_obs == 1,"posmam_hat"])
    cdr.hat[j,1]     <- mean(sim.data[,"posmam_hat"]*sim.data[,"seccancfu1yr_obs"])
    ppv.hat[j,1]     <- mean(sim.data[sim.data$posmam_hat == 1,"seccancfu1yr_obs"])
    recall.hat[j,1]  <- mean(sim.data[,"posmam_hat"])
    cancer.hat[j,1]  <- mean(temp[,"seccancfu1yr_obs"])
    fn.hat[j,1]      <- mean((temp[,"posmam_hat"]==0)*temp[,"seccancfu1yr_obs"])
    fp.hat[j,1]      <- mean(temp[temp$seccancfu1yr_obs == 0,"posmam_hat"])
    highrisk[j,1]    <- mean(temp[,"cancmod_prob"] > riskthresh)
    lowsens[j,1]     <- mean(temp[,"pos_prob"] < sensthresh)
  }  
  
  sens.hat <- data.frame(sens.hat)
  names(sens.hat) <- c("Overall",race_cat)
  cdr.hat <- data.frame(cdr.hat)
  names(cdr.hat) <- c("Overall",race_cat)
  ppv.hat <- data.frame(ppv.hat)
  names(ppv.hat) <- c("Overall",race_cat)
  recall.hat <- data.frame(recall.hat)
  names(recall.hat) <- c("Overall",race_cat)
  cancer.hat <- data.frame(cancer.hat)
  names(cancer.hat) <- c("Overall",race_cat)
  fn.hat <- data.frame(fn.hat)
  names(fn.hat) <- c("Overall",race_cat)
  fp.hat <- data.frame(fp.hat)
  names(fp.hat) <- c("Overall",race_cat)
  highrisk <- data.frame(highrisk)
  names(highrisk) <- c("Overall",race_cat)
  lowsens  <- data.frame(lowsens)
  names(lowsens) <- c("Overall",race_cat)
  
  return(list(SENS = sens.hat,CDR = cdr.hat,PPV = ppv.hat,RECALL = recall.hat,
              CANC = cancer.hat, FN = fn.hat, FP = fp.hat, HIGHRISK = highrisk,
              LOWSENS = lowsens))
}

## ----------------------------------------------------------------------------
## Simulate data for each race/ethnicity group
#   n = number of observations to simulate
#   i = race/ethnicity index
#   sdoh_flag = binary indicator of whether to simulate SDOH effects
#   cancer_flag = binary indicator of whether to simulate cancer/treatment effects
#   obscancprob = vector of probability of capturing second cancers by race

bc_simulation <- function(n, i, sdoh_flag, cancer_flag, obscancprob) { 

  data <- data.frame(dxage = rep(0,n))
  
  ##################### Patient clinical characteristics ####################
  
  # age at primary breast cancer diagnosis 
  if (cancer_flag == 0){ data$dxage<- data.frame(dxage = round(rtruncnorm(n, a=dxage_min, b=dxage_max, 
                                             mean=dxage_mean, sd=dxage_sd), 1)) 
  }else{
    minage  = dxage_min_race[[i]]
    maxage  = dxage_max_race[[i]]
    meanage = dxage_mean_race[[i]]
    sdage   = dxage_sd_race[[i]]
    data$dxage<- data.frame(dxage = round(rtruncnorm(n, a=minage, b=maxage, 
                                               mean=meanage, sd=sdage), 1)) 
  } 
  
  # breast density
  # 1: Almost entirely fat(<25% fibrglandular)
  # 2: Scttrd fibroglandular tiss(25%-50%)
  # 3: Heterogeneously dense(50%-75%)
  # 4: Extremely dense(>75%)
  if (cancer_flag == 0){ data$density = data.frame(t(rmultinom(n, size=1, 
                                      prob  = c(density_freq$n[-5])/sum(density_freq$n[-5])))) 
  }else{
    density_freq_race <- as.numeric(density_freq_strat[i,])
    data$density = data.frame(t(rmultinom(n, size=1, 
                   prob  = c(density_freq_race)/sum(density_freq_race)))) 
  } 
  
  # Time since breast cancer diagnosis in years
  data$yrsincdx = rexp(n, rate = 1/mean_time)

  ################ Patient social determinants of health ################  
  
  # BMI
  if (sdoh_flag == 0){ data$bmi = round(rtruncnorm(n, a=15, b=90, mean=bmi_mean, sd=bmi_sd), 1)  
  }else{
    meanbmi <- bmi_mean_race[i]
    sdbmi   <- bmi_sd_race[i]
    data$bmi = round(rtruncnorm(n, a=15, b=90, mean=meanbmi, sd=sdbmi), 1)  
  }
  
  # zip code median income
  if (sdoh_flag == 1){ 
    data$geo_income = rgamma(n,rate = geo_income_mean_race[i]/(geo_income_sd_race[i]^2),
                             shape = (geo_income_mean_race[i]^2)/(geo_income_sd_race[i]^2))
  } else{
    data$geo_income = rgamma(n,rate = geo_income_mean/(geo_income_sd^2),
                             shape = geo_income_mean^2/(geo_income_sd^2))
  }
  
  # Regular mammography utilization
  # Time since last mammogram (months): Regular: < 15, irregular >= 15 
  data$regscreen <- rbinom(n,1, prob = (sdoh_flag*c(as.numeric(regscreen_freq_race[i,2])/sum(regscreen_freq_race[i,]))+
                                          (1-sdoh_flag)*c(unlist(regscreen_freq[2,2])/sum(regscreen_freq[,2]))))
  
  # Education 
  if (sdoh_flag == 0){ data$educat = data.frame(t(rmultinom(n, size=1, 
                                   prob  = c(educat_freq$n[-5])/sum(educat_freq$n[-5])))) 
  }else{
    educat_freq_race <- as.numeric(educat_freq_strat[i,])
    data$educat = data.frame(t(rmultinom(n, size=1, 
                  prob  = c(educat_freq_race)/sum(educat_freq_race)))) 
  } 
  
  ########################## Cancer and Treatment Characteristics  ########################## 
  
  # Mode of detection of primary cancer 
  if (cancer_flag == 0){ data$mod1cat = data.frame(t(rmultinom(n, size=1, 
                                                          prob  = c(mod1cat_freq$n[-4])/sum(mod1cat_freq$n[-4])))) 
  }else{
    mod1cat_freq_race <- as.numeric(mod1cat_freq_strat[i,])
    data$mod1cat = data.frame(t(rmultinom(n, size=1, 
                                         prob  = c(mod1cat_freq_race)/sum(mod1cat_freq_race)))) 
  } 
  
  # ER/PR status, primary breast cancer - combine 1 and 2; 3 and 4
  # 1: ER-, PR-
  # 2: ER-, PR+
  # 3: ER+, PR-
  # 4: ER+, PR+
  data$erplus1 = rbinom(n, 1, (cancer_flag*c(as.numeric(erplus1_freq_race[i,2]))/sum(erplus1_freq_race[i,])+
                                 (1-cancer_flag)*c(unlist(erplus1_freq[2,2])/sum(erplus1_freq[,2]))))
  
  # HER2 status
  data$her1 = rbinom(n, 1, (cancer_flag*c(her1_freq_race[i,2])/sum(her1_freq_race[i,])+
                              (1-cancer_flag)*c(her1c_freq[2,2]/sum(her1c_freq[,2]))))

  # Primary breast cancer stage
  if (cancer_flag == 0){ data$stage = data.frame(t(rmultinom(n, size=1, 
                                                               prob  = c(stage_freq$n)/sum(stage_freq$n)))) 
  }else{
    stage_freq_race <- as.numeric(stage_freq_strat[i,])
    data$stage = data.frame(t(rmultinom(n, size=1, 
                                          prob  = c(stage_freq_race)/sum(stage_freq_race)))) 
  }   
  #Primary breast cancer surgery +/- radiation 
  # 1: BCS without RT
  # 2: BCS with RT
  # 3: Mastectomy
  if(cancer_flag == 1){
    data$srgrt1 = data.frame(t(rmultinom(n, size=1, prob  = c(as.numeric(srgrt1_freq_race[i,])/sum(srgrt1_freq_race[i,])))))
  } else {
    data$srgrt1 = data.frame(t(rmultinom(n, size=1, prob  = c(srgrt1_freq$n)/sum(srgrt1_freq$n))))
  }
  # Systemic adjuvant primary breast cancer treatment
  data$adjthrpy = rbinom(n, 1, (cancer_flag*c(adjthrpy_freq_race[i,2])/sum(adjthrpy_freq_race[i,])+
                                  (1-cancer_flag)*c(adjthrpy_freq[2,2]/sum(adjthrpy_freq[,2]))))
  
  ####################### Mammography Result and Second Cancer Dx ####################### 
  
  new_data <- as.matrix(cbind(1,data$dxage, data$yrsincdx, data$bmi,
                              data$density[,-1], data$geo_income, data$regscreen,
                              data$stage[,-1], data$erplus1,data$her1, data$srgrt1[,-1],
                              data$adjthrpy,data$educat[,-1],data$mod1cat[,-1]))
  
  co        <- c(cancmodcoef)$x
  cancmod_p <- new_data%*%co
  data$cancmod_prob <- exp(cancmod_p)/(1+exp((cancmod_p)))
  
  data$seccancfu1yr_cutoff_hat <- rbinom(n, 1, prob = data$cancmod_prob)
  
  # generate misclassified, second cancer indicator
  data$seccancfu1yr_obs <- 0 # assume observed cancer status has imperfect sensitivity but perfect specificity
  data$seccancfu1yr_obs <- ifelse(data$seccancfu1yr_cutoff_hat == 1, 
                                  rbinom(sum(data$seccancfu1yr_cutoff_hat),1,obscancprob),data$seccancfu1yr_obs)
  #Split data for mammography output
  data_posmamm <- data[data$seccancfu1yr_cutoff_hat == 1,]
  data_negmamm <- data[data$seccancfu1yr_cutoff_hat == 0,]
  
  new_negmamm <- as.matrix(cbind(1,data_negmamm$dxage, data_negmamm$yrsincdx, data_negmamm$bmi,
                                 data_negmamm$density[,-1], data_negmamm$geo_income, data_negmamm$regscreen,
                                 data_negmamm$stage[,-1], data_negmamm$erplus1,data_negmamm$her1, data_negmamm$srgrt1[,-1],
                                 data_negmamm$adjthrpy,data_negmamm$educat[,-1],data_negmamm$mod1cat[,-1]))
  if (nrow(data_posmamm)> 0){ # only run if non-zero number of cancers
    new_posmamm <- as.matrix(cbind(1,data_posmamm$dxage, data_posmamm$yrsincdx, data_posmamm$bmi,
                                   data_posmamm$density[,-1], data_posmamm$geo_income, data_posmamm$regscreen,
                                   data_posmamm$erplus1,data_posmamm$her1, data_posmamm$srgrt1[,-1],
                                   data_posmamm$adjthrpy,data_posmamm$educat[,-1],data_posmamm$mod1cat[,-1]))
    co <- c(sensmodcoef)$x
    sensmod_p <- c(co%*%t(new_posmamm))
    data_posmamm$pos_prob <- exp(sensmod_p)/(1+exp(sensmod_p)) 
    data_posmamm$posmam_hat <- rbinom(n = nrow(data_posmamm), 1, prob = data_posmamm$pos_prob)
  }

  # Sensitivity - positive mammogram & positive cancer
  #predicted values  
  co <- c(specmodcoef)$x
  specmod_p <- c(co%*%t(new_negmamm))
  data_negmamm$pos_prob <- exp(specmod_p)/(1+exp(specmod_p))
  
  # Generate positive mammography using FP prob if no cancer, TP prob if cancer
  data_negmamm$posmam_hat <- rbinom(n = nrow(data_negmamm), 1, prob = data_negmamm$pos_prob)
  
  if (nrow(data_posmamm)> 0) df <- rbind(data_negmamm, data_posmamm)
  else df <- data_negmamm
  
  # Generate predicted mammography sensitivity for everyone
  co <- c(sensmodcoef)$x
  allpred <- as.matrix(cbind(1,df$dxage, df$yrsincdx, df$bmi,
                             df$density[,-1], df$geo_income, df$regscreen,
                             df$erplus1,df$her1, df$srgrt1[,-1],
                             df$adjthrpy,df$educat[,-1],df$mod1cat[,-1]))
  sensmod_p <- c(co%*%t(allpred))
  df$pos_prob <- exp(sensmod_p)/(1+exp(sensmod_p)) 
#  print("predictions generated")
  
  return(df)
}    

