## Supporting code for "Identifying sources of disparities in surveillance mammography performance and personalized 
## recommendations for supplemental breast imaging: A simulation study" by Hubbard, Pujol, Alhajjar, Edoh, and Martin
##
## This file provides input values for simulations based on analysis of BCSC data
##
## Last updated: 3/15/23


## ----------------------------------------------------------------------------
## Patient characteristics
## ----------------------------------------------------------------------------

nptx      <- 65510 # total sample size
race_freq <- c(1453,5128,865,298,57702) # sample size by race/ethnicity

# Unstratified age distribution parameters
dxage_mean <- 58.46808
dxage_sd   <- 11.21753
dxage_min  <- 21
dxage_max  <- 90

# Stratified age distribution parameters
dxage_min_race <- read.delim(here("inputs","dxage_min_race"), sep = " ", header = F)
dxage_mean_race <- read.delim(here("inputs","dxage_mean_race"), sep = " ", header = F)
dxage_max_race <- read.delim(here("inputs","dxage_max_race"), sep = " ", header = F)
dxage_sd_race <- read.delim(here("inputs","dxage_sd_race"), sep = " ", header = F)

# Unstratified breast density 
density_freq <- read.csv(here("inputs","density_freq.csv"))

# Stratified breast density
density_freq_strat <- read.csv(here("inputs","density_freq_race.csv"))

# Time since diagnosis
mean_time  <- 6.731278

## ----------------------------------------------------------------------------
## Patient social determinants of health
## ----------------------------------------------------------------------------

# Unstratified BMI 
bmi_mean   <- 27.99185
bmi_sd     <- 6.375765

# Statified BMI
bmi_mean_race <- unlist(read.delim(here("inputs","bmi_mean_race"), sep = " ", header = F))
bmi_sd_race   <- unlist(read.delim(here("inputs","bmi_sd_race"), sep = " ", header = F))

# Unstratified income 
geo_income_mean <- 7.005789
geo_income_sd   <- 2.213084

# Stratified income
geo_income_mean_race <- unlist(read.delim(here("inputs","geo_income_mean_race"), sep = " ", header = F))
geo_income_sd_race   <- unlist(read.delim(here("inputs","geo_income_sd_race"), sep = " ", header = F))

# Unstratified regular mammography screening 
regscreen_freq <- read.csv(here("inputs","regscreen_freq.csv")) 

# Stratified regular mammography screening 
regscreen_freq_race <- read.csv(here("inputs","regscreen_freq_race.csv"))

# Unstratified education
educat_freq <- read.csv(here("inputs","educat_freq.csv"))

# Stratified education
educat_freq_strat <- read.csv(here("inputs","educat_freq_race.csv"))

## ----------------------------------------------------------------------------
## Primary cancer and treatment characteristics
## ----------------------------------------------------------------------------

# Unstratified mode of detection
mod1cat_freq <- read.csv(here("inputs","mod1cat_freq.csv"))

# Stratified  mode of detection
mod1cat_freq_strat <- read.csv(here("inputs","mod1cat_freq_race.csv"))

# Unstratified ER/PR status
erplus1_freq <- read.csv(here("inputs","erplus1_freq.csv"))

# Stratified ER/PR status
erplus1_freq_race <- read.csv(here("inputs","erplus1_freq_race.csv"))

# Unstratified HER2 status
her1c_freq <- read.csv(here("inputs","her1c_freq.csv")) 

# Stratified HER2 status
her1_freq_race <- read.csv(here("inputs","her1_freq_race.csv"))

# Unstratified stage
stage_freq <- read.csv(here("inputs","stage_freq.csv"))

# Stratified stage
stage_freq_strat <- read.csv(here("inputs","stage_freq_race.csv"))

# Unstratified surgery and RT
srgrt1_freq <- read.csv(here("inputs","srgrt1_freq.csv"))

# Stratified surgery and RT
srgrt1_freq_race <- read.csv(here("inputs","srgrt1_freq_race.csv"))

# Unstratified adjuvant therapy
adjthrpy_freq <- read.csv(here("inputs","adjthrpy_freq.csv")) 

# Stratified HER2 status
adjthrpy_freq_race <- read.csv(here("inputs","adjthrpy_freq_race.csv"))

## ----------------------------------------------------------------------------
## Read in regression coefficients relating patient characteristics, cancer characteristics,
## and SDOH to outcome measures
## ----------------------------------------------------------------------------

# Cancer risk model
cancmodcoef <- read.csv(here("inputs","cancmodcoef.csv"))

# Sensitivity - positive mammogram & cancer
sensmodcoef <- read.csv(here("inputs","sensmodcoef.csv"))

# Specificity - positive mammogram & no cancer 
specmodcoef <- read.csv(here("inputs","specmodcoef.csv"))


