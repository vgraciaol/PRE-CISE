#------------------------------------------------------------------------------#
# This script illustrates PRE-CISE and its implementation in a four-state      #
# Sick-Sicker cohort Markov model calibrated to survival probability,          #
# prevalence, and proportion of Sick using the IMIS algorithm.                 #
#                                                                              #
# Authors:                                                                     #
#     Fernando Alarid-Escudero, PhD, <falarid@stanford.edu>                    #
#     Valeria Gracia Olvera, MSc, <vgracia@stanford.edu>                       #
#                                                                              #
# Last update: 08/29/26                                                        #
#------------------------------------------------------------------------------#
rm(list = ls()) # clean environment

# Calibration Specifications ----------------------------------------------
# Model: 4-state Sick-Sicker Markov Model
# Inputs to be calibrated: 
#   p_S1S2 - probability to become sicker when sick
#   hr_S1  - hazard ratio of death in sick vs healthy
#   hr_S2  - hazard ratio of death in sicker vs healthy
# Targets: 
#   Surv     - survival probability
#   Prev     - prevalence (sick and sicker)
#   PropSick - proportion of Sick
# Search method: Random search using Latin-Hypercube Sampling
# Goodness-of-fit measure: Sum of Log-Likelihood

# Load libraries and functions --------------------------------------------
# calibration functionality
library(lhs)
library(matrixStats) # package used for summary statistics
library(FME)         # For identifiability and collinearity analysis
library(tidyverse)
options(dplyr.summarise.inform = FALSE) # do not show summarise info
library(mvtnorm)

# visualization
library(plotrix)
library(psych)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(scatterplot3d) # three inputs to estimate
library(ggcube)

# global sensitivity analysis
library(sensitivity)

# calibration
source("R/IMIS_function_SickSicker.R") # from IMIS package, adapted for this analysis
# Old version 
# devtools::install_version("IMIS", version = "0.1", repos = "http://cran.us.r-project.org")

# functions
source("R/helper_functions.R") # from dampack package
source("R/SickSicker_functions.R")

# Variables ---------------------------------------------------------------
# Number of simulations (coverage analysis)
n_sim <- 300

# Seed number
n_seed <- 111124

# Set seed
set.seed(n_seed)

# Vector of times (in years) at which the model will be calibrated
v_times <- c(10, 20, 30)

# Target data -------------------------------------------------------------
load("data/SickSicker_CalibTargets.RData")
l_targets <- SickSicker_targets 

# Gather targets in a data.frame
df_targets <- rbind(l_targets$Surv %>% 
                      mutate(type = "Target",
                             Outcome = "Survival",
                             abbrev_outcome = "Surv") %>% 
                      select(time, Outcome, abbrev_outcome, type, value, lb, ub, se),
                    l_targets$Prev %>% 
                      mutate(type = "Target",
                             Outcome = "Prevalence",
                             abbrev_outcome = "Prev") %>% 
                      select(time, Outcome, abbrev_outcome, type, value, lb, ub, se),
                    l_targets$PropSick %>% 
                      mutate(type = "Target",
                             Outcome = "Proportion of Sick",
                             abbrev_outcome = "PropSick") %>% 
                      select(time, Outcome, abbrev_outcome, type, value, lb, ub, se)) %>% 
  filter(time %in% v_times) # filter by the years to be calibrated

# Plot the targets
plot_targets(df_targets)

# Model -------------------------------------------------------------------
# -inputs are parameters to be estimated through calibration
# -outputs correspond to the target data

# Check that it works
v_params_test <- c(p_S1S2 = 0.105, hr_S1 = 3, hr_S2 = 10)
sick_sicker_cstm(v_params_test) # It works!

# Specify calibration parameters ------------------------------------------
# Names and number of input parameters to be calibrated
v_param_names <- c("p_S1S2", # probability to become sicker when sick
                   "hr_S1",  # hazard ratio of death in sick vs healthy
                   "hr_S2")  # hazard ratio of death in sicker vs healthy
n_param <- length(v_param_names)

# Range on input search space
v_lb <- c(p_S1S2 = 0.01, hr_S1 = 1.0, hr_S2 = 5.0)  # lower bound
v_ub <- c(p_S1S2 = 0.50, hr_S1 = 4.5, hr_S2 = 15.0) # upper bound

# Number of calibration targets
v_outcome_names <- c("Survival probability", "Prevalence", "Proportion of Sick")
names(v_outcome_names) <- c("Surv", "Prev", "PropSick")
n_outcomes <- length(v_outcome_names)

# Calibration functions (check that it works) -----------------------------
# View resulting parameter set samples
pairs.panels(sample_prior(n_samp = 1000, v_lb = v_lb, v_ub = v_ub))

# Prior
log_prior(v_params = v_params_test,
          v_lb     = v_lb, 
          v_ub     = v_ub)
log_prior(v_params = sample_prior(10, v_lb = v_lb, v_ub = v_ub),
          v_lb     = v_lb, 
          v_ub     = v_ub)

prior(v_params = v_params_test,
      v_lb     = v_lb, 
      v_ub     = v_ub)
prior(v_params = sample_prior(10, v_lb = v_lb, v_ub = v_ub),
      v_lb     = v_lb, 
      v_ub     = v_ub)

# Likelihood
log_lik(v_params = v_params_test)
log_lik(v_params = sample_prior(10, v_lb = v_lb, v_ub = v_ub))

likelihood(v_params = v_params_test)
likelihood(v_params = sample_prior(10, v_lb = v_lb, v_ub = v_ub))

# Posterior
log_post(v_params = v_params_test,
         v_lb     = v_lb, 
         v_ub     = v_ub)
log_post(v_params = sample_prior(10, v_lb = v_lb, v_ub = v_ub),
         v_lb     = v_lb, 
         v_ub     = v_ub)

posterior(v_params = v_params_test,
          v_lb     = v_lb, 
          v_ub     = v_ub)
posterior(v_params = sample_prior(10, v_lb = v_lb, v_ub = v_ub),
          v_lb     = v_lb, 
          v_ub     = v_ub)

# Define functions --------------------------------------------------------
## Coverage analysis ------------------------------------------------------
# Coverage function
# m_params: Matrix containing sample of parameter sets
run_coverage <- function(m_params){
  
  # Number of simulations
  n_sim <- nrow(m_params) 
  
  # Run model
  df_model <- data.frame(NULL) # empty data.frame to store results
  for(iter_i in 1:n_sim){ # iter_i = 1
    
    # Run model
    l_model_out <- calibration_out(v_params = m_params[iter_i, ])
    
    # Output of the model in a data.frame
    # Modify the next lines based on your model
    df_temp <- rbind(l_model_out$df_Surv,
                     l_model_out$df_Prev,
                     l_model_out$df_PropSick) %>% 
      mutate(iter = iter_i)
    
    # Gather results
    df_model <- rbind(df_model, df_temp) 
    
    # Print progress
    if(iter_i/100==round(iter_i/100,0)) { 
      cat('\r',paste(round(iter_i/n_sim*100,0),"% simulations done",sep=""))
    }
  }
  
  # Summarise
  df_model_summ <- df_model %>%
    group_by(time, Outcome, abbrev_outcome, type) %>%
    summarise(mean_val = mean(value),
              lb       = quantile(value,probs = 0.025),
              ub       = quantile(value, probs = 0.975),
              se       = sd(value)) %>%
    ungroup() %>%
    rename(value = mean_val)
  
  return(list(df_model        = df_model,
              df_model_summ   = df_model_summ)) 
}

## Sensitivity analysis --------------------------------------------------
# Cost function
# v_params: Vector of model parameters.
costFun <- function(v_params) {
  
  # Run model
  df_model_out <- sick_sicker_cstm(v_params)
  
  # Filter and select columns
  ## Surv
  df_model_out_Surv <- df_model_out %>% 
    filter(time %in% v_times) %>% 
    select(time, Surv)
  
  df_targets_Surv <- df_targets %>% 
    filter(abbrev_outcome == "Surv") %>% 
    select(time, value, se) %>% 
    rename(Surv = value)
  
  ## PropSick
  df_model_out_PropSick <- df_model_out %>% 
    filter(time %in% v_times) %>% 
    select(time, PropSick)
  
  df_targets_PropSick <- df_targets %>% 
    filter(abbrev_outcome == "PropSick") %>% 
    select(time, value, se) %>% 
    rename(PropSick = value)
  
  ## Prev
  df_model_out_Prev <- df_model_out %>% 
    filter(time %in% v_times) %>% 
    select(time, Prev)
  
  df_targets_Prev <- df_targets %>% 
    filter(abbrev_outcome == "Prev") %>% 
    select(time, value, se) %>% 
    rename(Prev = value)
  
  # Compute "costs" for each target
  cost <- modCost(model = df_model_out_Surv, 
                  obs   = df_targets_Surv, 
                  err   = "se")
  
  cost <- modCost(model = df_model_out_Prev,
                  obs   = df_targets_Prev,
                  err   = "se",
                  cost  = cost)
  
  cost <- modCost(model = df_model_out_PropSick,
                  obs   = df_targets_PropSick,
                  err   = "se",
                  cost  = cost)
  return(cost)
}

# Function to compute model outputs for global sensitivity analyses
# m_params:     Matrix containing sample of parameter sets
# outcome_name: Target outcome c("Survival probability","Prevalence",
#               "Proportion of sick")
# n_time:       Time (year) for which GSA is to be computed
model_out_GSA <- function(m_params, outcome_name, n_time){
  
  # Empty vector to store results
  v_res <- c()
  
  if(outcome_name == "Survival probability"){
    
    for(i in 1:nrow(m_params)){ # i = 1
      
      # Run model
      df_model_out <- sick_sicker_cstm(m_params[i,])
      
      # Data.frame for each calibration target
      # Survival
      df_Surv <- data.frame(time = df_model_out$time,
                            Outcome = "Survival",
                            abbrev_outcome = "Surv",
                            type = "Model",
                            value = df_model_out$Surv) %>%
        filter(time %in% n_time)
      
      # Store result
      v_res <- c(v_res, df_Surv$value)
      
    }
  }else if(outcome_name == "Prevalence"){
    
    for(i in 1:nrow(m_params)){ # i = 1
      
      # Run model
      df_model_out <- sick_sicker_cstm(m_params[i,])
      
      # Data.frame for each calibration target
      # Prevalence
      df_Prev <- data.frame(time = df_model_out$time,
                            Outcome = "Prevalence",
                            abbrev_outcome = "Prev",
                            type = "Model",
                            value = df_model_out$Prev) %>%
        filter(time %in% n_time)
      
      # Store result
      v_res <- c(v_res, df_Prev$value)
      
    }
  }else if(outcome_name == "Proportion of Sick"){
    
    for(i in 1:nrow(m_params)){ # i = 1
      
      # Run model
      df_model_out <- sick_sicker_cstm(m_params[i,])
      
      # Data.frame for each calibration target
      # Proportion who are Sick ("PropSick") among all those
      # afflicted (Sick+Sicker)
      df_PropSick <- data.frame(time = df_model_out$time,
                                Outcome = "Proportion of Sick",
                                abbrev_outcome = "PropSick",
                                type = "Model",
                                value = df_model_out$PropSick) %>%
        filter(time %in% n_time)
      
      # Store result
      v_res <- c(v_res, df_PropSick$value)
      
    }
  }
  
  # Return vector
  return(v_res)
}

# Global sensitivity analysis (GSA): this function runs soboljansen function 
# from sensitivity package.
# v_outcomes: Vector containing the name of the targets.
# v_times:    Vector containing the times (in years).
# X1:         First random sample (prior distribution).
# X2:         Second random sample (prior distribution).
# n_boot:     Number of bootstrap replicates
# n_conf:     Confidence level for bootstrap confidence intervals.
run_GSA <- function(v_outcomes, v_times, X1, X2, n_boot, n_conf){
  
  df_res_sobol <- data.frame(NULL)
  for(outcome_i in v_outcomes){ # outcome_i = "Prevalence"
    
    for(time_i in v_times){ # time_i = 20
      
      cat(outcome_i, time_i, "\n")
      
      # Compute Sobol Indices
      res_sobol <- soboljansen(model = model_out_GSA, X1 = X1, X2 = X2,
                               outcome_name = outcome_i, n_time = time_i,
                               nboot = n_boot, conf = n_conf)
      
      # Bind results
      df_res_sobol <- rbind(df_res_sobol,
                            data.frame(par = rownames(res_sobol$S),
                                       Outcome = outcome_i,
                                       type = "Main effect",
                                       x = time_i,
                                       res_sobol$S, row.names = NULL),
                            data.frame(par = rownames(res_sobol$T),
                                       Outcome = outcome_i,
                                       type = "Total effect",
                                       x = time_i,
                                       res_sobol$T, row.names = NULL))
      
    }
  }
  return(df_res_sobol)
}

# Compare sensitivity analysis --------------------------------------------
# Define input search space
# Range on input search space
v_lb_Sens <- c(p_S1S2 = 0.00, hr_S1 = 1.00, hr_S2 = 1.00) # lower bound
v_ub_Sens <- c(p_S1S2 = 0.50, hr_S1 = 8.00, hr_S2 = 25.0) # upper bound

# Get two sets of priors
set.seed(n_seed)
m_params_Sens_1 <- sample_prior(n_samp = n_sim,
                                v_lb   = v_lb_Sens,
                                v_ub   = v_ub_Sens)

m_params_Sens_2 <- sample_prior(n_samp = n_sim,
                                v_lb   = v_lb_Sens,
                                v_ub   = v_ub_Sens)

# Run GSA
df_res_sobol <- run_GSA(v_outcomes = as.character(v_outcome_names),
                        v_times    = v_times,
                        X1         = m_params_Sens_1,
                        X2         = m_params_Sens_2,
                        n_boot     = 100, 
                        n_conf     = 0.95)

# Bootstrap confidence interval bounds were manually restricted to 
# the plausible range [0,1]
df_res_sobol <- df_res_sobol %>% 
  mutate(min..c.i. = ifelse(min..c.i.<0, 0, min..c.i.),
         max..c.i. = ifelse(max..c.i.>1, 1, max..c.i.))

# Plot
plot_GSA(df_res_sobol)

# Local sensitivity
SensRes <- sensFun(func  = costFun, 
                   parms = colMeans(m_params_Sens_1) # mean values
)

# Plot
plot_Sfun(SensRes)

# Implement PRE-CISE (narrow) ---------------------------------------------
# 0. Define input search space and plot first coverage
v_lb_1 <- c(p_S1S2 = 0.10, hr_S1 = 1.0, hr_S2 = 1.0) # lower bound
v_ub_1 <- c(p_S1S2 = 0.20, hr_S1 = 1.7, hr_S2 = 2.0) # upper bound

# Run initial coverage analysis
set.seed(n_seed)  # set seed
m_params_1 <- sample_prior(n_samp = n_sim,
                           v_lb   = v_lb_1,
                           v_ub   = v_ub_1)

l_coverage_1  <- run_coverage(m_params = m_params_1)
df_coverage_1 <- l_coverage_1$df_model_summ %>% 
  mutate(type_coverage = "Initial model coverage")

plot_coverage(df_targets  = df_targets,
              df_coverage = df_coverage_1)

## Resize the prior distribution bounds -----------------------------------
### Adjust the prob. to become sicker when sick: p_S1S2 -------------------
# 1. Compute model-data residuals on initial coverage (mean values)
v_params_1 <- colMeans(m_params_1)
modCost1   <- costFun(v_params_1)

# Plot residuals
plot_costFun(Cfun = modCost1, wgt_res = T) # weighted residuals
plot_costFun(modCost1, wgt_res = F)        # unweighted residuals

# 2. Local sensitivity analysis
SensRes1 <- sensFun(func  = costFun, 
                    parms = v_params_1)

# Plot sensitivity
plot_Sfun(SensRes1)

# Save plot
plot_Sfun(SensRes1)

# Rank parameters according to their importance
m_Sens1 <- summary(SensRes1)

# value: value of the parameter 
# Mean: mean sensitivity
# L1: L1 norm
# L2: L2 norm
df_Sens1 <- as.data.frame(m_Sens1)

# Order by L2
df_Sens1[order(df_Sens1$L2, decreasing = T),] 

# 3. Compute the change in the prior distribution bounds based on the 
# sensitivity analysis
df_changes_1 <- data.frame(NULL)
for(target_i in c("PropSick","Prev")){ # target_i = "PropSick"
  
  # Filter elasticities
  df_epsi <- SensRes1 %>% 
    filter(var == target_i) %>% 
    select(x,p_S1S2) %>% 
    rename(epsilon = p_S1S2)
  
  # Change in y --> change in par
  df_temp <- modCost1$residuals %>% 
    filter(name == target_i) %>% 
    select(name,x,obs,mod) %>% 
    mutate(y_change = ifelse(mod != 0, (obs-mod)/mod, 0)) %>% 
    left_join(df_epsi,
              by = "x") %>% 
    mutate(par_change = y_change/epsilon,
           par_change_lb = ifelse((1+par_change)*v_lb_1["p_S1S2"]<0,
                                  NA,(1+par_change)*v_lb_1["p_S1S2"]))
  
  # Bind data
  df_changes_1 <- rbind(df_changes_1,
                        df_temp)
  
}

# Percent change for lower or upper bound
n_change_lb <- min(df_changes_1$par_change_lb, na.rm = T) 

# Resize bounds
v_lb_2 <- v_lb_1
v_ub_2 <- v_ub_1
v_lb_2["p_S1S2"] <- n_change_lb

# Run coverage
set.seed(n_seed)  # set seed
m_params_2 <- sample_prior(n_samp = n_sim,
                           v_lb = v_lb_2,
                           v_ub = v_ub_2)
l_coverage_2 <- run_coverage(m_params = m_params_2)
df_coverage_2 <- l_coverage_2$df_model_summ %>% 
  mutate(type_coverage = "A: Improved model coverage by resizing the probability to become Sicker when Sick")

plot_coverage(df_targets  = df_targets,
              df_coverage = rbind(df_coverage_1,df_coverage_2))

### Adjust the hazard ratio for death in Sicker: hr_S2 --------------------
# 1. Compute model-data residuals on previous coverage (mean values)
v_params_2 <- colMeans(m_params_2)        # mean values
modCost2   <- costFun(v_params_2)

# 2. Local sensitivity analysis
SensRes2 <- sensFun(func  = costFun, 
                    parms = v_params_2)

# Plot sensitivity
plot_Sfun(SensRes2)

# Rank parameters according to their importance
m_Sens2 <- summary(SensRes2)

# value: value of the parameter 
# Mean: mean sensitivity
# L1: L1 norm
# L2: L2 norm
df_Sens2 <- as.data.frame(m_Sens2)

# Order by L2
df_Sens2[order(df_Sens2$L2, decreasing = T),] 

# 3. Compute the change in the prior distribution bounds based on the 
# sensitivity analysis
df_changes_2 <- data.frame(NULL)
for(target_i in c("Surv","Prev")){ # target_i = "PropSick"
  
  # Filter elasticities
  df_epsi <- SensRes2 %>% 
    filter(var == target_i) %>% 
    select(x,hr_S2) %>% 
    rename(epsilon = hr_S2)
  
  # Change in y --> change in par
  df_temp <- modCost2$residuals %>% 
    filter(name == target_i) %>% 
    select(name,x,obs,mod) %>%  
    mutate(y_change = ifelse(mod != 0, (obs-mod)/mod, 0)) %>% 
    left_join(df_epsi,
              by = "x") %>% 
    mutate(par_change = y_change/epsilon,
           par_change_ub = ifelse((1+par_change)*v_ub_2["hr_S2"]<v_lb_2["hr_S2"], 
                                  NA, (1+par_change)*v_ub_2["hr_S2"]))
  
  # Bind data
  df_changes_2 <- rbind(df_changes_2,
                        df_temp)
  
}

# Percent change for lower or upper bound
n_change_ub <- max(df_changes_2$par_change_ub, na.rm = T)

# Resize bounds
v_lb_3 <- v_lb_2
v_ub_3 <- v_ub_2
v_ub_3["hr_S2"] <- n_change_ub

# Run coverage
set.seed(n_seed)  # set seed
m_params_3 <- sample_prior(n_samp = n_sim,
                           v_lb   = v_lb_3,
                           v_ub   = v_ub_3)
l_coverage_3 <- run_coverage(m_params = m_params_3)
df_coverage_3 <- l_coverage_3$df_model_summ %>% 
  mutate(type_coverage = "B: A + improved model coverage by resizing the hazard ratio for death in Sicker")

plot_coverage(df_targets  = df_targets,
              df_coverage = rbind(df_coverage_1,df_coverage_2,df_coverage_3))

### Adjust the hazard ratio for death in Sick: hr_S1 ----------------------
# 1. Compute model-data residuals on previous coverage (mean values)
v_params_3 <- colMeans(m_params_3)      # mean values
modCost3   <- costFun(v_params_3)

# 2. Local sensitivity analysis
SensRes3 <- sensFun(func  = costFun, 
                    parms = v_params_3)

# Plot sensitivity
plot_Sfun(SensRes3)

# 3. Compute the change in the prior distribution bounds based on the 
# sensitivity analysis 
df_changes_3 <- data.frame(NULL)
for(target_i in c("Surv")){ # target_i = "PropSick"
  
  # Filter elasticities
  df_epsi <- SensRes3 %>% 
    filter(var == target_i) %>% 
    select(x,hr_S1) %>% 
    rename(epsilon = hr_S1)
  
  # Change in y --> change in par
  df_temp <- modCost3$residuals %>% 
    filter(name == target_i) %>% 
    select(name,x,obs,mod) %>% 
    mutate(y_change = ifelse(mod != 0, (obs-mod)/mod, 0)) %>% 
    left_join(df_epsi,
              by = "x") %>% 
    mutate(par_change = y_change/epsilon,
           par_change_ub = ifelse((1+par_change)*v_ub_3["hr_S1"]<v_lb_3["hr_S1"], 
                                  NA, (1+par_change)*v_ub_3["hr_S1"]))
  
  # Bind data
  df_changes_3 <- rbind(df_changes_3,
                        df_temp)
  
}

# Percent change for lower or upper bound
n_change_ub <- max(df_changes_3$par_change_ub, na.rm = T)

# Resize bounds
v_lb_4 <- v_lb_3
v_ub_4 <- v_ub_3
v_ub_4["hr_S1"] <- n_change_ub

# Run coverage
set.seed(n_seed)  # set seed
m_params_4 <- sample_prior(n_samp = n_sim,
                           v_lb   = v_lb_4,
                           v_ub   = v_ub_4)
l_coverage_4 <- run_coverage(m_params = m_params_4)
df_coverage_4 <- l_coverage_4$df_model_summ %>% 
  mutate(type_coverage = "C: A + B + improved model coverage by resizing the hazard ratio for death in Sick")

plot_coverage(df_targets  = df_targets,
              df_coverage = rbind(df_coverage_1,
                                  df_coverage_2,
                                  df_coverage_3,
                                  df_coverage_4))

### Readjust the hazard ratio for death in Sick: hr_S1 --------------------
# 1. Compute model-data residuals on previous coverage (mean values)
v_params_4 <- colMeans(m_params_4)      # mean values
modCost4   <- costFun(v_params_4)

# 2. Local sensitivity analysis
SensRes4 <- sensFun(func  = costFun, 
                    parms = v_params_4)

# Plot sensitivity
plot_Sfun(SensRes4)

# 3. Compute the change in the prior distribution bounds based on the 
# sensitivity analysis 
df_changes_4 <- data.frame(NULL)
for(target_i in c("Surv")){ # target_i = "PropSick"
  
  # Filter elasticities
  df_epsi <- SensRes4 %>% 
    filter(var == target_i) %>% 
    select(x,hr_S1) %>% 
    rename(epsilon = hr_S1)
  
  # Change in y --> change in par
  df_temp <- modCost4$residuals %>% 
    filter(name == target_i) %>% 
    select(name,x,obs,mod) %>% 
    mutate(y_change = ifelse(mod != 0, (obs-mod)/mod, 0)) %>% 
    left_join(df_epsi,
              by = "x") %>% 
    mutate(par_change = y_change/epsilon,
           par_change_ub = ifelse((1+par_change)*v_ub_4["hr_S1"]<v_lb_4["hr_S1"],
                                  NA, (1+par_change)*v_ub_4["hr_S1"]))
  
  # Bind data
  df_changes_4 <- rbind(df_changes_4,
                        df_temp)
  
}

# Percent change for lower or upper bound
n_change_ub <- max(df_changes_4$par_change_ub, na.rm = T)

# Resize bounds
v_lb_5 <- v_lb_4
v_ub_5 <- v_ub_4
v_ub_5["hr_S1"] <- n_change_ub

# Run coverage
set.seed(n_seed)  # set seed
m_params_5 <- sample_prior(n_samp = n_sim,
                           v_lb   = v_lb_5,
                           v_ub   = v_ub_5)
l_coverage_5 <- run_coverage(m_params = m_params_5)
df_coverage_5 <- l_coverage_5$df_model_summ %>% 
  mutate(type_coverage = "D: A + B + C + improved model coverage by resizing the hazard ratio for death in Sick")

plot_coverage(df_targets  = df_targets,
              df_coverage = rbind(df_coverage_1,
                                  df_coverage_2,
                                  df_coverage_3,
                                  df_coverage_4,
                                  df_coverage_5))

### Readjust the hazard ratio for death in Sick: hr_S1 --------------------
# 1. Compute model-data residuals on previous coverage (mean values)
v_params_5 <- colMeans(m_params_5)      # mean values
modCost5   <- costFun(v_params_5)

# 2. Local sensitivity analysis
SensRes5 <- sensFun(func  = costFun, 
                    parms = v_params_5)

# Plot sensitivity
plot_Sfun(SensRes5)

# 3. Compute the change in the prior distribution bounds based on the 
# sensitivity analysis 
df_changes_5 <- data.frame(NULL)
for(target_i in c("Surv")){ # target_i = "PropSick"
  
  # Filter elasticities
  df_epsi <- SensRes5 %>% 
    filter(var == target_i) %>% 
    select(x,hr_S1) %>% 
    rename(epsilon = hr_S1)
  
  # Change in y --> change in par
  df_temp <- modCost5$residuals %>% 
    filter(name == target_i) %>% 
    select(name,x,obs,mod) %>% 
    mutate(y_change = ifelse(mod != 0, (obs-mod)/mod, 0)) %>% 
    left_join(df_epsi,
              by = "x") %>% 
    mutate(par_change = y_change/epsilon,
           par_change_ub = ifelse((1+par_change)*v_ub_5["hr_S1"]<v_lb_5["hr_S1"],
                                  NA, (1+par_change)*v_ub_5["hr_S1"]))
  
  # Bind data
  df_changes_5 <- rbind(df_changes_5,
                        df_temp)
  
}

# Percent change for lower or upper bound
n_change_ub <- max(df_changes_5$par_change_ub, na.rm = T)

# Resize bounds
v_lb_6 <- v_lb_5
v_ub_6 <- v_ub_5
v_ub_6["hr_S1"] <- n_change_ub

# Run coverage
set.seed(n_seed)  # set seed
m_params_6 <- sample_prior(n_samp = n_sim,
                           v_lb   = v_lb_6,
                           v_ub   = v_ub_6)
l_coverage_6 <- run_coverage(m_params = m_params_6)
df_coverage_6 <- l_coverage_6$df_model_summ %>% 
  mutate(type_coverage = "E: A + B + C + D + improved model coverage by resizing the hazard ratio for death in Sick")

plot_coverage(df_targets  = df_targets,
              df_coverage = rbind(df_coverage_1,
                                  df_coverage_2,
                                  df_coverage_3,
                                  df_coverage_4,
                                  df_coverage_5,
                                  df_coverage_6))

### Final coverage --------------------------------------------------------
# Range on input search space after resize bounds
v_lb_final <- v_lb_6
v_ub_final <- v_ub_6

set.seed(n_seed)  # set seed
m_params_final <- sample_prior(n_samp = n_sim,
                               v_lb   = v_lb_final,
                               v_ub   = v_ub_final)

l_coverage_final <- run_coverage(m_params = m_params_final)
df_coverage_final <- l_coverage_final$df_model_summ

# Plot final coverage
plot_model_out_vs_targets(df_targets = df_targets,
                          df_model   = df_coverage_final)

## Identifiability (Collinearity) analysis --------------------------------
# Run local sensitivity analysis
SensRes <- sensFun(func  = costFun, 
                   parms = colMeans(m_params_final) # mean values
)

# Combination of variables
m_vars <- matrix(1, ncol = n_outcomes, nrow = n_outcomes,
                 dimnames = list(names(v_outcome_names),
                                 names(v_outcome_names)))

m_vars[lower.tri(m_vars)] <- 0

# Empty data.frame to store results
df_collin <- data.frame(NULL)

# Run collinearity analysis on all combinations of outcomes
for(i in 1:n_outcomes){ # i = 1
  for(j in 1:n_outcomes){ # j = 2
    
    if(m_vars[i,j] != 1) next
    
    if(i == j){
      
      # Get variable names
      var_i <- rownames(m_vars)[i]
      outcome_i <- v_outcome_names[which(names(v_outcome_names) == var_i)] 
      
      # Compute collinearity analysis
      df_temp <- collin(SensRes, which = var_i) %>% 
        mutate(var = var_i,
               Outcome = outcome_i)
      
    }else{
      
      # Get variable names
      var_i <- rownames(m_vars)[i]
      outcome_i <- v_outcome_names[which(names(v_outcome_names) == var_i)] 
      var_j <- colnames(m_vars)[j]
      outcome_j <- v_outcome_names[which(names(v_outcome_names) == var_j)] 
      
      # Compute collinearity analysis
      df_temp <- collin(SensRes, which = c(var_i,var_j)) %>% 
        mutate(var = paste(var_i,var_j, sep = "_"),
               Outcome = paste(outcome_i,outcome_j, sep = " & ")) 
    }
    
    # Bind data.frames
    df_collin <- rbind(df_collin, df_temp)
  }
}

# All variables
df_collin <- rbind(df_collin, 
                   collin(SensRes) %>% 
                     mutate(var = "All",
                            Outcome = paste0(v_outcome_names,collapse = " & ")))

# Plot collinearity analysis
plot_CollIndx(df_collin = df_collin, log_y = T)

## Calibrate with PRE-CISE ---------------------------------------------
# Range on input search space
v_lb <- v_lb_final # lower bound
v_ub <- v_ub_final # upper bound

# Specify seed (for reproducible sequence of random numbers)
set.seed(n_seed)

# number of random samples
n_resamp <- 10000

# record start time of calibration
t_init <- Sys.time()

# Bayesian calibration using the IMIS algorithm (Raftery & Bao, 2010)
fit_imis <- IMIS(B        = 1000,      # the incremental sample size at each iteration of IMIS
                 B.re     = n_resamp,  # the desired posterior sample size
                 number_k = 30,        # the maximum number of iterations in IMIS
                 D        = 0,
                 v_lb     = v_lb,
                 v_ub     = v_ub) 

# Calculate computation time
comp_time <- Sys.time() - t_init
comp_time

# Obtain draws from posterior
m_calib_res <- fit_imis$resample

# Calculate log-likelihood (overall fit) and posterior probability of each sample
m_calib_res <- cbind(m_calib_res, 
                     "Overall_fit"    = log_lik(m_calib_res[,v_param_names]),
                     "Posterior_prob" = posterior(m_calib_res[,v_param_names],
                                                  v_lb = v_lb,
                                                  v_ub = v_ub))

# normalize posterior probability
m_calib_res[,"Posterior_prob"] <- m_calib_res[,"Posterior_prob"]/
  sum(m_calib_res[,"Posterior_prob"])

### Exploring best-fitting input sets -------------------------------------
# Plot the 10000 draws from the posterior
v_post_color <- scales::rescale(m_calib_res[,"Posterior_prob"])
s3d <- scatterplot3d(x = m_calib_res[, 1],
                     y = m_calib_res[, 2],
                     z = m_calib_res[, 3],
                     color = scales::alpha("black", v_post_color),
                     xlim = c(v_lb[1],v_ub[1]), 
                     ylim = c(v_lb[2],v_ub[2]), 
                     zlim = c(v_lb[3],v_ub[3]),
                     xlab = v_param_names[1], 
                     ylab = v_param_names[2], 
                     zlab = v_param_names[3])
# add center of Gaussian components
s3d$points3d(fit_imis$center, col = "red", pch = 8)

# Plot the 1000 draws from the posterior with marginal histograms
pairs.panels(m_calib_res[,v_param_names])

# Compute posterior mean
v_calib_post_mean <- colMeans(m_calib_res[,v_param_names])
v_calib_post_mean

# Compute posterior median and 95% credible interval
m_calib_res_95cr <- colQuantiles(m_calib_res[,v_param_names], 
                                 probs = c(0.025, 0.5, 0.975))
m_calib_res_95cr

### Plot maximum-a-posteriori (MAP) ---------------------------------------
# Compute MAP parameter set
v_calib_map <- m_calib_res[which.max(m_calib_res[,"Posterior_prob"]),]

# Run model with MAP
df_out_best <- sick_sicker_cstm(v_calib_map[v_param_names])

# TARGET 1: Survival ("Surv")
plotrix::plotCI(x = l_targets$Surv$time, y = l_targets$Surv$value, 
                ui = l_targets$Surv$ub,
                li = l_targets$Surv$lb,
                ylim = c(0, 1), 
                xlab = "Time", ylab = "Pr Survive")
points(x = df_out_best$time, 
       y = df_out_best$Surv, 
       pch = 8, col = "red")
legend("topright", 
       legend = c("Target", "Model-predicted output"),
       col = c("black", "red"), pch = c(1, 8))

# TARGET 2: "Prev"
plotrix::plotCI(x = l_targets$Prev$time, y = l_targets$Prev$value,
                ui = l_targets$Prev$ub,
                li = l_targets$Prev$lb,
                ylim = c(0, 1),
                xlab = "Time", ylab = "Prev")
points(x = df_out_best$time,
       y = df_out_best$Prev,
       pch = 8, col = "red")
legend("topright",
       legend = c("Target", "Model-predicted output"),
       col = c("black", "red"), pch = c(1, 8))

# TARGET 3: "PropSick"
plotrix::plotCI(x = l_targets$PropSick$time, y = l_targets$PropSick$value,
                ui = l_targets$PropSick$ub,
                li = l_targets$PropSick$lb,
                ylim = c(0, 1),
                xlab = "Time", ylab = "PropSick")
points(x = df_out_best$time,
       y = df_out_best$PropSick,
       pch = 8, col = "red")
legend("topright",
       legend = c("Target", "Model-predicted output"),
       col = c("black", "red"), pch = c(1, 8))

### Propagate calibrated parameter uncertainty ----------------------------
# Compute IMIS posterior predicted outputs
m_out_surv     <- matrix(NA, nrow = n_resamp, ncol = length(l_targets$Surv$value))  
m_out_prev     <- matrix(NA, nrow = n_resamp, ncol = length(l_targets$Prev$value))  
m_out_propsick <- matrix(NA, nrow = n_resamp, ncol = length(l_targets$PropSick$value))   

# Run model for each posterior parameter set
for(i in 1:n_resamp){ # i = 1
  df_model_res_temp <- sick_sicker_cstm(m_calib_res[i, ])
  m_out_surv[i, ] <- df_model_res_temp$Surv
  m_out_prev[i, ] <- df_model_res_temp$Prev
  m_out_propsick[i, ] <- df_model_res_temp$PropSick[df_model_res_temp$time %in% v_times]
  if(i/100==round(i/100,0)) { 
    cat('\r',paste(i/n_resamp*100,"% done",sep=""))
  }
}

# Posterior predicted mean
m_out_surv_postmean <- colMeans(m_out_surv)
m_out_prev_postmean <- colMeans(m_out_prev)
m_out_propsick_postmean <- colMeans(m_out_propsick)

# Implement PRE-CISE (aggressive approach) --------------------------------
# 0. Define input search space and plot first coverage
v_lb_1 <- c(p_S1S2 = 0.00, hr_S1 = 1.00, hr_S2 = 1.00) # lower bound
v_ub_1 <- c(p_S1S2 = 0.50, hr_S1 = 8.00, hr_S2 = 25.0) # upper bound

# Run initial coverage analysis
set.seed(n_seed)   # set seed
m_params_1 <- sample_prior(n_samp = n_sim,
                           v_lb   = v_lb_1,
                           v_ub   = v_ub_1)

l_coverage_1  <- run_coverage(m_params = m_params_1)
df_coverage_1 <- l_coverage_1$df_model_summ %>% 
  mutate(type_coverage = "Initial model coverage")

plot_coverage(df_targets  = df_targets,
              df_coverage = df_coverage_1)

## Resize the prior distribution bounds -----------------------------------
### Adjust the probability to become Sicker when Sick: p_S1S2 -------------
# 1. Compute the effect of the parameters on the variables of interest
# Compute model costs on initial coverage (mean values)
v_params_1 <- colMeans(m_params_1)
modCost1 <- costFun(v_params_1)

# Plot residuals
# plot_costFun(modCost1, wgt_res = T) # weighted residuals
# plot_costFun(modCost1, wgt_res = F) # unweighted residuals

# 2. Local sensitivity analysis
SensRes1 <- sensFun(func  = costFun, 
                    parms = v_params_1)

# Plot sensitivity*
plot_Sfun(SensRes1)

# Rank parameters according to their importance
m_Sens1 <- summary(SensRes1)

# value: value of the parameter 
# Mean: mean sensitivity
# L1: L1 norm
# L2: L2 norm
df_Sens1 <- as.data.frame(m_Sens1)

# Order by mean sensitivity
df_Sens1[order(abs(df_Sens1$L2), decreasing = T),] 

# 3. Resize the prior distribution bounds
# Start with the one the outcomes are most sensible to: p_S1S2
df_changes_1 <- data.frame(NULL)
for(target_i in c("PropSick","Prev")){ # target_i = "PropSick"
  
  # Filter elasticities
  df_epsi <- SensRes1 %>% 
    filter(var == target_i) %>% 
    select(x,p_S1S2) %>% 
    rename(epsilon = p_S1S2)
  
  # Change in y --> change in par
  df_temp <- modCost1$residuals %>% 
    filter(name == target_i) %>% 
    select(name,x,obs,mod) %>% 
    mutate(y_change = ifelse(mod != 0, (obs-mod)/mod, 0)) %>% 
    left_join(df_epsi,
              by = "x") %>% 
    mutate(par_change = y_change/epsilon,
           par_change_ub = ifelse((1+par_change)*v_ub_1["p_S1S2"] < v_lb_1["p_S1S2"],
                                  NA, (1+par_change)*v_ub_1["p_S1S2"]))
  
  # Bind data
  df_changes_1 <- rbind(df_changes_1,
                        df_temp)
  
}

# Percent change for lower or upper bound
n_change_ub <- max(df_changes_1$par_change_ub, na.rm = T)

# Resize bounds
v_lb_2 <- v_lb_1
v_ub_2 <- v_ub_1
v_ub_2["p_S1S2"] <- n_change_ub

# Run coverage
set.seed(n_seed)  # set seed
m_params_2 <- sample_prior(n_samp = n_sim,
                           v_lb = v_lb_2,
                           v_ub = v_ub_2)
l_coverage_2 <- run_coverage(m_params = m_params_2)
df_coverage_2 <- l_coverage_2$df_model_summ %>% 
  mutate(type_coverage = "A: Improved model coverage by resizing the probability to become Sicker when Sick")

plot_coverage(df_targets  = df_targets,
              df_coverage = rbind(df_coverage_1,df_coverage_2))

### Adjust the hazard ratio for death in Sicker: hr_S2 --------------------
# 1. Compute model-data residuals on previous coverage (mean values)
v_params_2 <- colMeans(m_params_2)
modCost2   <- costFun(v_params_2)

# 2. Local sensitivity analysis
SensRes2 <- sensFun(func  = costFun, 
                    parms = v_params_2)

# Plot sensitivity
# plot_Sfun(SensRes2)

# Rank parameters according to their importance
m_Sens2 <- summary(SensRes2)

# value: value of the parameter 
# Mean: mean sensitivity
# L1: L1 norm
# L2: L2 norm
df_Sens2 <- as.data.frame(m_Sens2)

# Order by L2
df_Sens2[order(df_Sens2$L2, decreasing = T),] 

# 3. Compute the change in the prior distribution bounds based on the 
# sensitivity analysis
df_changes_2 <- data.frame(NULL)
for(target_i in c("Surv","Prev")){ # target_i = "PropSick"
  
  # Filter elasticities
  df_epsi <- SensRes2 %>% 
    filter(var == target_i) %>% 
    select(x,hr_S2) %>% 
    rename(epsilon = hr_S2)
  
  # Change in y --> change in par
  df_temp <- modCost2$residuals %>% 
    filter(name == target_i) %>% 
    select(name,x,obs,mod) %>%  
    mutate(y_change = ifelse(mod != 0, (obs-mod)/mod, 0)) %>% 
    left_join(df_epsi,
              by = "x") %>% 
    mutate(par_change = y_change/epsilon,
           par_change_ub = ifelse((1+par_change)*v_ub_2["hr_S2"] < v_lb_2["hr_S2"],
                                  NA, (1+par_change)*v_ub_2["hr_S2"]))
  
  # Bind data
  df_changes_2 <- rbind(df_changes_2,
                        df_temp)
}

# Percent change for lower or upper bound
n_change_ub <- min(df_changes_2$par_change_ub, na.rm = T)

# Resize bounds
v_lb_3 <- v_lb_2
v_ub_3 <- v_ub_2
v_ub_3["hr_S2"] <- n_change_ub

# Run coverage
set.seed(n_seed)  # set seed
m_params_3 <- sample_prior(n_samp = n_sim,
                           v_lb   = v_lb_3,
                           v_ub   = v_ub_3)
l_coverage_3 <- run_coverage(m_params = m_params_3)
df_coverage_3 <- l_coverage_3$df_model_summ %>% 
  mutate(type_coverage = "B: A + improved model coverage by resizing the hazard ratio for death in Sicker")

plot_coverage(df_targets  = df_targets,
              df_coverage = rbind(df_coverage_1,df_coverage_2,df_coverage_3))

### Adjust the hazard ratio for death in Sick: hr_S1 ----------------------
# 1. Compute model-data residuals on previous coverage (mean values)
v_params_3 <- colMeans(m_params_3)      # mean values
modCost3   <- costFun(v_params_3)

# 2. Local sensitivity analysis
SensRes3 <- sensFun(func  = costFun, 
                    parms = v_params_3)

# Plot sensitivity
# plot_Sfun(SensRes3)

# 3. Compute the change in the prior distribution bounds based on the 
# sensitivity analysis 
df_changes_3 <- data.frame(NULL)
for(target_i in c("Surv")){ # target_i = "PropSick"
  
  # Filter elasticities
  df_epsi <- SensRes3 %>% 
    filter(var == target_i) %>% 
    select(x,hr_S1) %>% 
    rename(epsilon = hr_S1)
  
  # Change in y --> change in par
  df_temp <- modCost3$residuals %>% 
    filter(name == target_i) %>% 
    select(name,x,obs,mod) %>% 
    mutate(y_change = ifelse(mod != 0, (obs-mod)/mod, 0)) %>% 
    left_join(df_epsi,
              by = "x") %>% 
    mutate(par_change = y_change/epsilon,
           par_change_ub = ifelse((1+par_change)*v_ub_3["hr_S1"] < v_lb_3["hr_S1"],
                                  NA, (1+par_change)*v_ub_3["hr_S1"]))
  
  # Bind data
  df_changes_3 <- rbind(df_changes_3,
                        df_temp)
  
}

# Percent change for lower or upper bound, survival t = 10
n_change_ub <- df_changes_3$par_change_ub[df_changes_3$x == 10]

# Resize bounds
v_lb_4 <- v_lb_3
v_ub_4 <- v_ub_3
v_ub_4["hr_S1"] <- n_change_ub

# Run coverage
set.seed(n_seed)  # set seed
m_params_4 <- sample_prior(n_samp = n_sim,
                           v_lb   = v_lb_4,
                           v_ub   = v_ub_4)
l_coverage_4 <- run_coverage(m_params = m_params_4)
df_coverage_4 <- l_coverage_4$df_model_summ %>% 
  mutate(type_coverage = "C: A + B + improved model coverage by resizing the hazard ratio for death in Sick")

l_plot_coverage <- plot_coverage(df_targets  = df_targets,
                                 df_coverage = rbind(df_coverage_1,
                                                     df_coverage_2,
                                                     df_coverage_3,
                                                     df_coverage_4))

### Final coverage --------------------------------------------------------
# Range on input search space after resize bounds
v_lb_final <- v_lb_4
v_ub_final <- v_ub_4

set.seed(n_seed)   # set seed
m_params_final <- sample_prior(n_samp = n_sim,
                               v_lb   = v_lb_final,
                               v_ub   = v_ub_final)

l_coverage_final <- run_coverage(m_params = m_params_final)
df_coverage_final <- l_coverage_final$df_model_summ

plot_model_out_vs_targets(df_targets = df_targets,
                          df_model   = df_coverage_final)

# Implement PRE-CISE (conservative approach) ------------------------------
# 0. Define input search space and plot first coverage
v_lb_1 <- c(p_S1S2 = 0.00, hr_S1 = 1.00, hr_S2 = 1.00) # lower bound
v_ub_1 <- c(p_S1S2 = 0.50, hr_S1 = 8.00, hr_S2 = 25.0) # upper bound

# Run initial coverage analysis
set.seed(n_seed)   # set seed
m_params_1 <- sample_prior(n_samp = n_sim,
                           v_lb   = v_lb_1,
                           v_ub   = v_ub_1)

l_coverage_1  <- run_coverage(m_params = m_params_1)
df_coverage_1 <- l_coverage_1$df_model_summ %>% 
  mutate(type_coverage = "Initial model coverage")

l_plot_coverage_1 <- plot_coverage(df_targets  = df_targets,
                                   df_coverage = df_coverage_1)
l_plot_coverage_1

## Resize the prior distribution bounds -----------------------------------
### Adjust the probability to become Sicker when Sick: p_S1S2 -------------
# 1. Compute the effect of the parameters on the variables of interest
# Compute model costs on initial coverage (mean values)
v_params_1 <- colMeans(m_params_1)
modCost1 <- costFun(v_params_1)

# Plot residuals
plot_costFun(modCost1, wgt_res = T) # weighted residuals
plot_costFun(modCost1, wgt_res = F) # unweighted residuals

# 2. Local sensitivity analysis
SensRes1 <- sensFun(func  = costFun, 
                    parms = v_params_1)

# Plot sensitivity*
plot_Sfun(SensRes1)

# Rank parameters according to their importance
m_Sens1 <- summary(SensRes1)

# value: value of the parameter 
# Mean: mean sensitivity
# L1: L1 norm
# L2: L2 norm
df_Sens1 <- as.data.frame(m_Sens1)

# Order by mean sensitivity
df_Sens1[order(abs(df_Sens1$L2), decreasing = T),] 

# 3. Resize the prior distribution bounds
df_changes_1 <- data.frame(NULL)
for(target_i in c("Surv","PropSick","Prev")){ # target_i = "PropSick"
  
  # Filter elasticities
  df_epsi <- SensRes1 %>% 
    filter(var == target_i) %>% 
    select(x,p_S1S2) %>% 
    rename(epsilon = p_S1S2)
  
  # Change in y --> change in par
  df_temp <- modCost1$residuals %>% 
    filter(name == target_i) %>% 
    select(name,x,obs,mod) %>% 
    mutate(y_change = ifelse(mod != 0, (obs-mod)/mod, 0)) %>% 
    left_join(df_epsi,
              by = "x") %>% 
    mutate(par_change = y_change/epsilon,
           par_change_ub = ifelse((1+par_change)*v_ub_1["p_S1S2"]<v_lb_1["p_S1S2"] |
                                    (1+par_change)*v_ub_1["p_S1S2"]>1, 
                                  NA, (1+par_change)*v_ub_1["p_S1S2"]))
  
  # Bind data
  df_changes_1 <- rbind(df_changes_1,
                        df_temp)
  
}

# Percent change for lower or upper bound
n_change_ub <- max(df_changes_1$par_change_ub, na.rm = T)

# Resize bounds
v_lb_2 <- v_lb_1
v_ub_2 <- v_ub_1
v_ub_2["p_S1S2"] <- n_change_ub

# Run coverage
set.seed(n_seed)  # set seed
m_params_2 <- sample_prior(n_samp = n_sim,
                           v_lb = v_lb_2,
                           v_ub = v_ub_2)
l_coverage_2 <- run_coverage(m_params = m_params_2)
df_coverage_2 <- l_coverage_2$df_model_summ %>% 
  mutate(type_coverage = "A: Improved model coverage by resizing the probability to become Sicker when Sick")

plot_coverage(df_targets  = df_targets,
              df_coverage = rbind(df_coverage_1,df_coverage_2))

### Adjust the hazard ratio for death in Sicker: hr_S2 --------------------
# 1. Compute model-data residuals on previous coverage (mean values)
v_params_2 <- colMeans(m_params_2)
modCost2   <- costFun(v_params_2)

# 2. Local sensitivity analysis
SensRes2 <- sensFun(func  = costFun, 
                    parms = v_params_2)

# Plot sensitivity
plot_Sfun(SensRes2)

# Rank parameters according to their importance
m_Sens2 <- summary(SensRes2)

# value: value of the parameter 
# Mean: mean sensitivity
# L1: L1 norm
# L2: L2 norm
df_Sens2 <- as.data.frame(m_Sens2)

# Order by L2
df_Sens2[order(df_Sens2$L2, decreasing = T),] 

# 3. Compute the change in the prior distribution bounds based on the 
# sensitivity analysis
df_changes_2 <- data.frame(NULL)
for(target_i in c("Surv","Prev")){ # target_i = "PropSick"
  
  # Filter elasticities
  df_epsi <- SensRes2 %>% 
    filter(var == target_i) %>% 
    select(x,hr_S2) %>% 
    rename(epsilon = hr_S2)
  
  # Change in y --> change in par
  df_temp <- modCost2$residuals %>% 
    filter(name == target_i) %>% 
    select(name,x,obs,mod) %>%  
    mutate(y_change = ifelse(mod != 0, (obs-mod)/mod, 0)) %>% 
    left_join(df_epsi,
              by = "x") %>% 
    mutate(par_change = y_change/epsilon,
           par_change_ub = ifelse((1+par_change)*v_ub_2["hr_S2"]<v_lb_2["hr_S2"],
                                  NA,(1+par_change)*v_ub_2["hr_S2"]))
  
  # Bind data
  df_changes_2 <- rbind(df_changes_2,
                        df_temp)
}

# Percent change for lower or upper bound
n_change_ub <- mean(df_changes_2$par_change_ub, na.rm = T)

# Resize bounds
v_lb_3 <- v_lb_2
v_ub_3 <- v_ub_2
v_ub_3["hr_S2"] <- n_change_ub

# Run coverage
set.seed(n_seed)  # set seed
m_params_3 <- sample_prior(n_samp = n_sim,
                           v_lb   = v_lb_3,
                           v_ub   = v_ub_3)
l_coverage_3 <- run_coverage(m_params = m_params_3)
df_coverage_3 <- l_coverage_3$df_model_summ %>% 
  mutate(type_coverage = "B: A + improved model coverage by resizing the hazard ratio for death in Sicker")

l_plot_coverage <- plot_coverage(df_targets  = df_targets,
                                 df_coverage = rbind(df_coverage_1,
                                                     df_coverage_2,
                                                     df_coverage_3))

### Final coverage --------------------------------------------------------
# Range on input search space after resize bounds
v_lb_final <- v_lb_3
v_ub_final <- v_ub_3

set.seed(n_seed)   # set seed
m_params_final <- sample_prior(n_samp = n_sim,
                               v_lb   = v_lb_final,
                               v_ub   = v_ub_final)

l_coverage_final <- run_coverage(m_params = m_params_final)
df_coverage_final <- l_coverage_final$df_model_summ

plot_model_out_vs_targets(df_targets = df_targets,
                          df_model   = df_coverage_final)
