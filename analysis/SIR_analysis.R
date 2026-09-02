#------------------------------------------------------------------------------#
# This script illustrates the implementation of PRE-CISE in a deterministic    #
# Susceptible-Infectious-Recovered (SIR) transmission model calibrated to      #
# daily incident cases using the IMIS algorithm.                               #
#                                                                              #
# Authors:                                                                     #
#     Fernando Alarid-Escudero, PhD, <falarid@stanford.edu>                    #
#     Jeremy D. Goldhaber-Fiebert, PhD, <jeremygf@stanford.edu>                #
#     Valeria Gracia Olvera, MSc, <vgracia@stanford.edu>                       #
#                                                                              #
# Last update: 08/31/26                                                        #
#------------------------------------------------------------------------------#

rm(list = ls()) # clean environment

# Calibration specifications ----------------------------------------------
# Model: Susceptible-Infected-Recovered (SIR) model without demography
# Inputs to be calibrated: 
#   r_beta   - transmission rate
#   r_gamma  - recovery rate
# Targets: 
#   DXIncTot - total incident cases
# Search method: Random search using Latin-Hypercube Sampling
# Goodness-of-fit measure: Sum of Log-Likelihood

# Load libraries and functions --------------------------------------------
# calibration functionality
library(lhs)
library(matrixStats) # package used for summary statistics
library(FME)         # For identifiability and collinearity analysis
library(tidyverse)
options(dplyr.summarise.inform = FALSE) # do not show summarise info

# visualization
library(psych)
library(ggplot2)
library(ggthemes)
library(ggpubr)

# global sensitivity
library(sensitivity)

# calibration
library(mvtnorm)
source("R/IMIS_function_SIR.R") # from IMIS package, adapted for this analysis
# Old version 
# devtools::install_version("IMIS", version = "0.1", repos = "http://cran.us.r-project.org")

# load functions
source("R/helper_functions.R") # from dampack package
source("R/SIR_functions.R")

# Target data -------------------------------------------------------------
# Load target data: influenza outbreak at a British boarding school in 1978 
# (Anon 1978; El Attouga & El Khalifi, 2024 ; Keeling & Rohani, 2008) 
# Source:
# https://www.sciencedirect.com/science/article/pii/S0960077924005939#appendix
# https://www.bmj.com/content/bmj/1/6112/586.full.pdf
df_InfluenzaData_raw <- readxl::read_xlsx("data/Influenza1978Data.xlsx")
df_InfluenzaData_raw$Date <- as.Date(df_InfluenzaData_raw$Date)

df_targets <- df_InfluenzaData_raw %>% 
  mutate(Type = "Target",
         abbrev_outcome = "DXIncTot",
         Day = as.numeric(Date - Date[1] + 1),
         Date0 = Date - Date[1],
         lb = epitools::pois.exact(x = value, pt = population)$lower*population,
         ub = epitools::pois.exact(x = value, pt = population)$upper*population,
         se = (ub-lb)/(1.96*2)) %>% 
  select(Type,abbrev_outcome,Outcome,population,Day,Date,Date0,value,lb,ub,se) %>% 
  as.data.frame()

# Plot the targets
plot_targets(df_targets)

# Variables ---------------------------------------------------------------
# Calibrated dates
v_dates <- df_targets$Date

# First day of calibration
n_date_init <- first(v_dates)

# Last day of calibration
n_date_end <- last(v_dates)

# Number of days until infection become extinguished
n_t <- as.numeric(n_date_end - n_date_init)

# Number of model evaluations per day
time_step <- 1

# Total number of the boys in school
n_tot_pop <- unique(df_targets$population)

# Initial number of infected boys
n_inf_init <- df_targets$value[df_targets$Date == n_date_init]

# Number of simulations (coverage analysis)
n_sim <- 300

# Seed number
n_seed <- 111124

# Set seed
set.seed(n_seed)

# Number of calibration targets
v_outcome_names <- c("Number of infections")
names(v_outcome_names) <- c("DXIncTot")
n_outcomes <- length(v_outcome_names)

# Specify calibration parameters ------------------------------------------
# Names and number of input parameters to be calibrated
v_param_names <- c("r_beta",   # transmission rate
                   "r_gamma")  # recovery rate
n_param <- length(v_param_names)

# Range on input search space
v_lb <- c(r_beta = 0.50, r_gamma = 0.10)  # lower bound
v_ub <- c(r_beta = 3.00, r_gamma = 1.00)  # upper bound

# Model -------------------------------------------------------------------
# -inputs are parameters to be estimated through calibration
# -outputs correspond to the target data

# Set of parameter test
v_params_test <- c(r_beta = 1.66, r_gamma = 1/2.2)

# List of parameters
l_params_all <- load_params_all(n_t         = n_t,
                                time_step   = time_step,
                                n_date_init = n_date_init,
                                n_tot_pop   = n_tot_pop, 
                                n_inf_init  = n_inf_init,
                                r_beta      = v_params_test["r_beta"],
                                r_gamma     = v_params_test["r_gamma"])

# Run model, check it works
l_out_model <- SIR_model(l_params_all = l_params_all)
l_out_model$df_out_model

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
log_lik(v_params     = v_params_test, 
        l_params_all = l_params_all, 
        v_dates      = v_dates)
log_lik(v_params     = sample_prior(10, v_lb = v_lb, v_ub = v_ub),
        l_params_all = l_params_all, 
        v_dates      = v_dates)

likelihood(v_params     = v_params_test, 
           l_params_all = l_params_all, 
           v_dates      = v_dates)
likelihood(v_params     = sample_prior(10, v_lb = v_lb, v_ub = v_ub),
           l_params_all = l_params_all, 
           v_dates      = v_dates)

# Posterior
log_post(v_params     = v_params_test, 
         v_lb         = v_lb, 
         v_ub         = v_ub, 
         l_params_all = l_params_all, 
         v_dates      = v_dates)
log_post(v_params     = sample_prior(10, v_lb = v_lb, v_ub = v_ub),
         v_lb         = v_lb, 
         v_ub         = v_ub, 
         l_params_all = l_params_all, 
         v_dates      = v_dates)

posterior(v_params     = v_params_test,
          v_lb         = v_lb, 
          v_ub         = v_ub, 
          l_params_all = l_params_all,
          v_dates      = v_dates)
posterior(v_params     = sample_prior(10, v_lb = v_lb, v_ub = v_ub),
          v_lb         = v_lb, 
          v_ub         = v_ub,
          l_params_all = l_params_all, 
          v_dates      = v_dates)

# Define functions --------------------------------------------------------
## Coverage analysis ------------------------------------------------------
# Coverage function
# l_params_all: List containing all the parameters of decision model.
# m_params: Matrix containing sample of parameter sets
run_coverage <- function(l_params_all, m_params){ 
  
  # Number of simulations
  n_sim <- nrow(m_params)
  
  # Run model
  df_model <- data.frame(NULL) # empty data.frame to store results
  for(iter_i in 1:n_sim){ # iter_i = 1
    
    # Update list of parameters
    l_params_all$r_beta  <- m_params[iter_i, "r_beta"]
    l_params_all$r_gamma <- m_params[iter_i, "r_gamma"]
    
    # Run model
    l_out_model <- SIR_model(l_params_all = l_params_all)
    
    # Output of the model in a data.frame
    df_temp <- calc_DXIncTot(l_out_model  = l_out_model,
                             l_params_all = l_params_all)
    
    # Gather results
    df_model <- rbind(df_model, df_temp) 
    
    # Print progress
    if(iter_i/100==round(iter_i/100,0)) {
      cat('\r',paste(round(iter_i/n_sim*100,0),"% simulations done",sep=""))
    }
  }
  
  # Summarise
  df_model_summ <- df_model %>%
    group_by(Type, Outcome, Day, Date, Date0) %>%
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
# Cost function for one-at-a-time sensitivity analysis
# v_params:     Vector containing the model parameters.
# l_params_all: List containing all the parameters of decision model.
# v_dates:      Vector containing dates at which cost function is to be computed.
costFun <- function(v_params, l_params_all, v_dates) {
  
  # Update list of parameters
  l_params_all$r_beta  <- v_params["r_beta"]
  l_params_all$r_gamma <- v_params["r_gamma"]
  
  # Run model
  l_out_model <- SIR_model(l_params_all = l_params_all)
  
  # Output of the model in a data.frame
  df_DXIncTot <- calc_DXIncTot(l_out_model  = l_out_model, 
                               l_params_all = l_params_all)
  
  # Filter and select columns
  ## Number of infections
  df_model_DXIncTot <- df_DXIncTot %>% 
    filter(abbrev_outcome == "DXIncTot" & Date %in% v_dates) %>% 
    mutate(time = Day) %>%
    rename(DXIncTot = value) %>% 
    select(time, DXIncTot)
  
  df_targets_DXIncTot <- df_targets %>% 
    filter(abbrev_outcome == "DXIncTot" & Date %in% v_dates) %>% 
    mutate(time = Day) %>% 
    select(time, value, se) %>% 
    rename(DXIncTot = value)
  
  # Compute "costs" for each target
  cost <- modCost(model = df_model_DXIncTot, 
                  obs   = df_targets_DXIncTot, 
                  err   = "se")
  return(cost)
}

# Function to compute model outputs for global sensitivity analyses
# m_params:     Matrix containing sample of parameter sets.
# l_params_all: List containing all the parameters of decision model.
# n_Date:       Date at which GSA is to be computed.
model_out_GSA <- function(m_params, l_params_all, n_Date){
  
  v_res <- c()
  for(i in 1:nrow(m_params)){
    
    # Update list of parameters with calibrated parameter set
    l_params_all$r_beta  <- m_params[i,"r_beta"]
    l_params_all$r_gamma <- m_params[i,"r_gamma"]
    
    # Run model
    l_out_model  <- SIR_model(l_params_all)
    
    # Calculate total infections
    df_DXIncTot <- calc_DXIncTot(l_out_model  = l_out_model, 
                                 l_params_all = l_params_all) %>% 
      filter(Date == n_Date)
    
    # Bind results
    v_res <- c(v_res, df_DXIncTot$value)
    
  }
  
  # Return vector of values
  return(v_res)
}

# Global sensitivity analysis (GSA): this function runs soboljansen function 
# from sensitivity package.
# l_params_all: List containing all the parameters of decision model.
# v_dates:      Vector containing dates at which cost function is to be computed.
# X1:           First random sample (prior distribution).
# X2:           Second random sample (prior distribution).
# n_boot:       Number of bootstrap replicates
# n_conf:       Confidence level for bootstrap confidence intervals.
run_GSA <- function(l_params_all, v_dates, X1, X2, n_boot, n_conf){
  
  df_res_sobol <- data.frame(NULL)
  for(i in 2:length(v_dates)){ # i = 2
    
    cat("Day", i, "\n")
    
    # Compute Sobol Indices
    res_sobol <- soboljansen(model        = model_out_GSA, 
                             X1           = X1, 
                             X2           = X2, 
                             l_params_all = l_params_all, 
                             n_Date       = v_dates[i],
                             nboot        = n_boot, 
                             conf         = n_conf)
    
    # Bind results
    df_res_sobol <- rbind(df_res_sobol,
                          data.frame(par  = rownames(res_sobol$S),
                                     type = "Main effect",
                                     x    = i,
                                     res_sobol$S, 
                                     row.names = NULL),
                          data.frame(par  = rownames(res_sobol$T),
                                     type = "Total effect",
                                     x    = i,
                                     res_sobol$T,
                                     row.names = NULL))
    
  }
  return(df_res_sobol)
}

# Compare sensitivity analysis --------------------------------------------
# Define input search space
v_lb_Sens <- c(r_beta = 0.30, r_gamma = 0.10)  # lower bound
v_ub_Sens <- c(r_beta = 3.50, r_gamma = 1.20)  # upper bound

## Global sensitivity analysis
# Get two sets of priors
set.seed(n_seed)    # set seed
m_params_X1 <- sample_prior(n_samp = n_sim,
                            v_lb   = v_lb_Sens,
                            v_ub   = v_ub_Sens)

m_params_X2 <- sample_prior(n_samp = n_sim,
                            v_lb   = v_lb_Sens,
                            v_ub   = v_ub_Sens)

# Run GSA
df_res_sobol <- run_GSA(l_params_all = l_params_all,
                        v_dates      = v_dates,
                        X1           = m_params_X1,
                        X2           = m_params_X2, 
                        n_boot       = 100, 
                        n_conf       = 0.95)

# Bootstrap confidence interval bounds were manually restricted to 
# the plausible range [0,1]
df_res_sobol <- df_res_sobol %>% 
  mutate(min..c.i. = ifelse(min..c.i.<0, 0, min..c.i.),
         max..c.i. = ifelse(max..c.i.>1, 1, max..c.i.))

# Plot Sobol' indices
plot_GSA(df_sobol = df_res_sobol)

## Local sensitivity
# Compute model-data residuals on previous coverage (mean values) 
v_params <- colMeans(m_params_X1)        # mean values

# Sensitivity analysis
SensRes <- sensFun(func         = costFun,
                   parms        = v_params,
                   l_params_all = l_params_all,
                   v_dates = v_dates)

plot_Sfun(SensRes)

# Collinearity analysis on accumulated targets ----------------------------
# Define input search space
v_lb_0 <- c(r_beta = 0.40, r_gamma = 0.10)  # lower bound
v_ub_0 <- c(r_beta = 3.00, r_gamma = 1.00)  # upper bound

# Run initial coverage analysis
set.seed(n_seed)  # set seed 
m_params_1 <- sample_prior(n_samp = 1000,
                           v_lb   = v_lb_0,
                           v_ub   = v_ub_0)

# Dates index
v_dates_index  <- 1:(n_t+1)

# Empty data.frames to store results
df_collin_iter <- data.frame(NULL)
for(j in seq_along(v_dates_index)[-(n_t+1)]){ # j = 2
  
  # Select dates
  v_index_samp <- sort(v_dates_index[1:(1+j)]) 
  v_dates_samp <- v_dates[v_index_samp]
  
  # 1. Compute model-data residuals on previous coverage (mean values) 
  v_params <- colMeans(m_params_1)
  modCost   <- costFun(v_params     = v_params, 
                       l_params_all = l_params_all,
                       v_dates      = v_dates_samp)
  
  # 2. Local sensitivity analysis
  SensRes <- sensFun(func         = costFun,
                     parms        = v_params,
                     l_params_all = l_params_all,
                     v_dates      = v_dates_samp)
  
  m_Sens <- summary(SensRes)
  
  # 3. Collinearity analysis
  df_collin <- collin(SensRes) %>% 
    mutate(n_targets = j+1,
           var     = "All",
           Outcome = paste0(v_outcome_names,collapse = " & ")) %>% 
    relocate(n_targets, var)
  
  # Bind data.frames
  df_collin_iter <- rbind(df_collin_iter, df_collin)
  
}

# Implement PRE-CISE (wide) -----------------------------------------------
# 0. Define input search space and plot first coverage
v_lb_1 <- c(r_beta = 0.30, r_gamma = 0.10)  # lower bound
v_ub_1 <- c(r_beta = 3.50, r_gamma = 1.20)  # upper bound

# Run initial coverage analysis
set.seed(n_seed)  # set seed 
m_params_1 <- sample_prior(n_samp = n_sim,
                           v_lb   = v_lb_1,
                           v_ub   = v_ub_1)

l_coverage_1  <- run_coverage(l_params_all = l_params_all, m_params = m_params_1)
df_coverage_1 <- l_coverage_1$df_model_summ %>%
  mutate(type_coverage = "Initial model coverage")

l_plot_coverage_1 <- plot_coverage(df_targets  = df_targets,
                                   df_coverage = df_coverage_1)

l_plot_coverage_1

## Resize the prior distribution bounds -----------------------------------
### Adjust transmission rate: r_beta --------------------------------------
# 1. Compute model-data residuals on initial coverage (mean values) 
v_params_1 <- colMeans(m_params_1)
modCost1   <- costFun(v_params     = v_params_1, 
                      l_params_all = l_params_all,
                      v_dates      = v_dates)

# Plot residuals
# plot_costFun(Cfun = modCost1, wgt_res = T) # weighted residuals
# plot_costFun(modCost1, wgt_res = F)        # unweighted residuals

# 2. Local sensitivity analysis
SensRes1 <- sensFun(func         = costFun,
                    parms        = v_params_1,
                    l_params_all = l_params_all,
                    v_dates      = v_dates)

# Plot sensitivity
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
for(target_i in c("DXIncTot")){ # target_i = "DXIncTot"
  
  # Filter elasticities
  df_epsi <- SensRes1 %>% 
    filter(var == target_i) %>% 
    select(x,r_beta) %>% 
    rename(epsilon = r_beta)
  
  # Change in y --> change in par
  df_temp <- modCost1$residuals %>% 
    filter(name == target_i) %>% 
    select(name,x,obs,mod) %>% 
    mutate(y_change = ifelse(mod != 0, (obs-mod)/mod, 0)) %>% 
    left_join(df_epsi,
              by = "x") %>% 
    mutate(par_change = ifelse(epsilon == 0, NA, y_change/epsilon),
           par_change_lb = ifelse(par_change < 0,
                                  NA, (1+par_change)*v_lb_1["r_beta"]),
           par_change_ub = ifelse(par_change > 0 | 
                                    (1+par_change)*v_ub_1["r_beta"] < 0.3, 
                                  NA, (1+par_change)*v_ub_1["r_beta"]))
  
  # Bind data
  df_changes_1 <- rbind(df_changes_1,
                        df_temp)
  
}

# Change for lower or upper bound
n_change_lb <- mean(df_changes_1$par_change_lb, na.rm = T)
n_change_ub <- mean(df_changes_1$par_change_ub, na.rm = T)

# Resize bounds
v_lb_2 <- v_lb_1
v_ub_2 <- v_ub_1
v_lb_2["r_beta"] <- n_change_lb
v_ub_2["r_beta"] <- n_change_ub

# Run coverage
set.seed(n_seed)  # set seed 
m_params_2 <- sample_prior(n_samp = n_sim,
                           v_lb = v_lb_2,
                           v_ub = v_ub_2)
l_coverage_2 <- run_coverage(l_params_all = l_params_all, 
                             m_params     = m_params_2)
df_coverage_2 <- l_coverage_2$df_model_summ 
df_coverage_2$type_coverage <- 
  "A: Improved model coverage by resizing the transmission rate"

plot_coverage(df_targets  = df_targets,
              df_coverage = rbind(df_coverage_1,df_coverage_2))

### Adjust recovery rate: r_gamma -----------------------------------------
# 1. Compute model-data residuals on previous coverage (mean values)
v_params_2 <- colMeans(m_params_2)        # mean values
modCost2   <- costFun(v_params     = v_params_2, 
                      l_params_all = l_params_all,
                      v_dates      = v_dates)

# 2. Local sensitivity analysis
SensRes2 <- sensFun(func         = costFun,
                    parms        = v_params_2,
                    l_params_all = l_params_all,
                    v_dates      = v_dates)

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
for(target_i in c("DXIncTot")){ # target_i = "DXIncTot"
  
  # Filter elasticities
  df_epsi <- SensRes2 %>% 
    filter(var == target_i) %>% 
    select(x,r_gamma) %>% 
    rename(epsilon = r_gamma)
  
  # Change in y --> change in par
  df_temp <- modCost2$residuals %>% 
    filter(name == target_i) %>% 
    select(name,x,obs,mod) %>%  
    mutate(y_change = ifelse(mod != 0, (obs-mod)/mod, 0)) %>% 
    left_join(df_epsi,
              by = "x") %>% 
    mutate(par_change = ifelse(epsilon == 0, NA, y_change/epsilon),
           par_change_lb = ifelse(par_change < 0,
                                  NA, (1+par_change)*v_lb_2["r_gamma"]),
           par_change_ub = ifelse(par_change > 0 | 
                                    (1+par_change)*v_ub_2["r_gamma"] < 0.3, 
                                  NA, (1+par_change)*v_ub_2["r_gamma"]))
  
  # Bind data
  df_changes_2 <- rbind(df_changes_2,
                        df_temp)
  
}

# Change for lower or upper bound
n_change_lb <- mean(df_changes_2$par_change_lb, na.rm = T) 
n_change_ub <- mean(df_changes_2$par_change_ub, na.rm = T) 

# Resize bounds
v_lb_3 <- v_lb_2
v_ub_3 <- v_ub_2
v_lb_3["r_gamma"] <- n_change_lb
v_ub_3["r_gamma"] <- n_change_ub

# Run coverage
set.seed(n_seed)  # set seed 
m_params_3 <- sample_prior(n_samp = n_sim,
                           v_lb   = v_lb_3,
                           v_ub   = v_ub_3)
l_coverage_3 <- run_coverage(l_params_all = l_params_all, 
                             m_params     = m_params_3)
df_coverage_3 <- l_coverage_3$df_model_summ
df_coverage_3$type_coverage <- 
  "B: A + improved model coverage by resizing the recovery rate"

plot_coverage(df_targets  = df_targets,
              df_coverage = rbind(df_coverage_1,
                                  df_coverage_2,
                                  df_coverage_3))

### Readjust transmission rate: r_beta ------------------------------------
# 1. Compute model-data residuals on previous coverage (mean values) 
v_params_3 <- colMeans(m_params_3)
modCost3   <- costFun(v_params     = v_params_3, 
                      l_params_all = l_params_all,
                      v_dates      = v_dates)

# 2. Local sensitivity analysis
SensRes3 <- sensFun(func         = costFun,
                    parms        = v_params_3,
                    l_params_all = l_params_all,
                    v_dates      = v_dates)

# Plot sensitivity
plot_Sfun(SensRes3)

# Rank parameters according to their importance
m_Sens3 <- summary(SensRes3)

# value: value of the parameter 
# Mean: mean sensitivity
# L1: L1 norm
# L2: L2 norm
df_Sens3 <- as.data.frame(m_Sens3)

# Order by L2
df_Sens3[order(df_Sens3$L2, decreasing = T),] 

# 3. Compute the change in the prior distribution bounds based on the 
# sensitivity analysis
df_changes_3 <- data.frame(NULL)
for(target_i in c("DXIncTot")){ # target_i = "DXIncTot"
  
  # Filter elasticities
  df_epsi <- SensRes3 %>% 
    filter(var == target_i) %>% 
    select(x,r_beta) %>% 
    rename(epsilon = r_beta)
  
  # Change in y --> change in par
  df_temp <- modCost3$residuals %>% 
    filter(name == target_i) %>% 
    select(name,x,obs,mod) %>% 
    mutate(y_change = ifelse(mod != 0, (obs-mod)/mod, 0)) %>% 
    left_join(df_epsi,
              by = "x") %>% 
    mutate(par_change = ifelse(epsilon == 0, NA, y_change/epsilon),
           par_change_lb = ifelse(par_change < 0,
                                  NA, (1+par_change)*v_lb_3["r_beta"]),
           par_change_ub = ifelse(par_change > 0 |
                                    (1+par_change)*v_ub_3["r_beta"] < 0.3, 
                                  NA, (1+par_change)*v_ub_3["r_beta"]))
  
  # Bind data
  df_changes_3 <- rbind(df_changes_3,
                        df_temp)
  
}

# Change for lower or upper bound
n_change_lb <- mean(df_changes_3$par_change_lb, na.rm = T)
n_change_ub <- mean(df_changes_3$par_change_ub, na.rm = T)

# Resize bounds
v_lb_4 <- v_lb_3
v_ub_4 <- v_ub_3
v_lb_4["r_beta"] <- n_change_lb
v_ub_4["r_beta"] <- n_change_ub

# Run coverage
set.seed(n_seed)  # set seed 
m_params_4 <- sample_prior(n_samp = n_sim,
                           v_lb = v_lb_4,
                           v_ub = v_ub_4)
l_coverage_4 <- run_coverage(l_params_all = l_params_all, 
                             m_params     = m_params_4)
df_coverage_4 <- l_coverage_4$df_model_summ 
df_coverage_4$type_coverage <- 
  "C: A + B + improved model coverage by resizing the transmission rate"

plot_coverage(df_targets  = df_targets,
              df_coverage = rbind(df_coverage_1,
                                  df_coverage_2,
                                  df_coverage_3,
                                  df_coverage_4))

# Range on input search space after resize bounds
v_lb_final <- v_lb_4
v_ub_final <- v_ub_4

# Final coverage
set.seed(n_seed)  # set seed 
m_params_final <- sample_prior(n_samp = n_sim,
                               v_lb   = v_lb_final,
                               v_ub   = v_ub_final)

l_coverage_final <- run_coverage(l_params_all = l_params_all, 
                                 m_params     = m_params_final)
df_coverage_final <- l_coverage_final$df_model_summ
plot_model_out_vs_targets(df_targets = df_targets,
                          df_model   = df_coverage_final)

## Identifiability (Collinearity) analysis --------------------------------
# Compute model-data residuals on previous coverage (mean values) 
v_params <- colMeans(m_params_final)        # mean values

# Sensitivity analysis
SensRes <- sensFun(func         = costFun,
                   parms        = v_params,
                   l_params_all = l_params_all,
                   v_dates      = v_dates)

df_collin <- collin(SensRes) %>% 
  mutate(var     = "All",
         Outcome = paste0(v_outcome_names,collapse = " & "))

# Plot collinearity analysis
plot_CollIndx(df_collin = df_collin, flag_each = F, log_y = T)

## Calibrate --------------------------------------------------------------
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
fit_imis <- IMIS(B            = 1000,      # incremental sample size at each iteration of IMIS
                 B.re         = n_resamp,  # desired posterior sample size
                 number_k     = 20,        # maximum number of iterations in IMIS
                 D            = 0,
                 v_lb         = v_lb,
                 v_ub         = v_ub,
                 l_params_all = l_params_all,
                 v_dates      = v_dates) 

# Calculate computation time
comp_time <- Sys.time() - t_init
comp_time

# obtain draws from posterior
m_calib_res <- fit_imis$resample

# Calculate log-likelihood (overall fit) and posterior probability of each sample
m_calib_res <- cbind(m_calib_res, 
                     "Overall_fit"    = log_lik(m_calib_res[,v_param_names], 
                                                l_params_all = l_params_all,
                                                v_dates      = v_dates),
                     "Posterior_prob" = posterior(m_calib_res[,v_param_names],
                                                  v_lb         = v_lb, 
                                                  v_ub         = v_ub,
                                                  l_params_all = l_params_all,
                                                  v_dates      = v_dates))

# Normalize posterior probability
m_calib_res[,"Posterior_prob"] <- m_calib_res[,"Posterior_prob"]/
  sum(m_calib_res[,"Posterior_prob"])

### Exploring best-fitting input sets -------------------------------------
# Plot the 1000 draws from the posterior with marginal histograms
pairs.panels(m_calib_res[,v_param_names])

# Compute posterior mean
v_calib_post_mean <- colMeans(m_calib_res[,v_param_names])
v_calib_post_mean

# Compute posterior median and 95% credible interval
m_calib_res_95cr <- colQuantiles(m_calib_res[,v_param_names], 
                                 probs = c(0.025, 0.5, 0.975))
m_calib_res_95cr

### Plot maximum-a-posteriori ---------------------------------------------
# Compute MAP parameter set
v_calib_map <- m_calib_res[which.max(m_calib_res[,"Posterior_prob"]),]

# Run model with MAP
l_out_best <- SIR_calibration_out(v_params     = v_calib_map[v_param_names],
                                  l_params_all = l_params_all)

# Plot 
plotrix::plotCI(x = df_targets$Day, y = df_targets$value, 
                ui = df_targets$ub,
                li = df_targets$lb,
                ylim = c(0, 800), 
                xlab = "Time (days)", ylab = "Number of infections")
points(x = l_out_best$DXIncTot$Day, 
       y = l_out_best$DXIncTot$value, 
       pch = 8, col = "red")
legend("topright", 
       legend = c("Target", "Model-predicted output"),
       col = c("black", "red"), pch = c(1, 8))

### Propagate calibrated parameter uncertainty ----------------------------
# Compute IMIS posterior predicted outputs
m_out <- matrix(NA, nrow = n_resamp, ncol = nrow(df_targets))  

# Run model for each posterior parameter set
for(i in 1:n_resamp){ # i = 1
  l_out_model <- SIR_calibration_out(v_params     = m_calib_res[i, ],
                                     l_params_all = l_params_all)
  df_out_model <- l_out_model$DXIncTot
  m_out[i, ] <- df_out_model$value
  if(i/100==round(i/100,0)) { 
    cat('\r',paste(i/n_resamp*100,"% done",sep=""))
  }
}

# Posterior predicted mean
m_out_postmean <- colMeans(m_out)

# Implement PRE-CISE (narrow) ---------------------------------------------
# 0. Define input search space and plot first coverage
v_lb_1 <- c(r_beta = 0.30, r_gamma = 0.60)  # lower bound
v_ub_1 <- c(r_beta = 1.00, r_gamma = 0.90)  # upper bound

# Run initial coverage analysis
set.seed(n_seed) # set seed
m_params_1 <- sample_prior(n_samp = n_sim,
                           v_lb   = v_lb_1,
                           v_ub   = v_ub_1)

l_coverage_1  <- run_coverage(l_params_all = l_params_all, 
                              m_params     = m_params_1)
df_coverage_1 <- l_coverage_1$df_model_summ %>% 
  mutate(type_coverage = "Initial model coverage")

plot_coverage(df_targets  = df_targets,
              df_coverage = df_coverage_1)

## Resize the prior distribution bounds -----------------------------------
### Adjust transmission rate: r_beta --------------------------------------
# 1. Compute model-data residuals on initial coverage (mean values) 
v_params_1 <- colMeans(m_params_1)
modCost1   <- costFun(v_params     = v_params_1, 
                      l_params_all = l_params_all,
                      v_dates     = v_dates)

# Plot residuals
plot_costFun(Cfun = modCost1, wgt_res = T) # weighted residuals
plot_costFun(modCost1, wgt_res = F)        # unweighted residuals

# 2. Local sensitivity analysis
SensRes1 <- sensFun(func         = costFun,
                    parms        = v_params_1,
                    l_params_all = l_params_all,
                    v_dates      = v_dates)

# Plot sensitivity
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
for(target_i in c("DXIncTot")){ # target_i = "DXIncTot"
  
  # Filter elasticities
  df_epsi <- SensRes1 %>% 
    filter(var == target_i) %>% 
    select(x,r_beta) %>% 
    rename(epsilon = r_beta)
  
  # Change in y --> change in par
  df_temp <- modCost1$residuals %>% 
    filter(name == target_i) %>% 
    select(name,x,obs,mod) %>% 
    mutate(y_change = ifelse(mod != 0, (obs-mod)/mod, 0)) %>% 
    left_join(df_epsi,
              by = "x") %>% 
    mutate(par_change = ifelse(epsilon == 0, NA, y_change/epsilon),
           par_change_ub = ifelse(par_change < 0 | 
                                    (1+par_change)*v_ub_1["r_beta"]>5,
                                  NA, (1+par_change)*v_ub_1["r_beta"]))
  
  # Bind data
  df_changes_1 <- rbind(df_changes_1,
                        df_temp)
  
}

# Change for lower or upper bound
n_change_ub <- mean(df_changes_1$par_change_ub, na.rm = T)

# Resize bounds
v_lb_2 <- v_lb_1
v_ub_2 <- v_ub_1
v_ub_2["r_beta"] <- n_change_ub

# Run coverage
set.seed(n_seed) # set seed
m_params_2 <- sample_prior(n_samp = n_sim,
                           v_lb = v_lb_2,
                           v_ub = v_ub_2)
l_coverage_2 <- run_coverage(l_params_all = l_params_all, m_params = m_params_2)
df_coverage_2 <- l_coverage_2$df_model_summ 
df_coverage_2$type_coverage <- 
  "A: Improved model coverage by resizing the transmission rate"

plot_coverage(df_targets  = df_targets,
              df_coverage = rbind(df_coverage_1,
                                  df_coverage_2))

### Adjust recovery rate: r_gamma -----------------------------------------
# 1. Compute model-data residuals on previous coverage (mean values)
v_params_2 <- colMeans(m_params_2)        # mean values
modCost2   <- costFun(v_params     = v_params_2, 
                      l_params_all = l_params_all,
                      v_dates      = v_dates)

# 2. Local sensitivity analysis
SensRes2 <- sensFun(func         = costFun,
                    parms        = v_params_2,
                    l_params_all = l_params_all,
                    v_dates      = v_dates)

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
for(target_i in c("DXIncTot")){ # target_i = "DXIncTot"
  
  # Filter elasticities
  df_epsi <- SensRes2 %>% 
    filter(var == target_i) %>% 
    select(x,r_gamma) %>% 
    rename(epsilon = r_gamma)
  
  # Change in y --> change in par
  df_temp <- modCost2$residuals %>% 
    filter(name == target_i) %>% 
    select(name,x,obs,mod) %>%  
    mutate(y_change = ifelse(mod != 0, (obs-mod)/mod, 0)) %>% 
    left_join(df_epsi,
              by = "x") %>% 
    mutate(par_change = ifelse(epsilon == 0, NA, y_change/epsilon),
           par_change_lb = ifelse(par_change > 0 | 
                                    (1+par_change)*v_lb_2["r_gamma"] < 0, 
                                  NA, (1+par_change)*v_lb_2["r_gamma"]))
  
  # Bind data
  df_changes_2 <- rbind(df_changes_2,
                        df_temp)
  
}

# Change for lower or upper bound
n_change_lb <- min(df_changes_2$par_change_lb, na.rm = T) 

# Resize bounds
v_lb_3 <- v_lb_2
v_ub_3 <- v_ub_2
v_lb_3["r_gamma"] <- n_change_lb

# Run coverage
set.seed(n_seed) # set seed
m_params_3 <- sample_prior(n_samp = n_sim,
                           v_lb   = v_lb_3,
                           v_ub   = v_ub_3)
l_coverage_3 <- run_coverage(l_params_all = l_params_all, m_params = m_params_3)
df_coverage_3 <- l_coverage_3$df_model_summ
df_coverage_3$type_coverage <- 
  "B: A + improved model coverage by resizing the recovery rate"

plot_coverage(df_targets  = df_targets,
              df_coverage = rbind(df_coverage_1,
                                  df_coverage_2,
                                  df_coverage_3))

### Readjust transmission rate: r_beta ------------------------------------
# 1. Compute model-data residuals on previous coverage (mean values)
v_params_3 <- colMeans(m_params_3)
modCost3   <- costFun(v_params     = v_params_3,
                      l_params_all = l_params_all,
                      v_dates = v_dates)

# 2. Local sensitivity analysis
SensRes3 <- sensFun(func         = costFun,
                    parms        = v_params_3,
                    l_params_all = l_params_all,
                    v_dates = v_dates)

SensRes3

# Plot sensitivity
plot_Sfun(SensRes3)

# Rank parameters according to their importance
m_Sens3 <- summary(SensRes3)

# value: value of the parameter
# Mean: mean sensitivity
# L1: L1 norm
# L2: L2 norm
df_Sens3 <- as.data.frame(m_Sens3)

# Order by L2
df_Sens3[order(df_Sens3$L2, decreasing = T),]

# Transmission rate: r_beta
# 3. Compute the change in the prior distribution bounds based on the
# sensitivity analysis
df_changes_3 <- data.frame(NULL)
for(target_i in c("DXIncTot")){ # target_i = "DXIncTot"
  
  # Filter elasticities
  df_epsi <- SensRes3 %>%
    filter(var == target_i) %>%
    select(x,r_beta) %>%
    rename(epsilon = r_beta)
  
  # Change in y --> change in par
  df_temp <- modCost3$residuals %>%
    filter(name == target_i) %>%
    select(name,x,obs,mod) %>%
    mutate(y_change = ifelse(mod != 0, (obs-mod)/mod, 0)) %>%
    left_join(df_epsi,
              by = "x") %>%
    mutate(par_change = ifelse(epsilon == 0, NA, y_change/epsilon),
           par_change_lb = ifelse(par_change < 0,
                                  NA, (1+par_change)*v_lb_3["r_beta"]),
           par_change_ub = ifelse(par_change > 0 | 
                                    (1+par_change)*v_ub_3["r_beta"] < 0.3,
                                  NA, (1+par_change)*v_ub_3["r_beta"]))
  
  # Bind data
  df_changes_3 <- rbind(df_changes_3,
                        df_temp)
  
}

# Change for lower or upper bound
n_change_lb <- min(df_changes_3$par_change_lb, na.rm = T)
n_change_ub <- max(df_changes_3$par_change_ub, na.rm = T)

# Resize bounds
v_lb_4 <- v_lb_3
v_ub_4 <- v_ub_3
v_lb_4["r_beta"] <- n_change_lb
v_ub_4["r_beta"] <- n_change_ub

# Run coverage
set.seed(n_seed) # set seed
m_params_4 <- sample_prior(n_samp = n_sim,
                           v_lb = v_lb_4,
                           v_ub = v_ub_4)
l_coverage_4 <- run_coverage(l_params_all = l_params_all, m_params = m_params_4)
df_coverage_4 <- l_coverage_4$df_model_summ
df_coverage_4$type_coverage <-
  "C: A + B + improved model coverage by resizing the transmission rate"

plot_coverage(df_targets  = df_targets,
              df_coverage = rbind(df_coverage_1,
                                  df_coverage_2,
                                  df_coverage_3,
                                  df_coverage_4))

# Range on input search space after resize bounds
v_lb_final <- v_lb_4
v_ub_final <- v_ub_4

# Final coverage
set.seed(n_seed)
m_params_final <- sample_prior(n_samp = n_sim,
                               v_lb   = v_lb_final,
                               v_ub   = v_ub_final)

l_coverage_final <- run_coverage(l_params_all = l_params_all, 
                                 m_params     = m_params_final)
df_coverage_final <- l_coverage_final$df_model_summ
plot_model_out_vs_targets(df_targets = df_targets,
                          df_model   = df_coverage_final)

## Identifiability (Collinearity) analysis --------------------------------
# Compute model-data residuals on previous coverage (mean values) 
v_params <- colMeans(m_params_final)        # mean values

# Sensitivity analysis
SensRes <- sensFun(func         = costFun,
                   parms        = v_params,
                   l_params_all = l_params_all,
                   v_dates      = v_dates)

df_collin <- collin(SensRes) %>% 
  mutate(var     = "All",
         Outcome = paste0(v_outcome_names,collapse = " & "))

# Plot collinearity analysis
plot_CollIndx(df_collin = df_collin, flag_each = F, log_y = T)

