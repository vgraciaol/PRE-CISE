#------------------------------------------------------------------------------#
# This code was created by the DARTH workgroup (www.darthworkgroup.com).       #
# When using or modifying this code, please do so with attribution and         #
# cite our publications:                                                       #
#                                                                              #
# - Alarid-Escudero F, Maclehose RF, Peralta Y, Kuntz KM, Enns EA.             #
#   Non-identifiability in model calibration and implications for              #
#   medical decision making. Med Decis Making. 2018; 38(7):810-821.            #
#                                                                              #
# - Jalal H, Pechlivanoglou P, Krijkamp E, Alarid-Escudero F, Enns E,          #
#   Hunink MG. An Overview of R in Health Decision Sciences.                   #
#   Med Decis Making. 2017; 37(3): 735-746.                                    #
#                                                                              #
# A walkthrough of the code could be found in the follwing link:               #
# - https://darth-git.github.io/calibSMDM2018-materials/                       #
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# This script illustrate the uses of sensitivity and identifiability analyses. #
# for efficient calibration of health policy models by calibrating a Sick-     #
# Sicker 4-state Markov Model to survival probability, prevalence, and         #
# proportion of sick using the IMIS algorithm.                                 #
#                                                                              #
#   Authors:                                                                   #
#       Fernando Alarid-Escudero, PhD, <falarid@stanford.edu>                  #
#       Valeria Gracia Olvera, MSc, <vgracia@stanford.edu>                     #
#                                                                              #
#   Created: 11/14/24                                                          #
#   Last update: 11/14/24                                                      #
#------------------------------------------------------------------------------#

# Calibration Specifications ----------------------------------------------
# Model: Sick-Sicker 4-state Markov Model
# Inputs to be calibrated: p.S1S2, hr.S1, hr.S2
# Targets: Surv, Prev, PropSick
# Search method: Random search using Latin-Hypercube Sampling
# Goodness-of-fit measure: Sum of Log-Likelihood

# Load libraries and functions --------------------------------------------
# calibration functionality
library(lhs)
library(matrixStats) # package used for sumamry statistics
library(FME)         # For identifiability and collinearity analysis
library(tidyr)

# visualization
library(plotrix)
library(psych)
library(ggplot2)
library(ggpubr)
library(scatterplot3d) # three inputs to estimate

# helpful functions
source("R/helper_functions.R") # from dampack package
source("R/visualization_functions.R")

# calibration
source("R/IMIS_function.R") # from IMIS package
source("R/calibration_functions.R")

# Target data -------------------------------------------------------------
load("data/SickSicker_CalibTargets.RData")
l_targets <- SickSicker_targets

# Gather targets in a data.frame
df_targets <- rbind(l_targets$Surv %>% 
                      mutate(type = "Target",
                             Outcome = "Surv") %>% 
                      select(type, Outcome, time, value, lb, ub),
                    l_targets$Prev %>% 
                      mutate(type = "Target",
                             Outcome = "Prev") %>% 
                      select(type, Outcome, time, value, lb, ub),
                    l_targets$PropSick %>% 
                      mutate(type = "Target",
                             Outcome = "PropSick") %>% 
                      select(type, Outcome, time, value, lb, ub))

# Plot the targets
plot_targets(df_targets)

# Model -------------------------------------------------------------------
# -inputs are parameters to be estimated through calibration
# -outputs correspond to the target data

# Load model
source("R/sick-sicker-model.R") # creates the function sick_sicker_cstm()

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
v_target_names <- c("Surv", "Prev", "PropSick")
n_target <- length(v_target_names)

# Calibration functions (check that it works) -----------------------------
# View resulting parameter set samples
pairs.panels(sample_prior(1000))

# Prior
log_prior(v_params = v_params_test)
log_prior(v_params = sample_prior(10))

prior(v_params = v_params_test)
prior(v_params = sample_prior(10))

# Likelihood
log_lik(v_params = v_params_test)
log_lik(v_params = sample_prior(10))

likelihood(v_params = v_params_test)
likelihood(v_params = sample_prior(10))

# Posterior
log_post(v_params = v_params_test)
log_post(v_params = sample_prior(10))

posterior(v_params = v_params_test)
posterior(v_params = sample_prior(10))

# Coverage analysis -------------------------------------------------------
# Range on input search space
v_lb <- c(p_S1S2 = 0.10, hr_S1 = 1.0, hr_S2 = 1.0) # lower bound
v_ub <- c(p_S1S2 = 0.20, hr_S1 = 1.1, hr_S2 = 1.1) # upper bound

# Number and parameters set to run
n_sim <- 1000
m_params <- sample_prior(n_samp = n_sim)

df_model <- data.frame(NULL) # empty data.frame to store results
for(iter_i in 1:n_sim){ # iter_i = 1
  
  # Run model
  l_model_out <- calibration_out(v_params = m_params[iter_i, ])
  
  # Output of the model in a data.frame
  df_temp <- rbind(l_model_out$df_Surv,
                    l_model_out$df_Prev,
                    l_model_out$df_PropSick) %>% 
    mutate(iter = iter_i)
  
  # Gather results
  df_model <- rbind(df_model, df_temp) 
  
}

# Summarise
df_model_summ <- df_model %>%
  group_by(type, Outcome, time) %>%
  summarise(m_value = mean(value),
            lb_value = quantile(value,probs = 0.025),
            ub_value = quantile(value, probs = 0.975)) %>%
  ungroup() %>%
  rename(value = m_value,
         lb = lb_value,
         ub = ub_value)

# Plot model vs targets
plot_model_out_vs_targets(df_targets = df_targets,
                          df_model   = df_model_summ)

# Collinearity analysis ---------------------------------------------------
v_params_test <- colMeans(m_params)

# Cost function
costFun <- function(v_params) {
  
  model_out <- sick_sicker_cstm(v_params)
  
  # Get the data in an specific format to compute cost
  targetSurv <- l_targets$Surv[, c(1, 4, 5)]
  colnames(targetSurv) <- c("time","Surv","se")
  
  targetPrev <- l_targets$Prev[, c(1, 4, 5)]
  colnames(targetPrev) <- c("time","Prev","se")
  
  targetPropSick <- l_targets$PropSick[, c(1, 4, 5)]
  colnames(targetPropSick) <- c("time","PropSick","se")
  
  # Compute "costs" for each target
  cost <- modCost(model = model_out[[1]], 
                  obs = targetSurv, 
                  err = "se")
  
  cost <- modCost(model = model_out[[2]], 
                  obs = targetPrev, 
                  err="se", 
                  cost = cost)
  
  cost <- modCost(model = model_out[[3]], 
                  obs = targetPropSick, 
                  err="se", 
                  cost = cost)
  return(cost)
}

# Compute model costs
Cfun <- costFun(v_params_test)
Cfun

# Plot residuals
plot_costFun(Cfun, wgt_res = T) # weighted residuals

# Create Sensitivity function
Sfun <- sensFun(func = costFun, 
                parms = v_params_test, 
                senspar = 1:3)
summary(Sfun)

# Plot Sensitivity
plot_Sfun(Sfun)

# Run regression: model-predicted vs parameters (MISSING!)


# Identifiability (Collinearity) analysis ---------------------------------
ident <- collin(Sfun)
plot_CollIndx(ident)

plot_CollIndx(collin(Sfun, N = 2))
plot_CollIndx(collin(Sfun, N = 2), log_y = F)
plot_CollIndx(collin(Sfun, which = c("Surv")))
plot_CollIndx(collin(Sfun, which = c("Prev")))
plot_CollIndx(collin(Sfun, which = c("PropSick")))
plot_CollIndx(collin(Sfun, which = c("Surv", "PropSick")))
plot_CollIndx(collin(Sfun, which = c("Prev", "PropSick")))
plot_CollIndx(collin(Sfun, which = c("Surv", "Prev")))

# Calibrate! --------------------------------------------------------------
# Specify seed (for reproducible sequence of random numbers)
set.seed(111124)

# number of random samples
n_resamp <- 1000

# record start time of calibration
t_init <- Sys.time()

# Bayesian calibration using the IMIS algorithm (Raftery & Bao, 2010)
fit_imis <- IMIS(B        = 1000,      # the incremental sample size at each iteration of IMIS
                 B.re     = n_resamp,  # the desired posterior sample size
                 number_k = 10,        # the maximum number of iterations in IMIS
                 D        = 0) 

# obtain draws from posterior
m_calib_res <- fit_imis$resample

# Calculate log-likelihood (overall fit) and posterior probability of each sample
m_calib_res <- cbind(m_calib_res, 
                     "Overall_fit" = log_lik(m_calib_res[,v_param_names]),
                     "Posterior_prob" = posterior(m_calib_res[,v_param_names]))

# normalize posterior probability
m_calib_res[,"Posterior_prob"] <- m_calib_res[,"Posterior_prob"]/sum(m_calib_res[,"Posterior_prob"])

# Calculate computation time
comp_time <- Sys.time() - t_init
comp_time

# Exploring best-fitting input sets ---------------------------------------
# Plot the 1000 draws from the posterior
v_post_color <- scales::rescale(m_calib_res[,"Posterior_prob"])
s3d <- scatterplot3d(x = m_calib_res[, 1],
                     y = m_calib_res[, 2],
                     z = m_calib_res[, 3],
                     color = scales::alpha("black", v_post_color),
                     xlim = c(v_lb[1],v_ub[1]), ylim = c(v_lb[2],v_ub[2]), zlim = c(v_lb[3],v_ub[3]),
                     xlab = v_param_names[1], ylab = v_param_names[2], zlab = v_param_names[3])
# add center of Gaussian components
s3d$points3d(fit_imis$center, col = "red", pch = 8)

# Plot the 1000 draws from the posterior with marginal histograms
pairs.panels(m_calib_res[,v_param_names])

# Compute posterior mean
v_calib_post_mean <- colMeans(m_calib_res[,v_param_names])
v_calib_post_mean

# Compute posterior median and 95% credible interval
m_calib_res_95cr <- colQuantiles(m_calib_res[,v_param_names], probs = c(0.025, 0.5, 0.975))
m_calib_res_95cr

# Compute maximum-a-posteriori (MAP) parameter set
v_calib_map <- m_calib_res[which.max(m_calib_res[,"Posterior_prob"]),]

## Plot model-predicted output vs targets ---------------------------------

df_model_posterior <- data.frame(NULL) # empty data.frame to store results
for(iter_i in 1:n_resamp){ # iter_i = 1
  
  # Run model
  l_model_out <- calibration_out(v_params = m_calib_res[iter_i, ])
  
  # Output of the model in a data.frame
  df_temp <- rbind(l_model_out$df_Surv,
                   l_model_out$df_Prev,
                   l_model_out$df_PropSick) %>% 
    mutate(iter = iter_i)
  
  # Gather results
  df_model_posterior <- rbind(df_model_posterior, df_temp) 
  
}

# Summarise
df_model_posterior_summ <- df_model_posterior %>%
  group_by(type, Outcome, time) %>%
  summarise(m_value = mean(value),
            lb_value = quantile(value,probs = 0.025),
            ub_value = quantile(value, probs = 0.975)) %>%
  ungroup() %>%
  rename(value = m_value,
         lb = lb_value,
         ub = ub_value)

plot_model_out_vs_targets(df_targets = df_targets,
                          df_model   = df_model_posterior_summ)

l_out_best <- sick_sicker_cstm(v_calib_map[v_param_names])

# TARGET 1: Survival ("Surv")
plotrix::plotCI(x = l_targets$Surv$time, y = l_targets$Surv$value, 
                ui = l_targets$Surv$ub,
                li = l_targets$Surv$lb,
                ylim = c(0, 1), 
                xlab = "Time", ylab = "Pr Survive")
points(x = l_targets$Surv$time, 
       y = l_out_best[[1]]$Surv, 
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
points(x = l_targets$Prev$time,
       y = l_out_best[[2]]$Prev,
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
points(x = l_targets$PropSick$time,
       y = l_out_best[[3]]$PropSick,
       pch = 8, col = "red")
legend("topright",
       legend = c("Target", "Model-predicted output"),
       col = c("black", "red"), pch = c(1, 8))

# Propagate calibrated parameter uncertainty ------------------------------
# Compute IMIS posterior predicted outputs
m_out_surv     <- matrix(NA, nrow = n_resamp, ncol = length(l_targets$Surv$value))  
m_out_prev     <- matrix(NA, nrow = n_resamp, ncol = length(l_targets$Prev$value))  
m_out_propsick <- matrix(NA, nrow = n_resamp, ncol = length(l_targets$PropSick$value))   

# Run model for each posterior parameter set
for(i in 1:n_resamp){ # i = 1
  l_model_res_temp <- sick_sicker_cstm(m_calib_res[i, ])
  m_out_surv[i, ] <- l_model_res_temp[[1]]$Surv
  m_out_prev[i, ] <- l_model_res_temp[[2]]$Prev
  m_out_propsick[i, ] <- l_model_res_temp[[3]]$PropSick
  if(i/100==round(i/100,0)) { 
    cat('\r',paste(i/n_resamp*100,"% done",sep=""))
  }
}

# Posterior predicted mean
m_out_surv_postmean <- colMeans(m_out_surv)
m_out_prev_postmean <- colMeans(m_out_prev)
m_out_propsick_postmean <- colMeans(m_out_propsick)
