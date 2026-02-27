#' Generate model outputs for calibration from a parameter set
#'
#' \code{calibration_out} computes model outputs for the Sick-Sicker model
#' to be used for calibration routines.
#'
#' @param v_params Vector of parameters that need to be calibrated.
#' 
#' @return A data.frame with model-predicted outcomes for each target.
#' @export
calibration_out <- function(v_params){
  
  # Run model
  df_model_out <- sick_sicker_cstm(v_params)
  
  # Data.frame for each calibration target
  # Survival
  df_Surv <- data.frame(time = df_model_out$time,
                        Outcome = "Survival",
                        abbrev_outcome = "Surv",
                        type = "Model", 
                        value = df_model_out$Surv) %>% 
    filter(time %in% c(10,20,30))
  
  # Prevalence
  df_Prev <- data.frame(time = df_model_out$time,
                        Outcome = "Prevalence",
                        abbrev_outcome = "Prev",
                        type = "Model",
                        value = df_model_out$Prev) %>% 
    filter(time %in% c(10,20,30))
  
  # Proportion who are Sick ("PropSick")
  # among all those afflicted (Sick+Sicker)
  df_PropSick <- data.frame(time = df_model_out$time,
                            Outcome = "Proportion of sick",
                            abbrev_outcome = "PropSick",
                            type = "Model", 
                            value = df_model_out$PropSick) %>% 
    filter(time %in% c(10,20,30))
  
  # Return data.frame
  return(list(df_Surv     = df_Surv,
              df_Prev     = df_Prev,
              df_PropSick = df_PropSick))
}

#' Sample from prior distributions of calibrated parameters
#'
#' \code{sample_prior} generates a sample of calibrated parameters 
#' from their prior distribution.
#' 
#' @param n_samp Number of samples.
#' @param v_lb Vector containing the lower bounds of the calibrated 
#' parameters.
#' @param v_ub Vector containing the upper bounds of the calibrated 
#' parameters.
#'
#' @return A matrix with n number of columns as number of parameters to 
#' be calibrated and \code{n_samp} rows. Each row corresponds to a 
#' parameter set sampled from their prior distributions.
#' @export
sample_prior <- function(n_samp, v_lb, v_ub){
  
  m_lhs_unit   <- randomLHS(n = n_samp, k = n_param)
  m_param_samp <- matrix(nrow = n_samp, ncol = n_param)
  colnames(m_param_samp) <- v_param_names
  for (i in 1:n_param){
    m_param_samp[, i] <- qunif(m_lhs_unit[,i],
                               min = v_lb[i],
                               max = v_ub[i])
    # ALTERNATIVE prior using beta (or other) distributions
    # m_param_samp[, i] <- qbeta(m_lhs_unit[,i],
    #                            shape1 = 1,
    #                            shape2 = 1)
  }
  return(m_param_samp)
}

#' Evaluate log-prior of calibrated parameters
#'
#' \code{log_prior} computes a log-prior value for one (or multiple) parameter
#' set(s) based on their prior distributions.
#' 
#' @param v_params Vector (or matrix) of model parameters.
#' 
#' @return A scalar (or vector) with log-prior values
#' @export
log_prior <- function(v_params){
  
  if(is.null(dim(v_params))) { # If vector, change to matrix
    v_params <- t(v_params) 
  }
  n_samp <- nrow(v_params)
  colnames(v_params) <- v_param_names
  lprior <- rep(0, n_samp)
  for (i in 1:n_param){
    lprior <- lprior + dunif(v_params[, i],
                             min = v_lb[i],
                             max = v_ub[i], 
                             log = T)
    # ALTERNATIVE prior using beta distributions
    # lprior <- lprior + dbeta(v_params[, i],
    #                          shape1 = 1,
    #                          shape2 = 1, 
    #                          log = T)
  }
  return(lprior)
}

#' Evaluate prior of calibrated parameters
#'
#' \code{prior} computes a prior value for one (or multiple) parameter set(s).
#' 
#' @param v_params Vector (or matrix) of model parameters.
#'
#' @return A scalar (or vector) with prior values.
#' @export
prior <- function(v_params) { 
  return(exp(log_prior(v_params)))
}

#' Log likelihood normal distribution
#'
#' \code{log_lik} computes a log-likelihood value for one (or multiple)
#' parameter set(s).
#'
#' @param v_params_calib Vector (or matrix) of parameters that need to be
#' calibrated.
#'
#' @return A scalar (or vector) with log-likelihood values.
#' @export
log_lik <- function(v_params){
  
  # par_vector: a vector (or matrix) of model parameters 
  if(is.null(dim(v_params))) { # If vector, change to matrix
    v_params <- t(v_params) 
  }
  n_samp <- nrow(v_params)
  v_llik <- matrix(0, nrow = n_samp, ncol = n_outcomes) 
  llik_overall <- numeric(n_samp)
  for(j in 1:n_samp) { # j=1
    jj <- tryCatch( { 
      ###   Run model for parameter set "v_params" ###
      l_model_out <- calibration_out(v_params[j, ])
      
      ### Targets
      df_targets_Prev <- df_targets %>% 
        filter(abbrev_outcome == "Prev")
      df_targets_PropSick <- df_targets %>% 
        filter(abbrev_outcome == "PropSick")
      df_targets_Surv <- df_targets %>% 
        filter(abbrev_outcome == "Surv")
      
      ###  Calculate log-likelihood of model outputs to targets  ###
      # TARGET 1: Survival ("Surv")
      # log likelihood  
      v_llik[j, 1] <- sum(dnorm(x = df_targets_Surv$value,
                                mean = l_model_out$df_Surv$value,
                                sd = df_targets_Surv$se,
                                log = T))
      
      # TARGET 2: Prevalence ("Prev")
      # log likelihood
      v_llik[j, 2] <- sum(dnorm(x = df_targets_Prev$value,
                                mean = l_model_out$df_Prev$value,
                                sd = df_targets_Prev$se,
                                log = T))
      
      # TARGET 3: Proportion of sick in Sick state ("PropSick")
      # log likelihood
      v_llik[j, 3] <- sum(dnorm(x = df_targets_PropSick$value,
                                mean = l_model_out$df_PropSick$value,
                                sd = df_targets_PropSick$se,
                                log = T))
      
      # OVERALL 
      llik_overall[j] <- sum(v_llik[j, ])
    }, error = function(e) NA) 
    if (is.na(jj)) { llik_overall <- -Inf }
  } # End loop over sampled parameter sets
  # return LLIK
  return(llik_overall)
}

#' Likelihood
#'
#' \code{likelihood} computes a likelihood value for one (or multiple)
#' parameter set(s).
#'
#' @param v_params_calib Vector (or matrix) of model parameters.
#'
#' @return A scalar (or vector) with likelihood values.
#' @export
likelihood <- function(v_params){ 
  return(exp(log_lik(v_params))) 
}

#' Evaluate log-posterior of calibrated parameters
#'
#' \code{log_post} computes a log-posterior value for one (or multiple) 
#' parameter set(s) based on the simulation model, likelihood functions and 
#' prior distributions.
#' 
#' @param v_params_calib Vector (or matrix) of model parameters. 
#' 
#' @return A scalar (or vector) with log-posterior values.
#' @export
log_post <- function(v_params) { 
  lpost <- log_prior(v_params) + log_lik(v_params)
  return(lpost) 
}

#' Evaluate posterior of calibrated parameters
#'
#' \code{posterior} computes a posterior value for one (or multiple) 
#' parameterset(s).
#' 
#' @param v_params_calib Vector (or matrix) of model parameters.
#' 
#' @return A scalar (or vector) with posterior values.
#' @export
posterior <- function(v_params) { 
  exp(log_post(v_params)) 
}


