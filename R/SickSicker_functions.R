# Model-decision functions ------------------------------------------------

#' Sick Sicker Model
#' 
#' @param vparams Vector containing the model parameters. 
#' 
#' @return List containing the epidemiological outputs of the model
#' @export
sick_sicker_cstm <- function(v_params) {
  
  with(as.list(v_params), {
    
    ## Markov model parameters
    # age at baseline
    age     <- 25
    # maximum age of follow up
    max_age <- 55
    # time horizon, number of cycles
    n_t     <- max_age - age            
    # the 4 states of the model: Healthy (H), Sick (S1), Sicker (S2), Dead (D)
    v_n     <- c("H", "S1", "S2", "D")  
    # number of health states 
    n_s     <- length(v_n)              
    
    # Transition probabilities and hazard ratios
    # probability to die when healthy
    p_HD    <- 0.005
    # probability to become sick when healthy
    p_HS1   <- 0.15 
    # probability to become healthy when sick
    p_S1H   <- 0.5
    # CALIBRATED: probability to become sicker when sick
    # p_S1S2  <- 0.105 
    # CALIBRATED: hazard ratio of death in sick vs healthy
    # hr_S1   <- 3
    # CALIBRATED: hazard ratio of death in sicker vs healthy
    # hr_S2   <- 10          	        
    
    # Compute internal parameters as a function of external parameters
    # rate of death in healthy
    r_HD    <- - log(1 - p_HD)  
    # rate of death in sick
    r_S1D   <- hr_S1 * r_HD 	
    # rate of death in sicker
    r_S2D   <- hr_S2 * r_HD 
    # probability to die in sick
    p_S1D   <- 1 - exp(-r_S1D)    
    # probability to die in sicker
    p_S2D   <- 1 - exp(-r_S2D)           
    
    ### Initialization 
    # Create the cohort trace
    # (n_t + 1 because R doesn't understand  Cycle 0)
    m_M <- matrix(NA, nrow = n_t + 1 , 
                  ncol = n_s,
                  dimnames = list(0:n_t, v_n))     
    # Initialize Markov trace
    m_M[1, ] <- c(1, 0, 0, 0)                      
    
    # Create transition probability matrix for NO treatment
    m_P <- matrix(0,
                  nrow = n_s, 
                  ncol = n_s,
                  dimnames = list(v_n, v_n))
    
    # Fill in the transition probability array
    # From Healthy
    m_P["H", "H"]  <- 1 - (p_HS1 + p_HD)
    m_P["H", "S1"] <- p_HS1
    m_P["H", "D"]  <- p_HD
    
    # From Sick
    m_P["S1", "H"]  <- p_S1H
    m_P["S1", "S1"] <- 1 - (p_S1H + p_S1S2 + p_S1D)
    m_P["S1", "S2"] <- p_S1S2
    m_P["S1", "D"]  <- p_S1D
    
    # From Sicker
    m_P["S2", "S2"] <- 1 - p_S2D
    m_P["S2", "D"]  <- p_S2D
    
    # From Dead
    m_P["D", "D"] <- 1
    
    # Check rows add up to 1
    if (!isTRUE(all.equal(as.numeric(rowSums(m_P)), as.numeric(rep(1, n_s))))) {
      stop("This is not a valid transition Matrix")
    }
    
    ### Process
    for (t in 1:n_t){  # throughout the number of cycles
      # estimate the Markov trace for cycle the next cycle (t + 1)
      m_M[t + 1, ] <- m_M[t, ] %*% m_P             
    }
    
    ### Epidemiological outputs
    # Overall Survival (OS)
    # Calculate the overall survival (OS) probability for no treatment
    v_os <- 1 - m_M[, "D"]                
    
    # Disease prevalence
    v_prev <- rowSums(m_M[, c("S1", "S2")])/v_os
    
    # Proportion of sick in S1 state
    v_prop_S1 <- m_M[, "S1"] / v_prev
    
    # Return output
    out <- data.frame(time     = 1:30,
                      PropSick = v_prop_S1[-1],
                      Prev     = v_prev[-1],
                      Surv     = v_os[-1])
    
    return(out)
  }
  )
}


# Visualization functions -------------------------------------------------

#' Plot model outputs vs targets
#'
#' \code{plot_model_out_vs_targets} plots model outputs vs targets.
#' 
#' @param df_targets Data.frame containing the calibration targets.
#' 
#' @return A ggplot2 object.
#' 
#' @export
plot_targets <- function(df_targets){
  
  # Bind data.frames
  df_Surv <- df_targets %>% filter(abbrev_outcome == "Surv")
  df_Prev <- df_targets %>% filter(abbrev_outcome == "Prev")
  df_PropSick <- df_targets %>% filter(abbrev_outcome == "PropSick")
  
  # Plot Survival ("Surv")
  plotSurv <- ggplot(df_Surv, 
                     aes(x = time, y = value, ymin = lb, ymax = ub)) +
    geom_point(shape = 1) + geom_errorbar() +
    scale_y_continuous(breaks = number_ticks(6)) +
    scale_x_continuous(breaks = number_ticks(6)) +
    ylab("Survival probability") + xlab("Time (years)") +
    theme_bw(base_size = 18) +
    theme(legend.position = "top",
          strip.background = element_rect(color = "white", 
                                          fill = "white"),
          axis.text.x = element_text(angle = 0,
                                     vjust = 0.5,
                                     hjust = 0.5)) +
    coord_cartesian(ylim = c(0,1))
  
  # Plot Prevalence ("Prev")
  plotPrev <- ggplot(df_Prev, 
                     aes(x = time, y = value, ymin = lb, ymax = ub)) +
    geom_point(size = 1, shape = 1) + geom_errorbar() +
    scale_y_continuous(breaks = number_ticks(6)) +
    scale_x_continuous(breaks = number_ticks(6)) +
    ylab("Prevalence") + xlab("Time (years)") +
    theme_bw(base_size = 18) +
    theme(legend.position = "top",
          strip.background = element_rect(color = "white", 
                                          fill = "white"),
          axis.text.x = element_text(angle = 0,
                                     vjust = 0.5,
                                     hjust = 0.5)) +
    coord_cartesian(ylim = c(0,1))
  
  # Plot Proportion who are Sick ("PropSick")
  plotPropSick <- ggplot(df_PropSick, 
                         aes(x = time, y = value, ymin = lb, ymax = ub)) +
    geom_point(size = 1, shape = 1) + geom_errorbar() +
    scale_y_continuous(breaks = number_ticks(6)) +
    scale_x_continuous(breaks = number_ticks(6)) +
    ylab("Proportion of sick") + xlab("Time (years)") +
    theme_bw(base_size = 18) +
    theme(legend.position = "top",
          strip.background = element_rect(color = "white", 
                                          fill = "white"),
          axis.text.x = element_text(angle = 0,
                                     vjust = 0.5,
                                     hjust = 0.5)) +
    coord_cartesian(ylim = c(0,1))
  
  # Gather plots
  gg_targets <- ggarrange(plotSurv, NULL, 
                          plotPrev, NULL, 
                          plotPropSick,
                          heights = 1, 
                          widths = c(1,0.05,1,0.05,1),
                          ncol = 5, nrow = 1)
  
  return(gg_targets)
  
}

#' Plot model outputs vs targets
#'
#' \code{plot_model_out_vs_targets} plots model outputs vs targets.
#' 
#' @param df_targets Data.frame containing the calibration targets.
#' @param df_model Data.frame containing the model outputs.
#' 
#' @return A ggplot2 object.
#' 
#' @export
plot_model_out_vs_targets <- function(df_targets, df_model){
  
  # Times to be plotted
  v_times <- sort(unique(df_targets$time))
  
  # Bind data.frames
  df_Surv <- rbind(
    df_targets %>% select(time, Outcome, abbrev_outcome, type, value, lb, ub), 
    df_model %>% select(time, Outcome, abbrev_outcome, type, value, lb, ub)) %>% 
    filter(abbrev_outcome == "Surv")
  
  df_Prev <- rbind(
    df_targets %>% select(time, Outcome, abbrev_outcome, type, value, lb, ub), 
    df_model %>% select(time, Outcome, abbrev_outcome, type, value, lb, ub)) %>% 
    filter(abbrev_outcome == "Prev")
  
  df_PropSick <- rbind(
    df_targets %>% select(time, Outcome, abbrev_outcome, type, value, lb, ub), 
    df_model %>% select(time, Outcome, abbrev_outcome, type, value, lb, ub)) %>% 
    filter(abbrev_outcome == "PropSick")
  
  # Order
  df_Surv$type <- factor(df_Surv$type, levels = c("Target","Model"))
  df_Prev$type <- factor(df_Prev$type, levels = c("Target","Model"))
  df_PropSick$type <- factor(df_PropSick$type, levels = c("Target","Model"))
  
  # Plot Survival ("Surv")
  plotSurv <- ggplot(df_Surv, 
                     aes(x = factor(time), y = value, ymin = lb, ymax = ub,
                         color = type, shape = type)) +
    geom_point(size = 1.5, position = position_dodge(0.8)) +
    geom_errorbar(position = position_dodge(0.8), linewidth = 0.8) +
    scale_shape_manual(NULL, values = c(Target = 1, Model = 8)) +
    scale_color_manual(NULL, values = c(Target = "black", Model = "red")) +
    scale_y_continuous(breaks = number_ticks(6)) +
    scale_x_discrete(breaks = as.character(v_times)) +
    ylab("Survival probability") + xlab("") +
    theme_bw(base_size = 18) +
    theme(legend.position = "none",
          strip.background = element_rect(color = "white", 
                                          fill = "white"),
          axis.text.x = element_text(angle = 0,
                                     vjust = 0.5,
                                     hjust = 0.5)) +
    coord_cartesian(ylim = c(0,1))
  
  # Plot Prevalence ("Prev")
  plotPrev <- ggplot(df_Prev, 
                     aes(x = factor(time), y = value, ymin = lb, ymax = ub,
                         color = type, shape = type)) +
    geom_point(size = 1.5, position = position_dodge(0.8)) +
    geom_errorbar(position = position_dodge(0.8), linewidth = 0.8) +
    scale_shape_manual(NULL, values = c(Target = 1, Model = 8)) +
    scale_color_manual(NULL, values = c(Target = "black", Model = "red")) +
    scale_y_continuous(breaks = number_ticks(6)) +
    scale_x_discrete(breaks = as.character(v_times)) +
    ylab("Prevalence") + xlab("Time (years)") +
    theme_bw(base_size = 18) +
    theme(legend.position = "none",
          strip.background = element_rect(color = "white", 
                                          fill = "white"),
          axis.text.x = element_text(angle = 0,
                                     vjust = 0.5,
                                     hjust = 0.5)) +
    coord_cartesian(ylim = c(0,1))
  
  # Plot Proportion who are Sick ("PropSick")
  plotPropSick <- ggplot(df_PropSick, 
                         aes(x = factor(time), y = value, ymin = lb, ymax = ub,
                             color = type, shape = type)) +
    geom_point(size = 1.5, position = position_dodge(0.8)) +
    geom_errorbar(position = position_dodge(0.8), linewidth = 0.8) +
    scale_shape_manual(NULL, values = c(Target = 1, Model = 8)) +
    scale_color_manual(NULL, values = c(Target = "black", Model = "red")) +
    scale_y_continuous(breaks = number_ticks(6)) +
    scale_x_discrete(breaks = as.character(v_times)) +
    ylab("Proportion of Sick") + xlab("") +
    theme_bw(base_size = 18) +
    theme(legend.position = "inside",
          legend.position.inside = c(0.65,0.90),
          strip.background = element_rect(color = "white", 
                                          fill = "white"),
          axis.text.x = element_text(angle = 0,
                                     vjust = 0.5,
                                     hjust = 0.5)) +
    coord_cartesian(ylim = c(0,1))
  
  # Gather plots
  gg_model_out_vs_targets <- ggarrange(plotSurv, NULL, 
                                       plotPrev, NULL, 
                                       plotPropSick,
                                       heights = 1, 
                                       widths = c(1,0.05,1,0.05,1),
                                       ncol = 5, nrow = 1,
                                       common.legend = F)
  
  return(gg_model_out_vs_targets)
  
}

#' Plot coverage
#'
#' \code{plot_coverage} plots model outputs vs targets.
#' 
#' @param df_targets Data.frame containing the calibration targets.
#' @param df_coverage Data.frame containing the coverage analysis.
#' 
#' @return A ggplot2 object.
#' 
#' @export
plot_coverage <- function(df_targets, df_coverage){
  
  # Times to be plotted
  v_times <- sort(unique(df_targets$time))
  
  # Add type_coverage to targets
  df_targets <- df_targets %>% mutate(type_coverage = "Calibration target")
  
  # Get levels
  v_levels <- c("Calibration target", unique(df_coverage$type_coverage)) 
  
  # Bind data.frames
  df_model_target <- rbind(df_targets, df_coverage)
  df_model_target$type_coverage <- factor(df_model_target$type_coverage, 
                                          levels = v_levels)
  
  # Vector of shapes
  v_shape <- c(1,rep(NA,length(v_levels)-1))
  names(v_shape) <- v_levels
  
  # Vector of colors
  v_colors <- c("black","#FF8000","#56B4E9","#009E73","#CCBFFF","tan","firebrick")
  names(v_colors) <- c(
    "Calibration target",
    "Initial model coverage",                                                           
    "A: Improved model coverage by resizing the probability to become Sicker when Sick",
    "B: A + improved model coverage by resizing the hazard ratio for death in Sicker",  
    "C: A + B + improved model coverage by resizing the hazard ratio for death in Sick",
    "D: A + B + C + improved model coverage by resizing the hazard ratio for death in Sick",
    "E: A + B + C + D + improved model coverage by resizing the hazard ratio for death in Sick"
  )
  
  # Plot Survival ("Surv")
  plotSurv <- ggplot(subset(df_model_target, abbrev_outcome == "Surv"), 
                     aes(x = as.factor(time), 
                         y = value, ymin = lb, ymax = ub,
                         color = type_coverage, shape = type_coverage)) +
    geom_point(size = 1.5, position = position_dodge(0.8)) +
    geom_errorbar(position = position_dodge(0.8), linewidth = 0.8) +
    scale_color_manual("", values = v_colors) +
    scale_shape_manual("", values = v_shape) +
    scale_y_continuous(breaks = number_ticks(6)) +
    scale_x_discrete(breaks = as.character(v_times)) +
    ylab("Survival probability") + xlab("") +
    guides(color = guide_legend(ncol=1, byrow=TRUE)) +
    theme_bw(base_size = 18) +
    theme(legend.position = "bottom",
          strip.background = element_rect(color = "white", 
                                          fill = "white"),
          axis.text.x = element_text(angle = 0,
                                     vjust = 0.5,
                                     hjust = 0.5)) +
    coord_cartesian(ylim = c(0,1))
  
  # Plot Prevalence ("Prev")
  plotPrev <- ggplot(subset(df_model_target, abbrev_outcome == "Prev"), 
                     aes(x = as.factor(time), y = value, ymin = lb, ymax = ub,
                         color = type_coverage, shape = type_coverage)) +
    geom_point(size = 1.5, position = position_dodge(0.8)) +
    geom_errorbar(position = position_dodge(0.8), linewidth = 0.8) +
    scale_color_manual("", values = v_colors) +
    scale_shape_manual("", values = v_shape) +
    scale_y_continuous(breaks = number_ticks(6)) +
    scale_x_discrete(breaks = as.character(v_times)) +
    ylab("Prevalence") + xlab("Time (years)") +
    guides(color = guide_legend(ncol=1, byrow=TRUE)) +
    theme_bw(base_size = 18) +
    theme(legend.position = "bottom",
          strip.background = element_rect(color = "white", 
                                          fill = "white"),
          axis.text.x = element_text(angle = 0,
                                     vjust = 0.5,
                                     hjust = 0.5)) +
    coord_cartesian(ylim = c(0,1))
  
  # Plot Proportion who are Sick ("PropSick")
  plotPropSick <- ggplot(subset(df_model_target, abbrev_outcome == "PropSick"), 
                         aes(x = as.factor(time), y = value, ymin = lb, ymax = ub,
                             color = type_coverage, shape = type_coverage)) +
    geom_point(size = 1.5, position = position_dodge(0.8)) +
    geom_errorbar(position = position_dodge(0.8), linewidth = 0.8) +
    scale_color_manual("", values = v_colors) +
    scale_shape_manual("", values = v_shape) +
    scale_y_continuous(breaks = number_ticks(6)) +
    scale_x_discrete(breaks = as.character(v_times)) +
    ylab("Proportion of Sick") + xlab("") +
    guides(color = guide_legend(ncol=1, byrow=TRUE)) +
    theme_bw(base_size = 18) +
    theme(legend.position = "bottom",
          strip.background = element_rect(color = "white", 
                                          fill = "white"),
          axis.text.x = element_text(angle = 0,
                                     vjust = 0.5,
                                     hjust = 0.5)) +
    coord_cartesian(ylim = c(0,1))
  
  # Gather plots
  gg_model_out_vs_targets <- ggarrange(plotSurv, NULL, 
                                       plotPrev, NULL, 
                                       plotPropSick,
                                       heights = 1, 
                                       widths = c(1,0.05,1,0.05,1),
                                       ncol = 5, nrow = 1,
                                       common.legend = T,
                                       legend = "bottom")
  
  return(gg_model_out_vs_targets)
  
}


#' Plot residuals of a model vs observations
#'
#' @param Cfun modCost object obtained from \code{costFun}.
#' @param wgt_res Flag (default is TRUE) of whether weighted
#' residuals should be plotted.
#' 
#' @return A ggplot2 object.
#' 
#' @export
plot_costFun <- function(Cfun, wgt_res = TRUE){
  
  # Times to be plotted
  v_times <- sort(unique(df_targets$time))
  
  # Get residuals
  df_residuals <- Cfun$residuals
  
  # Add labels to outputs
  df_residuals$name_lab <- ""
  df_residuals$name_lab[df_residuals$name == "Surv"] <- "Survival probability"
  df_residuals$name_lab[df_residuals$name == "Prev"] <- "Prevalence"
  df_residuals$name_lab[df_residuals$name == "PropSick"] <- "Proportion of Sick"
  
  df_residuals$name_lab <- factor(df_residuals$name_lab,
                                  levels = c("Survival probability",
                                             "Prevalence",
                                             "Proportion of Sick"))
  
  if(wgt_res){
    
    plotCfun <- ggplot(df_residuals,
                       aes(x = x, y = res, color = name_lab, shape = name_lab)) +
      geom_point() +
      facet_wrap(~ name_lab, scales = "free_y") +
      scale_y_continuous(breaks = number_ticks(6)) +
      scale_x_continuous(breaks = v_times) +
      ylab("Weighted residuals") + xlab("Time (years)") +
      theme_bw(base_size = 18) +
      theme(legend.position = "none",
            strip.background = element_rect(color = "white", 
                                            fill = "white"),
            axis.text.x = element_text(angle = 0,
                                       vjust = 0.5,
                                       hjust = 0.5))
    
  }else{
    
    plotCfun <- ggplot(df_residuals,
                       aes(x = x, y = res.unweighted, color = name_lab, shape = name_lab)) +
      geom_point() +
      facet_wrap(~ name_lab, scales = "free_y") +
      scale_y_continuous(breaks = number_ticks(6)) +
      scale_x_continuous(breaks = v_times) +
      ylab("Residuals") + xlab("Time (years)") +
      theme_bw(base_size = 18) +
      theme(legend.position = "none",
            strip.background = element_rect(color = "white", 
                                            fill = "white"),
            axis.text.x = element_text(angle = 0,
                                       vjust = 0.5,
                                       hjust = 0.5))
  }
  
  return(plotCfun)
} 

#' Plot sensitivity function
#'
#' \code{plot_Sfun} plots the sensitivity function.
#'
#' @param Sfun Sensitivity function obtained from \code{sensFun}.
#' 
#' @return A ggplot2 object.
#' 
#' @export
plot_Sfun <- function(Sfun){
  
  # Times to be plotted
  v_times <- sort(unique(df_targets$time))
  
  # Get variables and target names
  v_abb_targets <- unique(Sfun$var)
  v_param_names <- colnames(Sfun)[3:length(Sfun)]
  n_params <- length(v_param_names)
  
  # Generate long data.frame for plots
  df_Sfun_long <- data.frame(NULL)
  for(abb_targ_i in v_abb_targets){
    
    df_temp <- subset(Sfun, var == abb_targ_i) %>%
      pivot_longer(cols = all_of(v_param_names),
                   names_to = "parm", 
                   values_to = "value")
    
    df_Sfun_long <- rbind(df_temp, df_Sfun_long)
    
  }
  df_plot <- df_Sfun_long %>% 
    rename(Parameter = parm)
  
  # Categorical variables
  # Parameters
  df_plot$Parameter <- factor(df_plot$Parameter,
                              levels = v_param_names,
                              labels = c("Probability to become Sicker when Sick",
                                         "Hazard ratio for death in Sick",
                                         "Hazard ratio for death in Sicker"),
                              ordered = TRUE)
  
  # Outcomes
  df_plot$var_label <- ""
  df_plot$var_label[df_plot$var == "Surv"] <- "Survival probability"
  df_plot$var_label[df_plot$var == "Prev"] <- "Prevalence"
  df_plot$var_label[df_plot$var == "PropSick"] <- "Proportion of Sick" 
  
  df_plot$var_label <- factor(df_plot$var_label,
                              levels = c("Survival probability",
                                         "Prevalence",
                                         "Proportion of Sick"))
  
  # Colors
  v_colors <- c("#56B4E9","#CCBFFF","#009E73")
  names(v_colors) <- c("Probability to become Sicker when Sick",
                       "Hazard ratio for death in Sick",
                       "Hazard ratio for death in Sicker")
  
  # Line type
  v_linetype <- c(1, 2, 3)
  names(v_linetype) <- c("Probability to become Sicker when Sick",
                         "Hazard ratio for death in Sick",
                         "Hazard ratio for death in Sicker")
  
  # Plot
  gg_Sfun <- ggplot(df_plot,
                    aes(x = x, y = value, color = Parameter, 
                        linetype = Parameter)) +
    geom_line(linewidth = 1.2) + #geom_point(size = 1.5) +
    facet_wrap( ~ var_label, nrow = 1, scales = "free_y") +
    scale_color_manual("", values = v_colors) +
    scale_linetype_manual("", values = v_linetype) +
    scale_y_continuous(breaks = number_ticks(6)) +
    scale_x_continuous(breaks = v_times) +
    ylab("Sensitivity") + xlab("Time (years)") +
    theme_bw(base_size = 18) +
    theme(legend.position = "bottom",
          legend.direction = "horizontal",
          strip.background = element_rect(color = "white", 
                                          fill = "white"),
          axis.text.x = element_text(angle = 0,
                                     vjust = 0.5,
                                     hjust = 0.5))
  
  return(gg_Sfun) 
}

#' Plot collinearity index 
#'
#' \code{plot_CollIndx} plots collinearity index
#'
#' @param df_collin Data.frame containing the collinearity indices 
#' from \code{collin}.
#' @param log_y Flag (default is TRUE) to whether transform y-axis to log.

#' @return A ggplot2 object.
#' 
#' @export
plot_CollIndx <- function(df_collin, log_y = TRUE){
  
  # Run coll
  collin_max <- max(df_collin$collinearity)
  
  # y-axis upper limit
  if(collin_max > 15){
    y_end <- collin_max + 3
  }else{
    y_end <- 18
  }
  
  # x-axis discrete values
  v_x_breaks <- unique(df_collin$N)
  
  # Order outcomes
  df_collin$Outcome <- factor(df_collin$Outcome,
                              levels = c("Survival probability",
                                         "Prevalence",
                                         "Proportion of Sick",
                                         "Survival probability & Prevalence",
                                         "Survival probability & Proportion of Sick",
                                         "Prevalence & Proportion of Sick",
                                         "Survival probability & Prevalence & Proportion of Sick"))
  
  if(log_y){
    
    plot_collindx <- ggplot(subset(df_collin, var != "All"), 
                            aes(x = factor(N), y = collinearity)) + 
      geom_point() +
      theme_bw() + 
      facet_wrap(~Outcome, ncol = 3) +
      scale_y_log10(n.breaks = 6) +
      scale_x_discrete(name = "Number of parameters",
                       breaks = paste0(v_x_breaks), 
                       labels  = paste0(v_x_breaks)) +
      ylab("Collinearity index") +
      geom_hline(yintercept = 15, 
                 col = "red",
                 linetype = "dashed") +
      theme_bw(base_size = 18) +
      theme(legend.position = "bottom",
            strip.background = element_rect(color = "white", 
                                            fill = "white"),
            strip.text = element_text(size = 14),
            axis.text.x = element_text(angle = 0,
                                       vjust = 0.5,
                                       hjust = 0.5),
            panel.grid.minor = element_line(color = "white"))  +
      coord_cartesian(ylim = c(NA, y_end))
    
    plot_collindx_all <- ggplot(subset(df_collin, var == "All"), 
                                aes(x = factor(N), y = collinearity)) + 
      geom_point() +
      theme_bw() + 
      facet_wrap(~Outcome, ncol = 1) +
      scale_y_log10(n.breaks = 6) +
      scale_x_discrete(name = "Number of parameters",
                       breaks = df_collin$N, 
                       labels  = df_collin$N) +
      ylab("Collinearity index") +
      geom_hline(yintercept = 15, 
                 col = "red",
                 linetype = "dashed") +
      theme_bw(base_size = 18) +
      theme(legend.position = "bottom",
            strip.background = element_rect(color = "white", 
                                            fill = "white"),
            strip.text = element_text(size = 14),
            axis.text.x = element_text(angle = 0,
                                       vjust = 0.5,
                                       hjust = 0.5),
            panel.grid.minor = element_line(color = "white"))  +
      coord_cartesian(ylim = c(NA, y_end))
    
  }else{
    
    plot_collindx <- ggplot(subset(df_collin, var != "All"), 
                            aes(x = factor(N), y = collinearity)) + 
      geom_point() +
      theme_bw() + 
      facet_wrap(~Outcome, ncol = 3) +
      scale_y_continuous(n.breaks = 7) +
      scale_x_discrete(name = "Number of parameters",
                       breaks = paste0(v_x_breaks), 
                       labels  = paste0(v_x_breaks)) +
      ylab("Collinearity index") +
      xlab("Number of parameters") +
      geom_hline(yintercept = 15, 
                 col = "red",
                 linetype = "dashed") +
      theme_bw(base_size = 18) +
      theme(legend.position = "bottom",
            strip.background = element_rect(color = "white", 
                                            fill = "white"),
            strip.text = element_text(size = 14),
            axis.text.x = element_text(angle = 0,
                                       vjust = 0.5,
                                       hjust = 0.5),
            panel.grid.minor = element_line(color = "white"))  +
      coord_cartesian(ylim = c(NA, y_end))
    
    plot_collindx_all <- ggplot(subset(df_collin, var == "All"), 
                                aes(x = factor(N), y = collinearity)) + 
      geom_point() +
      theme_bw() + 
      facet_wrap(~Outcome, ncol = 3) +
      scale_y_continuous(n.breaks = 7) +
      scale_x_discrete(name = "Number of parameters",
                       breaks = paste0(v_x_breaks), 
                       labels  = paste0(v_x_breaks)) +
      ylab("Collinearity index") +
      xlab("Number of parameters") +
      geom_hline(yintercept = 15, 
                 col = "red",
                 linetype = "dashed") +
      theme_bw(base_size = 18) +
      theme(legend.position = "bottom",
            strip.background = element_rect(color = "white", 
                                            fill = "white"),
            strip.text = element_text(size = 14),
            axis.text.x = element_text(angle = 0,
                                       vjust = 0.5,
                                       hjust = 0.5),
            panel.grid.minor = element_line(color = "white"))  +
      coord_cartesian(ylim = c(NA, y_end))
    
  }
  
  # Gather plots
  gg_collindx <- ggarrange(plot_collindx, NULL, 
                           plot_collindx_all,
                           heights = c(1,0.05,1), 
                           widths = 1,
                           ncol = 1, nrow = 3)
  
  return(gg_collindx)  
}

#' Plot global sensitivity analysis
#' 
#' \code{plot_GSA} plots Sobol' first-order and total sensitivity indices.
#'
#' @param df_sobol Data.frame containing the results from a GSA using function
#' \code{run_GSA}.
#'
#' @returns A ggplot object.
#' 
#' @export
plot_GSA <- function(df_sobol){
  
  # Outcomes
  v_outcomes <- unique(df_sobol$Outcome)
  n_outcomes <- length(v_outcomes)
  
  # Categorical variables
  # Parameters
  df_sobol$par <- factor(df_sobol$par,
                         levels = c("p_S1S2","hr_S1","hr_S2"),
                         labels = c("Probability to become Sicker when Sick",
                                    "Hazard ratio for death in Sick",
                                    "Hazard ratio for death in Sicker"),
                         ordered = TRUE)
  
  # Empty list to store plots
  l_GSA_plot <- vector(mode = "list", length = n_outcomes)
  i <- 1 # counter
  for(outcome_i in v_outcomes){ # outcome_i = "Prevalence"
    
    # Plot GSA by outcome
    l_GSA_plot[[i]] <- ggplot(subset(df_sobol, Outcome == outcome_i),
                              aes(x = as.factor(x), y = original, ymin = min..c.i., ymax = max..c.i.,
                                  color = type, shape = type)) +
      geom_point(position = position_dodge(width = 0.5), size = 3) +
      geom_errorbar(position = position_dodge(width = 0.5), 
                    linewidth = 1.05, 
                    width = 0.25) +
      scale_y_continuous(breaks = seq(0,1,0.25)) +
      ylab("Sobol index") + xlab("Time (years)") +
      scale_color_manual("", values = c("Main effect" = "#CD0BBC", 
                                        "Total effect" = "#999999")) +
      scale_shape_manual("", values = c("Main effect" = 1, 
                                        "Total effect" = 2)) +
      facet_wrap(~par, ncol = 1) +
      labs(subtitle = outcome_i) +
      theme_bw(base_size = 18) +
      theme(legend.position = "bottom",
            plot.subtitle = element_text(vjust = 0.5,
                                         hjust = 0.5),
            strip.text = element_text(size = 13),
            strip.background = element_rect(color = "white",
                                            fill = "white"),
            axis.text.x = element_text(angle = 0,
                                       vjust = 0.5,
                                       hjust = 0.5)) +
      coord_cartesian(ylim = c(0,1))
    
    # Increase counter
    i <- i+1
    
  }
  
  # Return plot
  gg_GSA <- ggarrange(l_GSA_plot[[1]] + xlab(""), NULL,
                      l_GSA_plot[[2]] + ylab(""), NULL,
                      l_GSA_plot[[3]] + xlab("") + ylab(""),
                      ncol = 5, nrow = 1,
                      widths = c(1,0.05,1,0.05,1),
                      heights = c(1), 
                      common.legend = T,
                      legend = "bottom",
                      label.y = "Sobol index",
                      label.x = "Time (years)")
  
  return(gg_GSA)
}

# Calibration functions ---------------------------------------------------

#' Generate model outputs for calibration from a parameter set
#'
#' \code{calibration_out} computes model outputs for the Sick-Sicker model
#' to be used for calibration routines.
#'
#' @param v_params Vector containing the calibrated parameters.
#' 
#' @return A data.frame containing the model-predicted outcomes.
#' 
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
    filter(time %in% v_times)
  
  # Prevalence
  df_Prev <- data.frame(time = df_model_out$time,
                        Outcome = "Prevalence",
                        abbrev_outcome = "Prev",
                        type = "Model",
                        value = df_model_out$Prev) %>% 
    filter(time %in% v_times)
  
  # Proportion who are Sick ("PropSick")
  # among all those afflicted (Sick+Sicker)
  df_PropSick <- data.frame(time = df_model_out$time,
                            Outcome = "Proportion of sick",
                            abbrev_outcome = "PropSick",
                            type = "Model", 
                            value = df_model_out$PropSick) %>% 
    filter(time %in% v_times)
  
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
#' 
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
#' @param v_params Vector (or matrix) containing the model parameters.
#' @param v_lb Vector containing the lower bounds for each parameter.
#' @param v_ub Vector containing the upper bounds for each parameter.
#' 
#' @return A scalar (or vector) containing the log-prior values.
#' 
#' @export
log_prior <- function(v_params, v_lb, v_ub){
  
  # Get param names
  v_param_names <- names(v_lb)
  
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
#' @param v_params Vector (or matrix) containing the model parameters.
#' @param v_lb Vector containing the with lower bounds for each parameter.
#' @param v_ub Vector containing the with upper bounds for each parameter.
#'
#' @return A scalar (or vector) containing the prior values.
#' 
#' @export
prior <- function(v_params, v_lb, v_ub) { 
  return(exp(log_prior(v_params, v_lb, v_ub)))
}

#' Log likelihood normal distribution
#'
#' \code{log_lik} computes a log-likelihood value for one (or multiple)
#' parameter set(s).
#'
#' @param v_params_calib Vector (or matrix) containing the calibrated parameters.
#' 
#' @return A scalar (or vector) containing the log-likelihood values.
#' 
#' @export
log_lik <- function(v_params){
  
  # par_vector: a vector (or matrix) of model parameters 
  if(is.null(dim(v_params))) { # If vector, change to matrix
    v_params <- t(v_params) 
  }
  n_samp <- nrow(v_params)
  v_abbrev_outcomes <- c("Surv", "Prev", "PropSick")
  n_outcomes <- length(v_abbrev_outcomes)
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
#' @param v_params_calib Vector (or matrix) containing the model parameters.
#'
#' @return A scalar (or vector) containing the likelihood values.
#' 
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
#' @param v_params_calib Vector (or matrix) containing the model parameters.
#' @param v_lb Vector containing the lower bounds for each parameter.
#' @param v_ub Vector containing the upper bounds for each parameter. 
#' 
#' @return A scalar (or vector) containing the log-posterior values.
#' 
#' @export
log_post <- function(v_params, v_lb, v_ub) { 
  lpost <- log_prior(v_params, v_lb, v_ub) + log_lik(v_params)
  return(lpost) 
}

#' Evaluate posterior of calibrated parameters
#'
#' \code{posterior} computes a posterior value for one (or multiple) 
#' parameter set(s).
#' 
#' @param v_params_calib Vector (or matrix) containing the model parameters.
#' @param v_lb Vector containing the lower bounds for each parameter.
#' @param v_ub Vector containing the upper bounds for each parameter.
#' 
#' @return A scalar (or vector) containing the posterior values.
#' 
#' @export
posterior <- function(v_params, v_lb, v_ub) { 
  exp(log_post(v_params, v_lb, v_ub)) 
}

