
# Model-decision functions ------------------------------------------------

#' Load all parameter sets
#'
#' \code{load_params_all} generates all the values of the SIR model parameters.
#'
#' @param n_t Time horizon in days.
#' @param time_step Model evaluations per day.
#' @param n_date_init Initial day on which the infections began.
#' @param n_inf_init Number of individuals initially infected.
#' @param n_tot_pop Total number of the population.
#' @param r_beta Transmission rate.
#' @param r_gamma Recovery rate.
#'
#' @returns List containing all the parameters of the decision model.
#' 
#' @export
load_params_all <- function(n_t,
                            time_step,
                            n_date_init,
                            n_inf_init,
                            n_tot_pop,
                            r_beta,
                            r_gamma){
  
  # Vector of times
  v_times <- seq(0,n_t,time_step)
  
  # Initial state
  S_0 <- (n_tot_pop - n_inf_init)/n_tot_pop
  I_0 <- n_inf_init/n_tot_pop 
  R_0 <- 0
  
  v_states_init <- c(S = S_0,
                     I = I_0,
                     R = R_0)
  v_states_names <- names(v_states_init)
  
  return(list(n_t            = n_t,
              time_step      = time_step,
              v_times        = v_times,
              n_date_init    = n_date_init,
              v_states_init  = v_states_init,
              v_states_names = v_states_names,
              r_beta         = r_beta,
              r_gamma        = r_gamma))
  
}

#' Derivatives of the SIR model
#' 
#' \code{SIR_dXdt} computes the derivatives of a Susceptible-Infectious-Recovered 
#' (SIR) transmission model.
#' 
#' @param time Time at which the derivatives will be computed.
#' @param v_pop Vector containing the population by health state
#' @param l_params_all List containing all the parameters of the decision model.
#' 
#' @return 
#' A list containing a vector of derivatives for each health state.
#' 
#' @export
SIR_dXdt <- function(time, v_pop, l_params_all){
  
  # Getting parameters
  r_beta  <- l_params_all$r_beta
  r_gamma <- l_params_all$r_gamma

  # Population
  S <- v_pop["S"]
  I <- v_pop["I"]
  R <- v_pop["R"]
  
  # ODE model
  dS <- -r_beta*S*I
  dI <- r_beta*S*I - r_gamma*I
  dR <- r_gamma*I
  
  return(list(c(dS,dI,dR)))
  
}

#' SIR model
#'
#' \code{SIR_model} implements a Susceptible-Infectious-Recovered (SIR) model.
#' 
#' @param l_params_all List containing all the parameters of the decision model.
#' 
#' @return 
#' A data.frame containing the population for each health state over time and 
#' a list containing all the parameters of decision model.
#' 
#' @export
SIR_model <- function(l_params_all){
  df_out_model <- lsoda(y     = l_params_all$v_states_init, 
                        times = l_params_all$v_times, 
                        func  = SIR_dXdt, 
                        parms = l_params_all)
  
  l_out_model <- list(df_out_model = as.data.frame(df_out_model),
                      l_params_all = l_params_all)
  
  return(l_out_model)
}

#' Total incident infections
#' 
#' \code{calc_dxinc_totals} calculate total number of incident diagnosed 
#' infections over time.
#' 
#' @param l_out_model List containing the model output.
#' @param l_params_all List containing all the parameters of the decision model.
#' 
#' @return 
#' A data.frame containing the total number of incident diagnosed infections 
#' over time.
#' 
#' @export
calc_DXIncTot <- function(l_out_model, l_params_all){
  
  # Get variables
  n_date_init <- l_params_all$n_date_init
  
  # Get model output
  df_out_model <- l_out_model$df_out_model
  
  df_DXIncTot <- data.frame(Type = "Model",
                            abbrev_outcome = "DXIncTot",
                            Outcome = "Incident confirmed cases",
                            time = df_out_model$time, 
                            Date = n_date_init + df_out_model$time,
                            value = df_out_model$I*n_tot_pop,
                            lb = NA,
                            ub = NA,
                            se = NA,
                            check.names = FALSE) %>% 
    mutate(Date0 = Date - Date[1],
           Day = time + 1) %>% 
    select(Type, abbrev_outcome, Outcome, Day, Date, Date0, value, lb, ub, se)
  
  return(df_DXIncTot)
}

# Visualization functions -------------------------------------------------

#' Plot model outputs vs targets
#'
#' \code{plot_model_out_vs_targets} plots model outputs vs targets.
#' 
#' @param df_targets Data.frame with calibration targets.
#' 
#' @return A ggplot2 object.
#' @export
plot_targets <- function(df_targets){
  
  # Plot targets
  gg_targets <- ggplot(df_targets, 
                     aes(x = Day, y = value, ymin = lb, ymax = ub)) +
    geom_point(shape = 1) + geom_errorbar() +
    scale_y_continuous(breaks = number_ticks(6)) +
    scale_x_continuous(breaks = number_ticks(10)) +
    # scale_x_date(date_labels = "%m/%d",
    #              breaks = number_ticks(10)) +
    ylab("Number of infections") + xlab("Time (days)") +
    theme_bw(base_size = 18) +
    theme(legend.position = "top",
          strip.background = element_rect(color = "white", 
                                          fill = "white"),
          axis.text.x = element_text(angle = 0,
                                     vjust = 0.5,
                                     hjust = 0.5)) + 
    coord_cartesian(ylim = c(0,800))
  
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
  
  # Bind data.frames
  df_plot <- rbind(
    df_targets %>% 
      select(Type, Outcome, Day, Date, Date0, value, lb, ub, se),
    df_model %>%
      select(Type, Outcome, Day, Date, Date0, value, lb, ub, se)
  )
  
  # Vector of dates
  v_dates <- as.Date(unique(df_plot$Date))
  
  # Order
  df_plot$Type <- factor(df_plot$Type, levels = c("Target","Model"))
  
  # Plot
  gg_model_out_vs_targets <- ggplot(df_plot, 
                                    aes(x = Day, y = value, ymin = lb, ymax = ub,
                                        color = Type, shape = Type)) +
    geom_point(size = 1.6) +
    geom_errorbar(linewidth = 0.8) +
    scale_shape_manual(NULL, values = c(Target = 1, 
                                        Model = 8)) +
    scale_color_manual(NULL, values = c(Target = "black", 
                                        Model = "red")) +
    scale_y_continuous(breaks = number_ticks(6)) +
    scale_x_continuous(breaks = number_ticks(10)) +
    # scale_x_date(date_labels = "%m/%d",
    #              breaks = number_ticks(10)) +
    ylab("Number of infections") + xlab("Time (days)") +
    theme_bw(base_size = 18) +
    theme(legend.position = "inside",
          legend.position.inside = c(0.15,0.85),
          strip.background = element_rect(color = "white", 
                                          fill = "white"),
          axis.text.x = element_text(angle = 0,
                                     vjust = 0.5,
                                     hjust = 0.5)) + 
    coord_cartesian(ylim = c(0,500))
   
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
  
  # Add type_coverage to targets
  df_targets <- df_targets %>% 
    mutate(type_coverage = case_when(Type == "Target" ~ "Calibration target"))
  
  # Get levels
  v_levels <- c("Calibration target", unique(df_coverage$type_coverage)) 
  
  # Bind data.frames
  df_model_target <- rbind(df_targets %>% 
                             select(Type, Outcome, type_coverage, Day, 
                                    Date, Date0, value, lb, ub, se),
                           df_coverage %>%
                             select(Type, Outcome, type_coverage, Day, 
                                    Date, Date0, value, lb, ub, se))
  df_model_target$type_coverage <- factor(df_model_target$type_coverage, 
                                          levels = v_levels)
  
  # Vector of shapes
  v_shape <- c(1,1,rep(NA,length(v_levels)-2))
  names(v_shape) <- v_levels
  
  # Vector of colors
  v_colors <- v_colors <- c("black","red","#FF8000","#56B4E9","#009E73","#CCBFFF","tan")
  names(v_colors) <- c(
    "Calibration target",
    "Calibration target sampled",
    "Initial model coverage",                                                           
    "A: Improved model coverage by resizing the transmission rate",
    "B: A + improved model coverage by resizing the recovery rate",  
    "C: A + B + improved model coverage by resizing the transmission rate",
    "D: A + B + C + improved model coverage by resizing the recovery rate"
  )
  
  # Plot Survival ("Surv")
  gg_coverage <- ggplot(df_model_target, 
                     aes(x = Day, 
                         y = value, ymin = lb, ymax = ub,
                         color = type_coverage, shape = type_coverage)) +
    geom_point(size = 1.5, position = position_dodge(0.8)) +
    geom_errorbar(position = position_dodge(0.8), linewidth = 0.8) +
    scale_color_manual(NULL, values = v_colors) +
    scale_shape_manual(NULL, values = v_shape) +
    scale_y_continuous(breaks = number_ticks(6)) +
    scale_x_continuous(breaks = number_ticks(10)) +
    # scale_x_date(date_labels = "%m/%d",
    #              breaks = number_ticks(10)) +
    ylab("Number of infections") + xlab("Time (days)") +
    guides(color = guide_legend(ncol=1, byrow=TRUE)) +
    theme_bw(base_size = 18) +
    theme(legend.position = "inside",
          legend.position.inside = c(0.7,0.8),
          strip.background = element_rect(color = "white", 
                                          fill = "white"),
          axis.text.x = element_text(angle = 0,
                                     vjust = 0.5,
                                     hjust = 0.5))  + 
    coord_cartesian(ylim = c(0,820))
  
  return(gg_coverage)
  
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
  
  df_residuals <- Cfun$residuals
  df_residuals$name_lab <- "Number of infections"
  
  if(wgt_res){
    
    plotCfun <- ggplot(df_residuals,
                       aes(x = x, y = res, color = name_lab, shape = name_lab)) +
      geom_point() +
      facet_wrap(~ name_lab, scales = "free_y") +
      scale_y_continuous(breaks = number_ticks(6)) +
      scale_x_continuous(breaks = number_ticks(6)) +
      ylab("Weighted residuals") + xlab("Time (days)") +
      theme_bw(base_size = 18) +
      theme(legend.position = "none",
            strip.background = element_rect(color = "white", 
                                            fill = "white"),
            axis.text.x = element_text(angle = 0,
                                       vjust = 0.5,
                                       hjust = 0.5))
    
  }else{
    
    plotCfun <- ggplot(df_residuals,
                       aes(x = x, y = res.unweighted, 
                           color = name_lab, shape = name_lab)) +
      geom_point() +
      facet_wrap(~ name_lab, scales = "free_y") +
      scale_y_continuous(breaks = number_ticks(6)) +
      scale_x_continuous(breaks = number_ticks(6)) +
      ylab("Residuals") + xlab("Time (days)") +
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
                              labels = c("Transmission rate",
                                         "Recovery rate"),
                              ordered = TRUE)
  
  # Outcomes
  df_plot$var_label <- "Number of infections"
  
  # Colors
  v_colors <- c("#56B4E9","#009E73")
  names(v_colors) <- c("Transmission rate",
                       "Recovery rate")
  
  # Line type
  v_linetype <- c(1, 3)
  names(v_linetype) <- c("Transmission rate",
                         "Recovery rate")
  
  # Plot
  gg_Sfun <- ggplot(df_plot,
                    aes(x = x, y = value, color = Parameter, 
                        linetype = Parameter)) +
    geom_line(linewidth = 1.2) +
    scale_color_manual("", values = v_colors) +
    scale_linetype_manual("", values = v_linetype) +
    scale_y_continuous(breaks = number_ticks(6)) +
    scale_x_continuous(breaks = number_ticks(6)) +
    ylab("Sensitivity") + xlab("Time (days)") +
    theme_bw(base_size = 18) +
    theme(legend.position = "inside",
          legend.position.inside = c(0.78,0.85),
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
#' @param df_collin Data.frame containing collinearity indices from \code{collin}.
#' @param flag_each Flag (default is FALSE) to whether individual outcome 
#' combinations is plotted.
#' @param log_y Flag (default is TRUE) to whether transform y-axis to log.

#' @return A ggplot2 object.
#' 
#' @export
plot_CollIndx <- function(df_collin, flag_each = FALSE, log_y = TRUE){
  
  # Run coll
  collin_max <- max(df_collin$collinearity)
  
  # y-axis upper limit
  if(collin_max > 15){
    y_end <- collin_max + 3
  }else{
    y_end <- 18
  }
  
  # Vars
  v_vars <- unique(df_collin$var)
  n_vars <- length(v_vars)
  
  # x-axis discrete values
  v_x_breaks <- unique(df_collin$N)
  
  if(log_y){
    if(flag_each){
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
              strip.text = element_text(size = 11),
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
              strip.text = element_text(size = 11),
              axis.text.x = element_text(angle = 0,
                                         vjust = 0.5,
                                         hjust = 0.5),
              panel.grid.minor = element_line(color = "white"))  +
        coord_cartesian(ylim = c(NA, y_end))
    }else{
      
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
              strip.text = element_text(size = 11),
              axis.text.x = element_text(angle = 0,
                                         vjust = 0.5,
                                         hjust = 0.5),
              panel.grid.minor = element_line(color = "white"))  +
        coord_cartesian(ylim = c(NA, y_end))
      
    } 
    
  }else{
    
    if(flag_each){
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
              strip.text = element_text(size = 11),
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
              strip.text = element_text(size = 11),
              axis.text.x = element_text(angle = 0,
                                         vjust = 0.5,
                                         hjust = 0.5),
              panel.grid.minor = element_line(color = "white"))  +
        coord_cartesian(ylim = c(NA, y_end))
    }else{
      
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
              strip.text = element_text(size = 11),
              axis.text.x = element_text(angle = 0,
                                         vjust = 0.5,
                                         hjust = 0.5),
              panel.grid.minor = element_line(color = "white"))  +
        coord_cartesian(ylim = c(NA, y_end))
      
    }
    
  }
  
  if(flag_each){
    
    # Gather plots
    gg_collindx <- ggarrange(plot_collindx, NULL, 
                             plot_collindx_all,
                             heights = c(1,0.05,1), 
                             widths = 1,
                             ncol = 1, nrow = 3)
    
  }else{
    gg_collindx <- plot_collindx_all
  }
  
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
  
  # Label variables
  df_sobol$par <- factor(df_sobol$par, 
                         levels = c("r_beta","r_gamma"),
                         labels = c("Transmission rate","Recovery rate"))
  
  # Plot
  ggGSA <- ggplot(df_sobol,
         aes(x = x, y = original, ymin = min..c.i., ymax = max..c.i.,
             color = type, shape = type)) +
    geom_point(position = position_dodge(width = 0.5), size = 3) +
    geom_errorbar(position = position_dodge(width = 0.5), 
                  linewidth = 1.05, 
                  width = 0.25) +
    ylab("Sobol index") + xlab("Time (days)") +
    scale_x_continuous(breaks = number_ticks(6)) +
    scale_y_continuous(breaks = seq(0,1,0.25)) +
    scale_color_manual("", values = c("Main effect" = "#CD0BBC", 
                                      "Total effect" = "#999999")) +
    scale_shape_manual("", values = c("Main effect" = 1, 
                                      "Total effect" = 2)) +
    facet_wrap(~par, ncol = 1) +
    theme_bw(base_size = 18) +
    theme(legend.position = "inside",
          legend.position.inside = c(0.83,0.95),
          legend.direction = "horizontal",
          strip.background = element_rect(color = "white", 
                                          fill = "white"),
          axis.text.x = element_text(angle = 0,
                                     vjust = 0.5,
                                     hjust = 0.5)) +
    coord_cartesian(ylim = c(0,1))
  
  return(ggGSA)
  
}

# Calibration functions ---------------------------------------------------
#' Generate model outputs for calibration from a parameter set
#'
#' \code{SIR_calibration_out} computes model outputs for the SIR transmission
#' model to be used for calibration routines.
#' 
#' @param v_params Vector (or matrix) containing the model parameters.
#' @param l_params_all List containing all the parameters of the decision model.
#' 
#' @return A data.frame containing the model-predicted outcomes.
#' 
#' @export
SIR_calibration_out <- function(v_params, l_params_all){
  
  # Update list of parameters with calibrated parameter set
  l_params_all$r_beta  <- v_params["r_beta"]
  l_params_all$r_gamma <- v_params["r_gamma"]
  
  # Run model
  l_out_model  <- SIR_model(l_params_all)
  
  # Calculate total infections
  df_DXIncTot <- calc_DXIncTot(l_out_model  = l_out_model, 
                               l_params_all = l_params_all)
  
  # Return data.frame
  return(list(DXIncTot = df_DXIncTot))
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
  
  # Number of parameters
  n_param <- length(v_lb) 
  
  # Create random Latin hypercube design
  m_lhs_unit   <- randomLHS(n = n_samp, k = n_param)
  
  # Sample
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
#' @param v_lb Vector containing the lower bounds for each parameter.
#' @param v_ub Vector containing the upper bounds for each parameter.
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
#' @param v_dates Vector containing the dates of the available data to perform 
#' calibration.
#' @param ... Further arguments to be passed to.
#'
#' @return A scalar (or vector) containing the log-likelihood values.
#' 
#' @export
log_lik <- function(v_params, v_dates, ...){
  
  # par_vector: a vector (or matrix) of model parameters 
  if(is.null(dim(v_params))) { # If vector, change to matrix
    v_params <- t(v_params) 
  }
  
  # Set variables
  n_samp         <- nrow(v_params)
  v_target_names <- c("DXIncTot")
  n_targets      <- length(v_target_names)
  v_llik         <- matrix(0, nrow = n_samp, ncol = n_targets) 
  v_llik_overall <- numeric(n_samp)
  
  for(j in 1:n_samp) { # j=1
    jj <- tryCatch( { 
      ###   Run model for parameter set "v_params" ###
      l_model_out <- SIR_calibration_out(v_params[j, ], 
                                         # l_params_all = l_params_all)
                                         ...)
      ## Model output
      df_model_DXIncTot <- l_model_out$DXIncTot %>% 
        filter(Date %in% v_dates)
      
      ## Targets
      df_targets_DXIncTot <- df_targets %>% 
        filter(abbrev_outcome == "DXIncTot" &
                 Date %in% v_dates)
      
      ###  Calculate log-likelihood of model outputs to targets  ###
      # TARGET 1: Number of infections (DXIncTot)
      ## Negative-Binomial log-likelihood
      v_llik[j, 1] <- sum(dnbinom(x = df_targets_DXIncTot$value,
                                  size = 5,#1,
                                  mu = df_model_DXIncTot$value,
                                  log = T), na.rm = T)
      # ## Poisson log-likelihood
      # v_llik[j, 1] <- sum(dpois(x = df_targets_DXIncTot$value,
      #                           lambda = df_model_DXIncTot$value,
      #                           log = T), na.rm = T)
      
      ## OVERALL ADD TO SIR
      ## can give different targets different weights (user must change this)
      v_weights <- rep(1, n_targets)
      ## weighted sum
      v_llik_overall[j] <- v_llik[j, ] %*% v_weights
      
    }, error = function(e) NA) 
    if (is.na(jj)) { v_llik_overall[j] <- -Inf }
  } # End loop over sampled parameter sets
  
  # return LLIK
  return(v_llik_overall)
}

#' Likelihood
#'
#' \code{likelihood} computes a likelihood value for one (or multiple)
#' parameter set(s).
#'
#' @param v_params_calib Vector (or matrix) containing the model parameters.
#' @param v_dates Vector containing the dates of the available data
#' to perform calibration.
#' @param ... Further arguments to be passed to.
#' 
#' @return A scalar (or vector) containing the likelihood values.
#' 
#' @export
likelihood <- function(v_params, v_dates, ...){ 
  return(exp(log_lik(v_params, v_dates, ...))) 
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
#' @param v_dates Vector containing the dates of the available data
#' to perform calibration.
#' @param ... Further arguments to be passed to.
#' 
#' @return A scalar (or vector) with log-posterior values.
#' 
#' @export
log_post <- function(v_params, v_lb, v_ub, v_dates, ...) { 
  lpost <- log_prior(v_params, v_lb, v_ub) + log_lik(v_params, v_dates, ...)
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
#' @param ... Further arguments to be passed to.
#' 
#' @return A scalar (or vector) containing the posterior values.
#' 
#' @export
posterior <- function(v_params, v_lb, v_ub, ...) { 
  exp(log_post(v_params, v_lb, v_ub, ...)) 
}

