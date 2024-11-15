#' Plot model outputs vs targets
#'
#' \code{plot_model_out_vs_targets} plots model outputs vs targets.
#' 
#' @param df_targets Data.frame with calibration targets.
#' 
#' @return A ggplot2 object.
#' @export
plot_targets <- function(df_targets){
  
  # Bind data.frames
  df_Surv <- df_targets %>% filter(Outcome == "Surv")
  df_Prev <- df_targets %>% filter(Outcome == "Prev")
  df_PropSick <- df_targets %>% filter(Outcome == "PropSick")
  
  # Plot Survival ("Surv")
  plotSurv <- ggplot(df_Surv, 
                     aes(x = time, y = value, ymin = lb, ymax = ub)) +
    geom_point(shape = 1) + geom_errorbar() +
    scale_y_continuous(breaks = number_ticks(6)) +
    scale_x_continuous(breaks = number_ticks(6)) +
    ylab("Survival Pr.") + xlab("Time") +
    theme_bw(base_size = 16) +
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
    ylab("Prevalence") + xlab("Time") +
    theme_bw(base_size = 16) +
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
    ylab("Prop. Sick") + xlab("Time") +
    theme_bw(base_size = 16) +
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
#' @param df_targets Data.frame with calibration targets.
#' @param df_model Data.frame with decision model outputs.
#' 
#' @return A ggplot2 object.
#' @export
plot_model_out_vs_targets <- function(df_targets, df_model){
  
  # Bind data.frames
  df_Surv <- rbind(df_targets %>% filter(Outcome == "Surv"),
                   df_model %>% filter(Outcome == "Surv"))
  
  df_Prev <- rbind(df_targets %>% filter(Outcome == "Prev"),
                   df_model %>% filter(Outcome == "Prev"))
  
  df_PropSick <- rbind(df_targets %>% filter(Outcome == "PropSick"),
                       df_model %>% filter(Outcome == "PropSick"))
  
  # Plot Survival ("Surv")
  plotSurv <- ggplot(df_Surv, 
         aes(x = time, y = value, ymin = lb, ymax = ub,
             color = type, shape = type)) +
    geom_point(size = 2) + geom_errorbar() +
    scale_shape_manual("", values = c(Target = 1, Model = 8)) +
    scale_color_manual("", values = c(Target = "red", Model = "black")) +
    scale_y_continuous(breaks = number_ticks(6)) +
    scale_x_continuous(breaks = number_ticks(6)) +
    ylab("Survival Pr.") + xlab("Time") +
    theme_bw(base_size = 16) +
    theme(legend.position = "top",
          strip.background = element_rect(color = "white", 
                                          fill = "white"),
          axis.text.x = element_text(angle = 0,
                                     vjust = 0.5,
                                     hjust = 0.5))
  
  # Plot Prevalence ("Prev")
  plotPrev <- ggplot(df_Prev, 
         aes(x = time, y = value, ymin = lb, ymax = ub,
             color = type, shape = type)) +
    geom_point(size = 2) + geom_errorbar() +
    scale_shape_manual("", values = c(Target = 1, Model = 8)) +
    scale_color_manual("", values = c(Target = "red", Model = "black")) +
    scale_y_continuous(breaks = number_ticks(6)) +
    scale_x_continuous(breaks = number_ticks(6)) +
    ylab("Prevalence") + xlab("Time") +
    theme_bw(base_size = 16) +
    theme(legend.position = "top",
          strip.background = element_rect(color = "white", 
                                          fill = "white"),
          axis.text.x = element_text(angle = 0,
                                     vjust = 0.5,
                                     hjust = 0.5))
  
  # Plot Proportion who are Sick ("PropSick")
  plotPropSick <- ggplot(df_PropSick, 
                      aes(x = time, y = value, ymin = lb, ymax = ub,
                          color = type, shape = type)) +
    geom_point(size = 2) + geom_errorbar() +
    scale_shape_manual("", values = c(Target = 1, Model = 8)) +
    scale_color_manual("", values = c(Target = "red", Model = "black")) +
    scale_y_continuous(breaks = number_ticks(6)) +
    scale_x_continuous(breaks = number_ticks(6)) +
    ylab("Prop. Sick") + xlab("Time") +
    theme_bw(base_size = 16) +
    theme(legend.position = "top",
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
                                       common.legend = T)
    
  return(gg_model_out_vs_targets)

}


#' Plot residuals of a model vs observations
#'
#' @param Cfun modCost object obtained from \code{costFun}.
#' @param wgt_res Flag (default is TRUE) of whether weighted
#' residuals should be plotted.
#' 
#' @return A ggplot2 object.
#' @export
plot_costFun <- function(Cfun, wgt_res = TRUE){
  
  df_residuals <- Cfun$residuals
  
  df_residuals$name <- factor(df_residuals$name,
                              levels = c("Surv","Prev","PropSick"))
  
  if(wgt_res){
    
    plotCfun <- ggplot(df_residuals,
           aes(x = x, y = weight, color = name, shape = name)) +
      geom_point() +
      geom_line() +
      facet_wrap(~ name, scales = "free_y") +
      scale_y_continuous(breaks = number_ticks(6)) +
      scale_x_continuous(breaks = number_ticks(6)) +
      ylab("Weighted residuals") + xlab("Time") +
      theme_bw(base_size = 16) +
      theme(legend.position = "none",
            strip.background = element_rect(color = "white", 
                                            fill = "white"),
            axis.text.x = element_text(angle = 0,
                                       vjust = 0.5,
                                       hjust = 0.5))
    
  }else{
    
    plotCfun <- ggplot(df_residuals,
           aes(x = x, y = res.unweighted, color = name, shape = name)) +
      geom_point() +
      geom_line() +
      facet_wrap(~ name, scales = "free_y") +
      scale_y_continuous(breaks = number_ticks(6)) +
      scale_x_continuous(breaks = number_ticks(6)) +
      ylab("Weighted residuals") + xlab("Time") +
      theme_bw(base_size = 16) +
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
  
  # Plot
  df_plot$Parameter <- factor(df_plot$Parameter,
                         levels = v_param_names,
                         labels = c("Prob. to become sicker when sick",
                                    "Hazard ratio of death in sick vs healthy",
                                    "Hazard ratio of death in sicker vs healthy"),
                         ordered = TRUE)
  
  df_plot$var <- factor(df_plot$var,
                        levels = c("Surv","Prev","PropSick"))
  
  gg_Sfun <- ggplot(df_plot,
                    aes(x = x, y = value, color = Parameter)) +
    geom_line(linewidth = 1.3) +
    facet_wrap( ~ var, nrow = 1, scales = "free_y") +
    scale_y_continuous(breaks = number_ticks(6)) +
    scale_x_continuous(breaks = number_ticks(6)) +
    ylab("Sensitivity") + xlab("Time") +
    labs(title = "") +
    guides(color = guide_legend(nrow=3, byrow=TRUE)) +
    theme_bw(base_size = 16) +
    theme(legend.position = "bottom",
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
#' @param ident Data.frame containing collinearity index from all combination
#' of parameters obtained from \code{collin}.
#' @param log_y Flag (default is TRUE) to whether transform y-axis to log.

#' @return A ggplot2 object.
#' @export
plot_CollIndx <- function(ident, log_y = TRUE){
  
  col_max <- max(ident$collinearity)
  
  if(col_max > 15){
    y_end <- col_max + 3
  }else{
    y_end <- 18
  }

  if(log_y){
    
    plot_collindx <- ggplot(ident, aes(x = N, y = collinearity)) + 
      geom_point() +
      theme_bw() +
      scale_y_log10(n.breaks = 6) +
      scale_x_continuous(name = "Number of parameters",
                          breaks = ident$N, 
                         labels  = ident$N) +
      ylab("Collinearity index") +
      geom_hline(yintercept = 15, 
                 col = "red",
                 linetype = "dashed") +
      theme_bw(base_size = 16) +
      theme(legend.position = "bottom",
            strip.background = element_rect(color = "white", 
                                            fill = "white"),
            axis.text.x = element_text(angle = 0,
                                       vjust = 0.5,
                                       hjust = 0.5))  +
      coord_cartesian(ylim = c(NA, y_end))
    
  }else{
    
    plot_collindx <- ggplot(ident, aes(x = N, y = collinearity)) + 
      geom_point() +
      theme_bw() +
      scale_y_continuous(n.breaks = 7) +
      scale_x_continuous(name = "Number of parameters",
                         breaks = ident$N, 
                         labels  = ident$N) +
      ylab("Collinearity index") +
      xlab("Number of parameters") +
      geom_hline(yintercept = 15, 
                 col = "red",
                 linetype = "dashed") +
      theme_bw(base_size = 16) +
      theme(legend.position = "bottom",
            strip.background = element_rect(color = "white", 
                                            fill = "white"),
            axis.text.x = element_text(angle = 0,
                                       vjust = 0.5,
                                       hjust = 0.5))  +
      coord_cartesian(ylim = c(NA, y_end))
    
  }
  
  return(plot_collindx)  
}