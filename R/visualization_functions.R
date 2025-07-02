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
  df_Surv <- df_targets %>% filter(abbrev_outcome == "Surv" & time %in% c(10,20,30))
  df_Prev <- df_targets %>% filter(abbrev_outcome == "Prev" & time %in% c(10,20,30))
  df_PropSick <- df_targets %>% filter(abbrev_outcome == "PropSick")
  
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
  df_Surv <- rbind(df_targets, df_model) %>% 
    filter(abbrev_outcome == "Surv" & time %in% c(10,20,30))
  
  df_Prev <- rbind(df_targets, df_model) %>% 
    filter(abbrev_outcome == "Prev" & time %in% c(10,20,30))
  
  df_PropSick <- rbind(df_targets, df_model) %>% 
    filter(abbrev_outcome == "PropSick" & time %in% c(10,20,30))
  
  # Plot Survival ("Surv")
  plotSurv <- ggplot(df_Surv, 
         aes(x = time, y = value, ymin = lb, ymax = ub,
             color = type, shape = type)) +
    geom_point(size = 2) + geom_errorbar() +
    scale_shape_manual("", values = c(Target = 1, Model = 8)) +
    scale_color_manual("", values = c(Target = "black", Model = "red")) +
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
         aes(x = time, y = value, ymin = lb, ymax = ub,
             color = type, shape = type)) +
    geom_point(size = 2) + geom_errorbar() +
    scale_shape_manual("", values = c(Target = 1, Model = 8)) +
    scale_color_manual("", values = c(Target = "black", Model = "red")) +
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
                      aes(x = time, y = value, ymin = lb, ymax = ub,
                          color = type, shape = type)) +
    geom_point(size = 2) + geom_errorbar() +
    scale_shape_manual("", values = c(Target = 1, Model = 8)) +
    scale_color_manual("", values = c(Target = "black", Model = "red")) +
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

#' Plot coverage
#'
#' \code{plot_coverage} plots model outputs vs targets.
#' 
#' @param df_targets Data.frame with calibration targets.
#' @param df_coverage Data.frame with coverage analysis.
#' 
#' @return A ggplot2 object.
#' @export
plot_coverage <- function(df_targets, df_coverage){
  
  # Add type_coverage to targets
  df_targets <- df_targets %>% mutate(type_coverage = "Target")
  
  # Get levels
  v_levels <- c("Target", unique(df_coverage$type_coverage)) 
  
  # Bind data.frames
  df_model_target <- rbind(df_targets, df_coverage)
  df_model_target$type_coverage <- factor(df_model_target$type_coverage, 
                                          levels = v_levels)
  
  # Plot Survival ("Surv")
  plotSurv <- ggplot(subset(df_model_target, abbrev_outcome == "Surv"), 
                     aes(x = as.factor(time), y = value, ymin = lb, ymax = ub,
                         color = type_coverage)) +
    # geom_point(size = 2) + 
    geom_errorbar(position = "dodge", linewidth = 0.8) +
    # scale_shape_manual("", values = c(Target = 1, Model = 8)) +
    ggthemes::scale_color_colorblind("") +
    scale_y_continuous(breaks = number_ticks(6)) +
    # scale_x_continuous(breaks = number_ticks(6)) +
    ylab("Prevalence") + xlab("Time") +
    guides(color = guide_legend(ncol=1, byrow=TRUE)) +
    theme_bw(base_size = 16) +
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
                         color = type_coverage)) +
    # geom_point(size = 2) + 
    geom_errorbar(position = "dodge", linewidth = 0.8) +
    # scale_shape_manual("", values = c(Target = 1, Model = 8)) +
    ggthemes::scale_color_colorblind("") +
    scale_y_continuous(breaks = number_ticks(6)) +
    # scale_x_continuous(breaks = number_ticks(6)) +
    ylab("Prevalence") + xlab("Time") +
    guides(color = guide_legend(ncol=1, byrow=TRUE)) +
    theme_bw(base_size = 16) +
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
                             color = type_coverage)) +
    # geom_point(size = 2) + 
    geom_errorbar(position = "dodge", linewidth = 0.8) +
    # scale_shape_manual("", values = c(Target = 1, Model = 8)) +
    ggthemes::scale_color_colorblind("") +
    scale_y_continuous(breaks = number_ticks(6)) +
    # scale_x_continuous(breaks = number_ticks(6)) +
    ylab("Prop. Sick") + xlab("Time") +
    guides(color = guide_legend(ncol=1, byrow=TRUE)) +
    theme_bw(base_size = 16) +
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
#' @export
plot_costFun <- function(Cfun, wgt_res = TRUE){
  
  df_residuals <- Cfun$residuals
  
  df_residuals$name <- factor(df_residuals$name,
                              levels = c("Surv","Prev","PropSick"))
  
  if(wgt_res){
    
    plotCfun <- ggplot(df_residuals,
           aes(x = x, y = res, color = name, shape = name)) +
      geom_point() +
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
      facet_wrap(~ name, scales = "free_y") +
      scale_y_continuous(breaks = number_ticks(6)) +
      scale_x_continuous(breaks = number_ticks(6)) +
      ylab("Residuals") + xlab("Time") +
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
#' @param df_collin Data.frame containing collinearity indices from \code{collin}.
#' @param log_y Flag (default is TRUE) to whether transform y-axis to log.

#' @return A ggplot2 object.
#' @export
plot_CollIndx <- function(df_collin, log_y = TRUE){
  
  # Run coll
  col_max <- max(df_collin$collinearity)
  
  if(col_max > 15){
    y_end <- col_max + 3
  }else{
    y_end <- 18
  }
  
  df_collin$var 

  if(log_y){
    
    plot_collindx <- ggplot(subset(df_collin, var != "All"), 
                            aes(x = N, y = collinearity)) + 
      geom_point() +
      theme_bw() + 
      facet_wrap(~Outcome, ncol = 3) +
      scale_y_log10(n.breaks = 6) +
      scale_x_continuous(name = "Number of parameters",
                         breaks = df_collin$N, 
                         labels  = df_collin$N) +
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
    
    plot_collindx_all <- ggplot(subset(df_collin, var == "All"), 
                            aes(x = N, y = collinearity)) + 
      geom_point() +
      theme_bw() + 
      facet_wrap(~Outcome, ncol = 1) +
      scale_y_log10(n.breaks = 6) +
      scale_x_continuous(name = "Number of parameters",
                         breaks = df_collin$N, 
                         labels  = df_collin$N) +
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
    
    plot_collindx <- ggplot(subset(df_collin, var != "All"), 
                            aes(x = N, y = collinearity)) + 
      geom_point() +
      theme_bw() + 
      facet_wrap(~Outcome, ncol = 3) +
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
    
    plot_collindx_all <- ggplot(subset(df_collin, var == "All"), 
                            aes(x = N, y = collinearity)) + 
      geom_point() +
      theme_bw() + 
      facet_wrap(~Outcome, ncol = 3) +
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
  
  # Gather plots
  gg_collindx <- ggarrange(plot_collindx, NULL, 
                           plot_collindx_all,
                           heights = c(1,0.05,1), 
                           widths = 1,
                           ncol = 1, nrow = 3)
  
  return(gg_collindx)  
}