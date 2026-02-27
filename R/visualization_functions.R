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
  df_Surv <- df_targets %>% 
    filter(abbrev_outcome == "Surv" & time %in% c(10,20,30))
  df_Prev <- df_targets %>% 
    filter(abbrev_outcome == "Prev" & time %in% c(10,20,30))
  df_PropSick <- df_targets %>% filter(abbrev_outcome == "PropSick")
  
  # Plot Survival ("Surv")
  plotSurv <- ggplot(df_Surv, 
                     aes(x = time, y = value, ymin = lb, ymax = ub)) +
    geom_point(shape = 1) + geom_errorbar() +
    scale_y_continuous(breaks = number_ticks(6)) +
    scale_x_continuous(breaks = number_ticks(6)) +
    ylab("Survival probability") + xlab("Time (years)") +
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
    ylab("Prevalence") + xlab("Time (years)") +
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
    ylab("Proportion of sick") + xlab("Time (years)") +
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
  df_Surv <- rbind(
    df_targets %>% 
      select(time, Outcome, abbrev_outcome, type, value, lb, ub), 
    df_model %>% 
      select(time, Outcome, abbrev_outcome, type, value, lb, ub)
    ) %>% 
    filter(abbrev_outcome == "Surv" & time %in% c(10,20,30))
  
  df_Prev <- rbind(
    df_targets %>% select(time, Outcome, abbrev_outcome, type, value, lb, ub), 
    df_model %>% select(time, Outcome, abbrev_outcome, type, value, lb, ub)
    ) %>% 
    filter(abbrev_outcome == "Prev" & time %in% c(10,20,30))
  
  df_PropSick <- rbind(
    df_targets %>% select(time, Outcome, abbrev_outcome, type, value, lb, ub), 
    df_model %>% select(time, Outcome, abbrev_outcome, type, value, lb, ub)
    ) %>% 
    filter(abbrev_outcome == "PropSick" & time %in% c(10,20,30))
  
  # Order
  df_Surv$type <- factor(df_Surv$type, levels = c("Target","Model"))
  df_Prev$type <- factor(df_Prev$type, levels = c("Target","Model"))
  df_PropSick$type <- factor(df_PropSick$type, levels = c("Target","Model"))
  
  # Plot Survival ("Surv")
  plotSurv <- ggplot(df_Surv, 
         aes(x = factor(time), y = value, ymin = lb, ymax = ub,
             color = type, shape = type)) +
    # geom_point(size = 2) + geom_errorbar() +
    geom_point(size = 1.5, position = position_dodge(0.8)) +
    geom_errorbar(position = position_dodge(0.8), linewidth = 0.8) +
    scale_shape_manual("", values = c(Target = 1, Model = 8)) +
    scale_color_manual("", values = c(Target = "black", Model = "red")) +
    scale_y_continuous(breaks = number_ticks(6)) +
    # scale_x_continuous(breaks = c(10,20,30)) +
    scale_x_discrete(breaks = c("10","20","30")) +
    ylab("Survival probability") + xlab("Time (years)") +
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
         aes(x = factor(time), y = value, ymin = lb, ymax = ub,
             color = type, shape = type)) +
    # geom_point(size = 2) + geom_errorbar() +
    geom_point(size = 1.5, position = position_dodge(0.8)) +
    geom_errorbar(position = position_dodge(0.8), linewidth = 0.8) +
    scale_shape_manual("", values = c(Target = 1, Model = 8)) +
    scale_color_manual("", values = c(Target = "black", Model = "red")) +
    scale_y_continuous(breaks = number_ticks(6)) +
    # scale_x_continuous(breaks = c(10,20,30)) +
    scale_x_discrete(breaks = c("10","20","30")) +
    ylab("Prevalence") + xlab("Time (years)") +
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
                      aes(x = factor(time), y = value, ymin = lb, ymax = ub,
                          color = type, shape = type)) +
    # geom_point(size = 2) + geom_errorbar() +
    geom_point(size = 1.5, position = position_dodge(0.8)) +
    geom_errorbar(position = position_dodge(0.8), linewidth = 0.8) +
    scale_shape_manual("", values = c(Target = 1, Model = 8)) +
    scale_color_manual("", values = c(Target = "black", Model = "red")) +
    scale_y_continuous(breaks = number_ticks(6)) +
    # scale_x_continuous(breaks = c(10,20,30)) +
    scale_x_discrete(breaks = c("10","20","30")) +
    ylab("Proportion of Sick") + xlab("Time (years)") +
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
  v_colors <- v_colors <- c("black","#FF8000","#56B4E9","#009E73","#CCBFFF")
  names(v_colors) <- c(
    "Calibration target",
    "Initial model coverage",                                                           
    "A: Improved model coverage by resizing the probability to become Sicker when Sick",
    "B: A + improved model coverage by resizing the hazard ratio for death in Sicker",  
    "C: A + B + improved model coverage by resizing the hazard ratio for death in Sick"
    )
  
  # Plot Survival ("Surv")
  plotSurv <- ggplot(subset(df_model_target, abbrev_outcome == "Surv"), 
                     aes(x = as.factor(time), 
                         y = value, ymin = lb, ymax = ub,
                         color = type_coverage, shape = type_coverage)) +
    geom_point(size = 1.5, position = position_dodge(0.8)) +
    geom_errorbar(position = position_dodge(0.8), linewidth = 0.8) +
    # ggthemes::scale_color_colorblind("") +
    scale_color_manual("", values = v_colors) +
    scale_shape_manual("", values = v_shape) +
    scale_y_continuous(breaks = number_ticks(6)) +
    scale_x_discrete(breaks = c("10","20","30")) +
    ylab("Survival probability") + xlab("Time (years)") +
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
                         color = type_coverage, shape = type_coverage)) +
    geom_point(size = 1.5, position = position_dodge(0.8)) +
    geom_errorbar(position = position_dodge(0.8), linewidth = 0.8) +
    # ggthemes::scale_color_colorblind("") +
    scale_color_manual("", values = v_colors) +
    scale_shape_manual("", values = v_shape) +
    scale_y_continuous(breaks = number_ticks(6)) +
    scale_x_discrete(breaks = c("10","20","30")) +
    ylab("Prevalence") + xlab("Time (years)") +
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
                             color = type_coverage, shape = type_coverage)) +
    geom_point(size = 1.5, position = position_dodge(0.8)) +
    geom_errorbar(position = position_dodge(0.8), linewidth = 0.8) +
    # ggthemes::scale_color_colorblind("") +
    scale_color_manual("", values = v_colors) +
    scale_shape_manual("", values = v_shape) +
    scale_y_continuous(breaks = number_ticks(6)) +
    scale_x_discrete(breaks = c("10","20","30")) +
    ylab("Proportion of Sick") + xlab("Time (years)") +
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
      scale_x_continuous(breaks = c(10,20,30)) +
      ylab("Weighted residuals") + xlab("Time (years)") +
      theme_bw(base_size = 16) +
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
      scale_x_continuous(breaks = c(10,20,30)) +
      ylab("Residuals") + xlab("Time (years)") +
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
    scale_x_continuous(breaks = c(10,20,30)) +
    ylab("Sensitivity") + xlab("Time (years)") +
    guides(color = guide_legend(ncol=1, byrow=TRUE)) +
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
  collin_max <- max(df_collin$collinearity)
  
  # y-axis upper limit
  if(collin_max > 15){
    y_end <- collin_max + 3
  }else{
    y_end <- 18
  }
  
  # x-axis discrete values
  v_x_breaks <- unique(df_collin$N)

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
      theme_bw(base_size = 16) +
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
      theme_bw(base_size = 16) +
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
      theme_bw(base_size = 16) +
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
      theme_bw(base_size = 16) +
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
  
  # Gather plots
  gg_collindx <- ggarrange(plot_collindx, NULL, 
                           plot_collindx_all,
                           heights = c(1,0.05,1), 
                           widths = 1,
                           ncol = 1, nrow = 3)
  
  return(gg_collindx)  
}