

produce_CI_plots <- function(plotting_data, plot_subdir) {
  
  require(tidyverse)
  require(lubridate)
  require(rhdf5)
  require(Cairo)
  require(cowplot)
  
  dt_forecasts <- plotting_data$hdf_files %>%
    pmap_dfr(function(file, state) {
      h5read(file = file,
             name = "/data/forecasts") %>%
        select(date, prob, ymin, ymax) %>%
        mutate(date = read_hdf_date(date),
               prob = factor(prob, levels = c(95, 50, 0)),
               state = state)
    })
  
  param_ordering <- c("sigma", "gamma", "gen_int",
                      "adjustment",
                      "R0_val", "R0_ix", "t0",
                      "mean_Ei", "mean_Et")
  
  param_ordering <- c(param_ordering,
                      "S_U", "S_V", "E1_U", "E1_V", "E2_U", "E2_V",
                      "I1_U", "I1_V", "I2_U", "I2_V", "R_U", "R_V")
  
  dt_parameter_CIs <- plotting_data$hdf_files %>%
    pmap_dfr(function(file, state) {
      h5read(file = file,
             name = "/data/model_cints") %>%
        select(date, prob, ymin, ymax, type, name) %>%
        filter(name != "R0") %>%
        mutate(date = read_hdf_date(date),
               prob = factor(prob, levels = c(100, 99, 95, 90, 50, 0)),
               name = factor(name, levels = param_ordering),
               state = state)
    })
  
  CairoPDF(paste0(plot_subdir, "/incidence.pdf"), 
           width = 12, height = 9)
  p <- make_obs_and_CI_plot(plotting_data$dt_notifications,
                            dt_forecasts,
                            plotting_data$min_date,
                            plotting_data$fs_date + 4,
                            use_columns = FALSE,
                            probs_plot = c(50,95)) +
    ggtitle("Backcast case counts (50% and 95% credible intervals)") +
    plotting_data$forecast_line
  print(p)

  p <- make_obs_and_CI_plot(plotting_data$dt_notifications,
                            dt_forecasts,
                            plotting_data$fs_date - 35,
                            plotting_data$fs_date + 4,
                            use_columns = TRUE,
                            probs_plot = c(50,95)) +
    ggtitle("Recent backcast (50% and 95% credible intervals)") +
    plotting_data$forecast_line
  print(p)
  
  p <- make_obs_and_CI_plot(plotting_data$dt_notifications,
                            dt_forecasts,
                            plotting_data$fs_date - 7,
                            plotting_data$fs_horizon,
                            probs_plot = c(50,95)) +
    ggtitle("Projected case counts (50% and 95% credible intervals)") +
    plotting_data$forecast_line
  print(p)
  
  p <- make_obs_and_CI_plot(plotting_data$dt_notifications,
                            dt_forecasts,
                            plotting_data$fs_date - 7,
                            plotting_data$fs_horizon,
                            probs_plot = c(50)) +
    ggtitle("Projected case counts (50% credible interval)") +
    plotting_data$forecast_line
  print(p)
  
  dev.off()
  
  
  CairoPDF(paste0(plot_subdir, "/parameter_estimates.pdf"), 
           width = 12, height = 9)
  
  for(i_name in param_ordering) {
    
    p1 <- make_parameter_CI_plot_by_name(dt_parameter_CIs,
                                         i_name,
                                         plotting_data$min_date,
                                         plotting_data$fs_date + 4) +
      ggtitle(paste0(i_name, " (50% and 95% credible intervals), Backcast")) +
      plotting_data$forecast_line
    
    p2 <- make_parameter_CI_plot_by_name(dt_parameter_CIs,
                                         i_name,
                                         plotting_data$fs_date - 7,
                                         plotting_data$fs_horizon) +
      ggtitle(paste0("Forecast")) +
      plotting_data$forecast_line
    
    cowplot::plot_grid(p1, p2,
                       align = 'h', axis = 'tb',
                       rel_widths = c(2,1)) %>%
      print()
  }
  
  dev.off()
  
  
  
  
}
