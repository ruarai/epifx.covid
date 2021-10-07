

produce_state_summaries  <- function(plotting_data, plot_subdir) {
  
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
  
  seeir_names <- c("S", "E1", "E2", "I1", "I2", "R",
                   "S_U", "E1_U", "E2_U", "I1_U", "I2_U", "R_U",
                   "S_V", "E1_V", "E2_V", "I1_V", "I2_V", "R_V")
  
  param_ordering <- c(seeir_names,
                      "sigma", "gamma", "gen_int",
                      "R_bias", "adjustment",
                      "R0_val", "R0_ix", "t0",
                      "mean_Ei", "mean_Et")
  
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
  
  dir.create(paste0(plot_subdir, "/state_summaries/"))
  
  for(i_state in plotting_data$states){
    p1_b <- make_obs_and_CI_plot(plotting_data$dt_notifications %>% filter(state == i_state),
                               dt_forecasts %>% filter(state == i_state),
                               plotting_data$min_date,
                               plotting_data$fs_date + 4,
                               use_columns = FALSE,
                               probs_plot = c(50,95)) +
      ggtitle("Incidence (50%, 95% credible intervals)") +
      plotting_data$forecast_line
    
    p2_b <- make_obs_and_CI_plot(plotting_data$dt_notifications %>% filter(state == i_state),
                               dt_forecasts %>% filter(state == i_state),
                               plotting_data$min_date,
                               plotting_data$fs_date + 4,
                               use_columns = FALSE,
                               probs_plot = c(50)) +
      ggtitle("Incidence (50% credible interval)") +
      plotting_data$forecast_line
    
    p3_b <- make_parameter_CI_plot(dt_parameter_CIs %>% filter(state == i_state),
                                 i_state,
                                 plotting_data$min_date,
                                 plotting_data$fs_date + 4) +
      ggtitle("Compartment and state values  (50%, 95% credible intervals)") +
      plotting_data$forecast_line
    
    
    p1_f <- make_obs_and_CI_plot(plotting_data$dt_notifications %>% filter(state == i_state),
                                 dt_forecasts %>% filter(state == i_state),
                                 plotting_data$fs_date - 7,
                                 plotting_data$fs_horizon,
                                 use_columns = TRUE,
                                 probs_plot = c(50,95)) +
      ggtitle(paste0("Incidence (50%, 95% credible intervals) ", i_state)) +
      plotting_data$forecast_line
    
    p2_f <- make_obs_and_CI_plot(plotting_data$dt_notifications %>% filter(state == i_state),
                                 dt_forecasts %>% filter(state == i_state),
                                 plotting_data$fs_date - 7,
                                 plotting_data$fs_horizon,
                                 use_columns = TRUE,
                                 probs_plot = c(50)) +
      ggtitle("Incidence (50% credible interval)") +
      plotting_data$forecast_line
    
    p3_f <- make_parameter_CI_plot(dt_parameter_CIs %>% filter(state == i_state),
                                   i_state,
                                   plotting_data$fs_date - 7,
                                   plotting_data$fs_horizon) +
      ggtitle("Compartment and state values (50%, 95% credible intervals)") +
      plotting_data$forecast_line
    
    CairoPDF(paste0(plot_subdir, "/state_summaries/state_highlight_", i_state, ".pdf"), 
             width = 12, height = 24)
    cowplot::plot_grid(p1_b, p2_b, p3_b,
                       ncol = 1,
                       align = 'v', axis = 'lr',
                       rel_heights = c(1,1,6)) %>%
      print()
    
    cowplot::plot_grid(p1_f, p2_f, p3_f,
                       ncol = 1,
                       align = 'v', axis = 'lr',
                       rel_heights = c(1,1,6)) %>%
      print()
    dev.off()
  }
  
  
}
  
