
library(tidyverse)
library(lubridate)
library(rhdf5)
library(Cairo)
library(furrr)


source("forecast_plotting/plotting_common.R")
plotting_data <- get_plotting_data("exps/vic_cf_vacc/")
# 
# 
# dir.create(plot_subdir, showWarnings = FALSE, recursive = TRUE)

plan(callr, workers = 8)

dt_statevec <- plotting_data$hdf_files %>%
  future_pmap(function(file, state) {
    h5read(file = file,
           name = "/data/seeiir")
  })

# %>%
#   pivot_longer(cols = -c(fs_date, date, ix, weight)) %>%
#   select(date, index = ix, weight, name, value) %>%
#   mutate(date = read_hdf_date(date),
#          state = state)

param_ordering <- c("sigma", "gamma",
                    "adjustment",
                    "R_eff", "t0",
                    "mean_Ei", "mean_Et")

param_ordering <- c(param_ordering,
                    "S_U", "S_V", "E1_U", "E1_V", "E2_U", "E2_V",
                    "I1_U", "I1_V", "I2_U", "I2_V", "R_U", "R_V")

dt_statevec_plot <- dt_statevec %>%
  filter(index %% 10 == 0) %>%
  filter(date >= plotting_data$min_date, date <= plotting_data$fs_horizon) %>%
  
  mutate(name = factor(name, levels = param_ordering))



for(i_state in plotting_data$states) {
  p <- ggplot() +
    geom_line(aes(x = date, y = value, group = index),
              size = 0.05, alpha = 0.5,
              dt_statevec_plot %>% filter(state == i_state)) +
    
    facet_wrap(~name,
               ncol = 1,
               scales = "free_y") +
    
    theme_minimal() +
    
    scale_y_log10()
  
  
  ggsave(paste0(plot_subdir, "/trace_params_", i_state, ".png"),
         plot = p, device = ragg::agg_png,
         width = 10, height = 12, dpi = 400,
         bg = 'white')
}





  
  
