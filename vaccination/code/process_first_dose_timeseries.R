



produce_first_dose_timeseries <- function(vacc_dose_data_file,
                                     run_name) {
  
  require(tidyverse)
  require(lubridate)
  
  diag_plot_dir <- paste0("results/", run_name, "/diagnostics")
  
  dir.create(diag_plot_dir, recursive = TRUE)
    
  AIR_ts <- read_csv(vacc_dose_data_file) %>%
    mutate(date = date + 6)
  
  AIR_first_dose_counts <-
    AIR_ts %>%
    group_by(state, date) %>%
    summarise(count = sum(effective_doses))
  
  
  AIR_first_dose_counts_interp <- expand_grid(
    date = seq(min(AIR_first_dose_counts$date), max(AIR_first_dose_counts$date), by = 'day'),
    state = unique(AIR_first_dose_counts$state)
  ) %>%
    left_join(AIR_first_dose_counts) %>%
    group_by(state) %>%
    mutate(count = zoo::na.approx(count)) %>%
    
    mutate(count_new = count - lag(count, default = 0))
  
  
  p1 <- ggplot() +
    geom_line(aes(x = date, y = count),
              AIR_first_dose_counts_interp) +
    
    facet_wrap(~state,
               scale = "free_y") +
    theme_minimal() +
    theme(legend.position = 'bottom')
  
  p2 <- ggplot() +
    geom_line(aes(x = date, y = count_new),
              AIR_first_dose_counts_interp) +
    
    facet_wrap(~state,
               scale = "free_y") +
    theme_minimal() +
    theme(legend.position = 'bottom')
  
  cowplot::plot_grid(p1,p2,ncol = 2)
  
  ggsave(paste0(diag_plot_dir,"/dose_timeseries.png"),
         device = ragg::agg_png, bg = 'white',
         width = 8, height = 4)

  AIR_first_dose_counts_interp

}