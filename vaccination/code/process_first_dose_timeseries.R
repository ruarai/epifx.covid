



produce_first_dose_timeseries <- function(vacc_dose_data_file,
                                     run_name) {
  
  require(tidyverse)
  require(lubridate)
  
  source("vaccination/code/functions_age_classes.R")
  
  diag_plot_dir <- paste0("results/", run_name, "/diagnostics")
  
  dir.create(diag_plot_dir, recursive = TRUE)
    
  AIR_ts <- read_csv(vacc_dose_data_file) %>%
    mutate(date = date + 6)
  
  AIR_first_dose_counts <-
    AIR_ts %>%
    group_by(state, date) %>%
    summarise(count = sum(effective_doses))
  
  
  state_populations <- read_csv("vaccination/data/2021-07-13-census-populations.csv")
  over12_pop <- state_populations %>% 
    filter(age_upper >= 12) %>% 
    group_by(state = ste_name16) %>% 
    summarise(population = sum(population)) %>%
    mutate(state = abbreviate_states(state)) %>%
    drop_na(state)
  
  
  

  AIR_first_dose_counts_interp <- expand_grid(
    date = seq(min(AIR_first_dose_counts$date), max(AIR_first_dose_counts$date), by = 'day'),
    state = unique(AIR_first_dose_counts$state)
  ) %>%
    left_join(AIR_first_dose_counts) %>%
    group_by(state) %>%
    mutate(count = zoo::na.approx(count, na.rm = FALSE))
  
  
  source("vaccination/code/simple_first_dose_projection.R")
  
  AIR_first_dose_counts_interp <- project_first_doses(AIR_first_dose_counts_interp,
                                                      over12_pop,
                                                      diag_plot_dir)
  
  p2 <- ggplot() +
    geom_line(aes(x = date, y = count_new),
              AIR_first_dose_counts_interp) +
    
    facet_wrap(~state,
               scale = "free_y") +
    theme_minimal() +
    theme(legend.position = 'bottom')
  
  ggsave(paste0(diag_plot_dir,"/dose_timeseries.png"),
         plot = p2,
         device = ragg::agg_png, bg = 'white',
         width = 8, height = 8)

  AIR_first_dose_counts_interp

}